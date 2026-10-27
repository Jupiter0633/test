#include "common/parallel.h"

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#if defined(_M_X64) || defined(__x86_64__)
#include <immintrin.h>
#endif

namespace hcw::parallel {
namespace {

thread_local bool t_inside = false;  // 이 스레드가 병렬 구간의 일을 하는 중

// fn 을 병렬 구간 표시와 함께 부른다
void call_inside(void (*fn)(void*, int), void* ctx, int tid) {
  t_inside = true;
  fn(ctx, tid);
  t_inside = false;
}

// 작업 스레드 묶음. 압력 해법은 한 반복에 여러 번 짧은 병렬 구간을 돌므로
// 스레드를 매번 만들지 않고, 잠깐은 바쁜 대기로 다음 일을 기다린다.
// 일이 한동안(kSpinTime) 없으면 잠들어 CPU 를 놓아 준다.
class Pool {
 public:
  explicit Pool(int n) : n_(n) {
    for (int tid = 1; tid < n_; ++tid) workers_.emplace_back([this, tid] { loop(tid); });
  }
  ~Pool() {
    stop_.store(true);
    generation_.fetch_add(1);
    {
      std::lock_guard<std::mutex> lk(m_);
    }
    cv_.notify_all();
    for (auto& w : workers_) w.join();
  }
  int size() const { return n_; }

  void run(void (*fn)(void*, int), void* ctx) {
    fn_ = fn;
    ctx_ = ctx;
    pending_.store(n_ - 1, std::memory_order_relaxed);
    generation_.fetch_add(1, std::memory_order_seq_cst);
    if (sleeping_.load(std::memory_order_seq_cst) > 0) {
      { std::lock_guard<std::mutex> lk(m_); }
      cv_.notify_all();
    }
    call_inside(fn, ctx, 0);
    while (pending_.load(std::memory_order_acquire) != 0) detail::pause();
  }

 private:
  static constexpr std::chrono::microseconds kSpinTime{2000};

  void loop(int tid) {
    std::uint64_t seen = 0;
    for (;;) {
      auto start = std::chrono::steady_clock::now();
      int spins = 0;
      while (generation_.load(std::memory_order_acquire) == seen) {
        detail::pause();
        if (++spins % 256 != 0) continue;
        if (std::chrono::steady_clock::now() - start < kSpinTime) continue;
        std::unique_lock<std::mutex> lk(m_);
        sleeping_.fetch_add(1, std::memory_order_seq_cst);
        while (generation_.load(std::memory_order_seq_cst) == seen) cv_.wait(lk);
        sleeping_.fetch_sub(1, std::memory_order_seq_cst);
      }
      seen = generation_.load(std::memory_order_acquire);
      if (stop_.load()) return;
      call_inside(fn_, ctx_, tid);
      pending_.fetch_sub(1, std::memory_order_release);
    }
  }

  int n_;
  std::vector<std::thread> workers_;
  std::atomic<std::uint64_t> generation_{0};
  std::atomic<int> pending_{0};
  std::atomic<int> sleeping_{0};
  std::atomic<bool> stop_{false};
  void (*fn_)(void*, int) = nullptr;
  void* ctx_ = nullptr;
  std::mutex m_;
  std::condition_variable cv_;
};

std::unique_ptr<Pool> g_pool;
int g_threads = 1;

}  // namespace

int threads() { return g_threads; }

int available() { return t_inside ? 1 : g_threads; }

void set_threads(int n) {
  if (n < 1) n = 1;
  if (n == g_threads) return;
  g_pool.reset();
  g_threads = n;
  if (n > 1) g_pool = std::make_unique<Pool>(n);
}

int threads_from_environment(std::string& error) {
  const char* v = std::getenv("HCW_THREADS");
  if (!v || !*v) return 1;
  char* end = nullptr;
  const long n = std::strtol(v, &end, 10);
  if (*end != '\0' || n < 1 || n > 256) {
    error = std::string("HCW_THREADS must be an integer from 1 to 256 (got '") + v + "')";
    return 0;
  }
  return static_cast<int>(n);
}

std::string describe() {
  if (g_threads <= 1) return "1 (single-threaded)";
  return std::to_string(g_threads) + " (results bitwise identical to 1 thread)";
}

namespace detail {

void run(void (*fn)(void*, int), void* ctx) {
  if (!g_pool) {
    fn(ctx, 0);
    return;
  }
  g_pool->run(fn, ctx);
}

void pause() {
#if defined(_M_X64) || defined(__x86_64__)
  _mm_pause();
#else
  std::this_thread::yield();
#endif
}

}  // namespace detail
}  // namespace hcw::parallel
