#pragma once

// 다중 코어 실행 (선택 기능, 기본은 1 스레드).
//
// 원칙: 스레드 수와 상관없이 결과가 비트 단위로 같아야 한다.
// 그래서 여기서 나누는 일은 "칸마다 따로 계산되는" 일뿐이다. 각 칸(행, 면)의 계산 순서와 식은
// 1 스레드일 때와 똑같고, 어느 스레드가 맡는지만 달라진다. 합(내적, 노름 등)처럼 더하는 순서가
// 결과를 바꾸는 계산은 이 모듈로 나누지 않고 한 스레드가 원래 순서대로 한다.
//
// 켜는 방법: 환경 변수 HCW_THREADS=<수> 또는 hcw_verify --threads <수>. 없으면 1.

#include <exception>
#include <string>
#include <utility>
#include <vector>

namespace hcw::parallel {

// 정한 스레드 수 (1 이면 모든 계산이 원래의 직렬 코드로 돈다)
int threads();

// 지금 나눠 쓸 수 있는 스레드 수. 병렬 구간 안(작업 스레드 위)에서 부르면 1 이다
// (병렬 구간 안에서 다시 나누지 않고 그 자리에서 직렬로 돈다).
int available();

// 스레드 수를 정한다 (1 이상). 작업 스레드를 새로 만든다. 계산 도중에 부르지 않는다.
void set_threads(int n);

// HCW_THREADS 를 읽는다. 없으면 1. 잘못된 값이면 0 을 돌려주고 error 에 이유를 쓴다.
int threads_from_environment(std::string& error);

// "1 thread" 또는 "N threads (results bitwise identical to 1 thread)" 같은 기록용 문구
std::string describe();

namespace detail {
// fn(ctx, tid) 를 tid = 0 .. threads()-1 에서 동시에 한 번씩 부르고 모두 끝날 때까지 기다린다.
// tid 0 은 부른 스레드 자신이다. fn 은 예외를 던지면 안 된다.
void run(void (*fn)(void*, int), void* ctx);
// 바쁜 대기용 (CPU 에 쉬라고 알림)
void pause();
}  // namespace detail

// f(tid) 를 모든 스레드에서 한 번씩 부른다.
// 병렬 구간 안에서 부르면 f(0) 만 부른다 (available() 참조).
template <class F>
void for_each_thread(F&& f) {
  if (available() <= 1) {
    f(0);
    return;
  }
  detail::run([](void* c, int tid) { (*static_cast<F*>(c))(tid); }, &f);
}

// 이보다 짧은 반복은 나누지 않는다 (나누는 비용이 더 큼). 결과에는 영향이 없다.
constexpr int kMinParallelRange = 2048;
// 격자 한 줄(열 i 의 모든 j)을 한 인덱스로 나눌 때의 최소 줄 수
constexpr int kRowGrain = 16;

// [begin, end) 를 스레드 수만큼 연속 구간으로 나눠 f(b, e) 를 부른다.
// f 안의 각 인덱스 계산은 다른 인덱스의 결과를 읽지 않아야 한다 (그래야 결과가 같다).
// min_range: 이보다 짧으면 나누지 않는다 (인덱스 하나의 일이 많으면, 예: 격자 한 줄, 작게 준다).
template <class F>
void parallel_for(int begin, int end, int min_range, F&& f) {
  const int t = available();
  if (t <= 1 || end - begin < min_range) {
    if (end > begin) f(begin, end);
    return;
  }
  const long long len = end - begin;
  // 예외는 구간마다 받아 두었다가 가장 앞 구간의 것을 다시 던진다. 구간 안은 원래 순서로 돌므로
  // 1 스레드일 때 처음 던져지는 예외와 같다.
  std::vector<std::exception_ptr> errors(static_cast<std::size_t>(t));
  for_each_thread([&](int tid) {
    const int b = begin + static_cast<int>(len * tid / t);
    const int e = begin + static_cast<int>(len * (tid + 1) / t);
    try {
      if (e > b) f(b, e);
    } catch (...) {
      errors[static_cast<std::size_t>(tid)] = std::current_exception();
    }
  });
  for (const auto& err : errors) {
    if (err) std::rethrow_exception(err);
  }
}

template <class F>
void parallel_for(int begin, int end, F&& f) {
  parallel_for(begin, end, kMinParallelRange, std::forward<F>(f));
}

}  // namespace hcw::parallel
