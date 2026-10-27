#pragma once

// 대칭 양의 (준)정부호 희소 행렬과 불완전 Cholesky(IC0) 전처리 켤레기울기법.
// P1 시험 격자(최대 256²)용. 큰 격자용 다중격자 해법은 이후 단계에서 교체한다.
//
// 다중 스레드 (common/parallel.h): 행렬–벡터 곱과 벡터 갱신은 행마다 나눈다. 내적·노름은 원래 순서
// 그대로 한 스레드가 더한다. 불완전 Cholesky 의 삼각 풀이는 행 순서(자연 순서)를 바꾸지 않고
// 파이프라인(파면, wavefront)으로 나눈다 (IncompleteCholesky::factor 참조). 그래서 모든 값이
// 스레드 수와 상관없이 1 스레드와 비트 단위로 같고, 반복 횟수도 같다.

#include <atomic>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

namespace hcw::flow {

struct CsrMatrix {
  int n = 0;
  std::vector<int> row_start;  // 크기 n+1
  std::vector<int> col;        // 행마다 열 번호 오름차순
  std::vector<double> val;

  void multiply(const std::vector<double>& x, std::vector<double>& y) const;
};

class IncompleteCholesky {
 public:
  // 실패(비양수 피벗) 시 std::runtime_error.
  // inner > 0: 미지수가 격자 순서 c = i·inner + j (j 가 빠른 쪽) 임을 알려 준다. 스레드가 2 개
  // 이상이면 삼각 풀이를 파이프라인으로 나눈다: j 를 구간(segment) s 로 나눠 스레드 s 가 맡고,
  // 열 i 의 구간 s 는 (i, s−1) 이 끝난 뒤 푼다 (뒤 풀이는 반대 방향). 각 행의 식과 더하는 순서는
  // 직렬 풀이와 똑같다. 행렬의 의존 관계가 이 순서로 다 지켜지는지 인수분해 때 검사하고,
  // 아니면 (또는 격자가 작으면) 직렬 풀이를 쓴다.
  void factor(const CsrMatrix& a, int inner = 0);
  void apply(const std::vector<double>& r, std::vector<double>& z) const;
  // side() 를 삼각 풀이와 동시에 (파이프라인이 쓰지 않는 스레드에서, 없으면 풀이 전에) 부른다.
  // side 는 r 을 읽기만 하고 z 와 이 객체는 건드리지 않아야 한다.
  void apply(const std::vector<double>& r, std::vector<double>& z,
             const std::function<void()>& side) const;
  // 파이프라인 구간 수 (0 이면 직렬 풀이)
  int pipeline_segments() const { return seg_; }

 private:
  void plan_pipeline(int inner);
  void run(const std::vector<double>& r, std::vector<double>& z,
           const std::function<void()>* side) const;
  void apply_serial(const std::vector<double>& r, std::vector<double>& z) const;
  void apply_segment(int s, std::uint32_t epoch, const std::vector<double>& r,
                     std::vector<double>& z) const;

  // 구간 완료 표시 (복사하면 새 0 배열)
  struct Flags {
    std::unique_ptr<std::atomic<std::uint32_t>[]> a;
    std::size_t n = 0;
    Flags() = default;
    Flags(const Flags& o) { reset(o.n); }
    Flags& operator=(const Flags& o) {
      if (this != &o) reset(o.n);
      return *this;
    }
    Flags(Flags&&) = default;
    Flags& operator=(Flags&&) = default;
    void reset(std::size_t size);
    // 표시 하나에 캐시 줄(64 바이트) 하나: 다른 스레드의 표시와 같은 줄을 쓰면 서로 방해한다
    static constexpr std::size_t kStride = 16;
    std::atomic<std::uint32_t>& at(std::size_t k) const { return a[k * kStride]; }
  };

  int n_ = 0;
  std::vector<int> row_start_;  // 하삼각(대각 포함) 패턴
  std::vector<int> col_;
  std::vector<double> val_;
  std::vector<double> diag_;
  mutable std::vector<double> work_;

  // 파이프라인 (seg_ > 0 일 때)
  int seg_ = 0, inner_ = 0, ncol_ = 0;
  std::vector<int> seg_begin_;            // 크기 seg_+1, 열 안의 j 구간 경계
  std::vector<int> lt_start_, lt_row_;    // Lᵀ: 열 j 의 (i > j) 를 i 내림차순으로 (뒤 풀이 모음 형태)
  std::vector<double> lt_val_;
  mutable Flags fwd_done_, bwd_done_;     // (열, 구간) 마다 끝난 apply 번호
  mutable std::uint32_t epoch_ = 0;
};

struct SolveResult {
  bool converged = false;
  int iterations = 0;          // 켤레기울기 반복 횟수 (보정 해법은 내부 반복의 합)
  double relative_residual = 0.0;
  int refinements = 0;         // 보정 해법의 바깥 반복 횟수
  bool stagnated = false;      // 잔차가 더 줄지 않아 조기 종료
};

// x 는 초기 추정값으로 쓰이고 해로 덮어쓴다. 허용치는 ||r|| / ||b||.
SolveResult pcg(const CsrMatrix& a, const IncompleteCholesky& m, const std::vector<double>& b,
                std::vector<double>& x, double tolerance, int max_iterations);

// 반복 보정 해법. 해를 두 배 정밀도(x = hi + lo)로 들고, 잔차 b − A(hi + lo) 를
// 오차 없는 변환(TwoSum, FMA)으로 계산한다. 보정량은 pcg 로 느슨하게(내부 허용치) 푼다.
// 배정밀도 하나로는 격자가 고울수록 잔차가 반올림 한계(≈ ε·|A||x|)에서 멈춰
// 1e-12 같은 허용치에 닿지 못하기 때문이다. 판정 잔차는 위 정밀 잔차의 ||·|| / ||b||.
// 잔차가 더 줄지 않으면 max_iterations 까지 기다리지 않고 실패(stagnated)로 끝낸다.
// hi, lo 는 초기 추정값으로 쓰이고 해로 덮어쓴다 (lo 가 비어 있으면 0).
SolveResult refined_solve(const CsrMatrix& a, const IncompleteCholesky& m,
                          const std::vector<double>& b, std::vector<double>& hi,
                          std::vector<double>& lo, double tolerance, int max_iterations);

}  // namespace hcw::flow
