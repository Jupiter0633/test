#pragma once

// 대칭 양의 (준)정부호 희소 행렬과 불완전 Cholesky(IC0) 전처리 켤레기울기법.
// P1 시험 격자(최대 256²)용. 큰 격자용 다중격자 해법은 이후 단계에서 교체한다.

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
  void factor(const CsrMatrix& a);
  void apply(const std::vector<double>& r, std::vector<double>& z) const;

 private:
  int n_ = 0;
  std::vector<int> row_start_;  // 하삼각(대각 포함) 패턴
  std::vector<int> col_;
  std::vector<double> val_;
  std::vector<double> diag_;
  mutable std::vector<double> work_;
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
