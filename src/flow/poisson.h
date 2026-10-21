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
  int iterations = 0;
  double relative_residual = 0.0;
};

// x 는 초기 추정값으로 쓰이고 해로 덮어쓴다. 허용치는 ||r|| / ||b||.
SolveResult pcg(const CsrMatrix& a, const IncompleteCholesky& m, const std::vector<double>& b,
                std::vector<double>& x, double tolerance, int max_iterations);

}  // namespace hcw::flow
