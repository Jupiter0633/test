#pragma once

// 비대칭 희소 연립방정식: ILU(0) 전처리 BiCGSTAB.
// P2 의 암시적 점성 운동량 식(물체 경계 보간 행 때문에 비대칭)에 쓴다.

#include <vector>

#include "flow/poisson.h"

namespace hcw::flow {

class Ilu0 {
 public:
  // 실패(0 피벗) 시 std::runtime_error. 행마다 열 번호는 오름차순이어야 한다.
  void factor(const CsrMatrix& a);
  void apply(const std::vector<double>& r, std::vector<double>& z) const;

 private:
  int n_ = 0;
  std::vector<int> row_start_, col_, diag_pos_;
  std::vector<double> val_;
};

// x 는 초기 추정값으로 쓰이고 해로 덮어쓴다. 허용치는 참 잔차 ||b − Ax|| / ||b||.
// 잔차가 더 줄지 않으면(정체) 일찍 끝내고 converged = false.
SolveResult bicgstab(const CsrMatrix& a, const Ilu0& m, const std::vector<double>& b,
                     std::vector<double>& x, double tolerance, int max_iterations);

}  // namespace hcw::flow
