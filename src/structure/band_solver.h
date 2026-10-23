#pragma once

// 대칭 스카이라인(열 높이 가변 띠) 행렬과 LDLᵀ 직접 해법.
// 열 j 는 행 top(j) .. j 의 값을 저장한다 (위 삼각). 분해 뒤 같은 자리에 L 의 열과 D 가 들어간다.

#include <vector>

namespace hcw::structure {

class SkylineMatrix {
 public:
  SkylineMatrix() = default;
  // top[j] ≤ j: 열 j 의 가장 위 0 이 아닌 행
  explicit SkylineMatrix(const std::vector<int>& top);

  int size() const { return static_cast<int>(top_.size()); }
  std::size_t stored() const { return values_.size(); }
  // i ≤ j 이고 i ≥ top(j) 인 자리에 더한다 (그 밖은 예외)
  void add(int i, int j, double v);
  double get(int i, int j) const;

  // 제자리 LDLᵀ 분해. 0 에 가까운 피벗(|d| ≤ 1e−14 |a_jj|)이면 예외.
  void factor();
  bool factored() const { return factored_; }
  // 분해된 행렬로 A x = b 를 푼다 (b 를 x 로 덮어씀)
  void solve(std::vector<double>& b) const;

 private:
  std::size_t at(int i, int j) const { return start_[j] + static_cast<std::size_t>(i - top_[j]); }
  std::vector<int> top_;
  std::vector<std::size_t> start_;
  std::vector<double> values_;
  bool factored_ = false;
};

}  // namespace hcw::structure
