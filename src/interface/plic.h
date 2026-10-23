#pragma once

// P3: 기하학적 VOF 의 경계면 재구성 (PLIC, 셀마다 직선 하나).
// 물 영역 = { x : n·x ≤ α } (n: 물에서 공기 쪽 단위 법선). 셀 체적은 축대칭이면 r 가중(∫∫ r dr dz, 2π 제외),
// 평면이면 면적. 볼록 다각형의 r 가중 면적은 꼭짓점 공식으로 정확히 구한다.
// 근거: docs/design/02_design_document.md §3.2, 03_verification_plan.md §1.3

#include <array>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>

#include "flow/grid.h"

namespace hcw::interface {

struct Point {
  double r = 0.0, z = 0.0;
};

// 볼록 다각형 (꼭짓점 최대 kMax 개). 직사각형을 반평면으로 한 번 자르면 최대 5 개다.
// 힙 메모리를 쓰지 않는다: solve_alpha 의 이분법이 셀마다 수백 번 자르므로, std::vector 로는
// 여러 스레드가 힙을 두고 다퉈 오히려 느려졌다 (2026-10-23 확인). 계산식은 그대로라 결과도 같다.
class Polygon {
 public:
  static constexpr std::size_t kMax = 8;
  Polygon() = default;
  Polygon(std::initializer_list<Point> pts) {
    for (const Point& p : pts) push_back(p);
  }
  void push_back(const Point& p) {
    if (n_ >= kMax) throw std::logic_error("PLIC polygon: too many vertices");
    p_[n_++] = p;
  }
  std::size_t size() const { return n_; }
  bool empty() const { return n_ == 0; }
  const Point& operator[](std::size_t k) const { return p_[k]; }
  const Point* begin() const { return p_.data(); }
  const Point* end() const { return p_.data() + n_; }

 private:
  std::array<Point, kMax> p_{};
  std::size_t n_ = 0;
};

// 직사각형 [r0, r1] × [z0, z1] (반시계 방향 꼭짓점)
Polygon rectangle(double r0, double r1, double z0, double z1);

// 볼록 다각형을 반평면 n·x ≤ α 로 자른다
Polygon clip(const Polygon& poly, std::array<double, 2> n, double alpha);

// 다각형 체적: 축대칭이면 ∫∫ r dA, 평면이면 ∫∫ dA
double polygon_volume(const Polygon& poly, flow::Geometry geometry);

// 직사각형 안 물 체적이 target 이 되는 α (이분법, 기계 정밀도까지)
double solve_alpha(const Polygon& rect, std::array<double, 2> n, double target,
                   flow::Geometry geometry);

// Youngs 방법 법선: n = −∇f / |∇f| (3×3 가중 차분). 기울기가 0 이면 (0, 1) (수면이 위).
// 영역 밖 이웃: 축(r = 0)은 거울 대칭, 나머지는 가장자리 값을 복제.
std::array<double, 2> youngs_normal(const flow::Grid& g, const flow::Array2D& f, int i, int j);

}  // namespace hcw::interface
