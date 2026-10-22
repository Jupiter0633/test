#pragma once

// P3: 기하학적 VOF 의 경계면 재구성 (PLIC, 셀마다 직선 하나).
// 물 영역 = { x : n·x ≤ α } (n: 물에서 공기 쪽 단위 법선). 셀 체적은 축대칭이면 r 가중(∫∫ r dr dz, 2π 제외),
// 평면이면 면적. 볼록 다각형의 r 가중 면적은 꼭짓점 공식으로 정확히 구한다.
// 근거: docs/design/02_design_document.md §3.2, 03_verification_plan.md §1.3

#include <array>
#include <vector>

#include "flow/grid.h"

namespace hcw::interface {

struct Point {
  double r = 0.0, z = 0.0;
};

using Polygon = std::vector<Point>;

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
