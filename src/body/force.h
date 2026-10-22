#pragma once

// P2: 물체에 작용하는 유체력 (표면 적분).
// 근거: docs/design/03_verification_plan.md §1.2 (힘 = 물체 표면 유체력의 축 방향 성분)

#include "body/cut_cell.h"
#include "flow/grid.h"

namespace hcw::body {

// 축 방향 압력 힘 F_z = −∮ p n_out dA (n_out: 물체 바깥 법선). 축대칭은 2π 를 곱한 3차원 힘.
// 경계 조각마다 셀 중심 압력을 주변 유체 셀의 최소제곱 기울기로 대표점까지 1차 외삽한다.
// 압력의 상수 차이는 닫힌 물체에서 힘에 영향을 주지 않는다 (Σ nz_area = 0).
// geom.active 가 거짓인 셀은 압력이 없는 것으로 보고, 그런 셀의 조각은 주변 유효 셀의 1차식으로 구한다.
double axial_pressure_force(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p);

// 축 방향 점성 힘 (no-slip 강체 표면). 벽에서 점성 견인력 = μ ∂u/∂n (n: 물체 바깥 법선) 이므로
// z 성분은 μ ∂u_z/∂n. 경계 조각 대표점 근처(3Δ 안)의 유체 쪽 z 면 속도로
// u_z − W ≈ d (c1 + c2 t + c3 d + c4 d² + c5 t²) (d: 표면까지 거리, t: 접선 거리) 를 최소제곱으로 맞추고
// 벽 기울기 c1 을 쓴다. 이 꼴은 실제 벽 위에서 u_z = W (no-slip) 를 만족한다. 물체 안의 면은 쓰지 않는다.
// 축대칭은 2π 를 곱한다.
// 경계 조각이 격자면과 겹치는 경우(격자면에 맞는 물체의 옆면)는 다루지 않는다.
double axial_viscous_force(const flow::Grid& g, const CutGeometry& geom, const Shape& shape,
                           double z_center, const flow::Array2D& w, double mu, double body_w = 0.0);

}  // namespace hcw::body
