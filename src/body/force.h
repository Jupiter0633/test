#pragma once

// P2: 물체에 작용하는 유체력 (표면 적분).
// 근거: docs/design/03_verification_plan.md §1.2 (힘 = 물체 표면 유체력의 축 방향 성분)

#include <functional>

#include "body/cut_cell.h"
#include "flow/grid.h"

namespace hcw::body {

// 축 방향 압력 힘 F_z = −∮ p n_out dA (n_out: 물체 바깥 법선). 축대칭은 2π 를 곱한 3차원 힘.
// 경계 조각마다 셀 중심 압력을 주변 유체 셀의 최소제곱 기울기로 대표점까지 1차 외삽한다.
// 압력의 상수 차이는 닫힌 물체에서 힘에 영향을 주지 않는다 (Σ nz_area = 0).
// geom.active 가 거짓인 셀은 압력이 없는 것으로 보고, 그런 셀의 조각은 주변 유효 셀의 1차식으로 구한다.
double axial_pressure_force(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p);

// 조각 하나의 대표점 압력 (axial_pressure_force 와 같은 외삽). 개발용 진단에서 조각별 기여를 볼 때 쓴다.
double piece_pressure(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p, const BoundaryPiece& piece);

// 축 방향 점성 힘 (no-slip 강체 표면). 벽에서 점성 견인력 = μ ∂u/∂n (n: 물체 바깥 법선) 이므로
// z 성분은 μ ∂u_z/∂n. 경계 조각 대표점 근처(3Δ 안)의 유체 쪽 z 면 속도로
// u_z − W ≈ d (c1 + c2 t + c3 d + c4 d² + c5 t²) (d: 표면까지 거리, t: 접선 거리) 를 최소제곱으로 맞추고
// 벽 기울기 c1 을 쓴다. 이 꼴은 실제 벽 위에서 u_z = W (no-slip) 를 만족한다. 물체 안의 면은 쓰지 않는다.
// 축대칭은 2π 를 곱한다.
// 경계 조각이 격자면과 겹치는 경우(격자면에 맞는 물체의 옆면)는 다루지 않는다.
double axial_viscous_force(const flow::Grid& g, const CutGeometry& geom, const Shape& shape,
                           double z_center, const flow::Array2D& w, double mu, double body_w = 0.0);

// 경계 조각 중 select 가 참인 것들의 압력 면적평균 Σ p n_z dA / Σ n_z dA (P4 바닥면·센서 압력).
// 압력은 axial_pressure_force 와 같은 1차 외삽. 고른 조각의 Σ n_z dA 가 0 이면 0 을 돌려준다.
double average_surface_pressure(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p,
                                const std::function<bool(const BoundaryPiece&)>& select);

// 위와 같되 조각 길이 × r 가중 평균 (법선이 수평인 옆벽 압력, 3-I)
double average_surface_pressure_by_length(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p,
                                          const std::function<bool(const BoundaryPiece&)>& select);

// 경계 조각 대표점 압력 (1차 외삽) 의 최댓값 (P4 한 점 최대 압력, 3-D)
double max_surface_pressure(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p);

// 제어면 운동량 수지로 구한 물체의 축 방향 힘 (정상 상태, 축대칭, 2π 포함).
// 닫힌 면: 옆면 r = rs, 아랫면 z = z_bottom, 윗면 z = z_top (모두 격자면에 맞아야 한다).
//   F_z = ∮ [ −p n_z + (τ·n)_z − ρ u_z (u·n) ] dA
// 윗·아랫면: p 는 이웃 두 셀 평균, τ_zz = 2μ ∂w/∂z (중심 차분). 옆면: τ_rz = μ(∂u/∂z + ∂w/∂r) (중심 차분).
// 면과 물체 사이에 다른 힘이 없으면 물체 표면의 유체력과 같다. 근거: 03_verification_plan.md §1.2 1-G.
double axial_force_control_surface(const flow::Grid& g, const flow::Array2D& u,
                                   const flow::Array2D& w, const flow::Array2D& p, double mu,
                                   double rho, double rs, double z_bottom, double z_top);

}  // namespace hcw::body
