#pragma once

// P2: 물체에 작용하는 유체력 (표면 적분).
// 근거: docs/design/03_verification_plan.md §1.2 (힘 = 물체 표면 유체력의 축 방향 성분)

#include "body/cut_cell.h"
#include "flow/grid.h"

namespace hcw::body {

// 축 방향 압력 힘 F_z = −∮ p n_out dA (n_out: 물체 바깥 법선). 축대칭은 2π 를 곱한 3차원 힘.
// 경계 조각마다 셀 중심 압력을 주변 유체 셀의 최소제곱 기울기로 대표점까지 1차 외삽한다.
// 압력의 상수 차이는 닫힌 물체에서 힘에 영향을 주지 않는다 (Σ nz_area = 0).
double axial_pressure_force(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p);

}  // namespace hcw::body
