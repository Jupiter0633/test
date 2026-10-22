#pragma once

// P2: 고정 격자 위에서 물체가 차지하는 부분 (컷셀 기하).
// 면마다 "열린 비율"(유체 쪽 면적 / 전체 면적, 축대칭은 r 가중 면적)을 형상에서 정확히 구한다.
// 압력(투영) 방정식은 이 비율을 곱한 면 플럭스로 쌓는다 → 물체 경계는 날카롭게(셀 안에서) 표현된다.
// 근거: docs/design/02_design_document.md §3.4 (날카로운 경계 IBM), 03_verification_plan.md §1.2

#include <vector>

#include "body/shape.h"
#include "flow/grid.h"

namespace hcw::body {

// 셀 안의 물체 경계 조각. 축 방향 힘과 강체 병진의 경계 플럭스에 쓴다.
struct BoundaryPiece {
  int i = 0, j = 0;
  // ∫ n_z dA (r 가중, 2π 제외). n 은 유체에서 물체 쪽(유체 밖)으로 향하는 법선.
  // 발산 정리로 = A_full (f_bottom − f_top). 강체가 속도 W 로 z 이동하면 이 셀의 경계 유출 = W · nz_area.
  double nz_area = 0.0;
  // 경계 조각의 대표점 (압력 재구성 위치)
  double r = 0.0, z = 0.0;
};

struct CutGeometry {
  flow::Array2D open_u;              // r 면 (nr+1, nz): 열린 비율 0..1
  flow::Array2D open_w;              // z 면 (nr, nz+1): 열린 비율 0..1 (r 가중)
  std::vector<char> active;          // nr*nz: 유체가 있는 셀 (면 중 하나라도 열림)
  std::vector<BoundaryPiece> pieces;  // nz_area ≠ 0 인 유체 셀
  bool grid_fitted = true;           // 모든 열린 비율이 정확히 0 또는 1 (물체가 격자면에 맞음)
};

// 물체가 없을 때: 모든 면 열림, 모든 셀 유체
CutGeometry open_geometry(const flow::Grid& g);

// 중심 높이 z_center 에 놓인 형상의 컷셀 기하. 축대칭 격자는 r0 = 0 이어야 한다.
CutGeometry cut_geometry(const flow::Grid& g, const Shape& shape, double z_center);

}  // namespace hcw::body
