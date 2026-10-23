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
  // 경계 조각의 길이 (r–z 평면, 현의 길이). 점성 응력 적분에 쓴다.
  double length = 0.0;
};

struct CutGeometry {
  flow::Array2D open_u;              // r 면 (nr+1, nz): 열린 비율 0..1
  flow::Array2D open_w;              // z 면 (nr, nz+1): 열린 비율 0..1 (r 가중)
  std::vector<char> active;          // nr*nz: 유체가 있는 셀 (면 중 하나라도 열림)
  std::vector<BoundaryPiece> pieces;  // 물체 경계가 지나가는 유체 셀 (nz_area = 0 일 수도 있음)
  bool grid_fitted = true;           // 모든 열린 비율이 정확히 0 또는 1 (물체가 격자면에 맞음)
};

// 물체가 없을 때: 모든 면 열림, 모든 셀 유체
CutGeometry open_geometry(const flow::Grid& g);

// 중심 높이 z_center 에 놓인 형상의 컷셀 기하. 축대칭 격자는 r0 = 0 이어야 한다.
CutGeometry cut_geometry(const flow::Grid& g, const Shape& shape, double z_center);

// 셀마다 유체가 차지하는 체적 비율 (0..1, 축대칭은 r 가중). P4: VOF 와 물체를 함께 쓸 때
// 물 체적분율의 상한과 공기 체적(압축성) 계산에 쓴다. 수평선마다 고체 반경 R(z) 를 써서
// ∫∫ r dr dz 를 z 방향 합성 가우스 적분 (형상의 위·아래 끝에서 구간을 나눔) 으로 구한다.
// 원기둥·쐐기는 반올림 수준, 회전타원체는 끝의 제곱근 특이성 때문에 약 1e−6 정도.
flow::Array2D fluid_fraction(const flow::Grid& g, const Shape& shape, double z_center);

// 운동량(속도 면) 분류 — 점성 흐름의 no-slip 처리 (Fadlun 형 직접 강제, 2차 보간).
// 고체 면: 면 중심이 물체 안 → 물체 속도.
// 경계 보간면: 면 중심은 유체이고 같은 성분의 이웃 면(±r 또는 ±z)이 고체 →
//   그 방향의 물체 경계점 b 와 반대쪽 바깥 면 o1, o2 를 지나는 2차 라그랑주 보간:
//   u_f = w_b u_b + w1 u_o1 + w2 u_o2  (d: 면에서 경계까지 거리, h: 격자 간격)
//   w_b = 2h²/((d+h)(d+2h)), w1 = 2d/(d+h), w2 = −d/(d+2h).
//   o2 가 없거나 고체면 1차 (직선) 보간 w1 = d/(d+h), w2 = 0.
// 2차 보간은 벽 근처 속도 3차, 벽 기울기(점성 견인력) 2차 정확도를 준다 (1-G 정상 Stokes 유동).
// 다만 음의 가중치(w2) 때문에, 물체가 움직이고 점성이 작은 비정상 계산(1-E(b), μ = 0.01)에서
// 불안정해짐을 확인했다 (2026-10-22, 속도가 지수적으로 커짐; 1차 보간은 안정). 그래서 차수를 고를 수 있게 한다.
struct FaceForcing {
  int comp = 0;    // 0: r 면(u), 1: z 면(w)
  int i = 0, j = 0;
  int oi = 0, oj = 0;    // 바깥 면 o1 (같은 성분)
  int oi2 = 0, oj2 = 0;  // 바깥 면 o2
  double weight = 0.0;   // w1 (0 이면 o1 을 쓰지 않음 → 물체 속도)
  double weight2 = 0.0;  // w2 (0 이면 o2 를 쓰지 않음)
};

struct MomentumClasses {
  std::vector<char> solid_u;  // (nr+1)*nz, 행 우선 (i * nz + j)
  std::vector<char> solid_w;  // nr*(nz+1), 행 우선 (i * (nz+1) + j)
  std::vector<FaceForcing> forcing;
};

// order: 경계 보간 차수 (1: 직선, 2: 2차)
MomentumClasses momentum_classes(const flow::Grid& g, const Shape& shape, double z_center, int order = 2);

}  // namespace hcw::body
