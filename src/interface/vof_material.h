#pragma once

// P3: VOF 체적분율에서 밀도·점성을 만드는 물성.
// 셀 값: ρ = f ρ_w + (1 − f) ρ_a.
// 면 값: 엇갈림 격자의 면 제어체적(이웃 두 셀의 마주 보는 반쪽)에서 PLIC 로 재구성한 물의 비율 f_face 로
//   ρ_face = f_face ρ_w + (1 − f_face) ρ_a. 셀 f 의 단순 평균을 쓰면 수면이 셀 안 낮은 곳에 있을 때
//   얇은 물막을 면 제어체적 절반에 퍼뜨려 공기 쪽에 가짜 수평 압력차가 생기고, 그 크기가 격자를
//   곱게 해도 줄지 않아 공기가 계속 가속된다 (2026-10-22 확인). 반쪽 셀의 실제 물 체적을 쓰면
//   수면이 셀 안에 있어도 정수압 단차가 맞는다.
// 꼭짓점·셀 점성은 이웃 셀 f 의 산술 평균.
// 두 시각 (t_a, 면 비율 a), (t_b, 면 비율 b) 사이를 시간에 대해 직선 보간한다.

#include "flow/material.h"

namespace hcw::interface {

// 위치별 물 비율 (셀, r 면, z 면, 꼭짓점)
struct Fractions {
  flow::Array2D cell;    // (nr, nz)
  flow::Array2D face_u;  // (nr+1, nz)
  flow::Array2D face_w;  // (nr, nz+1)
  flow::Array2D corner;  // (nr+1, nz+1)
  // 셀 전체 체적 대비 공기 비율 (P4, 물체가 있을 때만). 비어 있으면 1 − cell.
  flow::Array2D air;
};

// 셀 체적분율 f 에서 PLIC 재구성으로 면 제어체적의 물 비율을 구한다
Fractions fractions_from_vof(const flow::Grid& g, const flow::Array2D& f);

// 물체가 있을 때 (P4). f: 셀 전체 체적 대비 물 비율, phi: 유체 체적 비율 (body::fluid_fraction).
// 밀도·점성은 유체 부분 안의 물 비율 f/φ 로 만든다 (물에 잠긴 컷셀이 물 밀도를 갖도록; 고체 셀은 공기 값).
// 공기 비율(압축성)은 고체를 뺀 φ − f.
Fractions fractions_with_solid(const flow::Grid& g, const flow::Array2D& f, const flow::Array2D& phi);

class VofMaterial : public flow::Material {
 public:
  VofMaterial(double rho_water, double mu_water, double rho_air, double mu_air)
      : rho_w_(rho_water), mu_w_(mu_water), rho_a_(rho_air), mu_a_(mu_air) {}

  bool time_dependent() const override { return true; }
  void evaluate(const flow::Grid& g, double t, flow::MaterialArrays& out) const override;

  // 시각 t_a 에서 a, t_b 에서 b. t_a == t_b 이면 a 만 쓴다.
  void set(const Fractions& a, double t_a, const Fractions& b, double t_b);

  double rho_water() const { return rho_w_; }
  double rho_air() const { return rho_a_; }

 private:
  double rho_w_, mu_w_, rho_a_, mu_a_;
  Fractions a_, b_;
  double ta_ = 0.0, tb_ = 0.0;
};

}  // namespace hcw::interface
