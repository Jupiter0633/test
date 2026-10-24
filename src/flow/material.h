#pragma once

// 밀도·점성 분포. P1 에서는 주어진 함수(시간 의존 가능)이며,
// P3 에서 VOF 체적분율로부터 계산하는 구현이 추가된다.

#include <functional>

#include "flow/grid.h"

namespace hcw::flow {

// 해석기가 쓰는 위치별 물성.
struct MaterialArrays {
  Array2D rho_u;      // r 면 (nr+1, nz): u 운동량 제어체적의 밀도
  Array2D rho_w;      // z 면 (nr, nz+1): w 운동량 제어체적의 밀도
  Array2D mu_u;       // r 면 (nr+1, nz): 원주 응력 항
  Array2D mu_cell;    // 셀 중심 (nr, nz): 수직 응력
  Array2D mu_corner;  // 격자 꼭짓점 (nr+1, nz+1): 전단 응력
  Array2D air_fraction;  // 셀 중심 (nr, nz): 공기 체적분율 (압축성 공기용, 단상 물성은 0)
};

class Material {
 public:
  virtual ~Material() = default;
  virtual bool time_dependent() const = 0;
  virtual void evaluate(const Grid& g, double t, MaterialArrays& out) const = 0;
};

// 매끈한 함수로 주어진 물성. 각 위치의 점값을 쓴다 (2차 정확도).
class SmoothMaterial : public Material {
 public:
  using Fn = std::function<double(double r, double z, double t)>;
  SmoothMaterial(Fn rho, Fn mu, bool time_dependent)
      : rho_(std::move(rho)), mu_(std::move(mu)), time_dependent_(time_dependent) {}
  bool time_dependent() const override { return time_dependent_; }
  void evaluate(const Grid& g, double t, MaterialArrays& out) const override;

 private:
  Fn rho_, mu_;
  bool time_dependent_;
};

// z = z_interface 를 경계로 아래(below)·위(above) 두 층. 밀도·점성은 각
// 제어체적에 대한 정확한 평균을 쓴다 (정수압 평형이 이산적으로 정확히 성립).
class LayeredMaterial : public Material {
 public:
  LayeredMaterial(double z_interface, double rho_below, double mu_below, double rho_above,
                  double mu_above)
      : zi_(z_interface), rho_b_(rho_below), mu_b_(mu_below), rho_a_(rho_above),
        mu_a_(mu_above) {}
  bool time_dependent() const override { return false; }
  void evaluate(const Grid& g, double t, MaterialArrays& out) const override;

  // 구간 [za, zb] 에서 아래층이 차지하는 비율
  double below_fraction(double za, double zb) const;

 private:
  double zi_, rho_b_, mu_b_, rho_a_, mu_a_;
};

}  // namespace hcw::flow
