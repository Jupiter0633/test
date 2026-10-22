#pragma once

// P2: 축대칭 물체 형상 (r–z 평면). 물체는 축(r = 0)을 따라 z 방향으로만 움직인다.
// 형상 좌표는 물체 중심 기준이며, 실제 위치는 중심 높이 z_c 만큼 옮겨서 쓴다.
// 볼록하고 축대칭인 형상만 다룬다: 수직선과 만나는 고체 부분은 z 구간 하나,
// 수평선과 만나는 고체 부분은 r 구간 [0, R] 하나이다. 이 성질로 면의 열린 비율을 정확히 구한다.
// 근거: docs/design/02_design_document.md §3.4, 03_verification_plan.md §1.2

#include <array>

namespace hcw::body {

class Shape {
 public:
  virtual ~Shape() = default;

  // 고체 내부(경계 포함)이면 true
  virtual bool inside(double r, double z) const = 0;

  // 수직선 r = const (r ≥ 0) 위의 고체 구간 [lo, hi]. 만나지 않으면 false.
  virtual bool z_extent(double r, double& lo, double& hi) const = 0;

  // 수평선 z = const 위의 고체 반경 R (고체 구간은 [0, R]). 만나지 않으면 false.
  virtual bool r_extent(double z, double& radius) const = 0;

  // 축 방향 크기 (중심에서 위·아래 끝까지 거리)
  virtual double half_height() const = 0;

  // 표면 위 점 (r, z) 에서 물체 바깥쪽 단위 법선 (n_r, n_z)
  virtual std::array<double, 2> outward_normal(double r, double z) const = 0;

  // 표면까지의 부호 거리 (바깥 양수). 구·원기둥은 정확, 회전타원체는 표면 근처에서 2차 근사.
  virtual double distance(double r, double z) const = 0;
};

// 회전타원체 r²/a² + z²/c² ≤ 1 (a: 반경 방향 반축, c: 축 방향 반축). a = c 이면 구.
class Spheroid final : public Shape {
 public:
  Spheroid(double a, double c);
  bool inside(double r, double z) const override;
  bool z_extent(double r, double& lo, double& hi) const override;
  bool r_extent(double z, double& radius) const override;
  double half_height() const override { return c_; }
  std::array<double, 2> outward_normal(double r, double z) const override;
  double distance(double r, double z) const override;
  double a() const { return a_; }
  double c() const { return c_; }

 private:
  double a_, c_;
};

// 평판 원판(원기둥) r ≤ a, |z| ≤ h/2. 모서리는 직각.
class Cylinder final : public Shape {
 public:
  Cylinder(double a, double h);
  bool inside(double r, double z) const override;
  bool z_extent(double r, double& lo, double& hi) const override;
  bool r_extent(double z, double& radius) const override;
  double half_height() const override { return 0.5 * h_; }
  std::array<double, 2> outward_normal(double r, double z) const override;
  double distance(double r, double z) const override;
  double a() const { return a_; }
  double h() const { return h_; }

 private:
  double a_, h_;
};

}  // namespace hcw::body
