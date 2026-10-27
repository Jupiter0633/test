#pragma once

// 회전타원체가 축 방향으로 움직일 때의 퍼텐셜 유동 (무한 유체, 비점성).
// 근거: H. Lamb, Hydrodynamics, 6th ed. (1932), Art. 114 식 (6)–(12) [06_references.md R-21].
//   φ = C z ∫_λ^∞ dλ' / ((c²+λ') Δ),  Δ = sqrt((c²+λ')(a²+λ')²),  C = a²c U / (2 − α₀),
//   α₀ = a²c ∫_0^∞ dλ' / ((c²+λ') Δ),  부가질량 = α₀/(2−α₀) × ρ (4/3)π a² c.
// 속도는 Lamb 의 부호 약속대로 u = −∇φ. 구(a = c)는 Art. 92 의 이중 원천 해와 같다.
// 반경 방향 반축 a ≥ 축 방향 반축 c 인 경우(구, 납작한 회전타원체)만 다룬다. c → 0 은 얇은 원판 (Art. 102).
// 적분은 s = sqrt(c²+λ) 로 바꾸면 닫힌 꼴: ∫ = (2/e²)(1/s)(1 − atan(x)/x), e² = a² − c², x = e/s.

#include <array>

namespace hcw::verify {

class SpheroidPotential {
 public:
  SpheroidPotential(double a, double c);

  double alpha0() const { return alpha0_; }
  // 부가질량 계수 k = α₀/(2−α₀) (밀어낸 질량에 대한 비)
  double added_mass_coefficient() const;
  // 부가질량 (3차원). c = 0 (원판) 이면 (8/3)ρa³.
  double added_mass(double rho) const;

  // 물체 중심 기준 위치 (r, z) 의 유속 (u_r, u_z), 물체 속도 U = 1. 물체 바깥 점에서만 쓴다.
  std::array<double, 2> velocity(double r, double z) const;

  // ∫_λ^∞ dλ' / ((c²+λ') Δ)
  double integral(double lambda) const;
  // z²/(c²+λ) + r²/(a²+λ) = 1 의 양의 근 (안쪽이면 0)
  double lambda(double r, double z) const;

 private:
  double a_, c_, e_;
  double alpha0_, two_minus_alpha0_;
  double coef_;  // C / U
};

// 고정된 구(반지름 a)를 지나는 Stokes 유동, 무한 원방에서 균일 흐름 U (z 방향).
// Lamb (1932) Art. 337 식 (11) [R-30] 을 축대칭 좌표로 옮긴 것 (Lamb 의 x 축 = 여기 z 축):
//   u_z = U(1 − a/R) + ¼Ua(R² − a²)(1/R³ − 3z²/R⁵),  u_r = ¼Ua(R² − a²)(−3zr/R⁵),  R² = r² + z².
// 압력 p = −(3/2) μ a U z / R³ (상수 제외)는 Stokes 식 ∇p = μ∇²u 를 만족하도록 유도한 것이며
// 단위 시험에서 그 잔차를 수치로 확인한다. 항력 = 6πμaU (압력 2πμaU + 점성 4πμaU).
std::array<double, 2> stokes_sphere_velocity(double a, double U, double r, double z);
double stokes_sphere_pressure(double a, double U, double mu, double r, double z);

}  // namespace hcw::verify
