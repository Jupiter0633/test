#pragma once

// P5 검증 시험의 정확해 (03 §1.5, 기준 고정 커밋 9bf931b).
// 아래 식은 sympy 로 만든 것입니다 (tools/p5_reference_values.py 와 같은 유도, E = 1, ν = 0.3).
// - mms_*: 4-A(b)·4-G 의 제조해 U = (ε r cos(πr²/2) sin(πz), ε cos(πr²) cos(πz/2)), 그 응력과 응력 발산 div σ.
// - plate_*: 4-C 원판 굽힘의 다항식 정확해 (p = 1, a = 1, 두께 h).
// 응력 순서 (σ_rr, σ_θθ, σ_zz, σ_rz). div σ 의 성분: (∂σ_rr/∂r + ∂σ_rz/∂z + (σ_rr − σ_θθ)/r, ∂σ_rz/∂r + ∂σ_zz/∂z + σ_rz/r).

#include <cmath>

namespace hcw::verify::p5exact {

constexpr double kPi = 3.14159265358979323846;
using std::cos;
using std::pow;
using std::sin;

// mms_u
inline void mms_u(double r, double z, double eps, double* out) {
  out[0] = eps*r*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2));
  out[1] = eps*cos(kPi*pow(r, 2))*cos((1.0/2.0)*kPi*z);
}

// mms_stress
inline void mms_stress(double r, double z, double eps, double* out) {
  out[0] = (5.0/52.0)*eps*(-14*kPi*pow(r, 2)*sin((1.0/2.0)*kPi*pow(r, 2))*sin(kPi*z) - 3*kPi*sin((1.0/2.0)*kPi*z)*cos(kPi*pow(r, 2)) + 20*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2)));
  out[1] = (5.0/52.0)*eps*(-6*kPi*pow(r, 2)*sin((1.0/2.0)*kPi*pow(r, 2))*sin(kPi*z) - 3*kPi*sin((1.0/2.0)*kPi*z)*cos(kPi*pow(r, 2)) + 20*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2)));
  out[2] = (5.0/52.0)*eps*(-6*kPi*pow(r, 2)*sin((1.0/2.0)*kPi*pow(r, 2))*sin(kPi*z) - 7*kPi*sin((1.0/2.0)*kPi*z)*cos(kPi*pow(r, 2)) + 12*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2)));
  out[3] = (5.0/13.0)*kPi*eps*r*(-2*sin(kPi*pow(r, 2))*cos((1.0/2.0)*kPi*z) + cos((1.0/2.0)*kPi*pow(r, 2))*cos(kPi*z));
}

// mms_div_stress
inline void mms_div_stress(double r, double z, double eps, double* out) {
  out[0] = (5.0/26.0)*kPi*eps*r*(-7*kPi*pow(r, 2)*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2)) - 28*sin((1.0/2.0)*kPi*pow(r, 2))*sin(kPi*z) + 5*kPi*sin(kPi*pow(r, 2))*sin((1.0/2.0)*kPi*z) - 2*kPi*sin(kPi*z)*cos((1.0/2.0)*kPi*pow(r, 2)));
  out[1] = (5.0/104.0)*kPi*eps*(-20*kPi*pow(r, 2)*sin((1.0/2.0)*kPi*pow(r, 2))*cos(kPi*z) - 32*kPi*pow(r, 2)*cos(kPi*pow(r, 2))*cos((1.0/2.0)*kPi*z) - 32*sin(kPi*pow(r, 2))*cos((1.0/2.0)*kPi*z) + 40*cos((1.0/2.0)*kPi*pow(r, 2))*cos(kPi*z) - 7*kPi*cos(kPi*pow(r, 2))*cos((1.0/2.0)*kPi*z));
}

// plate_u
inline void plate_u(double r, double z, double h, double* out) {
  out[0] = (13.0/400.0)*r*z*(60*pow(h, 2) + 21*pow(r, 2) - 68*pow(z, 2) - 21)/pow(h, 3);
  out[1] = (13.0/11200.0)*(-320*pow(h, 3)*z - 147*pow(r, 4) + 42*pow(r, 2)*(7 - 12*pow(z, 2)) + 728*pow(z, 4) + 12*pow(z, 2)*(21 - 100*pow(h, 2)) - 147)/pow(h, 3);
}

// plate_stress
inline void plate_stress(double r, double z, double h, double* out) {
  out[0] = (1.0/280.0)*(-60*pow(h, 3) + 600*pow(h, 2)*z + 693*pow(r, 2)*z - 644*pow(z, 3) - 273*z)/pow(h, 3);
  out[1] = (1.0/280.0)*(-60*pow(h, 3) + 600*pow(h, 2)*z + 399*pow(r, 2)*z - 644*pow(z, 3) - 273*z)/pow(h, 3);
  out[2] = -1.0/2.0 - 3.0/2.0*z/h + 2*pow(z, 3)/pow(h, 3);
  out[3] = (3.0/4.0)*r*(pow(h, 2) - 4*pow(z, 2))/pow(h, 3);
}

}  // namespace hcw::verify::p5exact
