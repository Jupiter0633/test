#include "verify/potential.h"

#include <cmath>
#include <stdexcept>

namespace hcw::verify {
namespace {

constexpr double kPi = 3.14159265358979323846;

// g(x) = (1 − atan(x)/x) / x² . 작은 x 는 급수 (1/3 − x²/5 + x⁴/7 − …) 로 상쇄를 피한다.
double g_of(double x) {
  if (x < 0.1) {
    double sum = 0.0, term = 1.0;
    for (int k = 0; k < 12; ++k) {
      sum += term / (2 * k + 3);
      term *= -x * x;
    }
    return sum;
  }
  return (1.0 - std::atan(x) / x) / (x * x);
}

}  // namespace

SpheroidPotential::SpheroidPotential(double a, double c) : a_(a), c_(c) {
  if (!(a > 0.0 && c >= 0.0 && a >= c)) {
    throw std::invalid_argument("spheroid potential needs a >= c >= 0 (sphere or oblate)");
  }
  e_ = std::sqrt(a * a - c * c);
  if (c > 0.0) {
    // α₀ = a² c ∫_0^∞ = a² c (2/c³) g(e/c) = 2 (a/c)² g(e/c)
    const double x = e_ / c;
    alpha0_ = 2.0 * (a / c) * (a / c) * g_of(x);
    // 2 − α₀ 를 상쇄 없이: 2 − 2(1+x²)·(1 − atan(x)/x)/x² , y = atan(x)/x 일 때
    // = 2[y − (1 − y)/x²] ... x 가 작으면 직접 빼도 α₀ ≈ 2/3 이라 상쇄가 없다.
    if (x < 0.1) {
      two_minus_alpha0_ = 2.0 - alpha0_;
    } else {
      const double y = std::atan(x) / x;
      two_minus_alpha0_ = 2.0 * (y - (1.0 - y) / (x * x));
    }
    coef_ = a * a * c / two_minus_alpha0_;
  } else {
    // 얇은 원판 극한 (c → 0): α₀ → 2, C → a³/π (Lamb Art. 102 와 같은 해)
    alpha0_ = 2.0;
    two_minus_alpha0_ = 0.0;
    coef_ = a * a * a / kPi;
  }
}

double SpheroidPotential::added_mass_coefficient() const {
  if (c_ == 0.0) throw std::logic_error("added mass coefficient is undefined for a zero-thickness disk");
  return alpha0_ / two_minus_alpha0_;
}

double SpheroidPotential::added_mass(double rho) const {
  // k ρ (4/3)π a² c = α₀ ρ (4/3)π a² c / (2 − α₀) = α₀ (4/3)π ρ C
  return alpha0_ * (4.0 / 3.0) * kPi * rho * coef_;
}

double SpheroidPotential::integral(double lambda) const {
  const double s = std::sqrt(c_ * c_ + lambda);
  // (2/e²)(1/s)(1 − atan(x)/x) = (2/s³) g(x), x = e/s
  return 2.0 / (s * s * s) * g_of(e_ / s);
}

double SpheroidPotential::lambda(double r, double z) const {
  const double a2 = a_ * a_, c2 = c_ * c_;
  const double b = a2 + c2 - r * r - z * z;
  const double c0 = a2 * c2 - a2 * z * z - c2 * r * r;
  if (c0 >= 0.0 && (c2 > 0.0 ? (z * z / c2 + r * r / a2 <= 1.0) : (z == 0.0 && r <= a_))) return 0.0;
  const double disc = std::sqrt(std::max(0.0, b * b - 4.0 * c0));
  return b > 0.0 ? -2.0 * c0 / (b + disc) : 0.5 * (-b + disc);
}

std::array<double, 2> SpheroidPotential::velocity(double r, double z) const {
  const double l = lambda(r, z);
  const double a2l = a_ * a_ + l, c2l = c_ * c_ + l;
  const double ival = integral(l);
  const double di = -1.0 / (std::sqrt(c2l) * c2l * a2l);  // dI/dλ
  const double d = z * z / (c2l * c2l) + r * r / (a2l * a2l);
  double lr = 0.0, lz = 0.0;
  if (l > 0.0 && d > 0.0) {
    lr = 2.0 * r / a2l / d;
    lz = 2.0 * z / c2l / d;
  }
  // φ = C z I(λ),  u = −∇φ
  const double phi_r = coef_ * z * di * lr;
  const double phi_z = coef_ * (ival + z * di * lz);
  return {-phi_r, -phi_z};
}

}  // namespace hcw::verify
