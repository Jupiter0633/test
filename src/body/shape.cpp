#include "body/shape.h"

#include <cmath>
#include <stdexcept>

namespace hcw::body {

Spheroid::Spheroid(double a, double c) : a_(a), c_(c) {
  if (!(a > 0.0 && c > 0.0)) throw std::invalid_argument("spheroid semi-axes must be positive");
}

bool Spheroid::inside(double r, double z) const {
  return (r / a_) * (r / a_) + (z / c_) * (z / c_) <= 1.0;
}

bool Spheroid::z_extent(double r, double& lo, double& hi) const {
  const double s = 1.0 - (r / a_) * (r / a_);
  if (s <= 0.0) return false;
  hi = c_ * std::sqrt(s);
  lo = -hi;
  return true;
}

bool Spheroid::r_extent(double z, double& radius) const {
  const double s = 1.0 - (z / c_) * (z / c_);
  if (s <= 0.0) return false;
  radius = a_ * std::sqrt(s);
  return true;
}

std::array<double, 2> Spheroid::outward_normal(double r, double z) const {
  const double nr = r / (a_ * a_), nz = z / (c_ * c_);
  const double len = std::hypot(nr, nz);
  return {nr / len, nz / len};
}

Cylinder::Cylinder(double a, double h) : a_(a), h_(h) {
  if (!(a > 0.0 && h > 0.0)) throw std::invalid_argument("cylinder radius and height must be positive");
}

bool Cylinder::inside(double r, double z) const {
  return r <= a_ && std::abs(z) <= 0.5 * h_;
}

bool Cylinder::z_extent(double r, double& lo, double& hi) const {
  if (r > a_) return false;
  lo = -0.5 * h_;
  hi = 0.5 * h_;
  return true;
}

bool Cylinder::r_extent(double z, double& radius) const {
  if (std::abs(z) > 0.5 * h_) return false;
  radius = a_;
  return true;
}

std::array<double, 2> Cylinder::outward_normal(double r, double z) const {
  // 옆면과 윗면·아랫면 중 더 가까운 면의 법선 (모서리는 정의되지 않으므로 가까운 쪽)
  const double side = std::abs(a_ - r), cap = std::abs(0.5 * h_ - std::abs(z));
  if (side < cap) return {1.0, 0.0};
  return {0.0, z >= 0.0 ? 1.0 : -1.0};
}

}  // namespace hcw::body
