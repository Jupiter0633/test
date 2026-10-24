#pragma once

// 2변수 (r, z) 2차 절단 Taylor 전개. 제조해의 정확한 1·2차 도함수를
// 손 계산 없이 연쇄법칙으로 얻는다 (기호 미분 오류 방지).

#include <cmath>

namespace hcw::verify {

struct Jet {
  double v = 0, r = 0, z = 0, rr = 0, zz = 0, rz = 0;

  static Jet constant(double c) { return {c, 0, 0, 0, 0, 0}; }
  static Jet var_r(double r) { return {r, 1, 0, 0, 0, 0}; }
  static Jet var_z(double z) { return {z, 0, 1, 0, 0, 0}; }
};

inline Jet operator+(const Jet& a, const Jet& b) {
  return {a.v + b.v, a.r + b.r, a.z + b.z, a.rr + b.rr, a.zz + b.zz, a.rz + b.rz};
}
inline Jet operator-(const Jet& a, const Jet& b) {
  return {a.v - b.v, a.r - b.r, a.z - b.z, a.rr - b.rr, a.zz - b.zz, a.rz - b.rz};
}
inline Jet operator*(double s, const Jet& a) {
  return {s * a.v, s * a.r, s * a.z, s * a.rr, s * a.zz, s * a.rz};
}
inline Jet operator+(double s, const Jet& a) { return Jet::constant(s) + a; }
inline Jet operator-(double s, const Jet& a) { return Jet::constant(s) - a; }
inline Jet operator*(const Jet& a, const Jet& b) {
  return {a.v * b.v,
          a.r * b.v + a.v * b.r,
          a.z * b.v + a.v * b.z,
          a.rr * b.v + 2 * a.r * b.r + a.v * b.rr,
          a.zz * b.v + 2 * a.z * b.z + a.v * b.zz,
          a.rz * b.v + a.r * b.z + a.z * b.r + a.v * b.rz};
}

// f(g) 의 연쇄법칙: f0 = f(g), f1 = f'(g), f2 = f''(g)
inline Jet chain(const Jet& g, double f0, double f1, double f2) {
  return {f0,
          f1 * g.r,
          f1 * g.z,
          f2 * g.r * g.r + f1 * g.rr,
          f2 * g.z * g.z + f1 * g.zz,
          f2 * g.r * g.z + f1 * g.rz};
}
inline Jet sin(const Jet& g) { return chain(g, std::sin(g.v), std::cos(g.v), -std::sin(g.v)); }
inline Jet cos(const Jet& g) { return chain(g, std::cos(g.v), -std::sin(g.v), -std::cos(g.v)); }

}  // namespace hcw::verify
