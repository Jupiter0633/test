#include "body/force.h"

#include <algorithm>
#include <cmath>

namespace hcw::body {
namespace {

constexpr double kTwoPi = 6.28318530717958647692;

// 셀 (i, j) 중심에서 점 (r, z) 까지 압력을 1차 외삽. 이웃 유체 셀이 부족하면 셀 값 그대로.
double extrapolate(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p, int i,
                   int j, double r, double z) {
  const double rc = g.rc(i), zc = g.zc(j), pc = p(i, j);
  double arr = 0, arz = 0, azz = 0, br = 0, bz = 0;
  int count = 0;
  for (int di = -1; di <= 1; ++di) {
    for (int dj = -1; dj <= 1; ++dj) {
      const int ii = i + di, jj = j + dj;
      if ((di == 0 && dj == 0) || ii < 0 || jj < 0 || ii >= g.nr || jj >= g.nz) continue;
      if (!geom.active[static_cast<std::size_t>(ii) * g.nz + jj]) continue;
      const double dr = g.rc(ii) - rc, dz = g.zc(jj) - zc, dp = p(ii, jj) - pc;
      arr += dr * dr;
      arz += dr * dz;
      azz += dz * dz;
      br += dr * dp;
      bz += dz * dp;
      ++count;
    }
  }
  const double det = arr * azz - arz * arz;
  if (count < 3 || !(det > 1e-12 * (arr * azz))) return pc;
  const double gr = (azz * br - arz * bz) / det;
  const double gz = (arr * bz - arz * br) / det;
  return pc + gr * (r - rc) + gz * (z - zc);
}

}  // namespace

double axial_viscous_force(const flow::Grid& g, const CutGeometry& geom, const Shape& shape,
                           double zc, const MomentumClasses& mc, const flow::Array2D& u,
                           const flow::Array2D& w, double mu) {
  const int nz = g.nz;
  // 꼭짓점 (r_f(i), z_f(j)) 의 와도. 쓸 수 있으면 true.
  auto omega = [&](int i, int j, double& value) {
    if (i < 1 || j < 1 || i > g.nr - 1 || j > nz - 1) return false;
    if (shape.inside(g.rf(i), g.zf(j) - zc)) return false;
    if (mc.solid_u[static_cast<std::size_t>(i) * nz + j] || mc.solid_u[static_cast<std::size_t>(i) * nz + j - 1] ||
        mc.solid_w[static_cast<std::size_t>(i) * (nz + 1) + j] ||
        mc.solid_w[static_cast<std::size_t>(i - 1) * (nz + 1) + j]) {
      return false;
    }
    value = (u(i, j) - u(i, j - 1)) / g.dz - (w(i, j) - w(i - 1, j)) / g.dr;
    return true;
  };
  const double reach = 2.5 * std::max(g.dr, g.dz);
  double f = 0.0;
  for (const auto& piece : geom.pieces) {
    if (piece.length <= 0.0) continue;
    // 대표점 근처 꼭짓점들로 ω ≈ c0 + c1 (r − r_b) + c2 (z − z_b) 최소제곱
    double s[3][3] = {{0}}, t[3] = {0};
    int count = 0;
    const int i0 = static_cast<int>(std::floor((piece.r - g.r0) / g.dr));
    const int j0 = static_cast<int>(std::floor((piece.z - g.z0) / g.dz));
    for (int i = i0 - 2; i <= i0 + 3; ++i) {
      for (int j = j0 - 2; j <= j0 + 3; ++j) {
        double om;
        if (!omega(i, j, om)) continue;
        const double dr = g.rf(i) - piece.r, dz = g.zf(j) - piece.z;
        if (std::hypot(dr, dz) > reach) continue;
        const double v[3] = {1.0, dr, dz};
        for (int a = 0; a < 3; ++a) {
          for (int b = 0; b < 3; ++b) s[a][b] += v[a] * v[b];
          t[a] += v[a] * om;
        }
        ++count;
      }
    }
    if (count < 3) continue;
    // 3×3 대칭 연립방정식 (크래머)
    auto det3 = [](const double m[3][3]) {
      return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
             m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
             m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    };
    const double d = det3(s);
    if (!(std::abs(d) > 0.0)) continue;
    double m0[3][3];
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) m0[a][b] = b == 0 ? t[a] : s[a][b];
    }
    const double omega_b = det3(m0) / d;
    const auto n = shape.outward_normal(piece.r, piece.z - zc);
    f += -mu * omega_b * n[0] * g.metric(piece.r) * piece.length;
  }
  return g.geometry == flow::Geometry::Axisymmetric ? kTwoPi * f : f;
}

double axial_pressure_force(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p) {
  double f = 0.0;
  for (const auto& piece : geom.pieces) {
    // nz_area 의 법선은 유체에서 물체 쪽 → 물체가 받는 힘은 + p · nz_area
    f += extrapolate(g, geom, p, piece.i, piece.j, piece.r, piece.z) * piece.nz_area;
  }
  return g.geometry == flow::Geometry::Axisymmetric ? kTwoPi * f : f;
}

}  // namespace hcw::body
