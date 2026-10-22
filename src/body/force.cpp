#include "body/force.h"

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

double axial_pressure_force(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p) {
  double f = 0.0;
  for (const auto& piece : geom.pieces) {
    // nz_area 의 법선은 유체에서 물체 쪽 → 물체가 받는 힘은 + p · nz_area
    f += extrapolate(g, geom, p, piece.i, piece.j, piece.r, piece.z) * piece.nz_area;
  }
  return g.geometry == flow::Geometry::Axisymmetric ? kTwoPi * f : f;
}

}  // namespace hcw::body
