#include "interface/plic.h"

#include <algorithm>
#include <cmath>

namespace hcw::interface {

Polygon rectangle(double r0, double r1, double z0, double z1) {
  return {{r0, z0}, {r1, z0}, {r1, z1}, {r0, z1}};
}

Polygon clip(const Polygon& poly, std::array<double, 2> n, double alpha) {
  Polygon out;
  const std::size_t m = poly.size();
  if (m == 0) return out;
  auto side = [&](const Point& p) { return n[0] * p.r + n[1] * p.z - alpha; };
  for (std::size_t k = 0; k < m; ++k) {
    const Point& a = poly[k];
    const Point& b = poly[(k + 1) % m];
    const double sa = side(a), sb = side(b);
    if (sa <= 0.0) out.push_back(a);
    if ((sa < 0.0 && sb > 0.0) || (sa > 0.0 && sb < 0.0)) {
      const double t = sa / (sa - sb);
      out.push_back({a.r + t * (b.r - a.r), a.z + t * (b.z - a.z)});
    }
  }
  return out;
}

double polygon_volume(const Polygon& poly, flow::Geometry geometry) {
  const std::size_t m = poly.size();
  if (m < 3) return 0.0;
  double s = 0.0;
  for (std::size_t k = 0; k < m; ++k) {
    const Point& a = poly[k];
    const Point& b = poly[(k + 1) % m];
    const double cross = a.r * b.z - b.r * a.z;
    // ∫∫ r dA = (1/6) Σ cross (r_a + r_b),  ∫∫ dA = (1/2) Σ cross
    s += geometry == flow::Geometry::Axisymmetric ? cross * (a.r + b.r) / 6.0 : 0.5 * cross;
  }
  return std::abs(s);
}

double solve_alpha(const Polygon& rect, std::array<double, 2> n, double target,
                   flow::Geometry geometry) {
  double lo = 1e300, hi = -1e300;
  for (const auto& p : rect) {
    const double v = n[0] * p.r + n[1] * p.z;
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  // 물 체적은 α 에 대해 단조 증가
  for (int it = 0; it < 200 && hi - lo > 1e-15 * std::max(1.0, std::abs(hi)); ++it) {
    const double mid = 0.5 * (lo + hi);
    if (polygon_volume(clip(rect, n, mid), geometry) < target) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}

std::array<double, 2> youngs_normal(const flow::Grid& g, const flow::Array2D& f, int i, int j) {
  auto at = [&](int ii, int jj) {
    if (ii < 0) ii = (g.geometry == flow::Geometry::Axisymmetric && g.r0 == 0.0) ? -ii - 1 : 0;
    if (ii >= g.nr) ii = g.nr - 1;
    jj = std::max(0, std::min(g.nz - 1, jj));
    return f(ii, jj);
  };
  const double gr = ((at(i + 1, j + 1) + 2.0 * at(i + 1, j) + at(i + 1, j - 1)) -
                     (at(i - 1, j + 1) + 2.0 * at(i - 1, j) + at(i - 1, j - 1))) /
                    (4.0 * (g.drn(i) + g.drn(i + 1)));  // 균일 격자: 8Δr
  const double gz = ((at(i + 1, j + 1) + 2.0 * at(i, j + 1) + at(i - 1, j + 1)) -
                     (at(i + 1, j - 1) + 2.0 * at(i, j - 1) + at(i - 1, j - 1))) /
                    (4.0 * (g.dzn(j) + g.dzn(j + 1)));
  const double len = std::hypot(gr, gz);
  if (!(len > 0.0)) return {0.0, 1.0};
  return {-gr / len, -gz / len};
}

}  // namespace hcw::interface
