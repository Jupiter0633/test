#include "body/cut_cell.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hcw::body {
namespace {

// [x0, x1] 구간에서 고체가 차지하는 길이 비율 (고체 구간 [lo, hi])
double solid_fraction_1d(double x0, double x1, double lo, double hi) {
  const double a = std::max(x0, lo), b = std::min(x1, hi);
  return b > a ? (b - a) / (x1 - x0) : 0.0;
}

// z 면 (r ∈ [r0, r1]) 에서 고체 [0, R] 가 차지하는 비율. 축대칭은 r 가중.
double solid_fraction_z_face(const flow::Grid& g, double r0, double r1, double radius) {
  const double b = std::min(r1, radius);
  if (b <= r0) return 0.0;
  if (g.geometry == flow::Geometry::Axisymmetric) return (b * b - r0 * r0) / (r1 * r1 - r0 * r0);
  return (b - r0) / (r1 - r0);
}

}  // namespace

CutGeometry open_geometry(const flow::Grid& g) {
  CutGeometry c;
  c.open_u = flow::Array2D(g.nr + 1, g.nz, 1.0);
  c.open_w = flow::Array2D(g.nr, g.nz + 1, 1.0);
  c.active.assign(static_cast<std::size_t>(g.nr) * g.nz, 1);
  return c;
}

CutGeometry cut_geometry(const flow::Grid& g, const Shape& shape, double zc) {
  if (g.geometry == flow::Geometry::Axisymmetric && g.r0 != 0.0) {
    throw std::invalid_argument("bodies need an axisymmetric grid that starts at the axis");
  }
  CutGeometry c = open_geometry(g);

  // r 면: 수직선 r = r_f 의 고체 z 구간
  for (int i = 0; i <= g.nr; ++i) {
    double lo, hi;
    if (!shape.z_extent(g.rf(i), lo, hi)) continue;
    lo += zc;
    hi += zc;
    for (int j = 0; j < g.nz; ++j) {
      c.open_u(i, j) = 1.0 - solid_fraction_1d(g.zf(j), g.zf(j + 1), lo, hi);
    }
  }
  // z 면: 수평선 z = z_f 의 고체 반경
  for (int j = 0; j <= g.nz; ++j) {
    double radius;
    if (!shape.r_extent(g.zf(j) - zc, radius)) continue;
    for (int i = 0; i < g.nr; ++i) {
      c.open_w(i, j) = 1.0 - solid_fraction_z_face(g, g.rf(i), g.rf(i + 1), radius);
    }
  }
  // 반올림으로 생긴 아주 작은 비율은 0/1 로 맞춘다 (격자면에 놓인 경계)
  for (double& v : c.open_u.data()) {
    if (std::abs(v) < 1e-12) v = 0.0;
    if (std::abs(1.0 - v) < 1e-12) v = 1.0;
    if (v != 0.0 && v != 1.0) c.grid_fitted = false;
  }
  for (double& v : c.open_w.data()) {
    if (std::abs(v) < 1e-12) v = 0.0;
    if (std::abs(1.0 - v) < 1e-12) v = 1.0;
    if (v != 0.0 && v != 1.0) c.grid_fitted = false;
  }

  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      // 볼록 물체: 네 면이 모두 막힌 셀은 고체로 가득 차 있다.
      // 축(r = 0)의 r 면은 넓이가 0 이므로 판단에서 뺀다.
      const bool west_open = (g.geometry == flow::Geometry::Axisymmetric && i == 0)
                                 ? false
                                 : c.open_u(i, j) > 0.0;
      const bool any_open = west_open || c.open_u(i + 1, j) > 0.0 || c.open_w(i, j) > 0.0 ||
                            c.open_w(i, j + 1) > 0.0;
      c.active[static_cast<std::size_t>(i) * g.nz + j] = any_open ? 1 : 0;
      if (!any_open) continue;

      const double full = g.metric(g.rc(i)) * g.dr;
      const double nz_area = full * (c.open_w(i, j) - c.open_w(i, j + 1));
      if (nz_area == 0.0) continue;

      BoundaryPiece p;
      p.i = i;
      p.j = j;
      p.nz_area = nz_area;
      // 대표점: 셀 경계와 물체 경계가 만나는 두 점의 중점.
      const double r0 = g.rf(i), r1 = g.rf(i + 1), z0 = g.zf(j), z1 = g.zf(j + 1);
      std::vector<std::pair<double, double>> pts;
      for (double rr : {r0, r1}) {
        double lo, hi;
        if (!shape.z_extent(rr, lo, hi)) continue;
        for (double zz : {lo + zc, hi + zc}) {
          if (zz > z0 && zz < z1) pts.emplace_back(rr, zz);
        }
      }
      for (double zz : {z0, z1}) {
        double radius;
        if (!shape.r_extent(zz - zc, radius)) continue;
        if (radius > r0 && radius < r1) pts.emplace_back(radius, zz);
      }
      if (pts.size() == 2) {
        p.r = 0.5 * (pts[0].first + pts[1].first);
        p.z = 0.5 * (pts[0].second + pts[1].second);
      } else {
        // 경계가 셀 면과 겹치는 경우 (격자면에 맞는 물체): 고체가 더 많은 z 면 위의 고체 부분 중점
        const bool bottom = (1.0 - c.open_w(i, j)) >= (1.0 - c.open_w(i, j + 1));
        const double zf = bottom ? z0 : z1;
        double radius = r1;
        shape.r_extent(zf - zc, radius);
        p.r = 0.5 * (r0 + std::min(r1, std::max(r0, radius)));
        p.z = zf;
      }
      c.pieces.push_back(p);
    }
  }
  return c;
}

}  // namespace hcw::body
