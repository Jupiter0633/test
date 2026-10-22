#include "interface/vof.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "interface/plic.h"

namespace hcw::interface {
namespace {

// 이 값보다 0 이나 1 에 가까운 셀은 비었거나 찬 것으로 본다 (재구성 생략)
constexpr double kEmpty = 1e-12;

}  // namespace

Vof::Vof(const flow::Grid& g)
    : g_(g), f_(g.nr, g.nz, 0.0), nr_(g.nr, g.nz), nz_(g.nr, g.nz), alpha_(g.nr, g.nz) {}

void Vof::reconstruct() {
  for (int i = 0; i < g_.nr; ++i) {
    for (int j = 0; j < g_.nz; ++j) {
      const double fd = f_(i, j);
      if (fd <= kEmpty || fd >= 1.0 - kEmpty) continue;
      const auto n = youngs_normal(g_, f_, i, j);
      const Polygon cell = rectangle(g_.rf(i), g_.rf(i + 1), g_.zf(j), g_.zf(j + 1));
      nr_(i, j) = n[0];
      nz_(i, j) = n[1];
      alpha_(i, j) = solve_alpha(cell, n, fd * g_.cell_volume(i), g_.geometry);
    }
  }
}

double Vof::strip_water(int i, int j, double r0, double r1, double z0, double z1) const {
  const Polygon strip = rectangle(r0, r1, z0, z1);
  const double fd = f_(i, j);
  if (fd <= kEmpty) return 0.0;
  if (fd >= 1.0 - kEmpty) return polygon_volume(strip, g_.geometry);
  const std::array<double, 2> n{nr_(i, j), nz_(i, j)};
  return polygon_volume(clip(strip, n, alpha_(i, j)), g_.geometry);
}

void Vof::sweep_r(const flow::Array2D& qu, double dt, const std::vector<char>& c) {
  reconstruct();
  flow::Array2D fw(g_.nr + 1, g_.nz, 0.0);  // 물 플럭스 (Δt 동안의 체적, +r 방향)
  for (int i = 0; i <= g_.nr; ++i) {
    const double area = g_.metric(g_.rf(i)) * g_.dz;
    for (int j = 0; j < g_.nz; ++j) {
      const double q = qu(i, j);
      if (q == 0.0) continue;
      const int donor = q > 0.0 ? i - 1 : i;
      if (donor < 0 || donor >= g_.nr || !(area > 0.0)) {
        throw std::invalid_argument("VOF: nonzero flux through a domain boundary face");
      }
      if (std::abs(q) / area * dt > 0.5 * g_.dr * (1.0 + 1e-12)) {
        throw std::invalid_argument("VOF: CFL > 0.5 in r");
      }
      // 띠 폭 d 를 띠의 (r 가중) 체적이 정확히 |q|Δt 가 되도록 고른다. 축대칭에서는
      // r_f d ∓ d²/2 = |q|Δt/Δz (상류가 안쪽이면 −, 바깥쪽이면 +). 그러면 물 플럭스 = 띠 안 물 체적으로,
      // 비율로 다시 맞출 필요가 없어 0 ≤ f ≤ 1 이 더 잘 지켜진다.
      const double vs = std::abs(q) * dt / g_.dz;
      const double rf = g_.rf(i);
      double dist = vs;
      if (g_.geometry == flow::Geometry::Axisymmetric) {
        dist = q > 0.0 ? vs / (0.5 * (rf + std::sqrt(std::max(0.0, rf * rf - 2.0 * vs))))
                       : vs / (0.5 * (rf + std::sqrt(rf * rf + 2.0 * vs)));
      }
      const double r0 = q > 0.0 ? rf - dist : rf;
      const double r1 = q > 0.0 ? rf : rf + dist;
      fw(i, j) = (q > 0.0 ? 1.0 : -1.0) * strip_water(donor, j, r0, r1, g_.zf(j), g_.zf(j + 1));
    }
  }
  for (int i = 0; i < g_.nr; ++i) {
    const double vol = g_.cell_volume(i);
    for (int j = 0; j < g_.nz; ++j) {
      const double div = (qu(i + 1, j) - qu(i, j)) * dt;
      const double cc = c[static_cast<std::size_t>(i) * g_.nz + j] ? 1.0 : 0.0;
      f_(i, j) = (f_(i, j) * vol - (fw(i + 1, j) - fw(i, j)) + cc * div) / vol;
    }
  }
}

void Vof::sweep_z(const flow::Array2D& qw, double dt, const std::vector<char>& c) {
  reconstruct();
  flow::Array2D fw(g_.nr, g_.nz + 1, 0.0);
  for (int i = 0; i < g_.nr; ++i) {
    const double area = g_.metric(g_.rc(i)) * g_.dr;
    for (int j = 0; j <= g_.nz; ++j) {
      const double q = qw(i, j);
      if (q == 0.0) continue;
      const int donor = q > 0.0 ? j - 1 : j;
      if (donor < 0 || donor >= g_.nz) {
        throw std::invalid_argument("VOF: nonzero flux through a domain boundary face");
      }
      const double dist = std::abs(q) / area * dt;
      if (dist > 0.5 * g_.dz * (1.0 + 1e-12)) throw std::invalid_argument("VOF: CFL > 0.5 in z");
      const double z0 = q > 0.0 ? g_.zf(j) - dist : g_.zf(j);
      const double z1 = q > 0.0 ? g_.zf(j) : g_.zf(j) + dist;
      // z 방향 띠의 체적은 폭에 정확히 비례 (면적 × 폭 = |q|Δt)
      fw(i, j) = (q > 0.0 ? 1.0 : -1.0) * strip_water(i, donor, g_.rf(i), g_.rf(i + 1), z0, z1);
    }
  }
  for (int i = 0; i < g_.nr; ++i) {
    const double vol = g_.cell_volume(i);
    for (int j = 0; j < g_.nz; ++j) {
      const double div = (qw(i, j + 1) - qw(i, j)) * dt;
      const double cc = c[static_cast<std::size_t>(i) * g_.nz + j] ? 1.0 : 0.0;
      f_(i, j) = (f_(i, j) * vol - (fw(i, j + 1) - fw(i, j)) + cc * div) / vol;
    }
  }
}

void Vof::advect(const flow::Array2D& qu, const flow::Array2D& qw, double dt) {
  std::vector<char> c(static_cast<std::size_t>(g_.nr) * g_.nz);
  for (int i = 0; i < g_.nr; ++i) {
    for (int j = 0; j < g_.nz; ++j) c[static_cast<std::size_t>(i) * g_.nz + j] = f_(i, j) > 0.5;
  }
  if (steps_ % 2 == 0) {
    sweep_r(qu, dt, c);
    sweep_z(qw, dt, c);
  } else {
    sweep_z(qw, dt, c);
    sweep_r(qu, dt, c);
  }
  ++steps_;
}

double Vof::water_volume() const {
  double v = 0.0;
  for (int i = 0; i < g_.nr; ++i) {
    for (int j = 0; j < g_.nz; ++j) v += f_(i, j) * g_.cell_volume(i);
  }
  return v;
}

double Vof::bound_violation() const {
  double b = 0.0;
  for (double v : f_.data()) b = std::max(b, std::max(-v, v - 1.0));
  return b;
}

}  // namespace hcw::interface
