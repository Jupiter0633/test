#include "flow/material.h"

#include <algorithm>

namespace hcw::flow {
namespace {

void allocate(const Grid& g, MaterialArrays& m) {
  if (m.rho_u.ni() == g.nr + 1 && m.rho_u.nj() == g.nz) return;
  m.rho_u = Array2D(g.nr + 1, g.nz);
  m.rho_w = Array2D(g.nr, g.nz + 1);
  m.mu_u = Array2D(g.nr + 1, g.nz);
  m.mu_cell = Array2D(g.nr, g.nz);
  m.mu_corner = Array2D(g.nr + 1, g.nz + 1);
}

}  // namespace

void SmoothMaterial::evaluate(const Grid& g, double t, MaterialArrays& m) const {
  allocate(g, m);
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      m.rho_u(i, j) = rho_(g.rf(i), g.zc(j), t);
      m.mu_u(i, j) = mu_(g.rf(i), g.zc(j), t);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) m.rho_w(i, j) = rho_(g.rc(i), g.zf(j), t);
    for (int j = 0; j < g.nz; ++j) m.mu_cell(i, j) = mu_(g.rc(i), g.zc(j), t);
  }
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) m.mu_corner(i, j) = mu_(g.rf(i), g.zf(j), t);
  }
}

double LayeredMaterial::below_fraction(double za, double zb) const {
  if (zb <= za) return zi_ > za ? 1.0 : 0.0;
  return std::clamp((zi_ - za) / (zb - za), 0.0, 1.0);
}

void LayeredMaterial::evaluate(const Grid& g, double /*t*/, MaterialArrays& m) const {
  allocate(g, m);
  auto rho = [&](double za, double zb) {
    const double f = below_fraction(za, zb);
    return f * rho_b_ + (1.0 - f) * rho_a_;
  };
  auto mu = [&](double za, double zb) {
    const double f = below_fraction(za, zb);
    return f * mu_b_ + (1.0 - f) * mu_a_;
  };
  // 경계 밖으로 나가는 제어체적은 영역 안쪽 부분만 평균한다.
  auto clip_lo = [&](double z) { return std::max(z, g.z0); };
  auto clip_hi = [&](double z) { return std::min(z, g.z1); };

  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      m.rho_u(i, j) = rho(g.zf(j), g.zf(j + 1));
      m.mu_u(i, j) = mu(g.zf(j), g.zf(j + 1));
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) {
      m.rho_w(i, j) = rho(clip_lo(g.zc(j - 1)), clip_hi(g.zc(j)));
    }
    for (int j = 0; j < g.nz; ++j) m.mu_cell(i, j) = mu(g.zf(j), g.zf(j + 1));
  }
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) {
      m.mu_corner(i, j) = mu(clip_lo(g.zc(j - 1)), clip_hi(g.zc(j)));
    }
  }
}

}  // namespace hcw::flow
