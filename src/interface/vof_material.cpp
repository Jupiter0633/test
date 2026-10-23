#include "interface/vof_material.h"

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "common/parallel.h"
#include "interface/plic.h"

namespace hcw::interface {
namespace {

constexpr double kEmpty = 1e-12;

// 셀마다 PLIC 를 한 번 재구성해 두고, 반쪽 직사각형 안 물 체적을 구한다
struct HalfWater {
  const flow::Grid& g;
  const flow::Array2D& f;
  flow::Array2D nr, nz, alpha;
  HalfWater(const flow::Grid& grid, const flow::Array2D& frac)
      : g(grid), f(frac), nr(grid.nr, grid.nz), nz(grid.nr, grid.nz), alpha(grid.nr, grid.nz) {
    parallel::parallel_for(0, g.nr, parallel::kRowGrain, [&](int i_lo, int i_hi) {
      for (int i = i_lo; i < i_hi; ++i) {
        for (int j = 0; j < g.nz; ++j) {
          const double fc = f(i, j);
          if (fc <= kEmpty || fc >= 1.0 - kEmpty) continue;
          const auto n = youngs_normal(g, f, i, j);
          nr(i, j) = n[0];
          nz(i, j) = n[1];
          alpha(i, j) = solve_alpha(rectangle(g.rf(i), g.rf(i + 1), g.zf(j), g.zf(j + 1)), n,
                                    fc * g.cell_volume(i), g.geometry);
        }
      }
    });
  }
  double operator()(int i, int j, double r0, double r1, double z0, double z1) const {
    const Polygon half = rectangle(r0, r1, z0, z1);
    const double fc = f(i, j);
    if (fc <= kEmpty) return 0.0;
    if (fc >= 1.0 - kEmpty) return polygon_volume(half, g.geometry);
    return polygon_volume(clip(half, {nr(i, j), nz(i, j)}, alpha(i, j)), g.geometry);
  }
};

}  // namespace

Fractions fractions_from_vof(const flow::Grid& g, const flow::Array2D& f) {
  Fractions fr;
  fr.cell = f;
  fr.face_u = flow::Array2D(g.nr + 1, g.nz);
  fr.face_w = flow::Array2D(g.nr, g.nz + 1);
  fr.corner = flow::Array2D(g.nr + 1, g.nz + 1);
  const HalfWater water(g, f);
  auto vol = [&](double r0, double r1, double z0, double z1) {
    return polygon_volume(rectangle(r0, r1, z0, z1), g.geometry);
  };
  // r 면 i: 셀 i−1 의 오른쪽 반 + 셀 i 의 왼쪽 반
  parallel::parallel_for(0, g.nr + 1, parallel::kRowGrain, [&](int i_lo, int i_hi) {
    for (int i = i_lo; i < i_hi; ++i) {
      for (int j = 0; j < g.nz; ++j) {
        const double z0 = g.zf(j), z1 = g.zf(j + 1);
        double w = 0.0, v = 0.0;
        if (i > 0) {
          const double rm = g.rc(i - 1);
          w += water(i - 1, j, rm, g.rf(i), z0, z1);
          v += vol(rm, g.rf(i), z0, z1);
        }
        if (i < g.nr) {
          const double rm = g.rc(i);
          w += water(i, j, g.rf(i), rm, z0, z1);
          v += vol(g.rf(i), rm, z0, z1);
        }
        fr.face_u(i, j) = v > 0.0 ? std::clamp(w / v, 0.0, 1.0) : f(std::min(i, g.nr - 1), j);
      }
    }
  });
  // z 면 j: 셀 j−1 의 위쪽 반 + 셀 j 의 아래쪽 반
  parallel::parallel_for(0, g.nr, parallel::kRowGrain, [&](int i_lo, int i_hi) {
    for (int i = i_lo; i < i_hi; ++i) {
      const double r0 = g.rf(i), r1 = g.rf(i + 1);
      for (int j = 0; j <= g.nz; ++j) {
        double w = 0.0, v = 0.0;
        if (j > 0) {
          w += water(i, j - 1, r0, r1, g.zc(j - 1), g.zf(j));
          v += vol(r0, r1, g.zc(j - 1), g.zf(j));
        }
        if (j < g.nz) {
          w += water(i, j, r0, r1, g.zf(j), g.zc(j));
          v += vol(r0, r1, g.zf(j), g.zc(j));
        }
        fr.face_w(i, j) = std::clamp(w / v, 0.0, 1.0);
      }
    }
  });
  auto cf = [&](int i, int j) {
    return f(std::clamp(i, 0, g.nr - 1), std::clamp(j, 0, g.nz - 1));
  };
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) {
      fr.corner(i, j) = 0.25 * (cf(i - 1, j - 1) + cf(i, j - 1) + cf(i - 1, j) + cf(i, j));
    }
  }
  return fr;
}

Fractions fractions_with_solid(const flow::Grid& g, const flow::Array2D& f, const flow::Array2D& phi) {
  flow::Array2D fe(g.nr, g.nz);
  flow::Array2D air(g.nr, g.nz);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      const double ph = phi(i, j);
      fe(i, j) = ph > kEmpty ? std::clamp(f(i, j) / ph, 0.0, 1.0) : 0.0;
      air(i, j) = std::max(0.0, ph - f(i, j));
    }
  }
  // 고체로 가득 찬 셀은 이웃 유체 셀의 값으로 채운다 (두 겹). 물체 옆 면의 반쪽 제어체적이 고체 셀에
  // 걸칠 때 고체를 공기로 보아 물 쪽 면 밀도가 반으로 줄어드는 것을 막는다.
  std::vector<char> known(static_cast<std::size_t>(g.nr) * g.nz);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) known[static_cast<std::size_t>(i) * g.nz + j] = phi(i, j) > kEmpty;
  }
  for (int pass = 0; pass < 2; ++pass) {
    std::vector<char> next = known;
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j < g.nz; ++j) {
        if (known[static_cast<std::size_t>(i) * g.nz + j]) continue;
        double s = 0.0;
        int c = 0;
        for (int di = -1; di <= 1; ++di) {
          for (int dj = -1; dj <= 1; ++dj) {
            const int ii = i + di, jj = j + dj;
            if (ii < 0 || jj < 0 || ii >= g.nr || jj >= g.nz || !known[static_cast<std::size_t>(ii) * g.nz + jj]) continue;
            s += fe(ii, jj);
            ++c;
          }
        }
        if (c > 0) {
          fe(i, j) = s / c;
          next[static_cast<std::size_t>(i) * g.nz + j] = 1;
        }
      }
    }
    known.swap(next);
  }
  Fractions fr = fractions_from_vof(g, fe);
  fr.air = air;
  return fr;
}

void VofMaterial::set(const Fractions& a, double t_a, const Fractions& b, double t_b) {
  a_ = a;
  b_ = b;
  ta_ = t_a;
  tb_ = t_b;
}

void VofMaterial::evaluate(const flow::Grid& g, double t, flow::MaterialArrays& m) const {
  if (a_.cell.ni() != g.nr || a_.cell.nj() != g.nz) {
    throw std::logic_error("VofMaterial: volume fractions not set");
  }
  const double s = tb_ > ta_ ? std::clamp((t - ta_) / (tb_ - ta_), 0.0, 1.0) : 0.0;
  auto lerp = [&](const flow::Array2D& x, const flow::Array2D& y, int i, int j) {
    return (1.0 - s) * x(i, j) + s * y(i, j);
  };
  auto rho = [&](double ff) { return ff * rho_w_ + (1.0 - ff) * rho_a_; };
  auto mu = [&](double ff) { return ff * mu_w_ + (1.0 - ff) * mu_a_; };

  m.rho_u = flow::Array2D(g.nr + 1, g.nz);
  m.mu_u = flow::Array2D(g.nr + 1, g.nz);
  m.rho_w = flow::Array2D(g.nr, g.nz + 1);
  m.mu_cell = flow::Array2D(g.nr, g.nz);
  m.mu_corner = flow::Array2D(g.nr + 1, g.nz + 1);
  m.air_fraction = flow::Array2D(g.nr, g.nz);
  const bool has_air = a_.air.ni() == g.nr && b_.air.ni() == g.nr;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      m.air_fraction(i, j) = has_air ? lerp(a_.air, b_.air, i, j) : 1.0 - lerp(a_.cell, b_.cell, i, j);
    }
  }
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      const double ff = lerp(a_.face_u, b_.face_u, i, j);
      m.rho_u(i, j) = rho(ff);
      m.mu_u(i, j) = mu(ff);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) m.rho_w(i, j) = rho(lerp(a_.face_w, b_.face_w, i, j));
    for (int j = 0; j < g.nz; ++j) m.mu_cell(i, j) = mu(lerp(a_.cell, b_.cell, i, j));
  }
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) m.mu_corner(i, j) = mu(lerp(a_.corner, b_.corner, i, j));
  }
}

}  // namespace hcw::interface
