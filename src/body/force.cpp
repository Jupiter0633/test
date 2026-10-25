#include "body/force.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace hcw::body {
namespace {

constexpr double kTwoPi = 6.28318530717958647692;

// 셀 (i, j) 중심에서 점 (r, z) 까지 압력을 1차 외삽. 이웃 유체 셀이 부족하면 셀 값 그대로.
// 자기 셀에 압력이 없을 때: 주변 유효 셀(5×5)로 p ≈ c0 + c1 (r − r_b) + c2 (z − z_b) 최소제곱, c0 반환
double fit_at(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p, int i, int j,
              double r, double z) {
  double s[3][3] = {{0}}, t[3] = {0};
  int count = 0;
  for (int ii = i - 2; ii <= i + 2; ++ii) {
    for (int jj = j - 2; jj <= j + 2; ++jj) {
      if (ii < 0 || jj < 0 || ii >= g.nr || jj >= g.nz) continue;
      if (!geom.active[static_cast<std::size_t>(ii) * g.nz + jj]) continue;
      const double v[3] = {1.0, g.rc(ii) - r, g.zc(jj) - z};
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) s[a][b] += v[a] * v[b];
        t[a] += v[a] * p(ii, jj);
      }
      ++count;
    }
  }
  auto det3 = [](const double m[3][3]) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  };
  const double d = det3(s);
  if (count < 3 || !(std::abs(d) > 0.0)) return 0.0;
  double m0[3][3];
  for (int a = 0; a < 3; ++a) {
    for (int b = 0; b < 3; ++b) m0[a][b] = b == 0 ? t[a] : s[a][b];
  }
  return det3(m0) / d;
}

double extrapolate(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p, int i,
                   int j, double r, double z) {
  if (!geom.active[static_cast<std::size_t>(i) * g.nz + j]) return fit_at(g, geom, p, i, j, r, z);
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

// n×n 연립방정식 (부분 피벗 가우스 소거). 특이하면 false.
bool solve_small(int n, double* m, double* rhs) {
  for (int c = 0; c < n; ++c) {
    int piv = c;
    for (int r = c + 1; r < n; ++r) {
      if (std::abs(m[r * n + c]) > std::abs(m[piv * n + c])) piv = r;
    }
    if (!(std::abs(m[piv * n + c]) > 1e-300)) return false;
    if (piv != c) {
      for (int k = 0; k < n; ++k) std::swap(m[c * n + k], m[piv * n + k]);
      std::swap(rhs[c], rhs[piv]);
    }
    for (int r = c + 1; r < n; ++r) {
      const double f = m[r * n + c] / m[c * n + c];
      for (int k = c; k < n; ++k) m[r * n + k] -= f * m[c * n + k];
      rhs[r] -= f * rhs[c];
    }
  }
  for (int c = n - 1; c >= 0; --c) {
    for (int k = c + 1; k < n; ++k) rhs[c] -= m[c * n + k] * rhs[k];
    rhs[c] /= m[c * n + c];
  }
  return true;
}

// 벽 거리 최소제곱의 조건: 칸 크기 h 로 무차원화한 정규방정식을 부분 피벗 소거해
// 가장 작은 피벗 / 가장 큰 대각 성분을 돌려준다 (0 에 가까우면 이웃 점들이 거의 한 줄에 놓여 계수가 정해지지 않음).
double scaled_pivot_ratio(int n, const double* m, const int* power, double h) {
  double a[25];
  double diag = 0.0;
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < n; ++c) a[r * n + c] = m[r * n + c] / std::pow(h, power[r] + power[c]);
    diag = std::max(diag, std::abs(a[r * n + r]));
  }
  double worst = 1e300;
  for (int c = 0; c < n; ++c) {
    int piv = c;
    for (int r = c + 1; r < n; ++r) {
      if (std::abs(a[r * n + c]) > std::abs(a[piv * n + c])) piv = r;
    }
    if (piv != c) {
      for (int k = 0; k < n; ++k) std::swap(a[c * n + k], a[piv * n + k]);
    }
    worst = std::min(worst, std::abs(a[c * n + c]));
    if (!(std::abs(a[c * n + c]) > 0.0)) return 0.0;
    for (int r = c + 1; r < n; ++r) {
      const double f = a[r * n + c] / a[c * n + c];
      for (int k = c; k < n; ++k) a[r * n + k] -= f * a[c * n + k];
    }
  }
  return diag > 0.0 ? worst / diag : 0.0;
}

double axial_viscous_force(const flow::Grid& g, const CutGeometry& geom, const Shape& shape,
                           double zc, const flow::Array2D& w, double mu, double body_w) {
  double f = 0.0;
  for (const auto& piece : geom.pieces) {
    if (piece.length <= 0.0) continue;
    const auto n = shape.outward_normal(piece.r, piece.z - zc);
    // 국소 좌표: d = 표면까지의 거리(형상이 주는 값), t = 접선 방향 거리.
    //   u_z − W ≈ d (c1 + c2 t + c3 d + c4 d² + c5 t²) 를 근처 유체 쪽 z 면 값으로 최소제곱 → ∂u_z/∂n = c1
    // 모든 항이 d 를 곱하므로 실제 (곡면) 벽 위에서 u_z = W (no-slip) 를 정확히 만족한다.
    constexpr int kN = 5;
    double m[kN * kN] = {0}, rhs[kN] = {0};
    int count = 0;
    const int i0 = g.cell_r(piece.r);
    const int j0 = g.cell_z(piece.z);
    const double h = std::max(g.drc(std::clamp(i0, 0, g.nr - 1)), g.dzc(std::clamp(j0, 0, g.nz - 1)));
    const double reach = 3.0 * h;
    for (int i = std::max(0, i0 - 4); i <= std::min(g.nr - 1, i0 + 4); ++i) {
      for (int j = std::max(0, j0 - 4); j <= std::min(g.nz, j0 + 4); ++j) {
        const double r = g.rc(i), z = g.zf(j);
        if (shape.inside(r, z - zc)) continue;  // 물체 안의 면은 쓰지 않는다
        const double dr = r - piece.r, dz = z - piece.z;
        const double dn = shape.distance(r, z - zc), tn = -dr * n[1] + dz * n[0];
        if (dn <= 0.0 || std::hypot(dr, dz) > reach) continue;
        const double v[kN] = {dn, dn * tn, dn * dn, dn * dn * dn, dn * tn * tn};
        for (int a = 0; a < kN; ++a) {
          for (int b = 0; b < kN; ++b) m[a * kN + b] += v[a] * v[b];
          rhs[a] += v[a] * (w(i, j) - body_w);
        }
        ++count;
      }
    }
    if (count < 2 * kN) continue;
    // 늘어나는 격자에서 한 방향 칸만 얇으면 이웃 점이 거의 한 줄이 되어 5 항 맞춤이 정해지지 않는다
    // (힘이 한 단계 10^28 N 까지 튀었음). 그때는 1 항 맞춤 u_z − W ≈ c1 d 로 ∂u_z/∂n 을 구한다.
    constexpr int kPower[kN] = {1, 2, 2, 3, 3};
    double dudn = 0.0;
    if (scaled_pivot_ratio(kN, m, kPower, h) > 1e-9) {
      if (!solve_small(kN, m, rhs)) continue;
      dudn = rhs[0];
    } else {
      dudn = rhs[0] / m[0];  // Σ d (u − W) / Σ d²
    }
    f += mu * dudn * g.metric(piece.r) * piece.length;
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

double average_surface_pressure(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p,
                                const std::function<bool(const BoundaryPiece&)>& select) {
  double sp = 0.0, sa = 0.0;
  for (const auto& piece : geom.pieces) {
    if (!select(piece)) continue;
    const double w = piece.nz_area;  // 이미 r 가중
    sp += extrapolate(g, geom, p, piece.i, piece.j, piece.r, piece.z) * w;
    sa += w;
  }
  return sa != 0.0 ? sp / sa : 0.0;
}

double average_surface_pressure_by_length(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p,
                                          const std::function<bool(const BoundaryPiece&)>& select) {
  double sp = 0.0, sa = 0.0;
  for (const auto& piece : geom.pieces) {
    if (!select(piece) || !(piece.length > 0.0)) continue;
    const double w = piece.length * g.metric(piece.r);
    sp += extrapolate(g, geom, p, piece.i, piece.j, piece.r, piece.z) * w;
    sa += w;
  }
  return sa > 0.0 ? sp / sa : 0.0;
}

double max_surface_pressure(const flow::Grid& g, const CutGeometry& geom, const flow::Array2D& p) {
  double m = -1e300;
  for (const auto& piece : geom.pieces) {
    if (piece.length <= 0.0 && piece.nz_area == 0.0) continue;
    m = std::max(m, extrapolate(g, geom, p, piece.i, piece.j, piece.r, piece.z));
  }
  return m;
}

double axial_force_control_surface(const flow::Grid& g, const flow::Array2D& u,
                                   const flow::Array2D& w, const flow::Array2D& p, double mu,
                                   double rho, double rs, double z_bottom, double z_top) {
  auto face = [](int k, const char* what) {
    if (k < 0) throw std::invalid_argument(std::string("control surface is not on a grid face: ") + what);
    return k;
  };
  const int is = face(g.face_r(rs), "r");
  const int jb = face(g.face_z(z_bottom), "z bottom");
  const int jt = face(g.face_z(z_top), "z top");
  if (is < 1 || is > g.nr - 1 || jb < 1 || jt > g.nz - 1 || jb >= jt) {
    throw std::invalid_argument("control surface must lie strictly inside the grid");
  }
  double f = 0.0;
  // 윗면 (n_z = +1), 아랫면 (n_z = −1): sgn · (−p + τ_zz − ρ w²) · r dr
  for (int i = 0; i < is; ++i) {
    for (int sgn : {+1, -1}) {
      const int j = sgn > 0 ? jt : jb;
      const double pf = 0.5 * (p(i, j - 1) + p(i, j));
      const double tzz = 2.0 * mu * (w(i, j + 1) - w(i, j - 1)) / (g.dzc(j - 1) + g.dzc(j));
      f += sgn * (-pf + tzz - rho * w(i, j) * w(i, j)) * g.metric(g.rc(i)) * g.drc(i);
    }
  }
  // 옆면 (n_r = +1): τ_rz − ρ u_z u_r, r 면 위치, 셀 중심 높이
  for (int j = jb; j < jt; ++j) {
    const double dudz = (u(is, j + 1) - u(is, j - 1)) / (g.dzn(j) + g.dzn(j + 1));
    const double dwdr =
        0.5 * ((w(is, j) - w(is - 1, j)) + (w(is, j + 1) - w(is - 1, j + 1))) / g.drn(is);
    const double wc = 0.25 * (w(is - 1, j) + w(is, j) + w(is - 1, j + 1) + w(is, j + 1));
    f += (mu * (dudz + dwdr) - rho * wc * u(is, j)) * g.metric(g.rf(is)) * g.dzc(j);
  }
  return g.geometry == flow::Geometry::Axisymmetric ? kTwoPi * f : f;
}

}  // namespace hcw::body
