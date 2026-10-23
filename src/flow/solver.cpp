#include "flow/solver.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <stdexcept>
#include <string>
#include <utility>

#include "common/parallel.h"
#include "flow/linear.h"

namespace hcw::flow {
namespace {

std::array<double, 2> wall_velocity(const Boundary& b, double r, double z, double t) {
  if (!b.velocity) return {0.0, 0.0};
  return b.velocity(r, z, t);
}

// 압력 행렬 한 행 (열 번호 오름차순, 최대 kMax 개). 예전의 std::map<int, double> 과 같은 값을
// (같은 칸에 같은 순서로 더하므로 비트 단위로 같게) 힙 할당 없이 만든다.
struct SparseRow {
  static constexpr int kMax = 8;
  int n = 0;
  std::array<int, kMax> col{};
  std::array<double, kMax> val{};
  double& operator[](int c) {
    int k = 0;
    while (k < n && col[k] < c) ++k;
    if (k < n && col[k] == c) return val[k];
    if (n == kMax) throw std::logic_error("pressure matrix: too many entries in a row");
    for (int q = n; q > k; --q) {
      col[q] = col[q - 1];
      val[q] = val[q - 1];
    }
    col[k] = c;
    val[k] = 0.0;
    ++n;
    return val[k];
  }
};

void remove_mean(std::vector<double>& b) {
  double s = 0.0;
  for (double v : b) s += v;
  const double mean = s / static_cast<double>(b.size());
  for (double& v : b) v -= mean;
}

}  // namespace

FlowSolver::FlowSolver(Problem problem, SolverOptions options)
    : p_(std::move(problem)), opt_(options) {
  const Grid& g = p_.grid;
  const Boundaries& bc = p_.bc;
  if (!p_.material) throw std::invalid_argument("problem has no material");

  periodic_r_ = bc.r_min.kind == BoundaryKind::Periodic;
  periodic_z_ = bc.z_min.kind == BoundaryKind::Periodic;
  if (periodic_r_ != (bc.r_max.kind == BoundaryKind::Periodic) ||
      periodic_z_ != (bc.z_max.kind == BoundaryKind::Periodic)) {
    throw std::invalid_argument("periodic boundaries must be paired");
  }
  if (bc.r_max.kind == BoundaryKind::Axis || bc.z_min.kind == BoundaryKind::Axis ||
      bc.z_max.kind == BoundaryKind::Axis) {
    throw std::invalid_argument("axis boundary is only allowed at r_min");
  }
  const bool axisym = g.geometry == Geometry::Axisymmetric;
  // 평면 격자의 Axis 는 대칭면 (P4 3-A 쐐기): 법선 속도 0, 접선 속도는 짝함수 (미끄럼). 축대칭과 같은 처리.
  if (bc.r_min.kind == BoundaryKind::Axis && g.r0 != 0.0) {
    throw std::invalid_argument("axis (or planar symmetry) boundary requires a grid starting at r = 0");
  }
  if (axisym && g.r0 == 0.0 && bc.r_min.kind != BoundaryKind::Axis) {
    throw std::invalid_argument("axisymmetric grid starting at r = 0 needs an axis boundary");
  }
  if (axisym && periodic_r_) throw std::invalid_argument("radial periodicity is not allowed");
  if ((periodic_r_ && g.nr < 3) || (periodic_z_ && g.nz < 3)) {
    throw std::invalid_argument("periodic directions need at least 3 cells");
  }

  iu_lo_ = periodic_r_ ? 0 : 1;
  iu_hi_ = g.nr;
  jw_lo_ = periodic_z_ ? 0 : 1;
  jw_hi_ = g.nz;

  u_ = Array2D(g.nr + 1, g.nz);
  w_ = Array2D(g.nr, g.nz + 1);
  su_ = Array2D(g.nr + 1, g.nz);
  sw_ = Array2D(g.nr, g.nz + 1);
  for (auto& phi : phi_) {
    phi.hi.assign(static_cast<std::size_t>(g.nr) * g.nz, 0.0);
    phi.lo.assign(static_cast<std::size_t>(g.nr) * g.nz, 0.0);
  }

  p_inc_ = Array2D(g.nr, g.nz);
  geom_ = body::open_geometry(g);
  geom_valid_ = !p_.body;
  if (p_.body) {
    if (!p_.body->shape) throw std::invalid_argument("body has no shape");
    // 축대칭은 축에서 시작하는 격자만. 평면 격자도 허용한다 (P4 3-A 쐐기의 반쪽, x = 0 대칭면).
    if (axisym && g.r0 != 0.0) {
      throw std::invalid_argument("bodies on axisymmetric grids need a grid starting at the axis");
    }
    if (periodic_z_) throw std::invalid_argument("bodies are not supported with periodic z");
  }
}

void FlowSolver::update_geometry(double t) {
  if (!p_.body) return;
  if (geom_valid_ && (!p_.body->moving() || geom_time_ == t)) return;
  geom_ = body::cut_geometry(p_.grid, *p_.body->shape, p_.body->z(t));
  geom_time_ = t;
  geom_valid_ = true;
  // 압력 정칙화 셀(마지막 셀)은 유체여야 한다 (build_operator 참조)
  if (!geom_.active.back()) {
    throw std::invalid_argument("the last cell (pressure reference) must be in the fluid");
  }
}

int FlowSolver::wrap_r(int i) const {
  if (!periodic_r_) return i;
  const int n = p_.grid.nr;
  return ((i % n) + n) % n;
}

int FlowSolver::wrap_z(int j) const {
  if (!periodic_z_) return j;
  const int n = p_.grid.nz;
  return ((j % n) + n) % n;
}

double FlowSolver::U(const Array2D& u, int i, int j, double t) const {
  const Grid& g = p_.grid;
  if (periodic_r_) i = wrap_r(i);
  if (j >= 0 && j < g.nz) return u(i, j);
  if (periodic_z_) return u(i, wrap_z(j));
  // 벽: 접선 속도의 유령값 (벽에서 평균이 벽 속도가 되도록)
  if (j < 0) return 2.0 * wall_velocity(p_.bc.z_min, g.rf(i), g.z0, t)[0] - u(i, 0);
  return 2.0 * wall_velocity(p_.bc.z_max, g.rf(i), g.z1, t)[0] - u(i, g.nz - 1);
}

double FlowSolver::W(const Array2D& w, int i, int j, double t) const {
  const Grid& g = p_.grid;
  if (periodic_z_) j = wrap_z(j);
  if (i >= 0 && i < g.nr) return w(i, j);
  if (periodic_r_) return w(wrap_r(i), j);
  if (i < 0) {
    if (p_.bc.r_min.kind == BoundaryKind::Axis) return w(0, j);  // 축에 대해 짝함수
    return 2.0 * wall_velocity(p_.bc.r_min, g.r0, g.zf(j), t)[1] - w(0, j);
  }
  return 2.0 * wall_velocity(p_.bc.r_max, g.r1, g.zf(j), t)[1] - w(g.nr - 1, j);
}

void FlowSolver::apply_boundary(Array2D& u, Array2D& w, double t) const {
  const Grid& g = p_.grid;
  for (int j = 0; j < g.nz; ++j) {
    if (periodic_r_) {
      u(g.nr, j) = u(0, j);
    } else {
      u(0, j) = p_.bc.r_min.kind == BoundaryKind::Axis
                    ? 0.0
                    : wall_velocity(p_.bc.r_min, g.r0, g.zc(j), t)[0];
      u(g.nr, j) = wall_velocity(p_.bc.r_max, g.r1, g.zc(j), t)[0];
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    if (periodic_z_) {
      w(i, g.nz) = w(i, 0);
    } else {
      w(i, 0) = wall_velocity(p_.bc.z_min, g.rc(i), g.z0, t)[1];
      w(i, g.nz) = wall_velocity(p_.bc.z_max, g.rc(i), g.z1, t)[1];
    }
  }
}

void FlowSolver::material_changed() {
  mat_valid_ = false;
  op_.valid = false;
  op_mom_.valid = false;
  visc_valid_ = false;
}

const MaterialArrays& FlowSolver::material_at(double t) {
  if (!mat_valid_ || (p_.material->time_dependent() && t != mat_time_)) {
    p_.material->evaluate(p_.grid, t, mat_);
    mat_time_ = t;
    mat_valid_ = true;
  }
  return mat_;
}

FlowSolver::Operator& FlowSolver::operator_at(double t_rho, double t_geom) {
  update_geometry(t_geom);
  const bool moving = p_.body && p_.body->moving();
  const bool fresh = op_.valid && (!p_.material->time_dependent() || op_.time == t_rho) &&
                     (!moving || op_.geometry_time == t_geom);
  if (!fresh) {
    build_operator(material_at(t_rho), op_, geom_);
    op_.time = t_rho;
    op_.geometry_time = t_geom;
    op_.valid = true;
  }
  return op_;
}

void FlowSolver::build_operator(const MaterialArrays& m, Operator& op,
                                const body::CutGeometry& gm) const {
  const Grid& g = p_.grid;
  const int n = g.nr * g.nz;
  std::vector<SparseRow> rows(static_cast<std::size_t>(n));
  auto cell = [&](int i, int j) { return wrap_r(i) * g.nz + wrap_z(j); };
  auto couple = [&](int a, int b, double coef) {
    rows[a][a] += coef;
    rows[b][b] += coef;
    rows[a][b] -= coef;
    rows[b][a] -= coef;
  };

  // 컷셀: 면 계수에 열린 비율을 곱한다. 막힌 면(0)은 연결하지 않는다.
  // 물체가 없으면 비율이 정확히 1 이라 P1 과 같은 값이 된다.
  for (int i = iu_lo_; i < iu_hi_; ++i) {
    const double area = g.metric(g.rf(i)) * g.dz;
    for (int j = 0; j < g.nz; ++j) {
      const double open = gm.open_u(i, j);
      if (open > 0.0) couple(cell(i - 1, j), cell(i, j), open * (area / (m.rho_u(i, j) * g.dr)));
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    const double area = g.metric(g.rc(i)) * g.dr;
    for (int j = jw_lo_; j < jw_hi_; ++j) {
      const double open = gm.open_w(i, j);
      if (open > 0.0) couple(cell(i, j - 1), cell(i, j), open * (area / (m.rho_w(i, j) * g.dz)));
    }
  }
  // 고체로 가득 찬 셀은 다른 셀과 연결되지 않는 항등 행 (해 0)
  for (int c = 0; c < n; ++c) {
    if (!gm.active[static_cast<std::size_t>(c)]) rows[c][c] = 1.0;
  }
  // 압력은 상수만큼 미정 → 마지막 셀의 대각을 두 배로 해 정칙화한다.
  // 우변의 합이 0 이면 그 셀의 해가 0 이 되며 나머지 해는 바뀌지 않는다.
  // 마지막 셀을 고르는 이유: 불완전 Cholesky 의 마지막 피벗이 특이성 때문에
  // 0 에 가까워지는데, 그 피벗을 직접 보강해야 큰 밀도비(1000)에서도 무너지지 않는다.
  if (p_.air_kappa > 0.0) {
    // 압축성 공기: ∇·(∇φ/ρ) − ψ φ/Δt² = ∇·u*  →  대각에 ψ V/Δt² 를 더한다.
    // 이 항이 있으면 행렬이 정칙이라 정칙화(고정 셀)를 하지 않는다 (절대 압력).
    for (int c = 0; c < n && static_cast<std::size_t>(c) < psi_vol_.size(); ++c) {
      if (gm.active[static_cast<std::size_t>(c)]) rows[c][c] += psi_vol_[static_cast<std::size_t>(c)];
    }
  } else {
    rows[n - 1][n - 1] *= 2.0;
  }

  CsrMatrix& a = op.matrix;
  a.n = n;
  a.row_start.assign(static_cast<std::size_t>(n) + 1, 0);
  a.col.clear();
  a.val.clear();
  for (int r = 0; r < n; ++r) {
    for (int k = 0; k < rows[r].n; ++k) {
      a.col.push_back(rows[r].col[k]);
      a.val.push_back(rows[r].val[k]);
    }
    a.row_start[r + 1] = static_cast<int>(a.col.size());
  }
  op.preconditioner.factor(a, g.nz);  // 셀 순서 c = i·nz + j (다중 스레드 삼각 풀이용)
}

void FlowSolver::face_flux_rhs(const Array2D& fu, const Array2D& fw, double body_rate,
                               std::vector<double>& b) const {
  const Grid& g = p_.grid;
  const Array2D& ou = geom_.open_u;
  const Array2D& ow = geom_.open_w;
  b.assign(static_cast<std::size_t>(g.nr) * g.nz, 0.0);
  for (int i = 0; i < g.nr; ++i) {
    const double aw = g.metric(g.rf(i)) * g.dz;
    const double ae = g.metric(g.rf(i + 1)) * g.dz;
    const double az = g.metric(g.rc(i)) * g.dr;
    for (int j = 0; j < g.nz; ++j) {
      const std::size_t c = static_cast<std::size_t>(i) * g.nz + j;
      if (!geom_.active[c]) continue;
      const double out = ae * ou(i + 1, j) * fu(i + 1, j) - aw * ou(i, j) * fu(i, j) +
                         az * (ow(i, j + 1) * fw(i, j + 1) - ow(i, j) * fw(i, j));
      b[c] = -out;
    }
  }
  // 물체 경계 유출: 강체 z 병진 W 에 대해 W · ∫ n_z dA (body/cut_cell.h)
  for (const auto& piece : geom_.pieces) {
    b[static_cast<std::size_t>(piece.i) * g.nz + piece.j] -= body_rate * piece.nz_area;
  }
}

void FlowSolver::remove_active_mean(std::vector<double>& b) const {
  if (p_.air_kappa > 0.0) return;  // 압축성: 정칙 행렬이라 우변 평균을 빼지 않는다
  if (!p_.body) {
    remove_mean(b);
    return;
  }
  double s = 0.0;
  std::size_t count = 0;
  for (std::size_t c = 0; c < b.size(); ++c) {
    if (geom_.active[c]) {
      s += b[c];
      ++count;
    }
  }
  const double mean = s / static_cast<double>(count);
  for (std::size_t c = 0; c < b.size(); ++c) {
    if (geom_.active[c]) b[c] -= mean;
  }
}

void FlowSolver::record(const SolveResult& r, double tolerance) {
  stats_.poisson_solves += 1;
  stats_.poisson_iterations += r.iterations;
  stats_.poisson_refinements += r.refinements;
  stats_.max_poisson_iterations = std::max(stats_.max_poisson_iterations, r.iterations);
  stats_.max_relative_residual = std::max(stats_.max_relative_residual, r.relative_residual);
  if (!r.converged || r.relative_residual > tolerance) stats_.solver_failed = true;
}

void FlowSolver::project(Array2D& u, Array2D& w, double t_rho, double t_bc, Potential& phi,
                         const std::vector<double>* extra) {
  const Grid& g = p_.grid;
  Operator& op = operator_at(t_rho, t_bc);
  const MaterialArrays& m = material_at(t_rho);
  const double body_w = p_.body ? p_.body->w(t_bc) : 0.0;

  std::vector<double> b;
  face_flux_rhs(u, w, body_w, b);
  remove_active_mean(b);
  if (extra) {
    for (std::size_t c = 0; c < b.size(); ++c) b[c] += (*extra)[c];
  }
  record(refined_solve(op.matrix, op.preconditioner, b, phi.hi, phi.lo,
                       opt_.projection_tolerance, opt_.max_iterations),
         opt_.projection_tolerance);

  // 두 셀 사이 퍼텐셜 차. hi 와 lo 를 따로 빼야 두 배 정밀도가 유지된다.
  auto diff = [&](int i1, int j1, int i0, int j0) {
    const std::size_t a = static_cast<std::size_t>(wrap_r(i1)) * g.nz + wrap_z(j1);
    const std::size_t c = static_cast<std::size_t>(wrap_r(i0)) * g.nz + wrap_z(j0);
    return (phi.hi[a] - phi.hi[c]) + (phi.lo[a] - phi.lo[c]);
  };
  for (int i = iu_lo_; i < iu_hi_; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      u(i, j) -= diff(i, j, i - 1, j) / (m.rho_u(i, j) * g.dr);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = jw_lo_; j < jw_hi_; ++j) {
      w(i, j) -= diff(i, j, i, j - 1) / (m.rho_w(i, j) * g.dz);
    }
  }
  if (periodic_r_) {
    for (int j = 0; j < g.nz; ++j) u(g.nr, j) = u(0, j);
  }
  if (periodic_z_) {
    for (int i = 0; i < g.nr; ++i) w(i, g.nz) = w(i, 0);
  }
  // 물체 안에 완전히 묻힌 면은 물체 속도 (r 성분 0, z 성분 W)
  if (p_.body) {
    for (int i = iu_lo_; i < iu_hi_; ++i) {
      for (int j = 0; j < g.nz; ++j) {
        if (geom_.open_u(i, j) == 0.0) u(i, j) = 0.0;
      }
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = jw_lo_; j < jw_hi_; ++j) {
        if (geom_.open_w(i, j) == 0.0) w(i, j) = body_w;
      }
    }
  }
}

void FlowSolver::rhs(const Array2D& u, const Array2D& w, double t, Array2D& fu, Array2D& fw,
                     bool convection, bool viscous, bool forces) {
  const Grid& g = p_.grid;
  const MaterialArrays& m = material_at(t);
  const bool axisym = g.geometry == Geometry::Axisymmetric;
  const bool has_source = forces && static_cast<bool>(p_.source);
  if (has_source) p_.source(t, su_, sw_);
  if (fu.ni() != g.nr + 1) fu = Array2D(g.nr + 1, g.nz);
  if (fw.ni() != g.nr) fw = Array2D(g.nr, g.nz + 1);
  fu.fill(0.0);
  fw.fill(0.0);

  // 제한자 대류 (Convection::VanLeer): 면 값 = 상류 값 + ½·(van Leer 제한 기울기).
  // qll, ql: 상류 쪽 두 값, qr: 하류 값 (속도 a 의 방향 기준으로 이미 고른 순서)
  const bool consistent = opt_.convection == Convection::Consistent;
  const bool limited = opt_.convection == Convection::VanLeer || consistent;
  auto face_value = [](double qll, double ql, double qr) {
    const double d1 = ql - qll, d2 = qr - ql;
    const double s = d1 * d2 > 0.0 ? 2.0 * d1 * d2 / (d1 + d2) : 0.0;
    return ql + 0.5 * s;
  };
  // a 방향의 상류 쪽으로 치우친 면 값. 면은 q0–q1 사이, qm 은 q0 뒤, qp 는 q1 앞.
  // 필요한 상류 두 번째 값이 격자 밖이면 중심 평균.
  auto upwind = [&](double a, double qm, bool has_m, double q0, double q1, double qp,
                    bool has_p) {
    if (a >= 0.0) return has_m ? face_value(qm, q0, q1) : 0.5 * (q0 + q1);
    return has_p ? face_value(qp, q1, q0) : 0.5 * (q0 + q1);
  };
  // 넓은 스텐실의 접근 가능 여부 (U/W 는 유령 한 줄까지만 정의)
  auto u_r_ok = [&](int i) { return periodic_r_ || (i >= 0 && i <= g.nr); };
  auto u_z_ok = [&](int j) { return periodic_z_ || (j >= -1 && j <= g.nz); };
  auto w_r_ok = [&](int i) { return periodic_r_ || (i >= -1 && i <= g.nr); };
  auto w_z_ok = [&](int j) { return periodic_z_ || (j >= 0 && j <= g.nz); };
  // 운동량 제어체적 밀도 (유령 위치는 가장 가까운 안쪽 값, 주기면 감음)
  auto RU = [&](int i, int j) {
    i = periodic_r_ ? wrap_r(i) : std::clamp(i, 0, g.nr);
    j = periodic_z_ ? wrap_z(j) : std::clamp(j, 0, g.nz - 1);
    return m.rho_u(i, j);
  };
  auto RW = [&](int i, int j) {
    i = periodic_r_ ? wrap_r(i) : std::clamp(i, 0, g.nr - 1);
    j = periodic_z_ ? wrap_z(j) : std::clamp(j, 0, g.nz);
    return m.rho_w(i, j);
  };
  // 일관 수송: 네 면의 (계수, 이동 속도 a, 제한 밀도, 상류 밀도, 수송 속도) → 가속도
  struct Face {
    double coef, a, rho_lim, rho_up, q;
  };
  auto consistent_accel = [&](const Face (&f)[4], double rho_c, double q_c) {
    double num = 0.0, div = 0.0;
    bool heavy_in = false;  // 2 배 넘게 무거운 유체가 들어오는 면 (coef·a < 0 이 유입)
    for (const Face& x : f) {
      num += x.coef * x.rho_lim * x.a * (x.q - q_c);
      div += x.coef * x.rho_lim * x.a;
      heavy_in = heavy_in || (x.coef * x.a < 0.0 && x.rho_lim > 2.0 * rho_c);
    }
    // 밀도가 매끈한 곳은 ρ 로 나눈다 (반이산 질량 가중 이류 a·∇u, 시간 2차 유지).
    // 갱신 밀도 ρ* 로 나누면 시간 1차 오차 (Δt ∇·ṁ / ρ) 가 생기므로 무거운 유입이 있는 칸에만 쓴다.
    if (!heavy_in) return num / rho_c;
    double den = rho_c - conv_dt_ * div;
    if (den < 0.5 * rho_c) {
      num = div = 0.0;
      for (const Face& x : f) {
        num += x.coef * x.rho_up * x.a * (x.q - q_c);
        div += x.coef * x.rho_up * x.a;
      }
      den = rho_c - conv_dt_ * div;
    }
    return num / den;
  };

  // 꼭짓점 (r 면 fi, z 면 fj) 의 전단응력 τ_rz = μ (∂u/∂z + ∂w/∂r)
  auto tau_rz = [&](int fi, int fj) {
    return m.mu_corner(fi, fj) * ((U(u, fi, fj, t) - U(u, fi, fj - 1, t)) / g.dz +
                                  (W(w, fi, fj, t) - W(w, fi - 1, fj, t)) / g.dr);
  };

  // r 방향 운동량 (r 면)
  // 면마다 따로 계산하므로 r 열(i) 단위로 스레드에 나눠도 결과가 같다.
  parallel::parallel_for(iu_lo_, iu_hi_, parallel::kRowGrain, [&](int i_lo, int i_hi) {
    for (int i = i_lo; i < i_hi; ++i) {
      const double m_f = g.metric(g.rf(i));
      const double m_l = g.metric(g.rc(i - 1));
      const double m_r = g.metric(g.rc(i));
      for (int j = 0; j < g.nz; ++j) {
        const double u_l = 0.5 * (U(u, i - 1, j, t) + U(u, i, j, t));
        const double u_r = 0.5 * (U(u, i, j, t) + U(u, i + 1, j, t));
        double conv_r, conv_z;
        if (consistent) {
          // r 면 (셀 c 중심): 면 c 와 c+1 사이. z 면 (꼭짓점 fj): u 면 fj−1 과 fj 사이.
          auto r_face = [&](int c, double a, double coef) {
            const bool hm = u_r_ok(c - 1), hp = u_r_ok(c + 2);
            return Face{coef, a,
                        upwind(a, RU(c - 1, j), hm, RU(c, j), RU(c + 1, j), RU(c + 2, j), hp),
                        a >= 0.0 ? RU(c, j) : RU(c + 1, j),
                        upwind(a, hm ? U(u, c - 1, j, t) : 0.0, hm, U(u, c, j, t), U(u, c + 1, j, t),
                               hp ? U(u, c + 2, j, t) : 0.0, hp)};
          };
          auto z_face = [&](int fj, double coef) {
            const double a = 0.5 * (W(w, i - 1, fj, t) + W(w, i, fj, t));
            const bool hm = u_z_ok(fj - 2), hp = u_z_ok(fj + 1);
            return Face{coef, a,
                        upwind(a, RU(i, fj - 2), hm, RU(i, fj - 1), RU(i, fj), RU(i, fj + 1), hp),
                        a >= 0.0 ? RU(i, fj - 1) : RU(i, fj),
                        upwind(a, hm ? U(u, i, fj - 2, t) : 0.0, hm, U(u, i, fj - 1, t), U(u, i, fj, t),
                               hp ? U(u, i, fj + 1, t) : 0.0, hp)};
          };
          const Face f[4] = {r_face(i, u_r, m_r / (m_f * g.dr)), r_face(i - 1, u_l, -m_l / (m_f * g.dr)),
                             z_face(j + 1, 1.0 / g.dz), z_face(j, -1.0 / g.dz)};
          conv_r = consistent_accel(f, m.rho_u(i, j), u(i, j));
          conv_z = 0.0;
        } else if (!limited) {
          conv_r = (m_r * u_r * u_r - m_l * u_l * u_l) / (m_f * g.dr);
          auto flux_z = [&](int fj) {
            const double uc = 0.5 * (U(u, i, fj - 1, t) + U(u, i, fj, t));
            const double wc = 0.5 * (W(w, i - 1, fj, t) + W(w, i, fj, t));
            return uc * wc;
          };
          conv_z = (flux_z(j + 1) - flux_z(j)) / g.dz;
        } else {
          // 셀 c (면 c 와 c+1 사이) 의 r 플럭스: 이동 속도 = 중심 평균, 수송 u = 제한 보간
          auto flux_r = [&](int c, double a) {
            const double q = upwind(a, u_r_ok(c - 1) ? U(u, c - 1, j, t) : 0.0, u_r_ok(c - 1),
                                    U(u, c, j, t), U(u, c + 1, j, t),
                                    u_r_ok(c + 2) ? U(u, c + 2, j, t) : 0.0, u_r_ok(c + 2));
            return a * q;
          };
          conv_r = (m_r * flux_r(i, u_r) - m_l * flux_r(i - 1, u_l)) / (m_f * g.dr);
          auto flux_z = [&](int fj) {
            const double wc = 0.5 * (W(w, i - 1, fj, t) + W(w, i, fj, t));
            const double q = upwind(wc, u_z_ok(fj - 2) ? U(u, i, fj - 2, t) : 0.0, u_z_ok(fj - 2),
                                    U(u, i, fj - 1, t), U(u, i, fj, t),
                                    u_z_ok(fj + 1) ? U(u, i, fj + 1, t) : 0.0, u_z_ok(fj + 1));
            return q * wc;
          };
          conv_z = (flux_z(j + 1) - flux_z(j)) / g.dz;
        }

        auto tau_rr = [&](int c) {
          return 2.0 * m.mu_cell(wrap_r(c), j) * (U(u, c + 1, j, t) - U(u, c, j, t)) / g.dr;
        };
        const double visc_r = (m_r * tau_rr(i) - m_l * tau_rr(i - 1)) / (m_f * g.dr);
        const double visc_z = (tau_rz(i, j + 1) - tau_rz(i, j)) / g.dz;
        const double r_f = g.rf(i);
        const double hoop = axisym ? -2.0 * m.mu_u(i, j) * u(i, j) / (r_f * r_f) : 0.0;

        // 항을 뺄 때는 0 으로 바꿀 뿐이라 전부 넣을 때(P1)의 계산 순서는 그대로다.
        double f = (convection ? -(conv_r + conv_z) : 0.0) +
                   (viscous ? (visc_r + visc_z + hoop) / m.rho_u(i, j) : 0.0) +
                   (forces ? p_.gravity[0] : 0.0);
        if (has_source) f += su_(i, j);
        fu(i, j) = f;
      }
    }
  });

  // z 방향 운동량 (z 면)
  parallel::parallel_for(0, g.nr, parallel::kRowGrain, [&](int i_lo, int i_hi) {
    for (int i = i_lo; i < i_hi; ++i) {
      const double m_c = g.metric(g.rc(i));
      const double m_w = g.metric(g.rf(i));
      const double m_e = g.metric(g.rf(i + 1));
      for (int j = jw_lo_; j < jw_hi_; ++j) {
        const double w_s = 0.5 * (W(w, i, j - 1, t) + W(w, i, j, t));
        const double w_n = 0.5 * (W(w, i, j, t) + W(w, i, j + 1, t));
        double conv_r, conv_z;
        if (consistent) {
          // r 면 (꼭짓점 fi): w 면 fi−1 과 fi 사이. z 면 (셀 c 중심): w 면 c 와 c+1 사이.
          auto r_face = [&](int fi, double coef) {
            const double a = 0.5 * (U(u, fi, j - 1, t) + U(u, fi, j, t));
            const bool hm = w_r_ok(fi - 2), hp = w_r_ok(fi + 1);
            return Face{coef, a,
                        upwind(a, RW(fi - 2, j), hm, RW(fi - 1, j), RW(fi, j), RW(fi + 1, j), hp),
                        a >= 0.0 ? RW(fi - 1, j) : RW(fi, j),
                        upwind(a, hm ? W(w, fi - 2, j, t) : 0.0, hm, W(w, fi - 1, j, t), W(w, fi, j, t),
                               hp ? W(w, fi + 1, j, t) : 0.0, hp)};
          };
          auto z_face = [&](int c, double a, double coef) {
            const bool hm = w_z_ok(c - 1), hp = w_z_ok(c + 2);
            return Face{coef, a,
                        upwind(a, RW(i, c - 1), hm, RW(i, c), RW(i, c + 1), RW(i, c + 2), hp),
                        a >= 0.0 ? RW(i, c) : RW(i, c + 1),
                        upwind(a, hm ? W(w, i, c - 1, t) : 0.0, hm, W(w, i, c, t), W(w, i, c + 1, t),
                               hp ? W(w, i, c + 2, t) : 0.0, hp)};
          };
          const Face f[4] = {r_face(i + 1, m_e / (m_c * g.dr)), r_face(i, -m_w / (m_c * g.dr)),
                             z_face(j, w_n, 1.0 / g.dz), z_face(j - 1, w_s, -1.0 / g.dz)};
          conv_r = consistent_accel(f, m.rho_w(i, j), w(i, j));
          conv_z = 0.0;
        } else if (!limited) {
          auto flux_r = [&](int fi) {
            const double uc = 0.5 * (U(u, fi, j - 1, t) + U(u, fi, j, t));
            const double wc = 0.5 * (W(w, fi - 1, j, t) + W(w, fi, j, t));
            return uc * wc;
          };
          conv_r = (m_e * flux_r(i + 1) - m_w * flux_r(i)) / (m_c * g.dr);
          conv_z = (w_n * w_n - w_s * w_s) / g.dz;
        } else {
          auto flux_r = [&](int fi) {
            const double uc = 0.5 * (U(u, fi, j - 1, t) + U(u, fi, j, t));
            const double q = upwind(uc, w_r_ok(fi - 2) ? W(w, fi - 2, j, t) : 0.0, w_r_ok(fi - 2),
                                    W(w, fi - 1, j, t), W(w, fi, j, t),
                                    w_r_ok(fi + 1) ? W(w, fi + 1, j, t) : 0.0, w_r_ok(fi + 1));
            return uc * q;
          };
          conv_r = (m_e * flux_r(i + 1) - m_w * flux_r(i)) / (m_c * g.dr);
          // 셀 c (z 면 c 와 c+1 사이) 의 z 플럭스
          auto flux_z = [&](int c, double a) {
            const double q = upwind(a, w_z_ok(c - 1) ? W(w, i, c - 1, t) : 0.0, w_z_ok(c - 1),
                                    W(w, i, c, t), W(w, i, c + 1, t),
                                    w_z_ok(c + 2) ? W(w, i, c + 2, t) : 0.0, w_z_ok(c + 2));
            return a * q;
          };
          conv_z = (flux_z(j, w_n) - flux_z(j - 1, w_s)) / g.dz;
        }

        const double visc_r = (m_e * tau_rz(i + 1, j) - m_w * tau_rz(i, j)) / (m_c * g.dr);
        auto tau_zz = [&](int c) {
          return 2.0 * m.mu_cell(i, wrap_z(c)) * (W(w, i, c + 1, t) - W(w, i, c, t)) / g.dz;
        };
        const double visc_z = (tau_zz(j) - tau_zz(j - 1)) / g.dz;

        double f = (convection ? -(conv_r + conv_z) : 0.0) +
                   (viscous ? (visc_r + visc_z) / m.rho_w(i, j) : 0.0) +
                   (forces ? p_.gravity[1] : 0.0);
        if (has_source) f += sw_(i, j);
        fw(i, j) = f;
      }
    }
  });

  if (periodic_r_) {
    for (int j = 0; j < g.nz; ++j) fu(g.nr, j) = fu(0, j);
  }
  if (periodic_z_) {
    for (int i = 0; i < g.nr; ++i) fw(i, g.nz) = fw(i, 0);
  }
}

void FlowSolver::initialize(double t0, const std::function<double(double, double)>& u0,
                            const std::function<double(double, double)>& w0, bool do_project) {
  const Grid& g = p_.grid;
  t_ = t0;
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) u_(i, j) = u0(g.rf(i), g.zc(j));
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) w_(i, j) = w0(g.rc(i), g.zf(j));
  }
  apply_boundary(u_, w_, t0);
  hist_valid_ = false;
  if (do_project) project(u_, w_, t0, t0, phi_[4]);
}

void FlowSolver::step(double dt) {
  if (p_.air_kappa > 0.0) {
    throw std::logic_error("compressible air needs step_implicit (absolute pressure is carried there)");
  }
  if (p_.body) {
    throw std::logic_error(
        "time stepping with a body is not implemented yet (P2: near-body momentum pending)");
  }
  const Grid& g = p_.grid;
  const double t = t_;
  const Array2D u0 = u_, w0 = w_;
  Array2D fu, fw;
  conv_dt_ = dt;  // 세 단계 모두 Δt 짜리 전진 오일러의 볼록 결합

  auto combine = [&](Array2D& u, Array2D& w, double a0, double a1) {
    // u ← a0 u0 + a1 (u + dt f), 계산하는 면만
    for (int i = iu_lo_; i < iu_hi_; ++i) {
      for (int j = 0; j < g.nz; ++j) u(i, j) = a0 * u0(i, j) + a1 * (u(i, j) + dt * fu(i, j));
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = jw_lo_; j < jw_hi_; ++j) {
        w(i, j) = a0 * w0(i, j) + a1 * (w(i, j) + dt * fw(i, j));
      }
    }
  };

  // 1단: t → t+dt (밀도 t)
  rhs(u_, w_, t, fu, fw);
  combine(u_, w_, 0.0, 1.0);
  apply_boundary(u_, w_, t + dt);
  project(u_, w_, t, t + dt, phi_[0]);
  // 2단: → t+dt/2 (밀도 t+dt)
  rhs(u_, w_, t + dt, fu, fw);
  combine(u_, w_, 0.75, 0.25);
  apply_boundary(u_, w_, t + 0.5 * dt);
  project(u_, w_, t + dt, t + 0.5 * dt, phi_[1]);
  // 3단: → t+dt (밀도 t+dt/2)
  rhs(u_, w_, t + 0.5 * dt, fu, fw);
  combine(u_, w_, 1.0 / 3.0, 2.0 / 3.0);
  apply_boundary(u_, w_, t + dt);
  project(u_, w_, t + 0.5 * dt, t + dt, phi_[2]);

  conv_dt_ = 0.0;
  t_ = t + dt;
  stats_.steps += 1;
}

void FlowSolver::compute_pressure(Array2D& p) {
  const Grid& g = p_.grid;
  Array2D fu, fw;
  rhs(u_, w_, t_, fu, fw);

  // 경계에서 지정된 법선 속도의 시간 변화율 (중심 차분)
  const double delta = 1e-6 * (1.0 + std::abs(t_));
  Array2D up = u_, wp = w_, um = u_, wm = w_;
  apply_boundary(up, wp, t_ + delta);
  apply_boundary(um, wm, t_ - delta);
  for (int j = 0; j < g.nz; ++j) {
    for (int i : {0, g.nr}) {
      if (!solves_u(i)) fu(i, j) = (up(i, j) - um(i, j)) / (2.0 * delta);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j : {0, g.nz}) {
      if (!solves_w(j)) fw(i, j) = (wp(i, j) - wm(i, j)) / (2.0 * delta);
    }
  }
  if (periodic_r_) {
    for (int j = 0; j < g.nz; ++j) fu(g.nr, j) = fu(0, j);
  }
  if (periodic_z_) {
    for (int i = 0; i < g.nr; ++i) fw(i, g.nz) = fw(i, 0);
  }

  // 물체가 있으면 경계 플럭스의 시간 변화율 = 물체 가속도
  Operator& op = operator_at(t_, t_);
  std::vector<double> b;
  face_flux_rhs(fu, fw, p_.body ? p_.body->a(t_) : 0.0, b);
  remove_active_mean(b);
  Potential& x = phi_[3];
  record(refined_solve(op.matrix, op.preconditioner, b, x.hi, x.lo, opt_.pressure_tolerance,
                       opt_.max_iterations),
         opt_.pressure_tolerance);

  p = Array2D(g.nr, g.nz);
  double sum = 0.0, vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      const std::size_t k = static_cast<std::size_t>(i) * g.nz + j;
      p(i, j) = x.hi[k] + x.lo[k];
      if (!geom_.active[k]) continue;  // 고체 셀은 평균에서 뺀다
      sum += p(i, j) * g.cell_volume(i);
      vol += g.cell_volume(i);
    }
  }
  const double mean = sum / vol;
  for (std::size_t k = 0; k < p.data().size(); ++k) {
    p.data()[k] = geom_.active[k] ? p.data()[k] - mean : 0.0;
  }
}

double FlowSolver::max_divergence() const {
  const Grid& g = p_.grid;
  // face_flux_rhs 와 같은 이산 발산 (b = −유출) 을 셀 체적으로 나눈 값
  std::vector<double> b;
  face_flux_rhs(u_, w_, p_.body ? p_.body->w(t_) : 0.0, b);
  double worst = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      worst = std::max(worst, std::abs(b[static_cast<std::size_t>(i) * g.nz + j] / g.cell_volume(i)));
    }
  }
  return worst;
}

void FlowSolver::viscous_matrix(double t, CsrMatrix& a) {
  const Grid& g = p_.grid;
  const int nz = g.nz;
  const int nu = (g.nr + 1) * nz;
  const int n = nu + g.nr * (nz + 1);
  Array2D zu(g.nr + 1, nz), zw(g.nr, nz + 1), fu0, fw0, fu, fw;
  rhs(zu, zw, t, fu0, fw0, false, true, false);  // 상수 부분 c (벽 속도)

  std::vector<std::vector<std::pair<int, double>>> cols(static_cast<std::size_t>(n));
  // 점성 이산식은 한 면에서 각 방향 ±1 칸 안의 면만 쓴다 → 3×3 색칠로 탐침이 겹치지 않는다.
  for (int comp = 0; comp < 2; ++comp) {
    for (int ci = 0; ci < 3; ++ci) {
      for (int cj = 0; cj < 3; ++cj) {
        Array2D pu(g.nr + 1, nz), pw(g.nr, nz + 1);
        Array2D& pert = comp == 0 ? pu : pw;
        for (int i = ci; i < pert.ni(); i += 3) {
          for (int j = cj; j < pert.nj(); j += 3) pert(i, j) = 1.0;
        }
        rhs(pu, pw, t, fu, fw, false, true, false);
        // 값이 바뀐 행마다, 같은 색의 섭동 면 중 ±1 칸 안의 유일한 면이 그 열이다.
        auto owner = [&](int ri, int rj) {
          for (int i = ri - 1; i <= ri + 1; ++i) {
            for (int j = rj - 1; j <= rj + 1; ++j) {
              if (i < 0 || j < 0 || i >= pert.ni() || j >= pert.nj()) continue;
              if (i % 3 == ci && j % 3 == cj) return comp == 0 ? i * nz + j : nu + i * (nz + 1) + j;
            }
          }
          throw std::logic_error("viscous probe: no owner column");
        };
        for (int i = 0; i <= g.nr; ++i) {
          for (int j = 0; j < nz; ++j) {
            const double d = fu(i, j) - fu0(i, j);
            if (d != 0.0) cols[static_cast<std::size_t>(i * nz + j)].emplace_back(owner(i, j), d);
          }
        }
        for (int i = 0; i < g.nr; ++i) {
          for (int j = 0; j <= nz; ++j) {
            const double d = fw(i, j) - fw0(i, j);
            if (d != 0.0) {
              cols[static_cast<std::size_t>(nu + i * (nz + 1) + j)].emplace_back(owner(i, j), d);
            }
          }
        }
      }
    }
  }
  a.n = n;
  a.row_start.assign(static_cast<std::size_t>(n) + 1, 0);
  a.col.clear();
  a.val.clear();
  for (int r = 0; r < n; ++r) {
    auto& row = cols[static_cast<std::size_t>(r)];
    std::sort(row.begin(), row.end());
    for (const auto& [c, v] : row) {
      a.col.push_back(c);
      a.val.push_back(v);
    }
    a.row_start[r + 1] = static_cast<int>(a.col.size());
  }
}

void FlowSolver::step_implicit(double dt) {
  const Grid& g = p_.grid;
  const int nz = g.nz;
  const int nu = (g.nr + 1) * nz;
  const int n = nu + g.nr * (nz + 1);
  const double t1 = t_ + dt;
  // 시간 차수: BDF2 는 직전 단계 이력이 있을 때만 (03 §1.4)
  const bool bdf2 = opt_.time_order == 2 && hist_valid_;
  const double om = bdf2 ? dt / dt_prev_ : 0.0;
  const double c0 = bdf2 ? (1.0 + 2.0 * om) / (1.0 + om) : 1.0;
  const double c1 = bdf2 ? 1.0 + om : 1.0;
  const double c2 = bdf2 ? om * om / (1.0 + om) : 0.0;
  const double tau = dt / c0;  // 압력 증분 투영의 유효 간격

  // 명시적 부분 (시각 t): 대류 + 중력·원천
  Array2D fu, fw;
  conv_dt_ = dt;
  rhs(u_, w_, t_, fu, fw, true, false, true);
  conv_dt_ = 0.0;
  // BDF2: 명시 항을 t1 로 외삽. 이력용으로 이번 값은 따로 둔다.
  const Array2D fu_now = fu, fw_now = fw;
  if (bdf2) {
    for (std::size_t q = 0; q < fu.data().size(); ++q) fu.data()[q] = (1.0 + om) * fu.data()[q] - om * fu_prev_.data()[q];
    for (std::size_t q = 0; q < fw.data().size(); ++q) fw.data()[q] = (1.0 + om) * fw.data()[q] - om * fw_prev_.data()[q];
  }
  // 점성 부분 (시각 t1): V(u) = A u + c
  if (!visc_valid_ || (p_.material->time_dependent() && visc_time_ != t1)) {
    viscous_matrix(t1, visc_);
    visc_time_ = t1;
    visc_valid_ = true;
  }
  Array2D zu(g.nr + 1, nz), zw(g.nr, nz + 1), cu, cw;
  rhs(zu, zw, t1, cu, cw, false, true, false);
  const MaterialArrays& m = material_at(t1);

  // 경계 면의 값 (시각 t1)
  Array2D bu = u_, bw = w_;
  apply_boundary(bu, bw, t1);

  body::MomentumClasses mc;
  double ub_w = 0.0;
  if (p_.body) {
    mc = body::momentum_classes(g, *p_.body->shape, p_.body->z(t1), p_.body->interpolation_order);
    ub_w = p_.body->w(t1);
  }
  // 행 종류: 0 유체(운동량), 1 고정값, 2 경계 보간
  std::vector<char> kind(static_cast<std::size_t>(n), 0);
  std::vector<double> fixed(static_cast<std::size_t>(n), 0.0);
  std::vector<const body::FaceForcing*> forcing(static_cast<std::size_t>(n), nullptr);
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const int k = i * nz + j;
      if (!solves_u(i)) {
        kind[k] = 1;
        fixed[k] = bu(i, j);
      } else if (p_.body && mc.solid_u[static_cast<std::size_t>(k)]) {
        kind[k] = 1;
        fixed[k] = 0.0;
      }
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= nz; ++j) {
      const int k = nu + i * (nz + 1) + j;
      if (!solves_w(j)) {
        kind[k] = 1;
        fixed[k] = bw(i, j);
      } else if (p_.body && mc.solid_w[static_cast<std::size_t>(i) * (nz + 1) + j]) {
        kind[k] = 1;
        fixed[k] = ub_w;
      }
    }
  }
  for (const auto& f : mc.forcing) {
    const int k = f.comp == 0 ? f.i * nz + f.j : nu + f.i * (nz + 1) + f.j;
    if (kind[k] != 0) continue;  // 바깥 경계 면은 그대로
    kind[k] = 2;
    forcing[k] = &f;
  }

  // 운동량이 값을 정한 면(고체 면·경계 보간면)은 투영이 고치지 않고 알려진 플럭스로 둔다.
  // 그렇지 않으면 투영과 운동량이 같은 면을 번갈아 고쳐 정상 상태에 닿지 못한다.
  // 열린 면이 하나도 고칠 수 없는 셀은 압력 미지수에서 빠진다 (발산은 그 셀에 남는다).
  if (p_.body) {
    update_geometry(t1);
    const std::vector<char> was_active = pgeom_.active;
    pgeom_ = geom_;
    for (int i = iu_lo_; i < iu_hi_; ++i) {
      for (int j = 0; j < nz; ++j) {
        if (kind[static_cast<std::size_t>(i * nz + j)] != 0) pgeom_.open_u(i, j) = 0.0;
      }
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = jw_lo_; j < jw_hi_; ++j) {
        if (kind[static_cast<std::size_t>(nu + i * (nz + 1) + j)] != 0) pgeom_.open_w(i, j) = 0.0;
      }
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j < nz; ++j) {
        // 축 (축대칭) 과 대칭면 (평면, P4) 의 r = 0 면은 속도가 정해져 압력을 잇지 않는다.
        // 대칭면 면을 열린 면으로 세면 물체 옆에서 이웃과 끊긴 셀이 압력 미지수가 되어
        // 압축성 항만 남아 폭주한다 (2026-10-23, 쐐기 3-A 개발 실행).
        const bool west = !(i == 0 && p_.bc.r_min.kind == BoundaryKind::Axis) && pgeom_.open_u(i, j) > 0.0;
        const bool any = west || pgeom_.open_u(i + 1, j) > 0.0 || pgeom_.open_w(i, j) > 0.0 ||
                         pgeom_.open_w(i, j + 1) > 0.0;
        pgeom_.active[static_cast<std::size_t>(i) * nz + j] = any ? 1 : 0;
      }
    }
    if (!pgeom_.active.back()) {
      throw std::invalid_argument("the last cell (pressure reference) must be in the fluid");
    }
    // 물체가 움직여 새로 압력 미지수가 된 셀: 이전 압력이 없으므로(0) 그대로 두면 예측 단계의
    // 압력 기울기가 터진다. 이전에 유효했던 이웃 셀의 평균으로 채운다.
    if (was_active.size() == pgeom_.active.size()) {
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < nz; ++j) {
          const std::size_t k = static_cast<std::size_t>(i) * nz + j;
          if (!pgeom_.active[k] || was_active[k]) continue;
          double sum = 0.0;
          int count = 0;
          for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
              const int ii = i + di, jj = j + dj;
              if (ii < 0 || jj < 0 || ii >= g.nr || jj >= nz) continue;
              if (!was_active[static_cast<std::size_t>(ii) * nz + jj]) continue;
              sum += p_inc_(ii, jj);
              ++count;
            }
          }
          if (count > 0) p_inc_(i, j) = sum / count;
        }
      }
    }
  }

  // 행렬과 우변
  CsrMatrix a;
  a.n = n;
  a.row_start.assign(static_cast<std::size_t>(n) + 1, 0);
  std::vector<double> b(static_cast<std::size_t>(n), 0.0);
  std::vector<std::pair<int, double>> row;
  for (int k = 0; k < n; ++k) {
    row.clear();
    if (kind[k] == 1) {
      row.emplace_back(k, 1.0);
      b[k] = fixed[k];
    } else if (kind[k] == 2) {
      const auto& f = *forcing[k];
      const double ub = f.comp == 0 ? 0.0 : ub_w;
      row.emplace_back(k, 1.0);
      if (f.weight != 0.0) {
        const int ko = f.comp == 0 ? f.oi * nz + f.oj : nu + f.oi * (nz + 1) + f.oj;
        row.emplace_back(ko, -f.weight);
      }
      if (f.weight2 != 0.0) {
        const int ko = f.comp == 0 ? f.oi2 * nz + f.oj2 : nu + f.oi2 * (nz + 1) + f.oj2;
        row.emplace_back(ko, -f.weight2);
      }
      b[k] = (1.0 - f.weight - f.weight2) * ub;
    } else {
      // (1/Δt) u − A u = u^n/Δt + 명시적 항 + c − (1/ρ)∇p^n
      // BDF2: (c0/Δt) u − A u = (c1 u^n − c2 u^{n−1})/Δt + 외삽 명시 항 + c − (1/ρ)∇p^n
      row.emplace_back(k, bdf2 ? c0 / dt : 1.0 / dt);
      for (int q = visc_.row_start[k]; q < visc_.row_start[k + 1]; ++q) {
        row.emplace_back(visc_.col[q], -visc_.val[q]);
      }
      if (k < nu) {
        const int i = k / nz, j = k % nz;
        const double grad = (p_inc_(i, j) - p_inc_(i - 1, j)) / (m.rho_u(i, j) * g.dr);
        const double hist = bdf2 ? (c1 * u_(i, j) - c2 * u_prev_(i, j)) / dt : u_(i, j) / dt;
        b[k] = hist + fu(i, j) + cu(i, j) - grad;
      } else {
        const int i = (k - nu) / (nz + 1), j = (k - nu) % (nz + 1);
        const double grad = (p_inc_(i, j) - p_inc_(i, j - 1)) / (m.rho_w(i, j) * g.dz);
        const double hist = bdf2 ? (c1 * w_(i, j) - c2 * w_prev_(i, j)) / dt : w_(i, j) / dt;
        b[k] = hist + fw(i, j) + cw(i, j) - grad;
      }
    }
    std::sort(row.begin(), row.end());
    for (const auto& [c, v] : row) {  // 같은 열은 합친다
      if (static_cast<int>(a.col.size()) > a.row_start[k] && a.col.back() == c) {
        a.val.back() += v;
      } else {
        a.col.push_back(c);
        a.val.push_back(v);
      }
    }
    a.row_start[k + 1] = static_cast<int>(a.col.size());
  }
  Ilu0 ilu;
  ilu.factor(a);
  if (mom_x_.size() != static_cast<std::size_t>(n)) {
    mom_x_.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i <= g.nr; ++i) {
      for (int j = 0; j < nz; ++j) mom_x_[static_cast<std::size_t>(i * nz + j)] = u_(i, j);
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j <= nz; ++j) mom_x_[static_cast<std::size_t>(nu + i * (nz + 1) + j)] = w_(i, j);
    }
  }
  const SolveResult res = bicgstab(a, ilu, b, mom_x_, 1e-12, opt_.max_iterations);
  stats_.momentum_solves += 1;
  stats_.momentum_iterations += res.iterations;
  stats_.max_momentum_residual = std::max(stats_.max_momentum_residual, res.relative_residual);
  if (!res.converged) stats_.solver_failed = true;

  Array2D un = u_, wn = w_;
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < nz; ++j) un(i, j) = mom_x_[static_cast<std::size_t>(i * nz + j)];
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= nz; ++j) wn(i, j) = mom_x_[static_cast<std::size_t>(nu + i * (nz + 1) + j)];
  }
  apply_boundary(un, wn, t1);
  // 증분이므로 0 이 가장 가까운 초기 추정
  phi_[0].hi.assign(phi_[0].hi.size(), 0.0);
  phi_[0].lo.assign(phi_[0].lo.size(), 0.0);
  // 회전형 압력 보정용: 투영 전 셀 발산 (b = −유출)
  update_geometry(t1);
  std::vector<double> b_raw;
  face_flux_rhs(un, wn, ub_w, b_raw);
  // BDF2 압축성 이력 항: 우변 += ψ V c2 (p^n − p^{n−1}) / Δt (03 §1.4)
  std::vector<double> bdf_extra;
  if (p_.air_kappa > 0.0) {
    // ψ V / τ², ψ = (공기 체적분율)/(κ p^n) — 압력은 직전 값으로 선형화 (1 차는 τ = Δt)
    psi_vol_.assign(static_cast<std::size_t>(g.nr) * nz, 0.0);
    if (bdf2) bdf_extra.assign(psi_vol_.size(), 0.0);
    long floor_cells = 0;
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j < nz; ++j) {
        const double pa = p_inc_(i, j);
        const std::size_t kc = static_cast<std::size_t>(i) * nz + j;
        // 공기가 없는 셀 (물, 고체) 과 이번 단계 압력 미지수가 아닌 셀 (물체 옆 컷셀 중 모든 면을
        // 운동량이 정한 셀) 은 압축성 항이 없다. 그런 셀의 압력은 0 으로 남아 있을 수 있다.
        const bool no_pressure = p_.body && !pgeom_.active[static_cast<std::size_t>(i) * nz + j];
        if (m.air_fraction(i, j) <= 0.0 || no_pressure) {
          psi_vol_[kc] = 0.0;
          continue;
        }
        // VOF 잔여 공기 (물 셀의 1e−6 이하) 는 실제 공기가 아니다. 그런 셀에 물의 음압이 걸리면
        // 압축성 항을 빼고 넘어간다 (P4 에서 확인; 압력이 양수인 P3 계산은 그대로).
        if (!(pa > 0.0) && m.air_fraction(i, j) <= 1e-6 && opt_.air_pressure_floor <= 0.0) {
          psi_vol_[static_cast<std::size_t>(i) * nz + j] = 0.0;
          continue;
        }
        // 공기 압력 하한 (03 §1.4): 인장 구간에서 ψ 가 정의되도록 하한 압력을 쓴다
        double p_use = pa;
        if (opt_.air_pressure_floor > 0.0 && !(pa >= opt_.air_pressure_floor)) {
          p_use = opt_.air_pressure_floor;
          ++floor_cells;
        }
        if (!(p_use > 0.0)) {
          throw std::runtime_error("compressible air: non-positive absolute pressure at cell (" + std::to_string(i) +
                                   ", " + std::to_string(j) + "), p = " + std::to_string(pa) +
                                   ", air fraction " + std::to_string(m.air_fraction(i, j)));
        }
        const double psi_v = m.air_fraction(i, j) / (p_.air_kappa * p_use) * g.cell_volume(i);
        psi_vol_[kc] = psi_v / (tau * tau);
        // 직전 단계에 압력이 없던 셀 (물체 옆에서 새로 생김) 과 하한이 걸린 셀 (이번 또는 직전 압력이
        // 하한 아래) 은 이력 항을 쓰지 않는다 (그 셀만 1 차). 하한으로 키운 ψ 에 실제 압력 차를 곱하면
        // 폭주한다 (2026-10-23, a/Δ = 32 Δt × ¼ 개발 실행에서 6.3 ms 에 발산).
        const bool floored = bdf2 && opt_.air_pressure_floor > 0.0 &&
                             (!(pa >= opt_.air_pressure_floor) || !(p_prev_(i, j) >= opt_.air_pressure_floor));
        if (bdf2 && !floored && p_prev_(i, j) > 0.0) {
          bdf_extra[static_cast<std::size_t>(i) * nz + j] = psi_v * c2 * (pa - p_prev_(i, j)) / dt;
        }
      }
    }
    stats_.pressure_floor_cells += floor_cells;
    op_.valid = false;
    op_mom_.valid = false;
  }
  if (!p_.body) {
    project(un, wn, t1, t1, phi_[0], bdf_extra.empty() ? nullptr : &bdf_extra);
  } else {
    // 압력 기하 pgeom_ 는 운동량 조립 전에 만들었다
    if (!op_mom_.valid || p_.body->moving() || p_.material->time_dependent()) {
      build_operator(m, op_mom_, pgeom_);
      op_mom_.valid = true;
    }
    std::vector<double> div_b = b_raw;  // 원래 열린 비율로 잰 발산
    double s = 0.0;
    std::size_t count = 0;
    for (std::size_t c = 0; c < div_b.size(); ++c) {
      if (!pgeom_.active[c]) {
        div_b[c] = 0.0;
      } else {
        s += div_b[c];
        ++count;
      }
    }
    // 비압축이면 우변 평균을 빼서 풀이 가능하게 한다. 압축성 공기는 행렬이 정칙이라 빼지 않는다 (P4).
    if (p_.air_kappa <= 0.0) {
      for (std::size_t c = 0; c < div_b.size(); ++c) {
        if (pgeom_.active[c]) div_b[c] -= s / static_cast<double>(count);
      }
    }
    if (!bdf_extra.empty()) {
      for (std::size_t c = 0; c < div_b.size(); ++c) {
        if (pgeom_.active[c]) div_b[c] += bdf_extra[c];
      }
    }
    Potential& phi = phi_[0];
    record(refined_solve(op_mom_.matrix, op_mom_.preconditioner, div_b, phi.hi, phi.lo,
                         opt_.projection_tolerance, opt_.max_iterations),
           opt_.projection_tolerance);
    auto diff = [&](std::size_t k1, std::size_t k0) {
      return (phi.hi[k1] - phi.hi[k0]) + (phi.lo[k1] - phi.lo[k0]);
    };
    for (int i = iu_lo_; i < iu_hi_; ++i) {
      for (int j = 0; j < nz; ++j) {
        if (pgeom_.open_u(i, j) <= 0.0) continue;
        const std::size_t k1 = static_cast<std::size_t>(i) * nz + j;
        un(i, j) -= diff(k1, k1 - static_cast<std::size_t>(nz)) / (m.rho_u(i, j) * g.dr);
      }
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = jw_lo_; j < jw_hi_; ++j) {
        if (pgeom_.open_w(i, j) <= 0.0) continue;
        const std::size_t k1 = static_cast<std::size_t>(i) * nz + j;
        wn(i, j) -= diff(k1, k1 - 1) / (m.rho_w(i, j) * g.dz);
      }
    }
  }
  // 증분 압력, 회전형 (Timmermans 형): p += φ/τ − μ ∇·u*.  φ = τ · δp (투영 퍼텐셜, 1 차는 τ = Δt).
  // 둘째 항은 정상 상태에서 0 이므로 정상 해를 바꾸지 않으며, 큰 Δt 에서 압력 수렴을 빠르게 한다
  // (Stokes 에서 Uzawa 반복과 같은 꼴).
  const body::CutGeometry& pg = pressure_geometry();
  const Array2D p_n = p_inc_;  // BDF2 이력용 p^n
  double sum = 0.0, vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const std::size_t k = static_cast<std::size_t>(i) * nz + j;
      if (!pg.active[k]) continue;
      // 회전형 항은 ∇·u = 0 을 전제로 하므로 압축성 공기에서는 쓰지 않는다
      const double div = p_.air_kappa > 0.0 ? 0.0 : -b_raw[k] / g.cell_volume(i);
      p_inc_(i, j) += (phi_[0].hi[k] + phi_[0].lo[k]) / tau - m.mu_cell(i, j) * div;
      sum += p_inc_(i, j) * g.cell_volume(i);
      vol += g.cell_volume(i);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      // 압축성 공기에서는 절대 압력이므로 평균을 빼지 않는다
      const double shift = p_.air_kappa > 0.0 ? 0.0 : sum / vol;
      p_inc_(i, j) = pg.active[static_cast<std::size_t>(i) * nz + j] ? p_inc_(i, j) - shift : 0.0;
    }
  }
  // 물 압력 하한 (공동화 근사): 절대압력의 하한
  if (opt_.water_pressure_floor > 0.0 && p_.air_kappa > 0.0) {
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j < nz; ++j) {
        if (!pg.active[static_cast<std::size_t>(i) * nz + j]) continue;
        if (p_inc_(i, j) < opt_.water_pressure_floor) {
          p_inc_(i, j) = opt_.water_pressure_floor;
          ++stats_.water_floor_cells;
        }
      }
    }
  }

  double du = 0.0, umax = 0.0;
  for (std::size_t q = 0; q < un.data().size(); ++q) {
    du = std::max(du, std::abs(un.data()[q] - u_.data()[q]));
    umax = std::max(umax, std::abs(un.data()[q]));
  }
  for (std::size_t q = 0; q < wn.data().size(); ++q) {
    du = std::max(du, std::abs(wn.data()[q] - w_.data()[q]));
    umax = std::max(umax, std::abs(wn.data()[q]));
  }
  last_change_ = umax > 0.0 ? du / umax : 0.0;
  if (opt_.time_order == 2) {
    u_prev_ = u_;
    w_prev_ = w_;
    fu_prev_ = fu_now;
    fw_prev_ = fw_now;
    p_prev_ = p_n;
    dt_prev_ = dt;
    hist_valid_ = true;
  }
  u_ = un;
  w_ = wn;
  t_ = t1;
  stats_.steps += 1;
}

}  // namespace hcw::flow
