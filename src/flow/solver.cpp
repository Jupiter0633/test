#include "flow/solver.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

#include "flow/linear.h"

namespace hcw::flow {
namespace {

std::array<double, 2> wall_velocity(const Boundary& b, double r, double z, double t) {
  if (!b.velocity) return {0.0, 0.0};
  return b.velocity(r, z, t);
}

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
  if (bc.r_min.kind == BoundaryKind::Axis && !(axisym && g.r0 == 0.0)) {
    throw std::invalid_argument("axis boundary requires an axisymmetric grid starting at r = 0");
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
    if (!axisym || g.r0 != 0.0) {
      throw std::invalid_argument("bodies are supported on axisymmetric grids from the axis only");
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
  std::vector<std::map<int, double>> rows(static_cast<std::size_t>(n));
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
  rows[n - 1][n - 1] *= 2.0;

  CsrMatrix& a = op.matrix;
  a.n = n;
  a.row_start.assign(static_cast<std::size_t>(n) + 1, 0);
  a.col.clear();
  a.val.clear();
  for (int r = 0; r < n; ++r) {
    for (const auto& [c, v] : rows[r]) {
      a.col.push_back(c);
      a.val.push_back(v);
    }
    a.row_start[r + 1] = static_cast<int>(a.col.size());
  }
  op.preconditioner.factor(a);
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

void FlowSolver::project(Array2D& u, Array2D& w, double t_rho, double t_bc, Potential& phi) {
  const Grid& g = p_.grid;
  Operator& op = operator_at(t_rho, t_bc);
  const MaterialArrays& m = material_at(t_rho);
  const double body_w = p_.body ? p_.body->w(t_bc) : 0.0;

  std::vector<double> b;
  face_flux_rhs(u, w, body_w, b);
  remove_active_mean(b);
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

  // 꼭짓점 (r 면 fi, z 면 fj) 의 전단응력 τ_rz = μ (∂u/∂z + ∂w/∂r)
  auto tau_rz = [&](int fi, int fj) {
    return m.mu_corner(fi, fj) * ((U(u, fi, fj, t) - U(u, fi, fj - 1, t)) / g.dz +
                                  (W(w, fi, fj, t) - W(w, fi - 1, fj, t)) / g.dr);
  };

  // r 방향 운동량 (r 면)
  for (int i = iu_lo_; i < iu_hi_; ++i) {
    const double m_f = g.metric(g.rf(i));
    const double m_l = g.metric(g.rc(i - 1));
    const double m_r = g.metric(g.rc(i));
    for (int j = 0; j < g.nz; ++j) {
      const double u_l = 0.5 * (U(u, i - 1, j, t) + U(u, i, j, t));
      const double u_r = 0.5 * (U(u, i, j, t) + U(u, i + 1, j, t));
      const double conv_r = (m_r * u_r * u_r - m_l * u_l * u_l) / (m_f * g.dr);
      auto flux_z = [&](int fj) {
        const double uc = 0.5 * (U(u, i, fj - 1, t) + U(u, i, fj, t));
        const double wc = 0.5 * (W(w, i - 1, fj, t) + W(w, i, fj, t));
        return uc * wc;
      };
      const double conv_z = (flux_z(j + 1) - flux_z(j)) / g.dz;

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

  // z 방향 운동량 (z 면)
  for (int i = 0; i < g.nr; ++i) {
    const double m_c = g.metric(g.rc(i));
    const double m_w = g.metric(g.rf(i));
    const double m_e = g.metric(g.rf(i + 1));
    for (int j = jw_lo_; j < jw_hi_; ++j) {
      auto flux_r = [&](int fi) {
        const double uc = 0.5 * (U(u, fi, j - 1, t) + U(u, fi, j, t));
        const double wc = 0.5 * (W(w, fi - 1, j, t) + W(w, fi, j, t));
        return uc * wc;
      };
      const double conv_r = (m_e * flux_r(i + 1) - m_w * flux_r(i)) / (m_c * g.dr);
      const double w_s = 0.5 * (W(w, i, j - 1, t) + W(w, i, j, t));
      const double w_n = 0.5 * (W(w, i, j, t) + W(w, i, j + 1, t));
      const double conv_z = (w_n * w_n - w_s * w_s) / g.dz;

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

  if (periodic_r_) {
    for (int j = 0; j < g.nz; ++j) fu(g.nr, j) = fu(0, j);
  }
  if (periodic_z_) {
    for (int i = 0; i < g.nr; ++i) fw(i, g.nz) = fw(i, 0);
  }
}

void FlowSolver::initialize(double t0, const std::function<double(double, double)>& u0,
                            const std::function<double(double, double)>& w0) {
  const Grid& g = p_.grid;
  t_ = t0;
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) u_(i, j) = u0(g.rf(i), g.zc(j));
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) w_(i, j) = w0(g.rc(i), g.zf(j));
  }
  apply_boundary(u_, w_, t0);
  project(u_, w_, t0, t0, phi_[4]);
}

void FlowSolver::step(double dt) {
  if (p_.body) {
    throw std::logic_error(
        "time stepping with a body is not implemented yet (P2: near-body momentum pending)");
  }
  const Grid& g = p_.grid;
  const double t = t_;
  const Array2D u0 = u_, w0 = w_;
  Array2D fu, fw;

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

  // 명시적 부분 (시각 t): 대류 + 중력·원천
  Array2D fu, fw;
  rhs(u_, w_, t_, fu, fw, true, false, true);
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
        const bool west = !(g.geometry == Geometry::Axisymmetric && i == 0) && pgeom_.open_u(i, j) > 0.0;
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
      row.emplace_back(k, 1.0 / dt);
      for (int q = visc_.row_start[k]; q < visc_.row_start[k + 1]; ++q) {
        row.emplace_back(visc_.col[q], -visc_.val[q]);
      }
      if (k < nu) {
        const int i = k / nz, j = k % nz;
        const double grad = (p_inc_(i, j) - p_inc_(i - 1, j)) / (m.rho_u(i, j) * g.dr);
        b[k] = u_(i, j) / dt + fu(i, j) + cu(i, j) - grad;
      } else {
        const int i = (k - nu) / (nz + 1), j = (k - nu) % (nz + 1);
        const double grad = (p_inc_(i, j) - p_inc_(i, j - 1)) / (m.rho_w(i, j) * g.dz);
        b[k] = w_(i, j) / dt + fw(i, j) + cw(i, j) - grad;
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
  if (!p_.body) {
    project(un, wn, t1, t1, phi_[0]);
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
    for (std::size_t c = 0; c < div_b.size(); ++c) {
      if (pgeom_.active[c]) div_b[c] -= s / static_cast<double>(count);
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
  // 증분 압력, 회전형 (Timmermans 형): p += φ/Δt − μ ∇·u*.  φ = Δt · δp (투영 퍼텐셜).
  // 둘째 항은 정상 상태에서 0 이므로 정상 해를 바꾸지 않으며, 큰 Δt 에서 압력 수렴을 빠르게 한다
  // (Stokes 에서 Uzawa 반복과 같은 꼴).
  const body::CutGeometry& pg = pressure_geometry();
  double sum = 0.0, vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const std::size_t k = static_cast<std::size_t>(i) * nz + j;
      if (!pg.active[k]) continue;
      const double div = -b_raw[k] / g.cell_volume(i);
      p_inc_(i, j) += (phi_[0].hi[k] + phi_[0].lo[k]) / dt - m.mu_cell(i, j) * div;
      sum += p_inc_(i, j) * g.cell_volume(i);
      vol += g.cell_volume(i);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      p_inc_(i, j) = pg.active[static_cast<std::size_t>(i) * nz + j] ? p_inc_(i, j) - sum / vol : 0.0;
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
  u_ = un;
  w_ = wn;
  t_ = t1;
  stats_.steps += 1;
}

}  // namespace hcw::flow
