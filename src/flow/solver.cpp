#include "flow/solver.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

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
  for (auto& phi : phi_) phi.assign(static_cast<std::size_t>(g.nr) * g.nz, 0.0);
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

FlowSolver::Operator& FlowSolver::operator_at(double t) {
  const bool fresh = op_.valid && (!p_.material->time_dependent() || op_.time == t);
  if (!fresh) {
    build_operator(material_at(t), op_);
    op_.time = t;
    op_.valid = true;
  }
  return op_;
}

void FlowSolver::build_operator(const MaterialArrays& m, Operator& op) const {
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

  for (int i = iu_lo_; i < iu_hi_; ++i) {
    const double area = g.metric(g.rf(i)) * g.dz;
    for (int j = 0; j < g.nz; ++j) {
      couple(cell(i - 1, j), cell(i, j), area / (m.rho_u(i, j) * g.dr));
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    const double area = g.metric(g.rc(i)) * g.dr;
    for (int j = jw_lo_; j < jw_hi_; ++j) {
      couple(cell(i, j - 1), cell(i, j), area / (m.rho_w(i, j) * g.dz));
    }
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

void FlowSolver::face_flux_rhs(const Array2D& fu, const Array2D& fw,
                               std::vector<double>& b) const {
  const Grid& g = p_.grid;
  b.assign(static_cast<std::size_t>(g.nr) * g.nz, 0.0);
  for (int i = 0; i < g.nr; ++i) {
    const double aw = g.metric(g.rf(i)) * g.dz;
    const double ae = g.metric(g.rf(i + 1)) * g.dz;
    const double az = g.metric(g.rc(i)) * g.dr;
    for (int j = 0; j < g.nz; ++j) {
      const double out = ae * fu(i + 1, j) - aw * fu(i, j) + az * (fw(i, j + 1) - fw(i, j));
      b[static_cast<std::size_t>(i) * g.nz + j] = -out;
    }
  }
}

void FlowSolver::record(const SolveResult& r, double tolerance) {
  stats_.poisson_solves += 1;
  stats_.poisson_iterations += r.iterations;
  stats_.max_poisson_iterations = std::max(stats_.max_poisson_iterations, r.iterations);
  stats_.max_relative_residual = std::max(stats_.max_relative_residual, r.relative_residual);
  if (!r.converged || r.relative_residual > tolerance) stats_.solver_failed = true;
}

void FlowSolver::project(Array2D& u, Array2D& w, double t_rho, std::vector<double>& phi) {
  const Grid& g = p_.grid;
  Operator& op = operator_at(t_rho);
  const MaterialArrays& m = material_at(t_rho);

  std::vector<double> b;
  face_flux_rhs(u, w, b);
  remove_mean(b);
  record(pcg(op.matrix, op.preconditioner, b, phi, opt_.projection_tolerance,
             opt_.max_iterations),
         opt_.projection_tolerance);

  auto at = [&](int i, int j) {
    return phi[static_cast<std::size_t>(wrap_r(i)) * g.nz + wrap_z(j)];
  };
  for (int i = iu_lo_; i < iu_hi_; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      u(i, j) -= (at(i, j) - at(i - 1, j)) / (m.rho_u(i, j) * g.dr);
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = jw_lo_; j < jw_hi_; ++j) {
      w(i, j) -= (at(i, j) - at(i, j - 1)) / (m.rho_w(i, j) * g.dz);
    }
  }
  if (periodic_r_) {
    for (int j = 0; j < g.nz; ++j) u(g.nr, j) = u(0, j);
  }
  if (periodic_z_) {
    for (int i = 0; i < g.nr; ++i) w(i, g.nz) = w(i, 0);
  }
}

void FlowSolver::rhs(const Array2D& u, const Array2D& w, double t, Array2D& fu, Array2D& fw) {
  const Grid& g = p_.grid;
  const MaterialArrays& m = material_at(t);
  const bool axisym = g.geometry == Geometry::Axisymmetric;
  const bool has_source = static_cast<bool>(p_.source);
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

      double f = -(conv_r + conv_z) + (visc_r + visc_z + hoop) / m.rho_u(i, j) + p_.gravity[0];
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

      double f = -(conv_r + conv_z) + (visc_r + visc_z) / m.rho_w(i, j) + p_.gravity[1];
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
  project(u_, w_, t0, phi_[4]);
}

void FlowSolver::step(double dt) {
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
  project(u_, w_, t, phi_[0]);
  // 2단: → t+dt/2 (밀도 t+dt)
  rhs(u_, w_, t + dt, fu, fw);
  combine(u_, w_, 0.75, 0.25);
  apply_boundary(u_, w_, t + 0.5 * dt);
  project(u_, w_, t + dt, phi_[1]);
  // 3단: → t+dt (밀도 t+dt/2)
  rhs(u_, w_, t + 0.5 * dt, fu, fw);
  combine(u_, w_, 1.0 / 3.0, 2.0 / 3.0);
  apply_boundary(u_, w_, t + dt);
  project(u_, w_, t + 0.5 * dt, phi_[2]);

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

  std::vector<double> b;
  face_flux_rhs(fu, fw, b);
  remove_mean(b);
  Operator& op = operator_at(t_);
  std::vector<double>& x = phi_[3];
  record(pcg(op.matrix, op.preconditioner, b, x, opt_.pressure_tolerance, opt_.max_iterations),
         opt_.pressure_tolerance);

  p = Array2D(g.nr, g.nz);
  double sum = 0.0, vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      p(i, j) = x[static_cast<std::size_t>(i) * g.nz + j];
      sum += p(i, j) * g.cell_volume(i);
      vol += g.cell_volume(i);
    }
  }
  const double mean = sum / vol;
  for (double& v : p.data()) v -= mean;
}

double FlowSolver::max_divergence() const {
  const Grid& g = p_.grid;
  double worst = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    const double aw = g.metric(g.rf(i)) * g.dz;
    const double ae = g.metric(g.rf(i + 1)) * g.dz;
    const double az = g.metric(g.rc(i)) * g.dr;
    for (int j = 0; j < g.nz; ++j) {
      const double out = ae * u_(i + 1, j) - aw * u_(i, j) + az * (w_(i, j + 1) - w_(i, j));
      worst = std::max(worst, std::abs(out / g.cell_volume(i)));
    }
  }
  return worst;
}

}  // namespace hcw::flow
