#include "interface/two_phase.h"

namespace hcw::interface {

flow::Problem TwoPhaseSolver::with_material(flow::Problem p, std::shared_ptr<VofMaterial> m) {
  p.material = std::move(m);
  return p;
}

TwoPhaseSolver::TwoPhaseSolver(flow::Problem problem, std::shared_ptr<VofMaterial> material,
                               const flow::Array2D& f0, bool implicit,
                               flow::SolverOptions options)
    : mat_((material->set(fractions_from_vof(problem.grid, f0), 0.0,
                          fractions_from_vof(problem.grid, f0), 0.0),
            material)),
      solver_(with_material(std::move(problem), material), options),
      vof_(solver_.grid()),
      implicit_(implicit) {
  vof_.f() = f0;
}

void TwoPhaseSolver::initialize(double t0, const std::function<double(double, double)>& u0,
                                const std::function<double(double, double)>& w0, bool project) {
  const Fractions fr = fractions_from_vof(solver_.grid(), vof_.f());
  mat_->set(fr, t0, fr, t0);
  solver_.material_changed();
  solver_.initialize(t0, u0, w0, project);
}

void TwoPhaseSolver::face_fluxes(const flow::Array2D& u, const flow::Array2D& w, flow::Array2D& qu,
                                 flow::Array2D& qw) const {
  const flow::Grid& g = solver_.grid();
  qu = flow::Array2D(g.nr + 1, g.nz);
  qw = flow::Array2D(g.nr, g.nz + 1);
  for (int i = 1; i < g.nr; ++i) {
    const double area = g.metric(g.rf(i)) * g.dz;
    for (int j = 0; j < g.nz; ++j) qu(i, j) = u(i, j) * area;
  }
  for (int i = 0; i < g.nr; ++i) {
    const double area = g.metric(g.rc(i)) * g.dr;
    for (int j = 1; j < g.nz; ++j) qw(i, j) = w(i, j) * area;
  }
}

void TwoPhaseSolver::step(double dt) {
  const double t = solver_.time();
  const int order = steps_ % 2 == 0 ? 1 : 0;
  const flow::Array2D u0 = solver_.u(), w0 = solver_.w();
  const flow::Array2D f0 = vof_.f();
  flow::Array2D qu, qw;

  // 1) 예측
  face_fluxes(u0, w0, qu, qw);
  vof_.advect(qu, qw, dt, order);
  // 2) 유동 (밀도: f^n → f* 직선 보간)
  const Fractions fr0 = fractions_from_vof(solver_.grid(), f0);
  mat_->set(fr0, t, fractions_from_vof(solver_.grid(), vof_.f()), t + dt);
  solver_.material_changed();
  if (implicit_) {
    solver_.step_implicit(dt);
  } else {
    solver_.step(dt);
  }
  // 3) 수정
  flow::Array2D um = solver_.u(), wm = solver_.w();
  for (std::size_t k = 0; k < um.data().size(); ++k) um.data()[k] = 0.5 * (um.data()[k] + u0.data()[k]);
  for (std::size_t k = 0; k < wm.data().size(); ++k) wm.data()[k] = 0.5 * (wm.data()[k] + w0.data()[k]);
  face_fluxes(um, wm, qu, qw);
  vof_.f() = f0;
  vof_.advect(qu, qw, dt, order);
  const Fractions fr1 = fractions_from_vof(solver_.grid(), vof_.f());
  mat_->set(fr1, t + dt, fr1, t + dt);
  solver_.material_changed();
  ++steps_;
}

}  // namespace hcw::interface
