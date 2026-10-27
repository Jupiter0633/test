#include "interface/two_phase.h"

#include <algorithm>
#include <stdexcept>

#include "body/cut_cell.h"

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
  body_ = solver_.body();
  if (body_) {
    if (!implicit_) throw std::invalid_argument("two-phase flow with a body needs the semi-implicit step");
    const Fractions fr = fractions_at(vof_.f(), 0.0);
    mat_->set(fr, 0.0, fr, 0.0);
  }
}

flow::Array2D TwoPhaseSolver::fluid_fraction(double t) const {
  const flow::Grid& g = solver_.grid();
  if (!body_) return flow::Array2D(g.nr, g.nz, 1.0);
  return body::fluid_fraction(g, *body_->shape, body_->z(t));
}

Fractions TwoPhaseSolver::fractions_at(const flow::Array2D& f, double t) const {
  if (!body_) return fractions_from_vof(solver_.grid(), f);
  return fractions_with_solid(solver_.grid(), f, fluid_fraction(t));
}

void TwoPhaseSolver::initialize(double t0, const std::function<double(double, double)>& u0,
                                const std::function<double(double, double)>& w0, bool project) {
  const Fractions fr = fractions_at(vof_.f(), t0);
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
    const double mf = g.metric(g.rf(i));
    for (int j = 0; j < g.nz; ++j) qu(i, j) = u(i, j) * (mf * g.dzc(j));
  }
  for (int i = 0; i < g.nr; ++i) {
    const double area = g.metric(g.rc(i)) * g.drc(i);
    for (int j = 1; j < g.nz; ++j) qw(i, j) = w(i, j) * area;
  }
}

void TwoPhaseSolver::face_fluxes_with_body(const flow::Array2D& u, const flow::Array2D& w, double t,
                                           flow::Array2D& qu, flow::Array2D& qw) const {
  const flow::Grid& g = solver_.grid();
  const body::CutGeometry gm = body::cut_geometry(g, *body_->shape, body_->z(t));
  const double wb = body_->w(t);
  qu = flow::Array2D(g.nr + 1, g.nz);
  qw = flow::Array2D(g.nr, g.nz + 1);
  for (int i = 1; i < g.nr; ++i) {
    const double mf = g.metric(g.rf(i));
    for (int j = 0; j < g.nz; ++j) qu(i, j) = gm.open_u(i, j) * u(i, j) * (mf * g.dzc(j));  // 물체 r 속도 0
  }
  for (int i = 0; i < g.nr; ++i) {
    const double area = g.metric(g.rc(i)) * g.drc(i);
    for (int j = 1; j < g.nz; ++j) {
      const double open = gm.open_w(i, j);
      qw(i, j) = (open * w(i, j) + (1.0 - open) * wb) * area;
    }
  }
}

void TwoPhaseSolver::step(double dt) {
  if (body_) {
    step_with_body(dt);
    return;
  }
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

void TwoPhaseSolver::step_with_body(double dt) {
  const flow::Grid& g = solver_.grid();
  const double t = solver_.time();
  const int order = steps_ % 2 == 0 ? 1 : 0;
  const flow::Array2D f0 = vof_.f();
  flow::Array2D qu0, qw0, qu1, qw1;

  // 1) 예측: 시각 t 의 합친 플럭스
  face_fluxes_with_body(solver_.u(), solver_.w(), t, qu0, qw0);
  vof_.advect(qu0, qw0, dt, order);
  const flow::Array2D phi1 = fluid_fraction(t + dt);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) vof_.f()(i, j) = std::clamp(vof_.f()(i, j), 0.0, phi1(i, j));
  }
  // 2) 유동 (밀도: f^n → f* 직선 보간, 각 시각의 고체 부분 제외)
  mat_->set(fractions_at(f0, t), t, fractions_with_solid(g, vof_.f(), phi1), t + dt);
  solver_.material_changed();
  solver_.step_implicit(dt);
  // 3) 수정: 두 시각 플럭스의 평균 (각각 발산 0)
  face_fluxes_with_body(solver_.u(), solver_.w(), t + dt, qu1, qw1);
  for (std::size_t k = 0; k < qu1.data().size(); ++k) qu1.data()[k] = 0.5 * (qu1.data()[k] + qu0.data()[k]);
  for (std::size_t k = 0; k < qw1.data().size(); ++k) qw1.data()[k] = 0.5 * (qw1.data()[k] + qw0.data()[k]);
  vof_.f() = f0;
  vof_.advect(qu1, qw1, dt, order);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      const double v = vof_.f()(i, j), c = std::clamp(v, 0.0, phi1(i, j));
      clipped_ += (v - c) * g.cell_volume(i, j);
      vof_.f()(i, j) = c;
    }
  }
  const Fractions fr1 = fractions_with_solid(g, vof_.f(), phi1);
  mat_->set(fr1, t + dt, fr1, t + dt);
  solver_.material_changed();
  ++steps_;
}

}  // namespace hcw::interface
