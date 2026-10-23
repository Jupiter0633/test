#include "structure/time_integration.h"

#include <algorithm>
#include <stdexcept>

namespace hcw::structure {

HhtIntegrator::HhtIntegrator(const Model& model, double dt, double alpha, Rayleigh damping)
    : model_(model), dt_(dt), alpha_(alpha), damp_(damping) {
  if (alpha > 0.0 || alpha < -1.0 / 3.0) throw std::invalid_argument("HHT alpha must be in [-1/3, 0]");
  gamma_ = (1.0 - 2.0 * alpha) / 2.0;
  beta_ = (1.0 - alpha) * (1.0 - alpha) / 4.0;
  // M + (1+α)γΔt C + (1+α)βΔt² K, C = a₀M + a₁K
  const double s = 1.0 + alpha;
  cM_ = 1.0 + s * gamma_ * dt * damp_.a0;
  cK_ = s * gamma_ * dt * damp_.a1 + s * beta_ * dt * dt;
  A_ = model.assemble(cK_, cM_);
  A_.factor();
  Mff_ = model.assemble(0.0, 1.0);
  Mff_.factor();
}

void HhtIntegrator::apply_C(const Vec& x, Vec& y) const { model_.apply(damp_.a1, damp_.a0, x, y); }

void HhtIntegrator::initialize(const Vec& u0, const Vec& v0, const Vec& F0, const Vec& a_prescribed) {
  u_ = u0;
  v_ = v0;
  F_ = F0;
  a_.assign(model_.num_dofs(), 0.0);
  for (int d : model_.fixed_dofs()) a_[d] = a_prescribed[d];
  Vec Ku, Cv, Ma;
  model_.apply_K(u_, Ku);
  apply_C(v_, Cv);
  model_.apply_M(a_, Ma);  // 자유 자리는 0 이라 M a_p 만 남음
  Vec r(model_.num_dofs());
  for (std::size_t i = 0; i < r.size(); ++i) r[i] = F0[i] - Ku[i] - Cv[i] - Ma[i];
  double res = 0.0;
  const Vec af = solve_checked(model_, Mff_, 0.0, 1.0, model_.restrict_free(r), &res);
  max_res_ = std::max(max_res_, res);
  model_.scatter_free(af, a_);
  t_ = 0.0;
}

void HhtIntegrator::step(const Vec& F_next, const Vec& up, const Vec& vp, const Vec& ap) {
  const std::size_t n = u_.size();
  const double dt = dt_, s = 1.0 + alpha_;
  // 예측값 (a_{n+1} = 0 으로 둔 u, v), 지정 자리는 주어진 값
  Vec us(n), vs(n), as(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    us[i] = u_[i] + dt * v_[i] + dt * dt * (0.5 - beta_) * a_[i];
    vs[i] = v_[i] + dt * (1.0 - gamma_) * a_[i];
  }
  for (int d : model_.fixed_dofs()) {
    us[d] = up[d];
    vs[d] = vp[d];
    as[d] = ap[d];
  }
  Vec Ma, Cvs, Kus, Cvn, Kun;
  model_.apply_M(as, Ma);
  apply_C(vs, Cvs);
  model_.apply_K(us, Kus);
  apply_C(v_, Cvn);
  model_.apply_K(u_, Kun);
  Vec r(n);
  for (std::size_t i = 0; i < n; ++i) {
    r[i] = s * F_next[i] - alpha_ * F_[i] - Ma[i] - s * (Cvs[i] + Kus[i]) + alpha_ * (Cvn[i] + Kun[i]);
  }
  double res = 0.0;
  const Vec af = solve_checked(model_, A_, cK_, cM_, model_.restrict_free(r), &res);
  max_res_ = std::max(max_res_, res);
  Vec afull(n, 0.0);
  model_.scatter_free(af, afull);
  for (std::size_t i = 0; i < n; ++i) {
    if (model_.is_fixed(static_cast<int>(i))) {
      u_[i] = us[i];
      v_[i] = vs[i];
      a_[i] = as[i];
    } else {
      u_[i] = us[i] + beta_ * dt * dt * afull[i];
      v_[i] = vs[i] + gamma_ * dt * afull[i];
      a_[i] = afull[i];
    }
  }
  F_ = F_next;
  t_ += dt;
}

double HhtIntegrator::energy() const {
  Vec Ku, Mv;
  model_.apply_K(u_, Ku);
  model_.apply_M(v_, Mv);
  return 0.5 * model_.dot(v_, Mv) + 0.5 * model_.dot(u_, Ku);
}

}  // namespace hcw::structure
