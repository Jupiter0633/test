#pragma once

// 시간 적분: HHT-α (α = 0 이면 Newmark 평균가속도), Rayleigh 감쇠 C = a₀M + a₁K.
// 정의 (03 §1.5 4-H):
//   M a_{n+1} + (1+α)(C v_{n+1} + K u_{n+1}) − α (C v_n + K u_n) = (1+α) F_{n+1} − α F_n
//   u_{n+1} = u_n + Δt v_n + Δt² [(½ − β) a_n + β a_{n+1}]
//   v_{n+1} = v_n + Δt [(1 − γ) a_n + γ a_{n+1}],  γ = (1 − 2α)/2, β = (1 − α)²/4
// 지정 자유도의 값 (u, v, a) 은 매 단계 호출한 쪽이 준다 (시간에 따라 변하는 경계 변위 가능).

#include "structure/assembly.h"

namespace hcw::structure {

struct Rayleigh {
  double a0 = 0.0;  // 질량 비례
  double a1 = 0.0;  // 강성 비례
};

class HhtIntegrator {
 public:
  HhtIntegrator(const Model& model, double dt, double alpha, Rayleigh damping = {});

  // 처음 상태. a₀ 은 운동방정식에서 구한다: M_ff a_f = F_f − [K u₀ + C v₀ + M a_p]_f.
  // u0, v0, a_prescribed 는 전역 벡터 (a_prescribed 는 지정 자리만 씀).
  void initialize(const Vec& u0, const Vec& v0, const Vec& F0, const Vec& a_prescribed);

  // 한 단계. F_next = F(t_{n+1}), up/vp/ap = 지정 자유도의 t_{n+1} 값 (전역 벡터, 지정 자리만 씀).
  void step(const Vec& F_next, const Vec& up, const Vec& vp, const Vec& ap);

  const Vec& u() const { return u_; }
  const Vec& v() const { return v_; }
  const Vec& a() const { return a_; }
  double time() const { return t_; }
  double dt() const { return dt_; }
  // 이산 에너지 ½ vᵀMv + ½ uᵀKu
  double energy() const;
  // 지금까지 선형 풀이 확인값의 최댓값
  const SolveCheck& check() const { return check_; }

 private:
  // y = C x
  void apply_C(const Vec& x, Vec& y) const;

  const Model& model_;
  double dt_, alpha_, beta_, gamma_;
  Rayleigh damp_;
  double cK_, cM_;  // 유효 행렬 = cK K + cM M
  SkylineMatrix A_;
  SkylineMatrix Mff_;
  Vec u_, v_, a_, F_;
  double t_ = 0.0;
  SolveCheck check_;
};

}  // namespace hcw::structure
