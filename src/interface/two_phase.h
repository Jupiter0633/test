#pragma once

// P3: 물–공기 2상 진행 (기하 VOF + 변수밀도 투영 해석기).
// 한 단계 (예측–수정):
//   1) 예측: u^n 의 면 플럭스로 f^n 을 Δt 수송 → f*.
//   2) 물성을 [t, t+Δt] 에서 f^n → f* 직선 보간으로 두고 유동 한 단계 (P1 의 3단 RK, 각 단계 투영) → u^{n+1}.
//   3) 수정: (u^n + u^{n+1})/2 의 면 플럭스로 f^n 을 다시 수송 → f^{n+1}. 이후 물성은 f^{n+1}.
// 두 투영 속도장의 평균도 이산 발산이 (해법 허용치 안에서) 0 이므로 물 체적이 보존된다.
// 분할 순서는 단계마다 번갈아 쓰되, 한 단계의 예측·수정은 같은 순서를 쓴다.

#include <functional>
#include <memory>

#include "flow/solver.h"
#include "interface/vof.h"
#include "interface/vof_material.h"

namespace hcw::interface {

class TwoPhaseSolver {
 public:
  // problem.material 은 쓰지 않고 material 로 바꾼다. f0: 초기 물 체적분율.
  TwoPhaseSolver(flow::Problem problem, std::shared_ptr<VofMaterial> material,
                 const flow::Array2D& f0, flow::SolverOptions options = {});

  void initialize(double t0, const std::function<double(double, double)>& u0,
                  const std::function<double(double, double)>& w0);
  void step(double dt);

  double time() const { return solver_.time(); }
  const flow::FlowSolver& flow() const { return solver_; }
  flow::FlowSolver& flow() { return solver_; }
  const Vof& vof() const { return vof_; }

  // 속도 → 면 체적 플럭스 (단위 시간, 2π 제외). 영역 가장자리 면은 0.
  void face_fluxes(const flow::Array2D& u, const flow::Array2D& w, flow::Array2D& qu,
                   flow::Array2D& qw) const;

 private:
  static flow::Problem with_material(flow::Problem p, std::shared_ptr<VofMaterial> m);

  std::shared_ptr<VofMaterial> mat_;
  flow::FlowSolver solver_;
  Vof vof_;
  long steps_ = 0;
};

}  // namespace hcw::interface
