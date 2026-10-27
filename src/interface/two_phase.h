#pragma once

// P3: 물–공기 2상 진행 (기하 VOF + 변수밀도 투영 해석기).
// 한 단계 (예측–수정):
//   1) 예측: u^n 의 면 플럭스로 f^n 을 Δt 수송 → f*.
//   2) 물성을 [t, t+Δt] 에서 f^n → f* 직선 보간으로 두고 유동 한 단계 (P1 의 3단 RK, 각 단계 투영) → u^{n+1}.
//   3) 수정: (u^n + u^{n+1})/2 의 면 플럭스로 f^n 을 다시 수송 → f^{n+1}. 이후 물성은 f^{n+1}.
// 두 투영 속도장의 평균도 이산 발산이 (해법 허용치 안에서) 0 이므로 물 체적이 보존된다.
// 분할 순서는 단계마다 번갈아 쓰되, 한 단계의 예측·수정은 같은 순서를 쓴다.
//
// 물체가 있을 때 (P4, problem.body, implicit 필수):
//   - VOF 는 물체 속까지 격자 전체에서 푼다. 고체 부분은 "물이 아닌 상" 으로 물체와 함께 움직인다.
//   - 면 플럭스 = 열린 부분 × 유체 속도 + 막힌 부분 × 물체 속도 (r 성분 0, z 성분 W).
//     컷셀에서 유체 쪽 발산 (= −물체 경계 유출) 과 고체 쪽 발산 (= +물체 경계 유출) 이 정확히 상쇄되어
//     이 플럭스는 모든 셀에서 발산이 0 이다 → 물 체적 보존. 수정 단계는 두 시각 플럭스의 평균.
//   - 수송 뒤 f 를 [0, φ] 로 자른다 (φ: 유체 체적 비율). 잘린 물 체적은 clipped_water() 에 누적.
//   - 물성은 fractions_with_solid (유체 부분 안의 물 비율로 밀도, 공기 체적 φ − f 로 압축성).

#include <functional>
#include <memory>

#include "flow/solver.h"
#include "interface/vof.h"
#include "interface/vof_material.h"

namespace hcw::interface {

class TwoPhaseSolver {
 public:
  // problem.material 은 쓰지 않고 material 로 바꾼다. f0: 초기 물 체적분율.
  // implicit = true 이면 유동 단계에 step_implicit (점성 암시, 증분 압력; 압축성 공기에 필요),
  // false 이면 P1 의 3단 RK (step).
  TwoPhaseSolver(flow::Problem problem, std::shared_ptr<VofMaterial> material,
                 const flow::Array2D& f0, bool implicit = false, flow::SolverOptions options = {});

  // project = false 이면 초기 속도를 투영하지 않는다 (압축성 공기의 초기 압축 속도 등)
  void initialize(double t0, const std::function<double(double, double)>& u0,
                  const std::function<double(double, double)>& w0, bool project = true);
  void step(double dt);

  // 물체가 있을 때 f 를 [0, φ] 로 자르며 버린 물 체적의 누적 (부호 있음: + 는 버림)
  double clipped_water() const { return clipped_; }
  // 시각 t 의 유체 체적 비율 (물체가 없으면 모두 1)
  flow::Array2D fluid_fraction(double t) const;

  double time() const { return solver_.time(); }
  const flow::FlowSolver& flow() const { return solver_; }
  flow::FlowSolver& flow() { return solver_; }
  const Vof& vof() const { return vof_; }

  // 속도 → 면 체적 플럭스 (단위 시간, 2π 제외). 영역 가장자리 면은 0.
  void face_fluxes(const flow::Array2D& u, const flow::Array2D& w, flow::Array2D& qu,
                   flow::Array2D& qw) const;
  // 물체가 있을 때: 시각 t 의 컷셀 열린 비율로 유체·물체 플럭스를 합친다 (위 설명)
  void face_fluxes_with_body(const flow::Array2D& u, const flow::Array2D& w, double t,
                             flow::Array2D& qu, flow::Array2D& qw) const;

 private:
  static flow::Problem with_material(flow::Problem p, std::shared_ptr<VofMaterial> m);
  Fractions fractions_at(const flow::Array2D& f, double t) const;
  void step_with_body(double dt);

  std::shared_ptr<VofMaterial> mat_;
  flow::FlowSolver solver_;
  Vof vof_;
  bool implicit_ = false;
  long steps_ = 0;
  std::shared_ptr<const flow::Body> body_;
  double clipped_ = 0.0;
};

}  // namespace hcw::interface
