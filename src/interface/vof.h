#pragma once

// P3: 기하학적 VOF 수송 (물 체적분율 f, 셀마다 0..1).
// 방향 분할 (단계마다 r→z, z→r 순서를 번갈아). 면마다 상류 셀에서 체적이 정확히 |q|Δt 인 띠를 잡고,
// 그 띠 안의 PLIC 물 체적을 물 플럭스로 한다.
// 분할 중간 단계의 1차원 발산을 Weymouth–Yue 보정항 c·Σ(면 플럭스) (c = 단계 시작 f > 1/2) 로 메워,
// 전체 단계에서 이산 발산이 0 인 흐름이면 물 체적이 반올림 수준까지 보존된다.
// 근거: 02_design_document.md §3.2 (기하 VOF), 방법 문헌은 06_references.md R-32 (확인 수준 참조).

#include <vector>

#include "flow/grid.h"

namespace hcw::interface {

class Vof {
 public:
  explicit Vof(const flow::Grid& g);

  const flow::Grid& grid() const { return g_; }
  flow::Array2D& f() { return f_; }
  const flow::Array2D& f() const { return f_; }

  // 면 체적 플럭스 q (단위 시간, 2π 제외, 축대칭은 r 가중 면적 기준)로 Δt 만큼 수송.
  // qu: r 면 (nr+1, nz), qw: z 면 (nr, nz+1). 영역 가장자리 면의 플럭스는 0 이어야 한다.
  // 면마다 |q|Δt / 면적 이 그 방향 셀 폭의 절반을 넘으면 std::invalid_argument.
  void advect(const flow::Array2D& qu, const flow::Array2D& qw, double dt);

  double water_volume() const;     // Σ f · 셀 체적
  double bound_violation() const;  // max(−f, f − 1, 0) 전체 최댓값

 private:
  void sweep_r(const flow::Array2D& qu, double dt, const std::vector<char>& c);
  void sweep_z(const flow::Array2D& qw, double dt, const std::vector<char>& c);
  // 섞인 셀(0 < f < 1)의 PLIC 법선·α 를 현재 f 로 다시 구한다 (분할 단계마다)
  void reconstruct();
  // 상류 셀 (i, j) 의 PLIC 물 영역 중 직사각형 띠 안의 체적 (r 가중)
  double strip_water(int i, int j, double r0, double r1, double z0, double z1) const;

  flow::Grid g_;
  flow::Array2D f_;
  flow::Array2D nr_, nz_, alpha_;  // 재구성 결과 (섞인 셀만 의미 있음)
  long steps_ = 0;
};

}  // namespace hcw::interface
