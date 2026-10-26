#pragma once

// 2차원(평면/축대칭) 비압축 단상·변수밀도 Navier–Stokes 해석기 (P1).
// 엇갈림 격자 유한체적, 2차 중심차분, 3단 SSP-RK 시간적분, 각 단계 투영.
// 근거: docs/design/02_design_document.md §3.1, §3.5, §3.6

#include <array>
#include <functional>
#include <memory>

#include "body/cut_cell.h"
#include "body/shape.h"
#include "flow/grid.h"
#include "flow/material.h"
#include "flow/poisson.h"

namespace hcw::flow {

// 지정 운동 물체 (P2). 축을 따라 z 방향 병진만 한다. 함수가 비어 있으면 0 (정지, 중심 z = 0).
struct BodyMotion {
  std::function<double(double t)> z;  // 중심 높이
  std::function<double(double t)> w;  // 속도
  std::function<double(double t)> a;  // 가속도
};

struct Body {
  std::shared_ptr<const body::Shape> shape;
  BodyMotion motion;
  // step_implicit 의 경계 보간 차수 (body::momentum_classes 참조). 2: 정상 계산(벽 견인력 2차),
  // 1: 움직이는 물체·작은 점성의 비정상 계산 (2차는 불안정함을 확인)
  int interpolation_order = 2;
  bool moving() const { return static_cast<bool>(motion.z) || static_cast<bool>(motion.w); }
  double z(double t) const { return motion.z ? motion.z(t) : 0.0; }
  double w(double t) const { return motion.w ? motion.w(t) : 0.0; }
  double a(double t) const { return motion.a ? motion.a(t) : 0.0; }
};

enum class BoundaryKind { Wall, Periodic, Axis };

// 벽 속도 (r 성분, z 성분). 비어 있으면 정지 벽.
using WallVelocity = std::function<std::array<double, 2>(double r, double z, double t)>;

struct Boundary {
  BoundaryKind kind = BoundaryKind::Wall;
  WallVelocity velocity;
};

struct Boundaries {
  Boundary r_min, r_max, z_min, z_max;
};

// 단위 질량당 원천항을 r 면(su)·z 면(sw) 배열에 채우는 함수 (제조해 시험용).
using SourceFn = std::function<void(double t, Array2D& su, Array2D& sw)>;

struct Problem {
  Grid grid;
  Boundaries bc;
  std::shared_ptr<const Material> material;
  std::array<double, 2> gravity{0.0, 0.0};
  SourceFn source;
  // 물체 (없으면 nullptr). 현재 압력 방정식(투영, 압력 출력)만 물체를 반영한다.
  // 물체 근처 운동량(점성 no-slip, 이동 경계의 새 셀)은 아직 구현 전이라 step() 은 거부한다.
  std::shared_ptr<const Body> body;
  // 압축성 공기 (P3, D-03): 폴리트로프 지수 κ (> 0 이면 켬). 공기 압축률 ψ = (공기 체적분율)/(κ p),
  // p 는 절대 압력. step_implicit 에서만 쓴다 (step 은 거부). 0 이면 비압축.
  double air_kappa = 0.0;
};

// 운동량 대류의 면 값 보간 방식.
// Central: 2차 중심차분 (P1·P2 검증 방식).
// VanLeer: 수송되는 속도를 van Leer 제한자로 상류 쪽에 치우쳐 보간한다 (P3 2상용, 03 §1.3 2-H).
//   이동시키는 속도는 중심 평균 그대로이다. 매끈한 곳은 2차, 극값 근처는 1차 상류.
//   격자 끝에서 상류 두 번째 값이 없으면 그 면은 중심 평균을 쓴다.
// Consistent: VanLeer 에 질량–운동량 일관 수송을 더한 방식 (P3 방식 B, 03 §1.3).
//   운동량 제어체적의 각 면 질량 플럭스 ṁ = ρ_면 · a (ρ_면 도 같은 제한자로 보간) 로 운동량을 옮기고,
//   가속도 = −Σ ṁ (u_면 − u_중심) / (V ρ_분모)   (연속 극한에서 a·∇u 와 같음).
//   물–공기 경계에서 가벼운 공기 쪽 면은 작은 ṁ 만 가져 물 속도를 거의 바꾸지 못한다.
//   ρ_분모: 보통은 그 칸의 ρ. 2 배 넘게 무거운 유체가 들어오는 칸 (예: 물이 들어오는 공기 칸) 은
//   같은 ṁ 로 한 단계 갱신한 밀도 ρ* = ρ − Δt ∇·ṁ 로 나누어 속도가 질량 가중 평균으로 유지되게 한다
//   (ρ 로 나누면 유입률 ρ_물/ρ_공기 배로 굳어 불안정). ρ* 를 모든 칸에 쓰면 매끈한 변수 밀도에서
//   시간 1차 오차가 생겨 (2026-10-23 개발 확인) 이렇게 나눈다.
//   제한된 ṁ 로 ρ* < ρ/2 가 되는 칸은 상류 밀도(1차) ṁ 로 바꾼다 (양수 보장).
//   Δt 는 step / step_implicit 의 시간 간격 (압력 출력 compute_pressure 에서는 0).
enum class Convection { Central, VanLeer, Consistent };

struct SolverOptions {
  double projection_tolerance = 1e-12;
  double pressure_tolerance = 1e-13;
  int max_iterations = 50000;
  Convection convection = Convection::Central;
  // step_implicit 의 시간 차수. 1: 후방 오일러 (P2·P3 판정 방식). 2: 가변 간격 BDF2 (03 §1.4, P4).
  //   BDF2: (c0 u^{n+1} − c1 u^n + c2 u^{n−1})/Δt, ω = Δt/Δt_prev, c0 = (1+2ω)/(1+ω), c1 = 1+ω,
  //   c2 = ω²/(1+ω). 명시 항 (대류·중력·원천) 은 (1+ω)E^n − ω E^{n−1} 로 외삽.
  //   압력 증분은 τ = Δt/c0 로 투영하고, 압축성 공기의 압력 식에는 c2 (p^n − p^{n−1}) 항이 들어간다.
  //   첫 단계 (이전 값이 없을 때) 는 1 차.
  int time_order = 1;
  // 압축성 항 ψ = 공기/(κ p) 에 쓰는 압력의 하한 [Pa] (03 §1.4, P4 는 1e3). 0 이면 하한 없음 —
  // 공기가 있는 셀의 압력이 0 이하이면 예외 (P3 까지의 동작).
  double air_pressure_floor = 0.0;
  // 물 압력 하한 [Pa] (압축성 공기에서만, 0 이면 끔): 압력 갱신 뒤 절대압력이 이보다 낮은 셀을 이 값으로
  // 올린다 (공동화의 단순 근사). 작은 공기 셀이 한 단계에 사라질 때 물에 비물리적 인장 (−600 kPa) 과
  // 한 단계 힘 튐이 생기던 것을 막는다 (2026-10-23 비교, P4 는 2.3 kPa = 20 ℃ 증기압).
  double water_pressure_floor = 0.0;
  // 절단 셀 이동 경계 (docs/design/07_cut_cell_momentum_design.md, 개발 중, 기본 끔 = 기존 경계 보간 방식).
  // D2: 투영이 원래 열린 비율의 모든 면을 고치고, 유체가 있는 모든 셀이 압력 미지수가 된다
  //     (경계 보간 면을 닫고 그 셀의 발산을 버리던 처리를 하지 않는다). 운동량 행은 아직 기존 방식.
  bool cut_cell = false;
};

struct SolverStats {
  long steps = 0;
  long poisson_solves = 0;
  long poisson_iterations = 0;
  long poisson_refinements = 0;
  int max_poisson_iterations = 0;
  double max_relative_residual = 0.0;
  bool solver_failed = false;  // 허용치 미달 종료 발생
  // 암시적 운동량 풀이 (step_implicit)
  long momentum_solves = 0;
  long momentum_iterations = 0;
  double max_momentum_residual = 0.0;
  // 공기 압력 하한이 쓰인 셀 수의 누적 (단계마다 센 셀 수의 합)
  long pressure_floor_cells = 0;
  // water_pressure_floor 가 올린 셀 수의 누적
  long water_floor_cells = 0;
  // 절단 셀 (07 설계 D3): 고칠 수 있는 면이 없어 죽은 면을 다시 연 셀 수의 누적 (0 이어야 정상)
  long cut_cell_reopened = 0;
};

class FlowSolver {
 public:
  FlowSolver(Problem problem, SolverOptions options = {});

  const Grid& grid() const { return p_.grid; }
  const std::shared_ptr<const Body>& body() const { return p_.body; }
  const Problem& problem() const { return p_; }
  double time() const { return t_; }
  const Array2D& u() const { return u_; }
  const Array2D& w() const { return w_; }
  const SolverStats& stats() const { return stats_; }

  // 계산하는(경계에서 지정되지 않은) 면인지
  bool solves_u(int i) const { return i >= iu_lo_ && i < iu_hi_; }
  bool solves_w(int j) const { return j >= jw_lo_ && j < jw_hi_; }

  // 면 속도를 함수로 설정한 뒤 경계값을 넣고 투영한다.
  // project = false 이면 초기 속도를 투영하지 않고 그대로 둔다 (압축성 공기의 초기 압축 속도 등).
  void initialize(double t0, const std::function<double(double, double)>& u0,
                  const std::function<double(double, double)>& w0, bool project = true);

  // step_implicit 이 쓰는 압력장을 지정한다 (예: 정수압 + 대기압). 압축성 공기에서는 절대 압력.
  void set_pressure(const Array2D& p) {
    p_inc_ = p;
    hist_valid_ = false;  // BDF2 이력은 새 압력부터 다시
  }

  void step(double dt);

  // 반암시적 1차 시간 진행 (P2): 점성 항 후방 오일러(암시적), 대류·체적력 명시적,
  // 증분 압력 투영. 물체가 있으면 속도 면을 운동량 분류(body::momentum_classes)로 처리
  // (고체 면 = 물체 속도, 경계 보간면 = 경계점–바깥 면 직선 보간). 정상 상태에서는 투영 보정이
  // 0 이 되어 no-slip 이 정확히 성립한다. 정상 상태 계산용이며 시간 정확도는 1차이다.
  void step_implicit(double dt);

  // step_implicit 이 누적한 압력 (셀 중심, 유효 셀 평균 0)
  const Array2D& pressure() const { return p_inc_; }
  // step_implicit 의 압력이 정의된 셀과 투영이 고치는 면 (운동량이 값을 정한 면은 막힘).
  // 물체가 없으면 geometry() 와 같다.
  const body::CutGeometry& pressure_geometry() const { return p_.body ? pgeom_ : geom_; }

  // 마지막 step_implicit 의 속도 변화: max|Δu| / max|u|
  double last_change() const { return last_change_; }

  // 현재 시각의 속도로부터 압력을 구한다 (셀 중심, 상수 미정 → 평균 0).
  void compute_pressure(Array2D& p);

  // 셀 중심 이산 발산의 최댓값 (검사용). 물체가 있으면 유체 셀만, 경계 플럭스 포함.
  double max_divergence() const;

  // 물성(Material)이 같은 시각에서도 바뀌었을 때 (예: VOF 체적분율 갱신) 캐시한 물성·연산자를 버린다
  void material_changed();

  // 현재 시각의 컷셀 기하 (물체가 없으면 모든 면이 열림)
  const body::CutGeometry& geometry() const { return geom_; }

 private:
  // 압력(투영 퍼텐셜) 해. 두 배 정밀도 hi + lo (refined_solve 참조).
  struct Potential {
    std::vector<double> hi, lo;
  };

  struct Operator {
    double time = -1.0;
    double geometry_time = -1.0;
    bool valid = false;
    CsrMatrix matrix;
    IncompleteCholesky preconditioner;
  };

  const MaterialArrays& material_at(double t);
  // 밀도는 시각 t_rho, 물체 위치는 시각 t_geom 의 압력 연산자
  Operator& operator_at(double t_rho, double t_geom);
  void build_operator(const MaterialArrays& m, Operator& op, const body::CutGeometry& gm) const;
  void update_geometry(double t);

  void apply_boundary(Array2D& u, Array2D& w, double t) const;
  // 단위 질량당 가속도. convection/viscous/forces(중력·원천) 를 골라 넣을 수 있다 (기본: 전부).
  void rhs(const Array2D& u, const Array2D& w, double t, Array2D& fu, Array2D& fw,
           bool convection = true, bool viscous = true, bool forces = true);
  // Convection::Consistent 의 밀도 갱신에 쓰는 시간 간격 (step 중에만 0 이 아님)
  double conv_dt_ = 0.0;
  // 점성 항 V(u) = A u + c 의 A 를 기존 이산식(rhs 의 점성 부분)을 탐침해 희소 행렬로 만든다.
  // 미지수 순서: u 면 (i*nz + j), 그 뒤 w 면 (nu + i*(nz+1) + j).
  void viscous_matrix(double t, CsrMatrix& a);
  // t_rho: 밀도 시각, t_bc: 경계값(벽·물체 속도) 시각
  // extra: 우변에 더할 셀 값 (BDF2 압축성 이력 항, 없으면 nullptr)
  void project(Array2D& u, Array2D& w, double t_rho, double t_bc, Potential& phi,
               const std::vector<double>* extra = nullptr);
  // b = −(열린 면 유출 + 물체 경계 유출). body_rate 는 물체의 z 속도(투영) 또는 가속도(압력).
  void face_flux_rhs(const Array2D& fu, const Array2D& fw, double body_rate,
                     std::vector<double>& b) const;
  void remove_active_mean(std::vector<double>& b) const;
  void record(const SolveResult& r, double tolerance);

  // 경계·주기를 고려한 값 접근
  double U(const Array2D& u, int i, int j, double t) const;  // r 면, j 는 유령 행 허용
  double W(const Array2D& w, int i, int j, double t) const;  // z 면, i 는 유령 열 허용
  int wrap_r(int i) const;
  int wrap_z(int j) const;

  Problem p_;
  SolverOptions opt_;
  double t_ = 0.0;
  Array2D u_, w_;
  int iu_lo_ = 0, iu_hi_ = 0, jw_lo_ = 0, jw_hi_ = 0;
  bool periodic_r_ = false, periodic_z_ = false;

  body::CutGeometry geom_;
  double geom_time_ = -1.0;
  bool geom_valid_ = false;

  MaterialArrays mat_;
  double mat_time_ = -1.0;
  bool mat_valid_ = false;
  Operator op_;
  // 투영 해 보관 (다음 풀이의 초기 추정): 0–2 RK 단계, 3 출력 압력, 4 초기화
  std::array<Potential, 5> phi_;
  Array2D su_, sw_;
  SolverStats stats_;

  // step_implicit 상태
  Array2D p_inc_;
  double last_change_ = 0.0;
  CsrMatrix visc_;
  double visc_time_ = -1.0;
  bool visc_valid_ = false;
  std::vector<double> mom_x_;  // 직전 운동량 해 (초기 추정)
  std::vector<double> psi_vol_;  // 압축성 공기: 셀마다 ψ V / τ² (투영 연산자의 대각에 더함)
  // BDF2 이력 (time_order = 2): 직전 단계의 속도·명시 항·압력과 간격
  bool hist_valid_ = false;
  Array2D u_prev_, w_prev_, fu_prev_, fw_prev_, p_prev_;
  double dt_prev_ = 0.0;
  body::CutGeometry pgeom_;    // 운동량 모드의 압력 기하
  Operator op_mom_;            // 운동량 모드의 투영 연산자
};

}  // namespace hcw::flow
