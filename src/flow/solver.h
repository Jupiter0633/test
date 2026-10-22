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
};

struct SolverOptions {
  double projection_tolerance = 1e-12;
  double pressure_tolerance = 1e-13;
  int max_iterations = 50000;
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
};

class FlowSolver {
 public:
  FlowSolver(Problem problem, SolverOptions options = {});

  const Grid& grid() const { return p_.grid; }
  double time() const { return t_; }
  const Array2D& u() const { return u_; }
  const Array2D& w() const { return w_; }
  const SolverStats& stats() const { return stats_; }

  // 계산하는(경계에서 지정되지 않은) 면인지
  bool solves_u(int i) const { return i >= iu_lo_ && i < iu_hi_; }
  bool solves_w(int j) const { return j >= jw_lo_ && j < jw_hi_; }

  // 면 속도를 함수로 설정한 뒤 경계값을 넣고 투영한다.
  void initialize(double t0, const std::function<double(double, double)>& u0,
                  const std::function<double(double, double)>& w0);

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
  // 점성 항 V(u) = A u + c 의 A 를 기존 이산식(rhs 의 점성 부분)을 탐침해 희소 행렬로 만든다.
  // 미지수 순서: u 면 (i*nz + j), 그 뒤 w 면 (nu + i*(nz+1) + j).
  void viscous_matrix(double t, CsrMatrix& a);
  // t_rho: 밀도 시각, t_bc: 경계값(벽·물체 속도) 시각
  void project(Array2D& u, Array2D& w, double t_rho, double t_bc, Potential& phi);
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
  body::CutGeometry pgeom_;    // 운동량 모드의 압력 기하
  Operator op_mom_;            // 운동량 모드의 투영 연산자
};

}  // namespace hcw::flow
