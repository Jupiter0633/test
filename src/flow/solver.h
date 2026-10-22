#pragma once

// 2차원(평면/축대칭) 비압축 단상·변수밀도 Navier–Stokes 해석기 (P1).
// 엇갈림 격자 유한체적, 2차 중심차분, 3단 SSP-RK 시간적분, 각 단계 투영.
// 근거: docs/design/02_design_document.md §3.1, §3.5, §3.6

#include <array>
#include <functional>
#include <memory>

#include "flow/grid.h"
#include "flow/material.h"
#include "flow/poisson.h"

namespace hcw::flow {

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

  // 현재 시각의 속도로부터 압력을 구한다 (셀 중심, 상수 미정 → 평균 0).
  void compute_pressure(Array2D& p);

  // 셀 중심 이산 발산의 최댓값 (검사용)
  double max_divergence() const;

 private:
  // 압력(투영 퍼텐셜) 해. 두 배 정밀도 hi + lo (refined_solve 참조).
  struct Potential {
    std::vector<double> hi, lo;
  };

  struct Operator {
    double time = -1.0;
    bool valid = false;
    CsrMatrix matrix;
    IncompleteCholesky preconditioner;
  };

  const MaterialArrays& material_at(double t);
  Operator& operator_at(double t);
  void build_operator(const MaterialArrays& m, Operator& op) const;

  void apply_boundary(Array2D& u, Array2D& w, double t) const;
  void rhs(const Array2D& u, const Array2D& w, double t, Array2D& fu, Array2D& fw);
  void project(Array2D& u, Array2D& w, double t_rho, Potential& phi);
  void face_flux_rhs(const Array2D& fu, const Array2D& fw, std::vector<double>& b) const;
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

  MaterialArrays mat_;
  double mat_time_ = -1.0;
  bool mat_valid_ = false;
  Operator op_;
  // 투영 해 보관 (다음 풀이의 초기 추정): 0–2 RK 단계, 3 출력 압력, 4 초기화
  std::array<Potential, 5> phi_;
  Array2D su_, sw_;
  SolverStats stats_;
};

}  // namespace hcw::flow
