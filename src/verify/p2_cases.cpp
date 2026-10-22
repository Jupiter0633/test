#include "verify/p2_cases.h"

#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "body/force.h"
#include "body/shape.h"
#include "verify/potential.h"
#include "verify/report_util.h"

namespace hcw::verify {

using flow::Array2D;
using flow::BoundaryKind;
using flow::FlowSolver;
using flow::Geometry;
using flow::Grid;
using flow::Problem;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct ForceSetup {
  Grid grid;
  std::shared_ptr<const body::Shape> shape;
  // 바깥 경계 속도를 이 퍼텐셜 해로 지정 (없으면 정지 벽)
  std::shared_ptr<const SpheroidPotential> outer;
  // 1-H: 물체 고정 좌표계 (물체 정지, 벽이 −A 로 가속, 유체에 관성력 −A)
  bool body_fixed_frame = false;
  double z_offset = 0.0;  // 시각 0 의 물체 중심 높이 (IBM 모드만)
  double accel = 1.0;
  double rho = 1.0;
};

// 시각 0 의 축 방향 힘 (§1.2 "시각 0 의 힘"). 유체·물체 정지에서 물체가 가속도 A 로 출발.
double initial_force(const ForceSetup& su, GridRun& run) {
  const auto t0 = std::chrono::steady_clock::now();
  const double A = su.accel;

  auto body = std::make_shared<flow::Body>();
  body->shape = su.shape;
  if (!su.body_fixed_frame) {
    const double z0 = su.z_offset;
    body->motion.z = [A, z0](double t) { return z0 + 0.5 * A * t * t; };
    body->motion.w = [A](double t) { return A * t; };
    body->motion.a = [A](double) { return A; };
  }

  Problem pb;
  pb.grid = su.grid;
  pb.body = body;
  pb.bc.r_min.kind = BoundaryKind::Axis;
  pb.bc.r_max.kind = pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Wall;
  if (su.outer) {
    auto pot = su.outer;
    const double z0 = su.z_offset;
    flow::WallVelocity v = [pot, A, z0](double r, double z, double t) {
      const double u = A * t, zc = z0 + 0.5 * A * t * t;
      const auto e = pot->velocity(r, z - zc);
      return std::array<double, 2>{u * e[0], u * e[1]};
    };
    pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
  }
  if (su.body_fixed_frame) {
    // 물체와 함께 가속하는 좌표계: 벽 속도 −A t, 관성력 −A (단위 질량당)
    flow::WallVelocity v = [A](double, double, double t) {
      return std::array<double, 2>{0.0, -A * t};
    };
    pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
    pb.gravity = {0.0, -A};
  }
  const double rho = su.rho;
  pb.material = std::make_shared<flow::SmoothMaterial>(
      [rho](double, double, double) { return rho; }, [](double, double, double) { return 0.0; },
      false);

  FlowSolver s(pb);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
  if (su.body_fixed_frame && !s.geometry().grid_fitted) {
    throw std::invalid_argument("body-fixed frame mode needs a grid-fitted body");
  }
  Array2D p;
  s.compute_pressure(p);
  // 압력은 실험실·물체 좌표계에서 같은 물리량이므로 표면 압력 힘에 좌표계 보정은 0 이다
  // (물체 좌표계의 관성력은 ∂u'/∂t = ∂u/∂t − A 로 흡수되고 압력은 바뀌지 않는다).
  const double f = body::axial_pressure_force(s.grid(), s.geometry(), p);
  run.seconds = seconds_since(t0);
  finish_run(run, s);
  run.metrics["cut_pieces"] = static_cast<double>(s.geometry().pieces.size());
  return f;
}

Grid axisym_grid(int nr, int nz, double r1, double z0, double z1) {
  return Grid(Geometry::Axisymmetric, nr, nz, 0.0, r1, z0, z1);
}

// 상대오차 ≤ limit (가장 고운 격자)
Criterion finest_at_most(const std::vector<GridRun>& runs, const std::string& key, double limit) {
  Criterion c;
  c.text = key + " on the finest grid <= " + fmt(limit);
  if (runs.empty()) return c;
  const double v = runs.back().metrics.at(key);
  c.passed = std::isfinite(v) && v <= limit;
  c.detail = fmt(v);
  return c;
}

bool all_ok(const std::vector<GridRun>& runs) {
  for (const auto& r : runs) {
    if (!r.failure.empty()) return false;
  }
  return true;
}

// 공통 통과 조건 (★): 상대오차 단조 감소, 최세 격자 ≤ 2 %, p_obs ≥ 0.8
void add_star_criteria(TestReport& rep, const std::vector<GridRun>& runs) {
  rep.criteria.push_back(monotone(runs, "force_relative_error"));
  rep.criteria.push_back(finest_at_most(runs, "force_relative_error", 0.02));
  rep.criteria.push_back(order_at_least(runs, "force_relative_error", 0.8));
}

TestReport p2_report(const std::string& id, const std::string& title) {
  TestReport rep;
  rep.id = id;
  rep.title = title;
  rep.phase = "P2";
  rep.criteria_section = "1.2";
  rep.criteria_commit = "93cf47d";
  return rep;
}

GridRun guarded(const std::string& label, int n, const std::function<void(GridRun&)>& body) {
  GridRun run;
  run.label = label;
  run.n = n;
  try {
    body(run);
  } catch (const std::exception& e) {
    run.failure = std::string("FAILED_EXCEPTION ") + e.what();
  }
  return run;
}

}  // namespace

// ---------------------------------------------------------------- 1-E
TestReport run_1e(bool quick) {
  TestReport rep = p2_report("1-E", "Sphere accelerated from rest: force at t = 0 (added mass)");
  const double a = 1.0, rho = 1.0, A = 1.0;
  auto pot = std::make_shared<SpheroidPotential>(a, a);
  const double exact = -pot->added_mass(rho) * A;  // −(2/3)πρa³A
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{8, 16, 32, 64};
  for (int n : ns) {
    rep.runs.push_back(guarded("(a) a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      ForceSetup su;
      su.grid = axisym_grid(3 * n, 6 * n, 3.0, -3.0, 3.0);
      su.shape = std::make_shared<body::Spheroid>(a, a);
      su.outer = pot;
      su.accel = A;
      su.rho = rho;
      const double f = initial_force(su, run);
      run.metrics["force"] = f;
      run.metrics["force_exact"] = exact;
      run.metrics["force_relative_error"] = std::abs(f - exact) / std::abs(exact);
    }));
  }
  std::vector<GridRun> runs_a = rep.runs;

  // (b) 보고만: 물체가 격자 3칸을 지나는 동안의 힘 이력 (점성 작게, 반암시적 진행).
  // §1.2 는 격자·시간을 정하지 않았으므로 여기서 정한다: a/Δ = 8, 16, 32, μ = 0.01, Δt = 0.05 Δ.
  const double mu_b = 0.01;
  const std::vector<int> nb = quick ? std::vector<int>{4} : std::vector<int>{8, 16, 32};
  std::vector<GridRun> runs_b;
  for (int n : nb) {
    runs_b.push_back(guarded("(b) moving, a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      const auto t0 = std::chrono::steady_clock::now();
      const double dx = a / n;
      const double t_end = std::sqrt(2.0 * 3.0 * dx / A);  // 중심이 3칸 이동
      const double dt = 0.05 * dx;
      auto shape = std::make_shared<body::Spheroid>(a, a);
      auto b = std::make_shared<flow::Body>();
      b->shape = shape;
      b->motion.z = [A](double t) { return 0.5 * A * t * t; };
      b->motion.w = [A](double t) { return A * t; };
      b->motion.a = [A](double) { return A; };
      b->interpolation_order = 1;  // 2차 보간은 이 비정상·저점성 계산에서 불안정 (body/cut_cell.h)
      Problem pb;
      pb.grid = axisym_grid(3 * n, 6 * n, 3.0, -3.0, 3.0);
      pb.body = b;
      pb.bc.r_min.kind = BoundaryKind::Axis;
      flow::WallVelocity v = [pot, A](double r, double z, double t) {
        const auto e = pot->velocity(r, z - 0.5 * A * t * t);
        return std::array<double, 2>{A * t * e[0], A * t * e[1]};
      };
      pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
      pb.material = std::make_shared<flow::SmoothMaterial>(
          [rho](double, double, double) { return rho; }, [mu_b](double, double, double) { return mu_b; },
          false);
      FlowSolver s(pb);
      s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
      std::vector<double> ratio;
      while (s.time() < t_end - 1e-12) {
        s.step_implicit(std::min(dt, t_end - s.time()));
        const double zc = b->z(s.time());
        const double fp = body::axial_pressure_force(s.grid(), s.pressure_geometry(), s.pressure());
        const double fv = body::axial_viscous_force(s.grid(), s.geometry(), *shape, zc, s.w(), mu_b, b->w(s.time()));
        ratio.push_back((fp + fv) / exact);
        if (s.stats().solver_failed) break;
      }
      run.seconds = seconds_since(t0);
      finish_run(run, s);
      double mean = 0.0, jump = 0.0;
      for (std::size_t k = 0; k < ratio.size(); ++k) {
        mean += ratio[k] / static_cast<double>(ratio.size());
        if (k > 0) jump = std::max(jump, std::abs(ratio[k] - ratio[k - 1]));
      }
      run.metrics["force_over_added_mass_mean"] = mean;
      run.metrics["force_over_added_mass_first"] = ratio.empty() ? 0.0 : ratio.front();
      run.metrics["force_over_added_mass_last"] = ratio.empty() ? 0.0 : ratio.back();
      run.metrics["max_step_jump_over_added_mass"] = jump;
      run.metrics["steps"] = static_cast<double>(ratio.size());
      std::string hist;
      const std::size_t stride = std::max<std::size_t>(1, ratio.size() / 10);
      for (std::size_t k = 0; k < ratio.size(); k += stride) hist += (k ? ", " : "") + fmt(ratio[k]);
      rep.notes.push_back("(b) a/dx=" + std::to_string(n) + " F/(-m_a A) every " + std::to_string(stride) +
                          " steps: " + hist);
    }));
  }

  rep.runs = runs_a;
  rep.runs.insert(rep.runs.end(), runs_b.begin(), runs_b.end());
  if (all_ok(runs_a)) add_star_criteria(rep, runs_a);
  rep.notes.push_back("(a) exact force = -(2/3) pi rho a^3 A (Lamb 1932 Art. 92); outer boundary = exact dipole flow");
  rep.notes.push_back("(b) report only: sphere moves 3 cells, mu = 0.01, dt = 0.05 dx, semi-implicit step, "
                      "linear boundary interpolation (quadratic was unstable here); "
                      "F includes pressure and viscous parts; max_step_jump measures fresh-cell spikes");
  if (!all_ok(runs_b)) rep.notes.push_back("(b) a report-only run failed; see its label");
  // 판정은 (a) 만으로 한다 (§1.2)
  TestReport judge = rep;
  judge.runs = runs_a;
  conclude(judge, quick);
  rep.verdict = judge.verdict;
  for (const auto& n : judge.notes) {
    if (n.find("quick mode") != std::string::npos) rep.notes.push_back(n);
  }
  return rep;
}

// ---------------------------------------------------------------- 1-F
TestReport run_1f(bool quick) {
  TestReport rep = p2_report("1-F", "Disk-like bodies accelerated from rest: force at t = 0");
  const double a = 1.0, c = 0.2, rho = 1.0, A = 1.0;

  // (a) 판정: 납작한 회전타원체
  auto pot = std::make_shared<SpheroidPotential>(a, c);
  const double exact = -pot->added_mass(rho) * A;
  std::vector<GridRun> runs_a;
  const std::vector<int> na = quick ? std::vector<int>{8, 16} : std::vector<int>{16, 32, 64, 128};
  for (int n : na) {
    runs_a.push_back(guarded("(a) spheroid c/a=0.2, a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      ForceSetup su;
      su.grid = axisym_grid(2 * n, 2 * n, 2.0, -1.0, 1.0);
      su.shape = std::make_shared<body::Spheroid>(a, c);
      su.outer = pot;
      su.accel = A;
      su.rho = rho;
      const double f = initial_force(su, run);
      run.metrics["force"] = f;
      run.metrics["force_exact"] = exact;
      run.metrics["force_relative_error"] = std::abs(f - exact) / std::abs(exact);
    }));
  }

  // (b) 보고만: 두께 있는 평판 원판, 바깥 경계는 두께 0 원판 해 (Lamb Art. 102)
  auto disk = std::make_shared<SpheroidPotential>(a, 0.0);
  const double reference = -disk->added_mass(rho) * A;  // −(8/3)ρa³A
  std::vector<GridRun> runs_b;
  const std::vector<int> nb = quick ? std::vector<int>{20, 40} : std::vector<int>{40, 80, 160};
  for (int n : nb) {
    runs_b.push_back(guarded("(b) flat disk h/a=0.1, a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      ForceSetup su;
      su.grid = axisym_grid(3 * n, 6 * n, 3.0, -3.0, 3.0);
      su.shape = std::make_shared<body::Cylinder>(a, 0.1);
      su.outer = disk;
      su.accel = A;
      su.rho = rho;
      const double f = initial_force(su, run);
      run.metrics["force"] = f;
      run.metrics["force_over_thin_disk_value"] = f / reference;
    }));
  }

  rep.runs = runs_a;
  rep.runs.insert(rep.runs.end(), runs_b.begin(), runs_b.end());
  if (all_ok(runs_a)) add_star_criteria(rep, runs_a);
  // (b) 보고: 연속 격자 사이 변화
  for (std::size_t k = 0; k + 1 < runs_b.size(); ++k) {
    if (!runs_b[k].failure.empty() || !runs_b[k + 1].failure.empty()) continue;
    rep.notes.push_back("(b) change " + runs_b[k].label.substr(runs_b[k].label.rfind(',') + 2) + " -> " +
                        runs_b[k + 1].label.substr(runs_b[k + 1].label.rfind(',') + 2) + ": " +
                        fmt(runs_b[k + 1].metrics.at("force") - runs_b[k].metrics.at("force")));
  }
  rep.notes.push_back("(a) exact added mass from Lamb 1932 Art. 114 (closed form for an oblate spheroid) = " +
                      fmt(pot->added_mass(rho)) + " rho a^3");
  rep.notes.push_back("(b) is report only (no exact solution for a finite-thickness disk); "
                      "reference (8/3) rho a^3 is the zero-thickness disk (Art. 102), also used as outer boundary");
  if (!all_ok(runs_b)) rep.notes.push_back("(b) a report-only run failed; see its label");
  // 판정은 (a) 만으로 한다 (§1.2). (b) 의 실행 실패는 판정에 넣지 않도록 conclude 전에 (a) 로 판단한다.
  TestReport judge = rep;
  judge.runs = runs_a;
  conclude(judge, quick);
  rep.verdict = judge.verdict;
  for (const auto& n : judge.notes) {
    if (n.find("quick mode") != std::string::npos) rep.notes.push_back(n);
  }
  return rep;
}

// ---------------------------------------------------------------- 1-G
TestReport run_1g(bool quick) {
  TestReport rep = p2_report("1-G", "Sphere in Stokes flow (Ua/nu = 1e-3): steady drag");
  const double a = 1.0, rho = 1.0, U = 1.0, mu = 1000.0;
  const double exact = 6.0 * kPi * mu * a * U;
  // 후방 오일러 Δt (정상 해는 Δt 와 무관). 점성 시간 ρa²/μ 의 10 배.
  const double dt = 10.0 * rho * a * a / mu;
  const int max_steps = 20000;
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{8, 16, 32, 64};
  for (int n : ns) {
    rep.runs.push_back(guarded("a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      const auto t0 = std::chrono::steady_clock::now();
      auto shape = std::make_shared<body::Spheroid>(a, a);
      auto b = std::make_shared<flow::Body>();
      b->shape = shape;  // 정지
      Problem pb;
      pb.grid = axisym_grid(4 * n, 8 * n, 4.0, -4.0, 4.0);
      pb.body = b;
      pb.bc.r_min.kind = BoundaryKind::Axis;
      flow::WallVelocity v = [a, U](double r, double z, double) {
        return stokes_sphere_velocity(a, U, r, z);
      };
      pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
      pb.material = std::make_shared<flow::SmoothMaterial>(
          [rho](double, double, double) { return rho; }, [mu](double, double, double) { return mu; },
          false);
      FlowSolver s(pb);
      s.initialize(0.0, [](double, double) { return 0.0; }, [U](double, double) { return U; });
      int steps = 0;
      do {
        s.step_implicit(dt);
        ++steps;
      } while (s.last_change() > 1e-12 && steps < max_steps && !s.stats().solver_failed);
      // 압력은 운동량 모드의 압력 기하(압력이 정의된 셀)로 적분한다
      const double fp = body::axial_pressure_force(s.grid(), s.pressure_geometry(), s.pressure());
      const double fv = body::axial_viscous_force(s.grid(), s.geometry(), *shape, 0.0, s.w(), mu);
      run.seconds = seconds_since(t0);
      finish_run(run, s);
      if (run.failure.empty() && s.last_change() > 1e-12) {
        run.failure = "FAILED_NOT_STEADY (step change " + fmt(s.last_change()) + " after " +
                      std::to_string(steps) + " steps)";
      }
      run.metrics["force"] = fp + fv;
      run.metrics["force_pressure"] = fp;
      run.metrics["force_viscous"] = fv;
      run.metrics["force_exact"] = exact;
      run.metrics["force_relative_error"] = std::abs(fp + fv - exact) / exact;
      run.metrics["steps"] = steps;
      run.metrics["final_step_change"] = s.last_change();
      run.metrics["momentum_max_residual"] = s.stats().max_momentum_residual;
    }));
  }
  if (all_ok(rep.runs)) add_star_criteria(rep, rep.runs);
  rep.notes.push_back("exact drag 6 pi mu a U (Lamb 1932 Art. 337); outer boundary = exact Stokes flow (Art. 337 eq. 11)");
  rep.notes.push_back("steady state: max velocity change per step / max velocity <= 1e-12 (as 1-B); "
                      "backward-Euler viscous step, dt = 10 rho a^2/mu");
  rep.notes.push_back("exact split: pressure 2 pi mu a U = " + fmt(exact / 3.0) + ", viscous 4 pi mu a U = " +
                      fmt(2.0 * exact / 3.0) + " (reported, not a criterion)");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 1-H
TestReport run_1h(bool quick) {
  TestReport rep = p2_report("1-H", "Body-fixed frame vs IBM (lab frame): flat disk, force at t = 0");
  rep.criteria_commit = "d4c2574";  // 1-H 는 IBM 원판을 반 칸 어긋나게 두도록 변경 (03 §1.2 변경 기록)
  const double a = 1.0, h = 0.1, rho = 1.0, A = 1.0;
  const std::vector<int> ns = quick ? std::vector<int>{20, 40} : std::vector<int>{40, 80, 160};
  std::vector<double> f_ibm, f_bf;
  for (int n : ns) {
    rep.runs.push_back(guarded("a/dx=" + std::to_string(n), n, [&](GridRun& run) {
      ForceSetup su;
      su.grid = axisym_grid(3 * n, 6 * n, 3.0, -3.0, 3.0);
      su.shape = std::make_shared<body::Cylinder>(a, h);
      su.accel = A;
      su.rho = rho;
      GridRun r_ibm, r_bf;
      su.z_offset = 0.5 * (a / n);  // IBM 모드: 반 칸 위 → 컷셀이 생긴다
      const double fi = initial_force(su, r_ibm);
      run.metrics["ibm_cut_pieces"] = r_ibm.metrics.at("cut_pieces");
      su.z_offset = 0.0;
      su.body_fixed_frame = true;
      const double fb = initial_force(su, r_bf);
      run.stats = r_ibm.stats;
      run.seconds = r_ibm.seconds + r_bf.seconds;
      run.failure = !r_ibm.failure.empty() ? r_ibm.failure : r_bf.failure;
      run.metrics["force_ibm"] = fi;
      run.metrics["force_body_fixed"] = fb;
      run.metrics["relative_difference"] = std::abs(fi - fb) / std::abs(fb);
      f_ibm.push_back(fi);
      f_bf.push_back(fb);
    }));
  }
  if (all_ok(rep.runs) && rep.runs.size() >= 3) {
    // ① 두 모드 각각 단조 수렴: 연속 격자 사이 변화량이 줄어든다
    for (const auto& [name, f] : {std::make_pair(std::string("ibm"), f_ibm),
                                  std::make_pair(std::string("body_fixed"), f_bf)}) {
      Criterion c;
      c.text = "force_" + name + ": change between successive grids decreases";
      c.passed = true;
      std::string detail;
      for (std::size_t k = 0; k + 2 < f.size(); ++k) {
        const double d1 = std::abs(f[k + 1] - f[k]), d2 = std::abs(f[k + 2] - f[k + 1]);
        if (!(d2 < d1)) c.passed = false;
        detail += (k ? ", " : "") + fmt(d1);
      }
      detail += ", " + fmt(std::abs(f.back() - f[f.size() - 2]));
      c.detail = detail;
      rep.criteria.push_back(c);
    }
    // ② 가장 고운 격자에서 차이 ≤ 1 %
    rep.criteria.push_back(finest_at_most(rep.runs, "relative_difference", 0.01));
    // ③ 차이가 격자를 곱게 할수록 감소
    rep.criteria.push_back(monotone(rep.runs, "relative_difference"));
  }
  rep.notes.push_back("closed wall box r in [0,3], z in [-3,3]; body-fixed mode: disk centred at z = 0 "
                      "(faces on grid faces); IBM mode: disk centred at z = dx/2 (cut cells)");
  rep.notes.push_back("frame correction for the surface pressure force is zero (pressure is frame independent)");
  conclude(rep, quick);
  return rep;
}

}  // namespace hcw::verify
