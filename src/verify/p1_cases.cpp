#include "verify/p1_cases.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>

#include "verify/jet.h"
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

struct ErrorNorm {
  double sum_sq = 0.0, volume = 0.0, max_abs = 0.0;
  void add(double err, double vol) {
    sum_sq += err * err * vol;
    volume += vol;
    max_abs = std::max(max_abs, std::abs(err));
  }
  double l2() const { return volume > 0.0 ? std::sqrt(sum_sq / volume) : 0.0; }
};

using Fn2 = std::function<double(double, double)>;

// 계산하는 모든 면의 속도 오차 (r·z 성분을 한 집합으로)
ErrorNorm velocity_error(const FlowSolver& s, const Fn2& u_exact, const Fn2& w_exact) {
  const Grid& g = s.grid();
  ErrorNorm e;
  for (int i = 0; i <= g.nr; ++i) {
    if (!s.solves_u(i)) continue;
    for (int j = 0; j < g.nz; ++j) {
      e.add(s.u()(i, j) - u_exact(g.rf(i), g.zc(j)), g.metric(g.rf(i)) * g.drn(i) * g.dzc(j));
    }
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) {
      if (!s.solves_w(j)) continue;
      e.add(s.w()(i, j) - w_exact(g.rc(i), g.zf(j)), g.metric(g.rc(i)) * g.drc(i) * g.dzn(j));
    }
  }
  return e;
}

// 셀 중심 압력 오차. 두 값 모두 체적 가중 평균을 뺀 뒤 비교.
ErrorNorm pressure_error(const Grid& g, const Array2D& p, const Fn2& p_exact) {
  double sum = 0.0, vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      sum += p_exact(g.rc(i), g.zc(j)) * g.cell_volume(i, j);
      vol += g.cell_volume(i, j);
    }
  }
  const double mean_exact = sum / vol;
  double mean_num = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) mean_num += p(i, j) * g.cell_volume(i, j);
  }
  mean_num /= vol;
  ErrorNorm e;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      e.add((p(i, j) - mean_num) - (p_exact(g.rc(i), g.zc(j)) - mean_exact), g.cell_volume(i, j));
    }
  }
  return e;
}

// ---------------------------------------------------------------- 1-A
// 시간 계수 (§1.1): A = 1 + 0.5 sin t, B = 0.5 cos t, C = 1 + 0.5 sin 2t
struct MmsTime {
  double a, ad, b, bd, c;
  explicit MmsTime(double t)
      : a(1.0 + 0.5 * std::sin(t)), ad(0.5 * std::cos(t)), b(0.5 * std::cos(t)),
        bd(-0.5 * std::sin(t)), c(1.0 + 0.5 * std::sin(2.0 * t)) {}
};

// 공간 부분과 그 도함수
struct MmsSpace {
  Jet u1, w1, p1, r1, r2, m1;
  MmsSpace(double r, double z) {
    const Jet R = Jet::var_r(r), Z = Jet::var_z(z);
    const Jet one_m_r2 = 1.0 - R * R;
    const Jet two_pi_z = (2.0 * kPi) * Z;
    u1 = (-2.0 * kPi) * R * one_m_r2 * one_m_r2 * cos(two_pi_z);
    w1 = 2.0 * one_m_r2 * (1.0 - 3.0 * R * R) * sin(two_pi_z);
    p1 = cos(kPi * R) * cos(two_pi_z);
    const Jet c_r2 = cos(kPi * R * R);
    r1 = c_r2 * cos(two_pi_z);
    r2 = c_r2 * sin(two_pi_z);
    m1 = 0.01 * (1.0 + 0.5 * c_r2 * sin(two_pi_z));
  }
};

double mms_rho(const MmsSpace& s, double t, bool variable) {
  return variable ? 5.5 + 4.5 * (std::cos(t) * s.r1.v + std::sin(t) * s.r2.v) : 1.0;
}

// 단위 질량당 원천항 (축대칭 운동량 방정식에 해석해를 대입)
std::array<double, 2> mms_source(const MmsSpace& s, double r, double t, bool variable) {
  const MmsTime k(t);
  const Jet u = k.a * s.u1;
  const Jet w = k.a * s.w1 + Jet::constant(k.b);
  const Jet p = k.c * s.p1;
  const Jet mu = variable ? s.m1 : Jet::constant(0.01);
  const double rho = mms_rho(s, t, variable);

  const double u_t = k.ad * s.u1.v;
  const double w_t = k.ad * s.w1.v + k.bd;
  const double conv_u = u.v * u.r + w.v * u.z;
  const double conv_w = u.v * w.r + w.v * w.z;
  const double shear = u.z + w.r;
  const double div_tau_r = 2 * mu.r * u.r + 2 * mu.v * (u.rr + u.r / r) + mu.z * shear +
                           mu.v * (u.zz + w.rz) - 2 * mu.v * u.v / (r * r);
  const double div_tau_z = mu.r * shear + mu.v * (u.rz + w.rr) + mu.v * shear / r +
                           2 * mu.z * w.z + 2 * mu.v * w.zz;
  return {u_t + conv_u + (p.r - div_tau_r) / rho, w_t + conv_w + (p.z - div_tau_z) / rho};
}

GridRun run_1a_grid_impl(int n, bool variable, flow::Convection conv, const GridPlan* plan = nullptr) {
  GridRun run;
  run.n = n;
  run.label = "N=" + std::to_string(n);
  const auto t0 = std::chrono::steady_clock::now();

  Problem pb;
  pb.grid = make_grid(plan, n, Geometry::Axisymmetric, n, n, 0.0, 1.0, 0.0, 1.0);
  const Grid g = pb.grid;
  pb.bc.r_min.kind = BoundaryKind::Axis;
  pb.bc.r_max.kind = BoundaryKind::Wall;
  pb.bc.r_max.velocity = [](double, double, double t) {
    return std::array<double, 2>{0.0, MmsTime(t).b};  // r=1: u_r = 0, u_z = B(t)
  };
  pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Periodic;

  auto rho_fn = [variable](double r, double z, double t) {
    return variable ? 5.5 + 4.5 * std::cos(kPi * r * r) * std::cos(2 * kPi * z - t) : 1.0;
  };
  auto mu_fn = [variable](double r, double z, double) {
    return variable ? 0.01 * (1.0 + 0.5 * std::cos(kPi * r * r) * std::sin(2 * kPi * z)) : 0.01;
  };
  pb.material = std::make_shared<flow::SmoothMaterial>(rho_fn, mu_fn, variable);

  // 면 위치의 공간 부분을 미리 계산해 두고 시간 계수만 곱한다.
  auto u_space = std::make_shared<std::vector<MmsSpace>>();
  auto w_space = std::make_shared<std::vector<MmsSpace>>();
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) u_space->emplace_back(std::max(g.rf(i), 1e-300), g.zc(j));
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) w_space->emplace_back(g.rc(i), g.zf(j));
  }
  pb.source = [g, u_space, w_space, variable](double t, Array2D& su, Array2D& sw) {
    for (int i = 1; i <= g.nr; ++i) {  // r = 0 면은 계산하지 않음
      for (int j = 0; j < g.nz; ++j) {
        su(i, j) = mms_source((*u_space)[static_cast<std::size_t>(i) * g.nz + j], g.rf(i), t,
                              variable)[0];
      }
    }
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j <= g.nz; ++j) {
        sw(i, j) = mms_source((*w_space)[static_cast<std::size_t>(i) * (g.nz + 1) + j], g.rc(i),
                              t, variable)[1];
      }
    }
  };

  const double t_end = 0.1;
  auto u_exact = [](double t) {
    return [t](double r, double z) {
      return MmsTime(t).a * (-2 * kPi * r * (1 - r * r) * (1 - r * r) * std::cos(2 * kPi * z));
    };
  };
  auto w_exact = [](double t) {
    return [t](double r, double z) {
      const MmsTime k(t);
      return k.a * 2 * (1 - r * r) * (1 - 3 * r * r) * std::sin(2 * kPi * z) + k.b;
    };
  };
  auto p_exact = [t_end](double r, double z) {
    return MmsTime(t_end).c * std::cos(kPi * r) * std::cos(2 * kPi * z);
  };

  flow::SolverOptions opt;
  opt.convection = conv;
  FlowSolver s(pb, opt);
  s.initialize(0.0, u_exact(0.0), w_exact(0.0));
  advance_to(s, t_end, 0.064 * g.dr_min());

  const ErrorNorm ev = velocity_error(s, u_exact(t_end), w_exact(t_end));
  Array2D p;
  s.compute_pressure(p);
  const ErrorNorm ep = pressure_error(g, p, p_exact);
  run.metrics["velocity_L2"] = ev.l2();
  run.metrics["velocity_Linf"] = ev.max_abs;
  run.metrics["pressure_L2"] = ep.l2();
  run.metrics["pressure_Linf"] = ep.max_abs;
  run.metrics["max_divergence"] = s.max_divergence();
  run.seconds = seconds_since(t0);
  finish_run(run, s);
  return run;
}

}  // namespace

GridRun run_1a_grid(int n, bool variable, flow::Convection conv) {
  return run_1a_grid_impl(n, variable, conv);
}

GridRun run_1a_grid_plan(const GridPlan& plan, int level, bool variable) {
  return run_1a_grid_impl(level, variable, flow::Convection::Central, &plan);
}

double observed_order(double e_coarse, double e_fine) {
  return std::log(e_coarse / e_fine) / std::log(2.0);
}

TestReport run_1a(bool variable, bool quick) {
  TestReport rep;
  rep.id = variable ? "1-A(b)" : "1-A(a)";
  rep.title = variable ? "MMS axisymmetric, variable density (1-10) and viscosity"
                       : "MMS axisymmetric, constant density and viscosity";
  const std::vector<int> grids = quick ? std::vector<int>{16, 32} : std::vector<int>{32, 64, 128, 256};
  for (int n : grids) {
    try {
      rep.runs.push_back(run_1a_grid(n, variable, flow::Convection::Central));
    } catch (const std::exception& e) {
      GridRun r;
      r.n = n;
      r.label = "N=" + std::to_string(n);
      r.failure = std::string("FAILED_EXCEPTION ") + e.what();
      rep.runs.push_back(r);
    }
  }
  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    for (const char* k : {"velocity_L2", "velocity_Linf", "pressure_L2", "pressure_Linf"}) {
      rep.criteria.push_back(monotone(rep.runs, k));
    }
    rep.criteria.push_back(order_at_least(rep.runs, "velocity_L2", 1.7));
    rep.criteria.push_back(order_at_least(rep.runs, "velocity_Linf", 1.7));
    rep.criteria.push_back(order_at_least(rep.runs, "pressure_L2", 1.7));
    rep.criteria.push_back(order_at_least(rep.runs, "pressure_Linf", 1.5));
  }
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 1-B
TestReport run_1b(bool quick) {
  TestReport rep;
  rep.id = "1-B";
  rep.title = "Hagen-Poiseuille pipe flow (axisymmetric), steady state";
  const std::vector<int> grids = quick ? std::vector<int>{8, 16} : std::vector<int>{16, 32, 64, 128};
  const double q_exact = kPi / 2.0;

  for (int n : grids) {
    GridRun run;
    run.n = n;
    run.label = "Nr=" + std::to_string(n);
    const auto t0 = std::chrono::steady_clock::now();
    try {
      Problem pb;
      pb.grid = Grid(Geometry::Axisymmetric, n, 4, 0.0, 1.0, 0.0, 1.0);
      const Grid g = pb.grid;
      pb.bc.r_min.kind = BoundaryKind::Axis;
      pb.bc.r_max.kind = BoundaryKind::Wall;
      pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Periodic;
      pb.material = std::make_shared<flow::SmoothMaterial>(
          [](double, double, double) { return 1.0; }, [](double, double, double) { return 1.0; },
          false);
      pb.gravity = {0.0, 4.0};

      FlowSolver s(pb);
      s.initialize(0.0, [](double, double) { return 0.0; },
                   [](double r, double) { return 1.0 - r * r; });
      const double dt = 0.25 * g.dr_min() * g.dr_min();
      const long max_steps = static_cast<long>(20.0 / dt);
      bool steady = false;
      double change = 0.0;
      for (long k = 0; k < max_steps; ++k) {
        const Array2D before = s.w();
        s.step(dt);
        double dmax = 0.0, wmax = 0.0;
        for (std::size_t q = 0; q < before.data().size(); ++q) {
          dmax = std::max(dmax, std::abs(s.w().data()[q] - before.data()[q]));
          wmax = std::max(wmax, std::abs(s.w().data()[q]));
        }
        change = dmax / wmax;
        if (change <= 1e-12) {
          steady = true;
          break;
        }
      }
      auto zero = [](double, double) { return 0.0; };
      auto w_ex = [](double r, double) { return 1.0 - r * r; };
      const ErrorNorm ev = velocity_error(s, zero, w_ex);
      ErrorNorm norm;  // 해석해의 L2 크기 (상대오차용)
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) {
          norm.add(w_ex(g.rc(i), 0.0), g.metric(g.rc(i)) * g.drc(i) * g.dzc(j));
        }
      }
      double q = 0.0;
      for (int i = 0; i < g.nr; ++i) {
        double avg = 0.0;
        for (int j = 0; j < g.nz; ++j) avg += s.w()(i, j);
        q += 2.0 * kPi * g.rc(i) * g.drc(i) * avg / g.nz;
      }
      run.metrics["velocity_L2"] = ev.l2();
      run.metrics["velocity_L2_relative"] = ev.l2() / norm.l2();
      run.metrics["velocity_Linf"] = ev.max_abs;
      run.metrics["flow_rate_relative_error"] = std::abs(q - q_exact) / q_exact;
      run.metrics["steps"] = static_cast<double>(s.stats().steps);
      run.metrics["final_step_change"] = change;
      run.seconds = seconds_since(t0);
      finish_run(run, s);
      if (!steady && run.failure.empty()) run.failure = "FAILED_NOT_STEADY";
    } catch (const std::exception& e) {
      run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    }
    rep.runs.push_back(run);
  }

  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    // 각 지표: (단조 감소 그리고 p_obs >= 1.7) 또는 모든 격자 상대오차 <= 1e-10
    auto either = [&](const std::string& key, const std::string& rel_key) {
      Criterion mono = monotone(rep.runs, key);
      Criterion ord = order_at_least(rep.runs, key, 1.7);
      bool tiny = true;
      for (const auto& r : rep.runs) tiny = tiny && r.metrics.at(rel_key) <= 1e-10;
      Criterion c;
      c.text = key + ": (monotone and p_obs >= 1.7) or relative error <= 1e-10 on all grids";
      c.passed = (mono.passed && ord.passed) || tiny;
      c.detail = mono.detail + "; " + ord.detail + (tiny ? "; all relative <= 1e-10" : "");
      return c;
    };
    rep.criteria.push_back(either("velocity_L2", "velocity_L2_relative"));
    rep.criteria.push_back(either("flow_rate_relative_error", "flow_rate_relative_error"));
  }
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 1-C
// 1-C 한 격자 (2-H 에서 대류 방식만 바꿔 다시 쓴다)
GridRun run_1c_grid(int n, flow::Convection conv) {
  const double nu = 0.01, t_end = 1.0, L = 2.0 * kPi;
  GridRun run;
  run.n = n;
  run.label = "N=" + std::to_string(n);
  const auto t0 = std::chrono::steady_clock::now();
  try {
    Problem pb;
    pb.grid = Grid(Geometry::Planar, n, n, 0.0, L, 0.0, L);
    const Grid g = pb.grid;
    pb.bc.r_min.kind = pb.bc.r_max.kind = BoundaryKind::Periodic;
    pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Periodic;
    pb.material = std::make_shared<flow::SmoothMaterial>(
        [](double, double, double) { return 1.0; }, [nu](double, double, double) { return nu; },
        false);
    auto u_ex = [nu](double t) {
      return [=](double x, double y) { return std::sin(x) * std::cos(y) * std::exp(-2 * nu * t); };
    };
    auto v_ex = [nu](double t) {
      return [=](double x, double y) { return -std::cos(x) * std::sin(y) * std::exp(-2 * nu * t); };
    };
    auto p_ex = [nu, t_end](double x, double y) {
      return 0.25 * (std::cos(2 * x) + std::cos(2 * y)) * std::exp(-4 * nu * t_end);
    };
    flow::SolverOptions opt;
    opt.convection = conv;
    FlowSolver s(pb, opt);
    s.initialize(0.0, u_ex(0.0), v_ex(0.0));
    advance_to(s, t_end, 0.2 * g.dr_min());
    const ErrorNorm ev = velocity_error(s, u_ex(t_end), v_ex(t_end));
    Array2D p;
    s.compute_pressure(p);
    const ErrorNorm ep = pressure_error(g, p, p_ex);
    run.metrics["velocity_L2"] = ev.l2();
    run.metrics["velocity_Linf"] = ev.max_abs;
    run.metrics["pressure_L2"] = ep.l2();
    run.metrics["pressure_Linf"] = ep.max_abs;
    run.metrics["max_divergence"] = s.max_divergence();
    run.seconds = seconds_since(t0);
    finish_run(run, s);
  } catch (const std::exception& e) {
    run.failure = std::string("FAILED_EXCEPTION ") + e.what();
  }
  return run;
}

TestReport run_1c(bool quick) {
  TestReport rep;
  rep.id = "1-C";
  rep.title = "Taylor-Green vortex (planar, periodic)";
  const std::vector<int> grids = quick ? std::vector<int>{16, 32} : std::vector<int>{32, 64, 128, 256};
  for (int n : grids) rep.runs.push_back(run_1c_grid(n, flow::Convection::Central));
  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    for (const char* k : {"velocity_L2", "velocity_Linf", "pressure_L2", "pressure_Linf"}) {
      rep.criteria.push_back(monotone(rep.runs, k));
    }
    rep.criteria.push_back(order_at_least(rep.runs, "velocity_L2", 1.7));
    rep.criteria.push_back(order_at_least(rep.runs, "velocity_Linf", 1.7));
    rep.criteria.push_back(order_at_least(rep.runs, "pressure_L2", 1.7));
  }
  rep.notes.push_back("pressure_Linf order is reported only (not a criterion)");
  rep.notes.push_back("time step 0.2*dx; the last step is shortened to land exactly on T = 1");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 1-D
TestReport run_1d(bool quick) {
  TestReport rep;
  rep.id = "1-D";
  rep.title = "Hydrostatic two-layer rest state, density ratio 1000";
  const double g_acc = 9.81, rho_w = 1000.0, rho_a = 1.0, height = 1.0;
  const std::vector<int> grids = quick ? std::vector<int>{16} : std::vector<int>{32, 64};
  const int steps = quick ? 50 : 1000;

  for (Geometry geom : {Geometry::Planar, Geometry::Axisymmetric}) {
    for (int n : grids) {
      for (bool centered : {false, true}) {
        GridRun run;
        run.n = n;
        run.label = std::string(geom == Geometry::Planar ? "planar" : "axisym") + " N=" +
                    std::to_string(n) + (centered ? " interface@cell-center" : " interface@face");
        const auto t0 = std::chrono::steady_clock::now();
        try {
          Problem pb;
          pb.grid = Grid(geom, n, n, 0.0, 1.0, 0.0, height);
          const Grid g = pb.grid;
          const double zi = 0.5 + (centered ? 0.5 * g.dz_min() : 0.0);
          pb.bc.r_min.kind = geom == Geometry::Axisymmetric ? BoundaryKind::Axis : BoundaryKind::Wall;
          pb.bc.r_max.kind = BoundaryKind::Wall;
          pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Wall;
          pb.material = std::make_shared<flow::LayeredMaterial>(zi, rho_w, 1e-3, rho_a, 1.8e-5);
          pb.gravity = {0.0, -g_acc};

          FlowSolver s(pb);
          auto zero = [](double, double) { return 0.0; };
          s.initialize(0.0, zero, zero);
          double umax = 0.0;
          auto track = [&]() {
            for (double v : s.u().data()) umax = std::max(umax, std::abs(v));
            for (double v : s.w().data()) umax = std::max(umax, std::abs(v));
          };
          track();
          for (int k = 0; k < steps; ++k) {
            s.step(1e-3);
            track();
          }
          Array2D p;
          s.compute_pressure(p);
          auto p_ex = [=](double, double z) {
            return z >= zi ? rho_a * g_acc * (height - z)
                           : rho_a * g_acc * (height - zi) + rho_w * g_acc * (zi - z);
          };
          const ErrorNorm ep = pressure_error(g, p, p_ex);
          run.metrics["max_speed_over_run"] = umax;
          run.metrics["pressure_error_relative"] = ep.max_abs / (rho_w * g_acc * height);
          run.seconds = seconds_since(t0);
          finish_run(run, s);
        } catch (const std::exception& e) {
          run.failure = std::string("FAILED_EXCEPTION ") + e.what();
        }
        rep.runs.push_back(run);
      }
    }
  }

  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    Criterion speed{"max |u| over the whole run <= 1e-8 m/s (all runs)", true, ""};
    Criterion pres{"max |p - p_exact| / (rho_w g H) <= 1e-10 (all runs)", true, ""};
    double worst_u = 0.0, worst_p = 0.0;
    for (const auto& r : rep.runs) {
      worst_u = std::max(worst_u, r.metrics.at("max_speed_over_run"));
      worst_p = std::max(worst_p, r.metrics.at("pressure_error_relative"));
    }
    speed.passed = worst_u <= 1e-8;
    speed.detail = "worst = " + fmt(worst_u);
    pres.passed = worst_p <= 1e-10;
    pres.detail = "worst = " + fmt(worst_p);
    rep.criteria.push_back(speed);
    rep.criteria.push_back(pres);
  }
  conclude(rep, quick);
  return rep;
}

}  // namespace hcw::verify
