#include "verify/p3_cases.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "interface/two_phase.h"
#include "interface/vof.h"
#include "interface/vof_material.h"
#include "verify/p1_cases.h"
#include "verify/report_util.h"

namespace hcw::verify {

using flow::Array2D;
using flow::Geometry;
using flow::Grid;

namespace {

constexpr double kPi = 3.14159265358979323846;

TestReport p3_report(const std::string& id, const std::string& title) {
  TestReport rep;
  rep.id = id;
  rep.title = title;
  rep.phase = "P3";
  rep.criteria_section = "1.3";
  rep.criteria_commit = "d74d3e8";
  return rep;
}

// 새 대류 방식 (§1.3 "2상용 대류 방식과 그 검증", 기준 고정 8ffed34) 으로 다시 하는 판정 시험.
// 2026-10-23: 방식 A (제한자) 만으로는 2-C 가 불통과 (validation/records/WS-2026-10-23-limiter.md) 여서
// 방식 B (제한자 + 질량–운동량 일관 수송, flow::Convection::Consistent) 로 바꾼다. 기준은 그대로.
constexpr const char* kLimiterCommit = "8ffed34";
constexpr flow::Convection kNewScheme = flow::Convection::Consistent;

flow::SolverOptions convection_options(bool new_scheme) {
  flow::SolverOptions o;
  if (new_scheme) o.convection = kNewScheme;
  return o;
}

void mark_limiter(TestReport& rep, bool new_scheme) {
  if (!new_scheme) return;
  rep.id += "(B)";
  rep.criteria_commit = kLimiterCommit;
  rep.notes.push_back("re-run with scheme B: van Leer limiter + consistent mass-momentum transport "
                      "(03 section 1.3, locked 8ffed34); criteria unchanged");
}

Criterion all_at_most(const std::vector<GridRun>& runs, const std::string& key, double limit) {
  Criterion c;
  c.text = key + " <= " + fmt(limit) + " on every grid";
  c.passed = !runs.empty();
  std::string detail;
  for (const auto& r : runs) {
    const double v = r.metrics.at(key);
    if (!(v <= limit)) c.passed = false;
    detail += (detail.empty() ? "" : ", ") + fmt(v);
  }
  c.detail = detail;
  return c;
}

}  // namespace

double sphere_cell_volume(double r0, double r1, double z0, double z1, double zc, double radius) {
  const double rb = std::min(r1, radius);
  if (rb <= r0) return 0.0;
  // 구간 [r0, rb] 를 수직 겹침 길이가 꺾이는 점에서 나누고 각 조각을 8점 가우스 적분
  std::vector<double> cuts{r0, rb};
  for (double zz : {z0, z1}) {
    const double d = radius * radius - (zz - zc) * (zz - zc);
    if (d > 0.0) {
      const double rr = std::sqrt(d);
      if (rr > r0 && rr < rb) cuts.push_back(rr);
    }
  }
  std::sort(cuts.begin(), cuts.end());
  static const double x[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                              -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                              0.7966664774136267,  0.9602898564975363};
  static const double w[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                              0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                              0.2223810344533745, 0.1012285362903763};
  auto overlap = [&](double r) {
    const double s = std::sqrt(std::max(0.0, radius * radius - r * r));
    return std::max(0.0, std::min(z1, zc + s) - std::max(z0, zc - s));
  };
  double v = 0.0;
  for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
    const double a = cuts[k], b = cuts[k + 1];
    // 겹침 길이는 꺾이는 점(구 가장자리 r = R, 수직 겹침이 바뀌는 r)에서 √ 꼴로 변한다.
    // r = a + (b − a)(1 − (1 − s)²)·… 대신 양 끝 모두에 √ 특이성을 없애는 치환
    // r = a + (b − a)(3s² − 2s³) (dr/ds 가 양 끝에서 0) 을 쓰고 4 등분해 가우스 적분한다.
    for (int sub = 0; sub < 4; ++sub) {
      const double sa = sub / 4.0, sb = (sub + 1) / 4.0;
      for (int q = 0; q < 8; ++q) {
        const double s = 0.5 * (sa + sb) + 0.5 * (sb - sa) * x[q];
        const double r = a + (b - a) * (3.0 * s * s - 2.0 * s * s * s);
        const double drds = (b - a) * (6.0 * s - 6.0 * s * s);
        v += 0.5 * (sb - sa) * w[q] * drds * r * overlap(r);
      }
    }
  }
  return v;
}

// ---------------------------------------------------------------- 2-B
TestReport run_2b(bool quick) {
  TestReport rep = p3_report("2-B", "Reversible transport of a spherical drop (axisymmetric)");
  const double T = 2.0, radius = 0.15, zc = 0.35;
  // Ψ₀ = (1/π) r² sin²(2πr) sin²(πz),  Ψ = Ψ₀ cos(πt/T)
  auto psi0 = [](double r, double z) {
    const double s1 = std::sin(2.0 * kPi * r), s2 = std::sin(kPi * z);
    return r * r * s1 * s1 * s2 * s2 / kPi;
  };
  const std::vector<int> ns = quick ? std::vector<int>{16, 32} : std::vector<int>{32, 64, 128, 256};
  for (int n : ns) {
    GridRun run;
    run.n = n;
    run.label = "N=" + std::to_string(n) + " (" + std::to_string(n) + "x" + std::to_string(2 * n) + ")";
    try {
      const auto t0 = std::chrono::steady_clock::now();
      const Grid g(Geometry::Axisymmetric, n, 2 * n, 0.0, 0.5, 0.0, 1.0);
      interface::Vof vof(g);
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) {
          vof.f()(i, j) = sphere_cell_volume(g.rf(i), g.rf(i + 1), g.zf(j), g.zf(j + 1), zc, radius) /
                          g.cell_volume(i);
        }
      }
      const Array2D f0 = vof.f();
      const double v0 = vof.water_volume();

      // 공간 부분 면 플럭스 (꼭짓점 Ψ₀ 의 차 → 이산 발산이 정확히 0). 영역 가장자리는 0.
      Array2D qu0(g.nr + 1, g.nz), qw0(g.nr, g.nz + 1);
      for (int i = 1; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) qu0(i, j) = -(psi0(g.rf(i), g.zf(j + 1)) - psi0(g.rf(i), g.zf(j)));
      }
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 1; j < g.nz; ++j) qw0(i, j) = psi0(g.rf(i + 1), g.zf(j)) - psi0(g.rf(i), g.zf(j));
      }
      double umax = 0.0;
      for (int i = 1; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) umax = std::max(umax, std::abs(qu0(i, j)) / (g.rf(i) * g.dz));
      }
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 1; j < g.nz; ++j) umax = std::max(umax, std::abs(qw0(i, j)) / (g.rc(i) * g.dr));
      }
      const int steps = static_cast<int>(std::ceil(T / (0.25 * g.dr / umax)));
      const double dt = T / steps;

      double vol_change = 0.0, bounds = 0.0;
      Array2D qu = qu0, qw = qw0;
      for (int k = 0; k < steps; ++k) {
        const double s = std::cos(kPi * (k + 0.5) * dt / T);  // 단계 중간 시각
        for (std::size_t q = 0; q < qu.data().size(); ++q) qu.data()[q] = s * qu0.data()[q];
        for (std::size_t q = 0; q < qw.data().size(); ++q) qw.data()[q] = s * qw0.data()[q];
        vof.advect(qu, qw, dt);
        vol_change = std::max(vol_change, std::abs(vof.water_volume() - v0) / v0);
        bounds = std::max(bounds, vof.bound_violation());
      }
      double err = 0.0;
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) err += std::abs(vof.f()(i, j) - f0(i, j)) * g.cell_volume(i);
      }
      run.metrics["shape_error"] = err / v0;
      run.metrics["volume_relative_change"] = vol_change;
      run.metrics["bound_violation"] = bounds;
      run.metrics["steps"] = steps;
      run.stats.steps = steps;
      run.seconds = seconds_since(t0);
      if (!std::isfinite(err)) run.failure = "FAILED_NONFINITE (shape_error)";
    } catch (const std::exception& e) {
      run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    }
    rep.runs.push_back(run);
  }
  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    rep.criteria.push_back(monotone(rep.runs, "shape_error"));
    rep.criteria.push_back(order_at_least(rep.runs, "shape_error", 1.0));
    rep.criteria.push_back(all_at_most(rep.runs, "volume_relative_change", 1e-8));
    rep.criteria.push_back(all_at_most(rep.runs, "bound_violation", 1e-6));
  }
  rep.notes.push_back("prescribed flow from Psi = (1/pi) r^2 sin^2(2 pi r) sin^2(pi z) cos(pi t / T), T = 2; "
                      "drop radius 0.15 centred on the axis at z = 0.35; dt = 0.25 dx / max|u|");
  rep.notes.push_back("geometric VOF: Youngs normals, PLIC, direction split with Weymouth-Yue correction");
  conclude(rep, quick);
  return rep;
}


// ---------------------------------------------------------------- 2-C
namespace {

// J₁ 의 첫 영점 (= J₀′ 의 첫 영점), 뉴턴 반복
double first_zero_j1() {
  double x = 3.8317;
  for (int it = 0; it < 20; ++it) {
    const double j1 = std::cyl_bessel_j(1.0, x);
    const double dj1 = std::cyl_bessel_j(0.0, x) - j1 / x;
    x -= j1 / dj1;
  }
  return x;
}

// 셀 [r0, r1] × [z0, z1] 에서 수면 z = surface(r) 아래 물이 차지하는 체적 (r 가중이면 축대칭), r 방향 8점 가우스
double water_below_surface(double r0, double r1, double z0, double z1,
                           const std::function<double(double)>& surface, bool axisym) {
  static const double x[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                              -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                              0.7966664774136267,  0.9602898564975363};
  static const double w[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                              0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                              0.2223810344533745, 0.1012285362903763};
  double v = 0.0;
  for (int q = 0; q < 8; ++q) {
    const double r = 0.5 * (r0 + r1) + 0.5 * (r1 - r0) * x[q];
    const double len = std::clamp(surface(r) - z0, 0.0, z1 - z0);
    v += 0.5 * (r1 - r0) * w[q] * (axisym ? r : 1.0) * len;
  }
  return v;
}

}  // namespace

TestReport run_2c(bool quick, bool limiter) {
  TestReport rep = p3_report("2-C", "Axisymmetric standing wave in a circular tank with an air layer");
  mark_limiter(rep, limiter);
  const double R = 1.0, h = 0.5, hp = 0.5, eps = 0.002, g0 = 9.81, rho_w = 1000.0, rho_a = 1.2;
  const double k = first_zero_j1() / R;
  // Lamb (1932) Art. 231 (11) + Art. 191·257: σ² = gk(ρ − ρ′)/(ρ coth kh + ρ′ coth kh′)
  const double sigma = std::sqrt(g0 * k * (rho_w - rho_a) /
                                 (rho_w / std::tanh(k * h) + rho_a / std::tanh(k * hp)));
  const double period = 2.0 * kPi / sigma;
  const std::vector<int> ns = quick ? std::vector<int>{8, 16} : std::vector<int>{16, 32, 64, 128};
  for (int n : ns) {
    GridRun run;
    run.n = n;
    run.label = "R/dx=" + std::to_string(n);
    try {
      const auto t0 = std::chrono::steady_clock::now();
      flow::Problem pb;
      pb.grid = Grid(Geometry::Axisymmetric, n, n, 0.0, R, 0.0, h + hp);
      pb.bc.r_min.kind = flow::BoundaryKind::Axis;
      pb.gravity = {0.0, -g0};
      const Grid g = pb.grid;
      Array2D f0(g.nr, g.nz);
      auto surface = [&](double r) { return h + eps * std::cyl_bessel_j(0.0, k * r); };
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) {
          f0(i, j) = water_below_surface(g.rf(i), g.rf(i + 1), g.zf(j), g.zf(j + 1), surface, true) /
                     g.cell_volume(i);
        }
      }
      auto mat = std::make_shared<interface::VofMaterial>(rho_w, 0.0, rho_a, 0.0);
      interface::TwoPhaseSolver s(pb, mat, f0, false, convection_options(limiter));
      s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
      const double v0 = s.vof().water_volume();
      // Δt = 0.25 Δ / √(g h)  (§1.3 이 정하지 않아 여기서 정함)
      const double dt = 0.25 * g.dz / std::sqrt(g0 * h);
      auto eta = [&]() {  // 축 옆 첫 열의 물 높이 − h
        double col = 0.0;
        for (int j = 0; j < g.nz; ++j) col += s.vof().f()(0, j) * g.dz;
        return col - h;
      };
      std::vector<double> crossings;
      double prev = eta(), vol_change = 0.0, amp_last = 0.0;
      // §1.3: "3 주기 동안" 측정 → 0 교차 7 개 (반주기 6 개 = 정확히 3 주기)까지 진행한다.
      // 위·아래 반주기 길이가 진폭에 비례해 조금 달라서(약 ±0.4 %), 홀수 개의 반주기로 평균하면
      // 한 반주기 몫이 남아 진동수가 치우친다 (1차 정식 실행에서 3.25 주기로 끊어 일부 격자가
      // 2.5 주기만 재어 생긴 불통과, 2026-10-23 확인). 안전 상한은 4 주기.
      while (crossings.size() < 7 && s.time() < 4.0 * period) {
        s.step(dt);
        const double e = eta();
        if ((prev > 0.0) != (e > 0.0)) crossings.push_back(s.time() - dt * e / (e - prev));
        if (crossings.size() >= 5) amp_last = std::max(amp_last, std::abs(e));  // 마지막 한 주기
        prev = e;
        vol_change = std::max(vol_change, std::abs(s.vof().water_volume() - v0) / v0);
        if (s.flow().stats().solver_failed) break;
      }
      run.stats = s.flow().stats();
      if (run.stats.solver_failed) run.failure = "FAILED_SOLVER";
      if (crossings.size() < 7) throw std::runtime_error("fewer than seven zero crossings in 4 periods");
      const double t_meas = (crossings[6] - crossings[0]) / 3.0;
      const double omega = 2.0 * kPi / t_meas;
      run.metrics["frequency"] = omega;
      run.metrics["frequency_exact"] = sigma;
      run.metrics["frequency_relative_error"] = std::abs(omega - sigma) / sigma;
      run.metrics["volume_relative_change"] = vol_change;
      run.metrics["amplitude_last_period_over_initial"] = amp_last / eps;
      run.metrics["zero_crossings"] = static_cast<double>(crossings.size());
      run.seconds = seconds_since(t0);
    } catch (const std::exception& e) {
      run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    }
    rep.runs.push_back(run);
  }
  if (std::all_of(rep.runs.begin(), rep.runs.end(), [](const GridRun& r) { return r.failure.empty(); })) {
    rep.criteria.push_back(monotone(rep.runs, "frequency_relative_error"));
    Criterion fine;
    fine.text = "frequency_relative_error on the finest grid <= 0.005";
    fine.passed = rep.runs.back().metrics.at("frequency_relative_error") <= 0.005;
    fine.detail = fmt(rep.runs.back().metrics.at("frequency_relative_error"));
    rep.criteria.push_back(fine);
    rep.criteria.push_back(order_at_least(rep.runs, "frequency_relative_error", 0.8));
    rep.criteria.push_back(all_at_most(rep.runs, "volume_relative_change", 1e-8));
  }
  rep.notes.push_back("exact sigma = " + fmt(sigma) + " rad/s (Lamb 1932 Art. 231 eq. 11 with k R = " +
                      fmt(k * R) + " from Art. 191; the combination is derived in 03 section 1.3)");
  rep.notes.push_back("inviscid, incompressible air; eps = 0.002 m; dt = 0.25 dx / sqrt(g h); frequency from "
                      "7 zero crossings (exactly 3 periods) of the water height in the axis column; amplitude (last period) is report only");
  conclude(rep, quick);
  return rep;
}


// ---------------------------------------------------------------- 공통 (2-A, 2-F, 2-D)
namespace {

constexpr double kRhoW = 1000.0, kRhoA = 1.2, kMuW = 1e-3, kMuA = 1.8e-5, kG = 9.81;

// 셀마다 z < zs (수평 수면) 아래 물 비율
Array2D flat_fraction(const Grid& g, double zs) {
  Array2D f(g.nr, g.nz);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) f(i, j) = std::clamp((zs - g.zf(j)) / g.dz, 0.0, 1.0);
  }
  return f;
}

// 같은 문제를 비압축으로 풀어 얻은 이산 정수압 (평균 0)
Array2D hydrostatic_pressure(const flow::Problem& pb, const Array2D& f) {
  flow::Problem q = pb;
  q.air_kappa = 0.0;
  auto mat = std::make_shared<interface::VofMaterial>(kRhoW, 0.0, kRhoA, 0.0);
  interface::TwoPhaseSolver s(q, mat, f);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
  Array2D p;
  s.flow().compute_pressure(p);
  return p;
}

double max_speed(const flow::FlowSolver& s) {
  double m = 0.0;
  for (double v : s.u().data()) m = std::max(m, std::abs(v));
  for (double v : s.w().data()) m = std::max(m, std::abs(v));
  return m;
}

bool runs_ok(const std::vector<GridRun>& runs) {
  return std::all_of(runs.begin(), runs.end(), [](const GridRun& r) { return r.failure.empty(); });
}

}  // namespace

// ---------------------------------------------------------------- 2-A
TestReport run_2a(bool quick, bool limiter) {
  TestReport rep = p3_report("2-A", "Water-air rest state (interface inside a cell), 1000 steps");
  mark_limiter(rep, limiter);
  const std::vector<int> ns = quick ? std::vector<int>{8, 16} : std::vector<int>{32, 64};
  const int steps = quick ? 100 : 1000;
  for (Geometry geom : {Geometry::Planar, Geometry::Axisymmetric}) {
    for (int n : ns) {
      for (bool compressible : {false, true}) {
        GridRun run;
        run.n = n;
        run.label = std::string(geom == Geometry::Planar ? "planar" : "axisym") + " N=" + std::to_string(n) +
                    (compressible ? " compressible air" : " incompressible air");
        try {
          const auto t0 = std::chrono::steady_clock::now();
          flow::Problem pb;
          pb.grid = Grid(geom, n, n, 0.0, 1.0, 0.0, 1.0);
          pb.bc.r_min.kind = geom == Geometry::Axisymmetric ? flow::BoundaryKind::Axis : flow::BoundaryKind::Wall;
          pb.gravity = {0.0, -kG};
          const Array2D f0 = flat_fraction(pb.grid, 0.5 + 0.3 * pb.grid.dz);
          Array2D p0;
          if (compressible) {
            p0 = hydrostatic_pressure(pb, f0);
            for (double& v : p0.data()) v += 1e5;
            pb.air_kappa = 1.4;
          }
          auto mat = std::make_shared<interface::VofMaterial>(kRhoW, kMuW, kRhoA, kMuA);
          interface::TwoPhaseSolver s(pb, mat, f0, compressible, convection_options(limiter));
          s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; },
                       !compressible);
          if (compressible) s.flow().set_pressure(p0);
          const double v0 = s.vof().water_volume();
          double umax = 0.0, vol = 0.0;
          for (int k = 0; k < steps; ++k) {
            s.step(1e-3);
            umax = std::max(umax, max_speed(s.flow()));
            vol = std::max(vol, std::abs(s.vof().water_volume() - v0) / v0);
          }
          run.metrics["max_speed_over_run"] = umax;
          run.metrics["volume_relative_change"] = vol;
          run.stats = s.flow().stats();
          if (run.stats.solver_failed) run.failure = "FAILED_SOLVER";
          run.seconds = seconds_since(t0);
        } catch (const std::exception& e) {
          run.failure = std::string("FAILED_EXCEPTION ") + e.what();
        }
        rep.runs.push_back(run);
      }
    }
  }
  if (runs_ok(rep.runs)) {
    rep.criteria.push_back(all_at_most(rep.runs, "max_speed_over_run", 1e-8));
    rep.criteria.push_back(all_at_most(rep.runs, "volume_relative_change", 1e-12));
  }
  rep.notes.push_back("interface at z = 0.5 + 0.3 dx; incompressible runs use the 3-stage RK step, compressible "
                      "runs (kappa = 1.4) the semi-implicit step started from the discrete hydrostatic pressure + 1e5 Pa");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 2-F
TestReport run_2f(bool quick, bool limiter) {
  TestReport rep = p3_report("2-F", "Air springs: water column between two closed air pockets");
  rep.criteria_commit = "a48d360";  // 해석해에 공기층 관성 포함 (03 §1.3 2-F 변경 기록)
  mark_limiter(rep, limiter);
  const double la = 0.2, lw = 1.0, p0 = 1e5, v0 = 0.01;
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{8, 16, 32, 64};
  std::vector<std::vector<GridRun>> groups;
  for (double kappa : {1.4, 1.0}) {
    // 유효 질량에 두 공기층의 관성 ρ_a L_a/3 씩을 넣는다 (§1.3 2-F, 2026-10-23 변경)
    const double mass = kRhoW * lw + 2.0 * kRhoA * la / 3.0;
    const double omega = std::sqrt(kappa * p0 * (2.0 / la) / mass);
    const double period = 2.0 * kPi / omega;
    const double dp_exact = kappa * p0 * (v0 / omega) / la;
    std::vector<GridRun> runs;
    for (int n : ns) {
      GridRun run;
      run.n = n;
      run.label = "kappa=" + fmt(kappa) + " L_a/dz=" + std::to_string(n);
      try {
        const auto t0 = std::chrono::steady_clock::now();
        const double dz = la / n;
        flow::Problem pb;
        pb.grid = Grid(Geometry::Planar, 2, static_cast<int>(std::lround((2 * la + lw) / dz)), 0.0, 0.1, 0.0,
                       2 * la + lw);
        pb.bc.r_min.kind = flow::BoundaryKind::Wall;
        pb.air_kappa = kappa;
        const Grid& g = pb.grid;
        Array2D f0(g.nr, g.nz);
        for (int i = 0; i < g.nr; ++i) {
          for (int j = 0; j < g.nz; ++j) {
            f0(i, j) = std::clamp((std::min(la + lw, g.zf(j + 1)) - std::max(la, g.zf(j))) / g.dz, 0.0, 1.0);
          }
        }
        auto mat = std::make_shared<interface::VofMaterial>(kRhoW, 0.0, kRhoA, 0.0);
        interface::TwoPhaseSolver s(pb, mat, f0, true, convection_options(limiter));
        auto w_init = [&](double, double z) {
          if (z <= la) return v0 * z / la;
          if (z >= la + lw) return v0 * (2 * la + lw - z) / la;
          return v0;
        };
        s.initialize(0.0, [](double, double) { return 0.0; }, w_init, false);
        s.flow().set_pressure(Array2D(g.nr, g.nz, p0));
        const double vol0 = s.vof().water_volume();
        // Δt ∝ Δ: 1 주기에 50·(L_a/Δz) 단계 (후방 오일러 감쇠가 격자와 함께 줄도록; §1.3 이 정하지 않아 여기서 정함)
        const double dt = period / (50.0 * n);
        auto pocket = [&]() {  // 아래 공기층 평균 압력 − p₀
          double sp = 0.0, sv = 0.0;
          for (int i = 0; i < g.nr; ++i) {
            for (int j = 0; j < g.nz && g.zf(j + 1) <= la + 1e-12; ++j) {
              sp += s.flow().pressure()(i, j) * g.cell_volume(i);
              sv += g.cell_volume(i);
            }
          }
          return sp / sv - p0;
        };
        std::vector<double> crossings;
        double prev = 0.0, amp = 0.0, vol = 0.0;
        bool first = true;
        // 2-C 와 같이 0 교차 7 개 (정확히 3 주기) 까지 진행한다. 안전 상한 4.5 주기.
        while (crossings.size() < 7 && s.time() < 4.5 * period) {
          s.step(dt);
          const double pp = pocket();
          if (s.time() <= period) amp = std::max(amp, std::abs(pp));
          if (!first && (prev > 0.0) != (pp > 0.0)) crossings.push_back(s.time() - dt * pp / (pp - prev));
          prev = pp;
          first = false;
          vol = std::max(vol, std::abs(s.vof().water_volume() - vol0) / vol0);
          if (s.flow().stats().solver_failed) break;
        }
        run.stats = s.flow().stats();
        if (run.stats.solver_failed) run.failure = "FAILED_SOLVER";
        if (crossings.size() < 7) throw std::runtime_error("fewer than seven zero crossings in 4.5 periods");
        const double om = 2.0 * kPi / ((crossings[6] - crossings[0]) / 3.0);
        run.metrics["frequency"] = om;
        run.metrics["frequency_exact"] = omega;
        run.metrics["frequency_relative_error"] = std::abs(om - omega) / omega;
        run.metrics["pressure_amplitude"] = amp;
        run.metrics["pressure_amplitude_exact"] = dp_exact;
        run.metrics["pressure_amplitude_relative_error"] = std::abs(amp - dp_exact) / dp_exact;
        run.metrics["volume_relative_change"] = vol;
        run.seconds = seconds_since(t0);
      } catch (const std::exception& e) {
        run.failure = std::string("FAILED_EXCEPTION ") + e.what();
      }
      runs.push_back(run);
    }
    if (runs_ok(runs)) {
      const std::string tag = "kappa=" + fmt(kappa) + ": ";
      auto add = [&](Criterion c) {
        c.text = tag + c.text;
        rep.criteria.push_back(c);
      };
      add(monotone(runs, "frequency_relative_error"));
      Criterion fine;
      fine.text = "frequency_relative_error on the finest grid <= 0.005";
      fine.passed = runs.back().metrics.at("frequency_relative_error") <= 0.005;
      fine.detail = fmt(runs.back().metrics.at("frequency_relative_error"));
      add(fine);
      add(order_at_least(runs, "frequency_relative_error", 0.8));
      Criterion pa;
      pa.text = "pressure_amplitude_relative_error on the finest grid <= 0.01";
      pa.passed = runs.back().metrics.at("pressure_amplitude_relative_error") <= 0.01;
      pa.detail = fmt(runs.back().metrics.at("pressure_amplitude_relative_error"));
      add(pa);
      add(all_at_most(runs, "volume_relative_change", 1e-8));
    }
    rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());
    rep.notes.push_back("kappa=" + fmt(kappa) + ": exact omega = " + fmt(omega) + " rad/s, pressure amplitude " +
                        fmt(dp_exact) + " Pa (derived with air-layer inertia, 03 section 1.3)");
  }
  rep.notes.push_back("planar tube 0.1 m wide (2 cells), no gravity, inviscid; semi-implicit step (backward Euler, "
                      "incremental pressure) with dt = period / (50 L_a/dz); pressure = bottom pocket mean");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 2-D
// 2026-10-23 보고만 → 판정 시험으로 변경 (03 §1.3, 기준 고정 8ffed34). 대류는 새 방식 (B).
// 이전 보고만 실행의 기록은 validation/records/WS-2026-10-22.md 에 있다.
TestReport run_2d(bool quick) {
  TestReport rep = p3_report("2-D", "Axisymmetric water column collapse (stability)");
  rep.criteria_commit = kLimiterCommit;
  const double a = 0.05;  // 기둥 반지름, 높이 2a
  const std::vector<int> ns = quick ? std::vector<int>{8} : std::vector<int>{16, 32};
  const double t_star_end = 3.0, t_scale = std::sqrt(a / kG);
  const double v_limit = 3.0 * std::sqrt(2.0 * kG * 2.0 * a);  // 기둥 높이 자유낙하 속도의 3 배
  for (int n : ns) {
    GridRun run;
    run.n = n;
    run.label = "a/dx=" + std::to_string(n);
    try {
      const auto t0 = std::chrono::steady_clock::now();
      flow::Problem pb;
      pb.grid = Grid(Geometry::Axisymmetric, 5 * n, 3 * n, 0.0, 5 * a, 0.0, 3 * a);
      pb.bc.r_min.kind = flow::BoundaryKind::Axis;
      pb.gravity = {0.0, -kG};
      const Grid& g = pb.grid;
      Array2D f0(g.nr, g.nz);
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) f0(i, j) = (g.rc(i) < a && g.zc(j) < 2 * a) ? 1.0 : 0.0;
      }
      auto mat = std::make_shared<interface::VofMaterial>(kRhoW, kMuW, kRhoA, kMuA);
      interface::TwoPhaseSolver s(pb, mat, f0, false, convection_options(true));
      s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
      const double v0 = s.vof().water_volume();
      // 적응 시간 간격: 고정 간격은 퍼지는 앞머리 속도를 따라가지 못해 VOF CFL 한계(0.5)를 넘었다
      // (2026-10-22 워크스테이션 1차 실행). 한 단계에 최대 1/4 칸만 움직이도록 줄인다.
      const double dt_max = 0.1 * g.dr / std::sqrt(kG * a);
      std::string hist;
      double next = 0.0, vol = 0.0, dt_min = dt_max, vmax_run = 0.0;
      while (s.time() < t_star_end * t_scale) {
        const double vmax = max_speed(s.flow());
        const double dt = vmax > 0.0 ? std::min(dt_max, 0.25 * std::min(g.dr, g.dz) / vmax) : dt_max;
        dt_min = std::min(dt_min, dt);
        s.step(dt);
        vol = std::max(vol, std::abs(s.vof().water_volume() - v0) / v0);
        vmax_run = std::max(vmax_run, max_speed(s.flow()));
        if (s.time() >= next * t_scale) {
          double front = 0.0;  // 바닥 행에서 물이 닿은 반경 (물 비율 합 × Δr)
          for (int i = 0; i < g.nr; ++i) front += s.vof().f()(i, 0) * g.dr;
          hist += (hist.empty() ? "" : ", ") + fmt(s.time() / t_scale) + ":" + fmt(front / a);
          next += 0.5;
        }
        // 폭주하면 끝까지 갈 필요가 없다 (속도 기준의 10 배에서 중단)
        if (s.flow().stats().solver_failed || !(vmax_run <= 10.0 * v_limit)) break;
      }
      run.stats = s.flow().stats();
      if (run.stats.solver_failed) run.failure = "FAILED_SOLVER";
      if (s.time() < t_star_end * t_scale - 1e-12 && run.failure.empty()) run.failure = "FAILED_UNSTABLE";
      run.metrics["volume_relative_change"] = vol;
      run.metrics["end_time_over_scale"] = s.time() / t_scale;
      run.metrics["smallest_dt_over_initial"] = dt_min / dt_max;
      run.metrics["max_speed_at_end"] = max_speed(s.flow());
      run.metrics["max_speed_over_run"] = vmax_run;
      rep.notes.push_back("a/dx=" + std::to_string(n) + " front r/a at t sqrt(g/a): " + hist);
      run.seconds = seconds_since(t0);
    } catch (const std::exception& e) {
      run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    }
    rep.runs.push_back(run);
  }
  if (runs_ok(rep.runs)) {
    Criterion done;
    done.text = "every grid reaches t sqrt(g/a) = 3";
    done.passed = true;
    for (const auto& r : rep.runs) {
      done.passed = done.passed && r.metrics.at("end_time_over_scale") >= t_star_end - 1e-9;
      done.detail += (done.detail.empty() ? "" : ", ") + fmt(r.metrics.at("end_time_over_scale"));
    }
    rep.criteria.push_back(done);
    rep.criteria.push_back(all_at_most(rep.runs, "max_speed_over_run", v_limit));
    rep.criteria.push_back(all_at_most(rep.runs, "volume_relative_change", 1e-8));
  }
  rep.notes.push_back("speed limit = 3 sqrt(2 g 2a) = " + fmt(v_limit) + " m/s (3 x free fall from the column height)");
  rep.notes.push_back("front position vs Martin & Moyce [R-17] is report only (data not verified)");
  rep.notes.push_back("scheme B convection (limiter + consistent mass-momentum transport); adaptive dt = min(0.1 dx / sqrt(g a), 0.25 dx / max|u|)");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 2-H
// 새 대류 방식의 정확도: 1-C 와 1-A(b) 를 조건 그대로, 대류만 새 방식 (B) 으로 (03 §1.3, 기준 고정 8ffed34)
TestReport run_2h(bool quick) {
  TestReport rep = p3_report("2-H", "New convection scheme accuracy: 1-C and 1-A(b) re-solved");
  rep.criteria_commit = kLimiterCommit;
  const std::vector<int> grids = quick ? std::vector<int>{16, 32} : std::vector<int>{32, 64, 128, 256};
  for (const char* which : {"1-C", "1-A(b)"}) {
    std::vector<GridRun> runs;
    for (int n : grids) {
      GridRun run;
      try {
        run = std::string(which) == "1-C" ? run_1c_grid(n, kNewScheme)
                                          : run_1a_grid(n, true, kNewScheme);
      } catch (const std::exception& e) {
        run.n = n;
        run.failure = std::string("FAILED_EXCEPTION ") + e.what();
      }
      run.label = std::string(which) + " N=" + std::to_string(n);
      runs.push_back(run);
    }
    if (runs_ok(runs)) {
      for (const char* k : {"velocity_L2", "pressure_L2"}) {
        Criterion c = monotone(runs, k);
        c.text = std::string(which) + ": " + c.text;
        rep.criteria.push_back(c);
        c = order_at_least(runs, k, 1.5);
        c.text = std::string(which) + ": " + c.text;
        rep.criteria.push_back(c);
      }
      for (const char* k : {"velocity_Linf", "pressure_Linf"}) {
        const auto& a = runs[runs.size() - 2].metrics.at(k);
        const auto& b = runs.back().metrics.at(k);
        rep.notes.push_back(std::string(which) + ": " + k + " order (two finest, report only) = " +
                            fmt(observed_order(a, b)));
      }
    }
    rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());
  }
  rep.notes.push_back("conditions as in 03 sections 1.1 1-C and 1-A(b); only the momentum convection changes "
                      "(scheme B: van Leer limited velocity and density on the faces, centred advecting velocity, "
                      "consistent mass-momentum update)");
  conclude(rep, quick);
  return rep;
}

}  // namespace hcw::verify
