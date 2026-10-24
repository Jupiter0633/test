#include "verify/report_util.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace hcw::verify {

using flow::FlowSolver;

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// 종료 시각까지 Δt 로 진행하고 마지막 단계만 남은 시간으로 맞춘다.
void advance_to(FlowSolver& s, double t_end, double dt) {
  while (s.time() < t_end - 1e-12 * std::max(1.0, t_end)) {
    s.step(std::min(dt, t_end - s.time()));
  }
}

void finish_run(GridRun& run, const FlowSolver& s) {
  run.stats = s.stats();
  if (run.stats.solver_failed) run.failure = "FAILED_SOLVER";
  for (const auto& [k, v] : run.metrics) {
    if (!std::isfinite(v)) run.failure = "FAILED_NONFINITE (" + k + ")";
  }
}

std::string fmt(double v) {
  std::ostringstream o;
  o.precision(4);
  o << v;
  return o.str();
}

// 오차가 격자 세분화마다 줄어드는지
Criterion monotone(const std::vector<GridRun>& runs, const std::string& key) {
  Criterion c;
  c.text = key + ": error decreases on every refinement";
  c.passed = runs.size() >= 2;
  std::string detail;
  for (std::size_t k = 0; k + 1 < runs.size(); ++k) {
    const double a = runs[k].metrics.at(key), b = runs[k + 1].metrics.at(key);
    if (!(b < a)) c.passed = false;
    detail += (k ? ", " : "") + fmt(a);
  }
  if (!runs.empty()) detail += ", " + fmt(runs.back().metrics.at(key));
  c.detail = detail;
  return c;
}

Criterion order_at_least(const std::vector<GridRun>& runs, const std::string& key, double min) {
  Criterion c;
  c.text = key + ": observed order (two finest grids) >= " + fmt(min);
  if (runs.size() < 2) {
    c.detail = "needs two grids";
    return c;
  }
  const double p = observed_order(runs[runs.size() - 2].metrics.at(key),
                                  runs.back().metrics.at(key));
  c.passed = std::isfinite(p) && p >= min;
  c.detail = "p_obs = " + fmt(p);
  return c;
}

flow::Grid make_grid(const GridPlan* plan, int level, flow::Geometry g, int nr, int nz, double r0, double r1, double z0,
                     double z1) {
  if (!plan) return flow::Grid(g, nr, nz, r0, r1, z0, z1);
  const double h = plan->h0;
  return flow::Grid::stretched(g, flow::stretched_faces(r0, r1, plan->rc0, plan->rc1, h, 1.05, 4.0 * h, level),
                               flow::stretched_faces(z0, z1, plan->zc0, plan->zc1, h, 1.05, 4.0 * h, level));
}

void conclude(TestReport& rep, bool quick) {
  for (const auto& run : rep.runs) {
    if (!run.failure.empty()) {
      rep.verdict = run.failure.substr(0, run.failure.find(' '));
      if (quick) rep.verdict = "QUICK_CHECK_ONLY (" + rep.verdict + ")";
      return;
    }
  }
  const bool all = std::all_of(rep.criteria.begin(), rep.criteria.end(),
                               [](const Criterion& c) { return c.passed; });
  if (quick) {
    rep.verdict = "QUICK_CHECK_ONLY";
    rep.notes.push_back("quick mode uses reduced grids; criteria shown are not a verdict");
  } else {
    rep.verdict = all ? "VERIFIED" : "FAILED_CRITERIA";
  }
}

}  // namespace hcw::verify
