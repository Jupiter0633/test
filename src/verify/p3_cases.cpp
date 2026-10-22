#include "verify/p3_cases.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "interface/vof.h"
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

}  // namespace hcw::verify
