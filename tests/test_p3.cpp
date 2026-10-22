// P3 단위 시험: PLIC 기하, VOF 수송의 보존과 범위.

#include <cmath>
#include <iostream>
#include <string>

#include "interface/plic.h"
#include "interface/vof.h"
#include "verify/p3_cases.h"

using namespace hcw;
using flow::Geometry;
using flow::Grid;

namespace {

int g_failures = 0, g_checks = 0;

void expect(bool ok, const std::string& what) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::cerr << "FAIL: " << what << "\n";
  }
}

void test_polygon() {
  const auto rect = interface::rectangle(0.2, 0.5, 1.0, 1.3);
  expect(std::abs(interface::polygon_volume(rect, Geometry::Planar) - 0.09) < 1e-15, "planar rectangle area");
  // ∫∫ r dr dz = (r1² − r0²)/2 · Δz
  expect(std::abs(interface::polygon_volume(rect, Geometry::Axisymmetric) - 0.5 * (0.25 - 0.04) * 0.3) < 1e-15,
         "axisymmetric rectangle volume");
  // 대각선으로 반 자르기: 평면은 정확히 절반
  const std::array<double, 2> n{std::sqrt(0.5), std::sqrt(0.5)};
  const double mid = n[0] * 0.35 + n[1] * 1.15;
  expect(std::abs(interface::polygon_volume(interface::clip(rect, n, mid), Geometry::Planar) - 0.045) < 1e-15,
         "diagonal clip halves a square");
  // α 풀이 왕복: 여러 법선·목표 비율
  for (double ang : {0.3, 1.1, 2.0, 2.9, 4.0, 5.5}) {
    for (double frac : {0.01, 0.3, 0.5, 0.97}) {
      const std::array<double, 2> m{std::cos(ang), std::sin(ang)};
      const double target = frac * interface::polygon_volume(rect, Geometry::Axisymmetric);
      const double a = interface::solve_alpha(rect, m, target, Geometry::Axisymmetric);
      const double got = interface::polygon_volume(interface::clip(rect, m, a), Geometry::Axisymmetric);
      expect(std::abs(got - target) < 1e-13, "alpha round trip ang=" + std::to_string(ang) +
                                                 " frac=" + std::to_string(frac));
    }
  }
}

void test_sphere_volume() {
  // 셀 겹침 체적의 합 = 구 체적 / 2π = (2/3) R³
  const Grid g(Geometry::Axisymmetric, 20, 40, 0.0, 0.5, 0.0, 1.0);
  double v = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      v += verify::sphere_cell_volume(g.rf(i), g.rf(i + 1), g.zf(j), g.zf(j + 1), 0.35, 0.15);
    }
  }
  expect(std::abs(v / (2.0 / 3.0 * 0.15 * 0.15 * 0.15) - 1.0) < 1e-9,
         "sphere overlap volumes sum to the sphere volume (" + std::to_string(v) + ")");
}

// 흐름함수 Ψ = sin²(2πr) sin²(πz) r² / π (축대칭) 또는 sin²(2πx) sin²(πy) / π (평면) 의 회전 흐름.
// 면 플럭스를 꼭짓점 Ψ 차로 만들면 이산 발산이 정확히 0 → 물 체적이 반올림 수준까지 보존돼야 한다.
void test_rotation(Geometry geom) {
  const Grid g(geom, 16, 32, 0.0, 0.5, 0.0, 1.0);
  const double kPi = 3.14159265358979323846;
  auto psi = [&](double r, double z) {
    const double s1 = std::sin(2 * kPi * r), s2 = std::sin(kPi * z);
    return (geom == Geometry::Axisymmetric ? r * r : 1.0) * s1 * s1 * s2 * s2 / kPi;
  };
  interface::Vof vof(g);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      vof.f()(i, j) =
          geom == Geometry::Axisymmetric
              ? verify::sphere_cell_volume(g.rf(i), g.rf(i + 1), g.zf(j), g.zf(j + 1), 0.35, 0.15) / g.cell_volume(i)
              : (g.zc(j) < 0.4 + 0.3 * g.rc(i) ? 1.0 : 0.0);
    }
  }
  if (geom == Geometry::Planar) vof.f()(3, 13) = 0.37;  // 섞인 셀 하나
  const double v0 = vof.water_volume();
  flow::Array2D qu(g.nr + 1, g.nz), qw(g.nr, g.nz + 1);
  for (int i = 1; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) qu(i, j) = -(psi(g.rf(i), g.zf(j + 1)) - psi(g.rf(i), g.zf(j)));
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 1; j < g.nz; ++j) qw(i, j) = psi(g.rf(i + 1), g.zf(j)) - psi(g.rf(i), g.zf(j));
  }
  double worst_bound = 0.0;
  for (int k = 0; k < 40; ++k) {
    vof.advect(qu, qw, 0.004);
    worst_bound = std::max(worst_bound, vof.bound_violation());
  }
  const std::string tag = geom == Geometry::Axisymmetric ? "axisymmetric" : "planar";
  expect(std::abs(vof.water_volume() - v0) <= 1e-13 * v0,
         "divergence-free flow conserves water volume (" + tag + ", " +
             std::to_string(std::abs(vof.water_volume() - v0) / v0) + ")");
  expect(worst_bound <= 1e-10, "f stays within [0, 1] (" + tag + ", " + std::to_string(worst_bound) + ")");
}

void test_refuses_bad_input() {
  const Grid g(Geometry::Planar, 4, 4, 0.0, 1.0, 0.0, 1.0);
  interface::Vof vof(g);
  flow::Array2D qu(5, 4), qw(4, 5);
  qw(1, 0) = 1.0;  // 가장자리 면
  bool threw = false;
  try {
    vof.advect(qu, qw, 0.01);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "flux through a domain boundary is refused");
  qw(1, 0) = 0.0;
  qw(1, 2) = 1.0;
  threw = false;
  try {
    vof.advect(qu, qw, 1.0);  // CFL 4
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "CFL above 0.5 is refused");
}

}  // namespace

int main() {
  test_polygon();
  test_sphere_volume();
  test_rotation(Geometry::Planar);
  test_rotation(Geometry::Axisymmetric);
  test_refuses_bad_input();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
