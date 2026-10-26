// P4 단위 시험: 셀 유체 비율, 쐐기 형상, 물체와 함께 쓰는 VOF 플럭스의 발산, GCI 계산.

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "body/cut_cell.h"
#include "body/shape.h"
#include "interface/two_phase.h"
#include "verify/p4_cases.h"

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

// 원기둥 (반지름 0.3, 높이 0.25, 격자면에 맞지 않음) 의 고체 체적 = Σ (1 − φ) V 가 π R² H / (2π) 와 같은지
void test_fluid_fraction_cylinder() {
  for (Geometry geom : {Geometry::Axisymmetric, Geometry::Planar}) {
    const Grid g(geom, 16, 20, 0.0, 1.0, -0.6, 0.6);
    const body::Cylinder c(0.3, 0.25);
    const double zc = 0.037;
    const flow::Array2D phi = body::fluid_fraction(g, c, zc);
    double solid = 0.0;
    for (int i = 0; i < g.nr; ++i) {
      for (int j = 0; j < g.nz; ++j) solid += (1.0 - phi(i, j)) * g.cell_volume(i, j);
    }
    const double exact = geom == Geometry::Axisymmetric ? 0.5 * 0.3 * 0.3 * 0.25 : 0.3 * 0.25;
    expect(std::abs(solid - exact) < 1e-12 * exact,
           std::string(geom == Geometry::Axisymmetric ? "axisym" : "planar") + " cylinder solid volume");
  }
}

void test_wedge_shape() {
  const double a = 0.05, beta = 20.0 * 3.14159265358979323846 / 180.0;
  const body::Wedge w(a, beta);
  const double h = a * std::tan(beta);
  expect(std::abs(w.half_height() - 0.5 * h) < 1e-15, "wedge half height");
  expect(w.inside(0.0, -0.5 * h + 1e-12) && !w.inside(0.5 * a, -0.5 * h), "wedge apex and slant");
  double lo, hi, r;
  expect(w.z_extent(0.5 * a, lo, hi) && std::abs(lo - (-0.5 * h + 0.5 * a * std::tan(beta))) < 1e-15, "wedge z extent");
  expect(w.r_extent(0.0, r) && std::abs(r - 0.5 * a) < 1e-12, "wedge r extent at mid height");
  // 비스듬한 면 위의 점에서 부호 거리 0, 바깥 법선 방향으로 d 만큼 가면 거리 d
  const double rs = 0.4 * a, zs = -0.5 * h + rs * std::tan(beta);
  const auto n = w.outward_normal(rs, zs);
  expect(std::abs(w.distance(rs, zs)) < 1e-15, "wedge distance on the surface");
  expect(std::abs(w.distance(rs + 1e-3 * n[0], zs + 1e-3 * n[1]) - 1e-3) < 1e-12, "wedge distance outside");
  // 셀 고체 체적 합 = 삼각형 넓이 a h / 2 (평면)
  const Grid g(Geometry::Planar, 24, 30, 0.0, 0.08, -0.03, 0.03);
  const flow::Array2D phi = body::fluid_fraction(g, w, 0.0013);
  double solid = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) solid += (1.0 - phi(i, j)) * g.cell_volume(i, j);
  }
  expect(std::abs(solid - 0.5 * a * h) < 1e-9 * a * h, "wedge solid area");
}

// 물체와 함께 쓰는 면 플럭스 (열린 부분 유체 속도 + 막힌 부분 물체 속도) 는
// 투영된 속도장에서 모든 셀의 발산이 0 이어야 한다 (VOF 물 체적 보존의 전제).
void test_combined_flux_divergence_free() {
  const double a = 0.1;
  flow::Problem pb;
  pb.grid = Grid(Geometry::Axisymmetric, 24, 36, 0.0, 3 * a, -2 * a, 2.5 * a);
  pb.bc.r_min.kind = flow::BoundaryKind::Axis;
  const Grid g = pb.grid;
  auto bd = std::make_shared<flow::Body>();
  bd->shape = std::make_shared<body::Cylinder>(a, 0.37 * a);
  bd->motion.z = [a](double t) { return 0.913 * a - 1.5 * t; };
  bd->motion.w = [](double) { return -1.5; };
  bd->motion.a = [](double) { return 0.0; };
  bd->interpolation_order = 1;
  pb.body = bd;
  flow::Array2D f0(g.nr, g.nz);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) f0(i, j) = g.zc(j) < 0.0 ? 1.0 : 0.0;
  }
  auto mat = std::make_shared<interface::VofMaterial>(1000.0, 1e-3, 1.2, 1.8e-5);
  flow::SolverOptions opt;
  opt.convection = flow::Convection::Consistent;
  // P4 와 같이 압축성 공기 (비압축이면 벽 옆 미투영 셀의 발산 때문에 우변 평균을 빼야 해서
  // 모든 셀에 고른 발산이 남는다)
  pb.air_kappa = 1.4;
  interface::TwoPhaseSolver s(pb, mat, f0, true, opt);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; }, false);
  s.flow().set_pressure(flow::Array2D(g.nr, g.nz, 1e5));
  const double dt = 1e-4;
  s.step(dt);
  s.step(dt);
  flow::Array2D qu, qw;
  s.face_fluxes_with_body(s.flow().u(), s.flow().w(), s.time(), qu, qw);
  // 경계 보간 1차 방식에서는 모든 면을 운동량이 정한 벽 옆 셀이 압력 미지수에서 빠져 발산이 남는다
  // (그 셀의 물은 [0, φ] 자르기로 처리). 투영된 물 셀 (pressure_geometry().active, f = 1) 에서만 0 이어야 한다.
  const auto& active = s.flow().pressure_geometry().active;
  double worst = 0.0, scale = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) {
      scale = std::max(scale, std::abs(qw(i, j)));
      if (!active[static_cast<std::size_t>(i) * g.nz + j]) continue;
      if (s.vof().f()(i, j) < 1.0 - 1e-9) continue;  // 공기가 있으면 압축으로 실제 발산이 있다
      const double div = qu(i + 1, j) - qu(i, j) + qw(i, j + 1) - qw(i, j);
      worst = std::max(worst, std::abs(div));
      scale = std::max(scale, std::abs(qw(i, j)));
    }
  }
  expect(worst <= 1e-9 * scale, "combined body/fluid face flux is divergence free (" + std::to_string(worst / scale) + ")");
  expect(std::abs(s.clipped_water()) < 1e-12, "no water clipped in two steps");
}

void test_gci() {
  double gci = 0.0, p = 0.0;
  // φ = 1 + h² : 1/64², 1/32², 1/16² → p = 2, GCI = 1.25 |ε21/φ1| / 3
  const double p1 = 1.0 + 1.0 / 4096.0, p2 = 1.0 + 1.0 / 1024.0, p3 = 1.0 + 1.0 / 256.0;
  expect(verify::gci_fine(p1, p2, p3, gci, p) && std::abs(p - 2.0) < 1e-9, "GCI observed order 2");
  expect(std::abs(gci - 1.25 * (p2 - p1) / p1 / 3.0) < 1e-12, "GCI value");
  expect(!verify::gci_fine(1.0, 1.1, 1.05, gci, p), "oscillatory convergence is refused");
}

// 늘어나는 격자 칸 경계: 끝점 정확, 단조 증가, core 안 등간격 h, 이웃 크기 비 ≤ 1.05 (가장 거친 격자),
// 세분하면 core 칸이 정확히 절반, 원래 경계를 그대로 포함
void test_stretched_faces() {
  const double h = 1.0 / 32.0;
  const auto f0 = flow::stretched_faces(0.0, 1.0, 0.3, 0.6, h, 1.05, 4.0 * h, 0);
  const auto f1 = flow::stretched_faces(0.0, 1.0, 0.3, 0.6, h, 1.05, 4.0 * h, 1);
  expect(f0.front() == 0.0 && f0.back() == 1.0 && f1.front() == 0.0 && f1.back() == 1.0, "stretched endpoints exact");
  expect(f1.size() == 2 * (f0.size() - 1) + 1, "refinement doubles the cells");
  bool mono = true, core = true, ratio = true, nested = true;
  double worst = 1.0;
  for (std::size_t k = 1; k < f0.size(); ++k) {
    const double d = f0[k] - f0[k - 1];
    mono = mono && d > 0.0;
    const double mid = 0.5 * (f0[k] + f0[k - 1]);
    if (mid > 0.3 && mid < 0.6) core = core && std::abs(d - h) < 1e-12;  // core 는 [0.284375, 0.615625] 로 넓어짐
    if (k > 1) worst = std::max(worst, std::max(d / (f0[k - 1] - f0[k - 2]), (f0[k - 1] - f0[k - 2]) / d));
    nested = nested && std::abs(f1[2 * k] - f0[k]) < 1e-14;
  }
  ratio = worst <= 1.05 + 1e-9;
  for (std::size_t k = 1; k < f1.size(); ++k) {
    const double mid = 0.5 * (f1[k] + f1[k - 1]);
    // core 끝 칸은 PCHIP 기울기가 바깥 칸에 섞이므로 core 안쪽만 본다
    if (mid > 0.3 + h && mid < 0.6 - h) core = core && std::abs(f1[k] - f1[k - 1] - 0.5 * h) < 1e-12;
  }
  expect(mono, "stretched faces increase");
  expect(core, "core spacing exact on base and refined grids");
  expect(ratio, "neighbour ratio <= 1.05 on the base grid (" + std::to_string(worst) + ")");
  expect(nested, "refined grid contains the base faces");
  // 고정점: 경계가 정확히 지나고 비 ≤ 1.05 유지
  const auto fp = flow::stretched_faces(-3.0, 1.5, -0.4, 0.2, 1.0 / 64.0, 1.05, 1.0 / 16.0, 0, {-2.0, 1.25});
  bool hit_b = false, hit_t = false;
  double wp = 1.0;
  for (std::size_t k = 0; k < fp.size(); ++k) {
    hit_b = hit_b || fp[k] == -2.0;
    hit_t = hit_t || fp[k] == 1.25;
    if (k > 1) {
      const double d1 = fp[k] - fp[k - 1], d0 = fp[k - 1] - fp[k - 2];
      wp = std::max(wp, std::max(d1 / d0, d0 / d1));
    }
  }
  expect(hit_b && hit_t, "pinned faces are exact");
  expect(wp <= 1.05 + 1e-9, "neighbour ratio <= 1.05 with pins (" + std::to_string(wp) + ")");
  // 축대칭 셀 체적 합 = π (r1² − r0²)(z1 − z0) / (2π)
  const auto g = Grid::stretched(Geometry::Axisymmetric, f0, flow::stretched_faces(-1.0, 0.5, -0.2, 0.1, h, 1.05, 4 * h, 1));
  double vol = 0.0;
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) vol += g.cell_volume(i, j);
  }
  expect(std::abs(vol - 0.5 * 1.5) < 1e-12, "stretched axisymmetric cell volumes sum to the domain");
}

}  // namespace

// 07 설계 D2: 절단 셀 스위치를 켜도 정지 수중 원통 (3-G 조건, 빠른 격자, 격자면에서 벗어난 위치) 의
// 정수압 부력과 정지 상태가 그대로인지
void test_cut_cell_static_buoyancy() {
  // 원통을 격자면에서 칸의 0.37 만큼 옮겨 위·아래 면에 절단 셀이 생기게 한다 (S-5 core 칸 = a/16)
  const double shift = 0.37 / 16.0;
  const verify::TestReport off = verify::dev_3g(true, false, shift);
  const verify::TestReport on = verify::dev_3g(true, true, shift);
  expect(off.runs.size() == 1 && on.runs.size() == 1 && off.runs[0].failure.empty() && on.runs[0].failure.empty(),
         "3-G quick runs finish with and without the cut-cell switch");
  if (on.runs.empty() || !on.runs[0].failure.empty()) return;
  const auto& m = on.runs[0].metrics;
  const double err = m.at("force_relative_error"), v = m.at("max_speed_over_run");
  std::cout << "  3-G quick: force error off " << off.runs[0].metrics.at("force_relative_error") << ", on " << err
            << "; max speed off " << off.runs[0].metrics.at("max_speed_over_run") << ", on " << v << " m/s"
            << "; cut cells " << m.at("cut_cells") << ", pressure unknowns off "
            << off.runs[0].metrics.at("pressure_unknowns") << ", on " << m.at("pressure_unknowns") << "\n";
  expect(m.at("cut_cells") > 0.0 && m.at("pressure_unknowns") > off.runs[0].metrics.at("pressure_unknowns"),
         "cut-cell switch adds the thin cut cells to the pressure unknowns");
  expect(err < 1e-6, "cut-cell switch: buoyancy force error " + std::to_string(err));
  expect(v < 1e-8, "cut-cell switch: still water stays still (max speed " + std::to_string(v) + " m/s)");
}

int main() {
  test_stretched_faces();
  test_cut_cell_static_buoyancy();
  test_fluid_fraction_cylinder();
  test_wedge_shape();
  test_combined_flux_divergence_free();
  test_gci();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
