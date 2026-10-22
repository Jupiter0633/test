// P2 단위 시험: 형상, 컷셀 기하, 퍼텐셜 해, 물체가 있는 압력 해법과 표면 힘.

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

#include "body/cut_cell.h"
#include "body/force.h"
#include "body/shape.h"
#include "flow/solver.h"
#include "verify/potential.h"

using namespace hcw;
using flow::Geometry;
using flow::Grid;

namespace {

int g_failures = 0, g_checks = 0;
constexpr double kPi = 3.14159265358979323846;

void expect(bool ok, const std::string& what) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::cerr << "FAIL: " << what << "\n";
  }
}

void test_shapes() {
  body::Spheroid s(1.0, 0.5);
  double lo, hi, radius;
  expect(s.z_extent(0.6, lo, hi) && std::abs(hi - 0.4) < 1e-15 && lo == -hi, "spheroid z extent");
  expect(s.r_extent(0.3, radius) && std::abs(radius - 0.8) < 1e-15, "spheroid r extent");
  expect(!s.z_extent(1.01, lo, hi) && !s.r_extent(0.51, radius), "spheroid misses outside lines");
  body::Cylinder c(1.0, 0.1);
  expect(c.z_extent(1.0, lo, hi) && lo == -0.05 && hi == 0.05, "cylinder side is part of the body");
  expect(c.inside(0.5, 0.0) && !c.inside(0.5, 0.06), "cylinder inside test");
}

void test_cut_geometry() {
  // 구: 닫힌 물체이므로 Σ nz_area = 0, 격자면에 맞지 않음
  const Grid g(Geometry::Axisymmetric, 24, 48, 0.0, 3.0, -3.0, 3.0);
  const auto sphere = body::cut_geometry(g, body::Spheroid(1.0, 1.0), 0.1);
  double sum = 0.0, top = 0.0;
  for (const auto& p : sphere.pieces) sum += p.nz_area;
  expect(std::abs(sum) < 1e-13, "closed body: sum of axial boundary areas is zero");
  expect(!sphere.grid_fitted, "sphere is not grid-fitted");
  int solid = 0;
  for (char a : sphere.active) solid += a ? 0 : 1;
  expect(solid > 0, "cells inside the sphere are inactive");

  // 격자면에 맞춘 원판: 열린 비율이 0 또는 1, 위쪽 면적 합 = ∫_0^a r dr = a²/2
  const Grid gd(Geometry::Axisymmetric, 30, 60, 0.0, 3.0, -3.0, 3.0);  // Δ = 0.1
  const auto disk = body::cut_geometry(gd, body::Cylinder(1.0, 0.2), 0.0);
  expect(disk.grid_fitted, "aligned disk is grid-fitted");
  sum = 0.0;
  for (const auto& p : disk.pieces) {
    sum += p.nz_area;
    if (p.nz_area < 0) top += p.nz_area;  // 물체 위쪽(법선 −z) 조각
  }
  expect(std::abs(sum) < 1e-13, "aligned disk: sum of axial boundary areas is zero");
  expect(std::abs(-top - 0.5) < 1e-12, "aligned disk: top area equals a^2/2 (2pi omitted)");

  // 물체가 없으면 모든 면이 열림
  const auto open = body::open_geometry(g);
  expect(open.pieces.empty() && open.grid_fitted && open.open_u(3, 4) == 1.0, "open geometry");
}

void test_potential() {
  verify::SpheroidPotential sphere(1.0, 1.0);
  expect(std::abs(sphere.added_mass_coefficient() - 0.5) < 1e-14, "sphere added mass coefficient 1/2");
  expect(std::abs(sphere.added_mass(1.0) - 2.0 * kPi / 3.0) < 1e-13, "sphere added mass (2/3) pi");
  // 축 위 이중 원천: u_z = a³/z³
  expect(std::abs(sphere.velocity(0.0, 2.0)[1] - 0.125) < 1e-13, "sphere axis velocity a^3/z^3");

  verify::SpheroidPotential oblate(1.0, 0.2);
  // 같은 Lamb Art. 114 공식을 독립 수치 적분한 값 (2026-10-21, P2 기준안 작성 시): k = 3.007758
  expect(std::abs(oblate.added_mass_coefficient() - 3.007758) < 2e-6, "oblate spheroid k = 3.007758");
  verify::SpheroidPotential disk(1.0, 0.0);
  expect(std::abs(disk.added_mass(1.0) - 8.0 / 3.0) < 1e-13, "thin disk added mass (8/3) rho a^3");
  // 아주 얇은 회전타원체는 원판 값에 가까워야 한다
  verify::SpheroidPotential thin(1.0, 1e-4);
  expect(std::abs(thin.added_mass(1.0) / (8.0 / 3.0) - 1.0) < 1e-3, "thin spheroid tends to the disk");

  // 표면 바로 바깥에서 법선 속도 = 물체 속도의 법선 성분 (U = 1, 부호 확인)
  for (double theta : {0.1, 0.5, 1.0, 1.4}) {
    const double a = 1.0, c = 0.2, s = 1.0 + 1e-7;
    const double r = s * a * std::sin(theta), z = s * c * std::cos(theta);
    double nr = r / (a * a), nz = z / (c * c);
    const double len = std::hypot(nr, nz);
    nr /= len;
    nz /= len;
    const auto v = oblate.velocity(r, z);
    expect(std::abs(v[0] * nr + v[1] * nz - nz) < 1e-5,
           "oblate: normal velocity on the surface, theta=" + std::to_string(theta));
  }
  // 축대칭 발산 (1/r) ∂(r u_r)/∂r + ∂u_z/∂z = 0
  const double r = 0.7, z = 0.45, h = 1e-5;
  const double div = ((r + h) * oblate.velocity(r + h, z)[0] - (r - h) * oblate.velocity(r - h, z)[0]) /
                         (2 * h * r) +
                     (oblate.velocity(r, z + h)[1] - oblate.velocity(r, z - h)[1]) / (2 * h);
  expect(std::abs(div) < 1e-7, "oblate potential flow is divergence-free");
}

// 시각 0 의 힘. outer: 바깥 경계 퍼텐셜 (nullptr 이면 정지 벽), body_fixed: 물체 좌표계
double initial_force(const Grid& g, std::shared_ptr<const body::Shape> shape,
                     std::shared_ptr<const verify::SpheroidPotential> outer, bool body_fixed,
                     bool* failed = nullptr) {
  auto b = std::make_shared<flow::Body>();
  b->shape = std::move(shape);
  if (!body_fixed) {
    b->motion.z = [](double t) { return 0.5 * t * t; };
    b->motion.w = [](double t) { return t; };
    b->motion.a = [](double) { return 1.0; };
  }
  flow::Problem pb;
  pb.grid = g;
  pb.body = b;
  pb.bc.r_min.kind = flow::BoundaryKind::Axis;
  if (outer) {
    flow::WallVelocity v = [outer](double rr, double zz, double t) {
      const auto e = outer->velocity(rr, zz - 0.5 * t * t);
      return std::array<double, 2>{t * e[0], t * e[1]};
    };
    pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
  }
  if (body_fixed) {
    flow::WallVelocity v = [](double, double, double t) { return std::array<double, 2>{0.0, -t}; };
    pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
    pb.gravity = {0.0, -1.0};
  }
  pb.material = std::make_shared<flow::SmoothMaterial>(
      [](double, double, double) { return 1.0; }, [](double, double, double) { return 0.0; }, false);
  flow::FlowSolver s(pb);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
  flow::Array2D p;
  s.compute_pressure(p);
  if (failed) *failed = s.stats().solver_failed;
  return body::axial_pressure_force(s.grid(), s.geometry(), p);
}

void test_body_pressure() {
  // 구, 거친 격자 (a/Δ = 8): 부가질량 힘 −(2/3)π 에 가까워야 한다
  const Grid g(Geometry::Axisymmetric, 24, 48, 0.0, 3.0, -3.0, 3.0);
  auto pot = std::make_shared<verify::SpheroidPotential>(1.0, 1.0);
  bool failed = true;
  const double f = initial_force(g, std::make_shared<body::Spheroid>(1.0, 1.0), pot, false, &failed);
  const double exact = -2.0 * kPi / 3.0;
  expect(!failed, "pressure solve with a sphere converges");
  expect(std::abs(f - exact) / std::abs(exact) < 0.1, "sphere force at t=0 within 10% on a coarse grid (" +
                                                        std::to_string(f / exact) + ")");

  // 상수 압력은 닫힌 물체에 힘을 주지 않는다
  const auto geom = body::cut_geometry(g, body::Spheroid(1.0, 1.0), 0.0);
  flow::Array2D pc(g.nr, g.nz, 3.7);
  expect(std::abs(body::axial_pressure_force(g, geom, pc)) < 1e-12, "constant pressure gives no force");

  // 격자면에 맞춘 원판: 실험실(IBM) 모드와 물체 좌표계 모드의 이산식이 같아 힘이 반올림 수준까지 같다
  const Grid gd(Geometry::Axisymmetric, 30, 60, 0.0, 3.0, -3.0, 3.0);
  auto disk = std::make_shared<body::Cylinder>(1.0, 0.2);
  const double f_lab = initial_force(gd, disk, nullptr, false);
  const double f_bf = initial_force(gd, disk, nullptr, true);
  expect(std::abs(f_lab - f_bf) <= 1e-10 * std::abs(f_bf), "grid-fitted disk: lab and body frames agree");

  // 물체가 있으면 시간 진행은 아직 거부한다
  bool threw = false;
  try {
    flow::Problem pb;
    pb.grid = g;
    auto b = std::make_shared<flow::Body>();
    b->shape = std::make_shared<body::Spheroid>(1.0, 1.0);
    pb.body = b;
    pb.bc.r_min.kind = flow::BoundaryKind::Axis;
    pb.material = std::make_shared<flow::SmoothMaterial>(
        [](double, double, double) { return 1.0; }, [](double, double, double) { return 0.0; }, false);
    flow::FlowSolver s(pb);
    s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
    s.step(1e-3);
  } catch (const std::logic_error&) {
    threw = true;
  }
  expect(threw, "time stepping with a body is refused until near-body momentum exists");
}

}  // namespace

int main() {
  test_shapes();
  test_cut_geometry();
  test_potential();
  test_body_pressure();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
