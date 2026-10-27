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

// 절단 셀 운동량 기하 (07 설계 D1): 제어체적이 영역을 채우는지, 발산 정리, 해석값, 연속성
void test_momentum_geometry() {
  auto check_grid = [](const Grid& g, const body::Shape& shape, double zc, double tol, const std::string& name) {
    const auto m = body::momentum_geometry(g, shape, zc);
    const auto phi = body::fluid_fraction(g, shape, zc);
    const bool axisym = g.geometry == Geometry::Axisymmetric;
    auto area = [&](double a, double b) { return axisym ? 0.5 * (b * b - a * a) : b - a; };
    double worst_u = 0.0, worst_w = 0.0, nz_sum = 0.0, scale = 0.0;
    for (int j = 0; j < g.nz; ++j) {  // u 제어체적은 각 줄을 채운다
      double a = 0.0, b = 0.0;
      for (int i = 0; i <= g.nr; ++i) {
        const auto& c = m.u[static_cast<std::size_t>(i) * g.nz + j];
        a += c.theta * area(c.r0, c.r1) * (c.z1 - c.z0);
      }
      for (int i = 0; i < g.nr; ++i) b += phi(i, j) * area(g.rf(i), g.rf(i + 1)) * g.dzc(j);
      worst_u = std::max(worst_u, std::abs(a - b) / b);
    }
    for (int i = 0; i < g.nr; ++i) {  // w 제어체적은 각 열을 채운다
      double a = 0.0, b = 0.0;
      for (int j = 0; j <= g.nz; ++j) {
        const auto& c = m.w[static_cast<std::size_t>(i) * (g.nz + 1) + j];
        a += c.theta * area(c.r0, c.r1) * (c.z1 - c.z0);
        nz_sum += c.nz_area;
        scale += std::abs(c.nz_area);
      }
      for (int j = 0; j < g.nz; ++j) b += phi(i, j) * area(g.rf(i), g.rf(i + 1)) * g.dzc(j);
      worst_w = std::max(worst_w, std::abs(a - b) / b);
    }
    expect(worst_u < tol && worst_w < tol, name + ": momentum volumes tile the fluid volume (u " +
                                               std::to_string(worst_u) + ", w " + std::to_string(worst_w) + ")");
    expect(scale > 0.0 && std::abs(nz_sum) < 1e-12 * scale, name + ": closed body, sum of w-cell axial areas is zero");
  };
  const body::Cylinder disk(1.0, 0.2);
  const Grid gu(Geometry::Axisymmetric, 30, 60, 0.0, 3.0, -3.0, 3.0);  // Δ = 0.1
  check_grid(gu, disk, 0.0137, 1e-12, "disk, uniform grid");
  // 구는 극의 제곱근 특이성 때문에 가우스 적분 오차가 남는다 (구간 나눔이 셀과 달라 약 4e−5)
  check_grid(gu, body::Spheroid(1.0, 1.0), 0.0137, 1e-4, "sphere, uniform grid");
  std::vector<double> rf, zf;
  for (int i = 0; i <= 40; ++i) rf.push_back(3.0 * std::pow(i / 40.0, 1.3));
  for (int j = 0; j <= 70; ++j) zf.push_back(-3.0 + 6.0 * (j / 70.0) + 0.2 * std::sin(kPi * j / 70.0));
  check_grid(Grid::stretched(Geometry::Axisymmetric, rf, zf), disk, -0.0211, 1e-12, "disk, stretched grid");

  // 해석값: 원판 (반경 1, 두께 0.2, 중심 zc) 의 옆벽을 가로지르는 u 제어체적
  const double zc = 0.0137;
  const auto m = body::momentum_geometry(gu, disk, zc);
  int checked = 0;
  for (int i = 1; i < gu.nr; ++i) {
    for (int j = 0; j < gu.nz; ++j) {
      const auto& c = m.u[static_cast<std::size_t>(i) * gu.nz + j];
      if (!(c.r0 < 1.0 && c.r1 > 1.0)) continue;
      const double h = std::max(0.0, std::min(c.z1, zc + 0.1) - std::max(c.z0, zc - 0.1));
      if (h <= 0.0) continue;
      const double v = 0.5 * (c.r1 * c.r1 - c.r0 * c.r0) * (c.z1 - c.z0);
      const double theta = 1.0 - 0.5 * (1.0 - c.r0 * c.r0) * h / v;
      // 옆벽 (r = 1) 의 법선은 유체 → 물체 = −r 방향, r 가중 넓이 1 · h
      expect(std::abs(c.theta - theta) < 1e-12 && std::abs(c.nr_area + h) < 1e-12,
             "disk side wall u-cell: theta and n_r area match the exact values");
      ++checked;
    }
  }
  expect(checked > 0, "disk side wall u-cells were found");

  // 연속성: 물체를 칸의 1e-4 만큼 옮기면 θ 변화도 그 정도
  const auto m2 = body::momentum_geometry(gu, disk, zc + 1e-5);
  double jump = 0.0;
  for (std::size_t k = 0; k < m.w.size(); ++k) jump = std::max(jump, std::abs(m.w[k].theta - m2.w[k].theta));
  for (std::size_t k = 0; k < m.u.size(); ++k) jump = std::max(jump, std::abs(m.u[k].theta - m2.u[k].theta));
  expect(jump > 0.0 && jump < 3e-4, "momentum volume fractions change continuously with the body position (" +
                                        std::to_string(jump) + ")");
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

// Lamb Art. 337 Stokes 해: 운동량 잔차, 그리고 격자에 넣었을 때 표면 힘 (압력 2πμaU + 점성 4πμaU)
void test_stokes_exact() {
  const double a = 1.0, U = 1.0, mu = 1.0;
  auto vel = [&](double r, double z) { return verify::stokes_sphere_velocity(a, U, r, z); };
  auto pre = [&](double r, double z) { return verify::stokes_sphere_pressure(a, U, mu, r, z); };
  const double r = 0.8, z = 1.1, h = 1e-3;
  auto lap = [&](int c, double rr, double zz) {
    const double f0 = vel(rr, zz)[c];
    const double frr = (vel(rr + h, zz)[c] - 2 * f0 + vel(rr - h, zz)[c]) / (h * h);
    const double fr = (vel(rr + h, zz)[c] - vel(rr - h, zz)[c]) / (2 * h);
    const double fzz = (vel(rr, zz + h)[c] - 2 * f0 + vel(rr, zz - h)[c]) / (h * h);
    return frr + fr / rr + fzz - (c == 0 ? f0 / (rr * rr) : 0.0);
  };
  const double pr = (pre(r + h, z) - pre(r - h, z)) / (2 * h);
  const double pz = (pre(r, z + h) - pre(r, z - h)) / (2 * h);
  expect(std::abs(pr - mu * lap(0, r, z)) < 1e-5 && std::abs(pz - mu * lap(1, r, z)) < 1e-5,
         "Stokes sphere: grad p = mu laplacian u");
  const double div = ((r + h) * vel(r + h, z)[0] - (r - h) * vel(r - h, z)[0]) / (2 * h * r) +
                     (vel(r, z + h)[1] - vel(r, z - h)[1]) / (2 * h);
  // 중심 차분 절단 오차 O(h²) ≈ 3.5e−7 (h = 1e−3) 만 남는다 (h 를 줄이면 h² 로 줄어듦을 확인)
  expect(std::abs(div) < 1e-6, "Stokes sphere: divergence-free");
  expect(std::abs(vel(a * std::sin(0.7), a * std::cos(0.7))[1]) < 1e-14, "Stokes sphere: no-slip");

  // 정확한 장을 격자에 넣고 표면 힘 적분 (a/Δ = 16)
  const Grid g(Geometry::Axisymmetric, 64, 128, 0.0, 4.0, -4.0, 4.0);
  body::Spheroid sphere(a, a);
  const auto geom = body::cut_geometry(g, sphere, 0.0);
  const auto mc = body::momentum_classes(g, sphere, 0.0);
  flow::Array2D u(g.nr + 1, g.nz), w(g.nr, g.nz + 1), p(g.nr, g.nz);
  for (int i = 0; i <= g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) u(i, j) = vel(g.rf(i), g.zc(j))[0];
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j <= g.nz; ++j) w(i, j) = vel(g.rc(i), g.zf(j))[1];
  }
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) p(i, j) = pre(g.rc(i), g.zc(j));
  }
  const double fp = body::axial_pressure_force(g, geom, p);
  const double fv = body::axial_viscous_force(g, geom, sphere, 0.0, w, mu);
  expect(std::abs(fp / (2 * kPi * mu * a * U) - 1.0) < 0.03,
         "Stokes sphere: pressure drag 2 pi mu a U (" + std::to_string(fp / (2 * kPi)) + ")");
  expect(std::abs(fv / (4 * kPi * mu * a * U) - 1.0) < 0.03,
         "Stokes sphere: viscous drag 4 pi mu a U (" + std::to_string(fv / (4 * kPi)) + ")");
  expect(!mc.forcing.empty(), "faces next to the sphere are forcing faces");
  // 제어면 운동량 수지 (Stokes 장이므로 대류 항 없이 ρ = 0): 6πμaU
  const double fcs = body::axial_force_control_surface(g, u, w, p, mu, 0.0, 2.5, -2.5, 2.5);
  expect(std::abs(fcs / (6 * kPi * mu * a * U) - 1.0) < 2e-3,
         "Stokes sphere: control-surface drag 6 pi mu a U (" + std::to_string(fcs / (6 * kPi)) + ")");
  bool threw = false;
  try {
    body::axial_force_control_surface(g, u, w, p, mu, 0.0, 2.51, -2.5, 2.5);
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  expect(threw, "control surface off the grid faces is refused");
}

// 반암시적 진행(step_implicit)의 정상 해 = 검증된 명시적 진행(step)의 정상 해 (물체 없음)
void test_implicit_matches_explicit() {
  auto make = [] {
    flow::Problem pb;
    pb.grid = Grid(Geometry::Axisymmetric, 8, 8, 0.0, 1.0, 0.0, 1.0);
    pb.bc.r_min.kind = flow::BoundaryKind::Axis;
    pb.bc.z_max.velocity = [](double r, double, double) {
      return std::array<double, 2>{r * (1.0 - r), 0.0};  // 윗벽이 반경 방향으로 끌림
    };
    pb.material = std::make_shared<flow::SmoothMaterial>(
        [](double, double, double) { return 1.0; }, [](double, double, double) { return 1.0; }, false);
    return pb;
  };
  auto zero = [](double, double) { return 0.0; };
  flow::FlowSolver ex(make()), im(make());
  ex.initialize(0.0, zero, zero);
  im.initialize(0.0, zero, zero);
  const double dx = 1.0 / 8.0;
  for (int k = 0; k < 4000; ++k) ex.step(0.2 * dx * dx);
  int steps = 0;
  do {
    im.step_implicit(0.1);
  } while (im.last_change() > 1e-13 && ++steps < 2000);
  double diff = 0.0, umax = 0.0;
  for (std::size_t q = 0; q < ex.u().data().size(); ++q) {
    diff = std::max(diff, std::abs(ex.u().data()[q] - im.u().data()[q]));
    umax = std::max(umax, std::abs(ex.u().data()[q]));
  }
  for (std::size_t q = 0; q < ex.w().data().size(); ++q) {
    diff = std::max(diff, std::abs(ex.w().data()[q] - im.w().data()[q]));
  }
  expect(steps < 2000, "implicit stepping reaches a steady state");
  expect(umax > 1e-3 && diff < 1e-9 * std::max(1.0, umax),
         "implicit and explicit steady states agree (diff " + std::to_string(diff) + ")");
  expect(!im.stats().solver_failed && !ex.stats().solver_failed, "solvers converged in both runs");
}

// 회귀 시험: 움직이는 구(저점성)가 격자 몇 칸을 지나는 동안 속도가 유한하게 머문다
// (2026-10-22: 2차 경계 보간에서 지수적으로 커졌음 → 이 계산은 1차 보간을 쓴다)
void test_moving_body_stable() {
  const int n = 8;
  auto pot = std::make_shared<verify::SpheroidPotential>(1.0, 1.0);
  auto b = std::make_shared<flow::Body>();
  b->shape = std::make_shared<body::Spheroid>(1.0, 1.0);
  b->motion.z = [](double t) { return 0.5 * t * t; };
  b->motion.w = [](double t) { return t; };
  b->motion.a = [](double) { return 1.0; };
  b->interpolation_order = 1;
  flow::Problem pb;
  pb.grid = Grid(Geometry::Axisymmetric, 3 * n, 6 * n, 0.0, 3.0, -3.0, 3.0);
  pb.body = b;
  pb.bc.r_min.kind = flow::BoundaryKind::Axis;
  flow::WallVelocity v = [pot](double r, double z, double t) {
    const auto e = pot->velocity(r, z - 0.5 * t * t);
    return std::array<double, 2>{t * e[0], t * e[1]};
  };
  pb.bc.r_max.velocity = pb.bc.z_min.velocity = pb.bc.z_max.velocity = v;
  pb.material = std::make_shared<flow::SmoothMaterial>(
      [](double, double, double) { return 1.0; }, [](double, double, double) { return 0.01; }, false);
  flow::FlowSolver s(pb);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
  const double dx = 1.0 / n, t_end = std::sqrt(6.0 * dx);
  double worst = 0.0;
  while (s.time() < t_end - 1e-12) {
    s.step_implicit(std::min(0.05 * dx, t_end - s.time()));
    double umax = 0.0;
    for (double x : s.w().data()) umax = std::max(umax, std::abs(x));
    worst = std::max(worst, umax / std::max(s.time(), 1e-12));  // 물체 속도 = t
  }
  expect(std::isfinite(worst) && worst < 3.0, "moving sphere stays bounded (max |w| / body speed = " +
                                                  std::to_string(worst) + ")");
  expect(!s.stats().solver_failed, "moving sphere: solvers converged");
}

}  // namespace

int main() {
  test_shapes();
  test_cut_geometry();
  test_momentum_geometry();
  test_potential();
  test_body_pressure();
  test_stokes_exact();
  test_implicit_matches_explicit();
  test_moving_body_stable();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
