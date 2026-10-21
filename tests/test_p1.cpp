// P1 단위 시험: Jet 도함수, 켤레기울기 해법, 투영 후 발산, 경계 처리.

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "flow/poisson.h"
#include "flow/solver.h"
#include "verify/jet.h"

using namespace hcw::flow;
using hcw::verify::Jet;

namespace {

int g_failures = 0, g_checks = 0;

void expect(bool ok, const std::string& what) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::cerr << "FAIL: " << what << "\n";
  }
}

void test_jet() {
  // f(r, z) = sin(r z) cos(r) + r² z, 도함수를 중심 차분과 비교
  auto f = [](double r, double z) { return std::sin(r * z) * std::cos(r) + r * r * z; };
  auto fj = [](double r, double z) {
    const Jet R = Jet::var_r(r), Z = Jet::var_z(z);
    return sin(R * Z) * cos(R) + R * R * Z;
  };
  const double r = 0.7, z = -0.4, h = 1e-4;
  const Jet j = fj(r, z);
  const double fr = (f(r + h, z) - f(r - h, z)) / (2 * h);
  const double fz = (f(r, z + h) - f(r, z - h)) / (2 * h);
  const double frr = (f(r + h, z) - 2 * f(r, z) + f(r - h, z)) / (h * h);
  const double fzz = (f(r, z + h) - 2 * f(r, z) + f(r, z - h)) / (h * h);
  const double frz = (f(r + h, z + h) - f(r + h, z - h) - f(r - h, z + h) + f(r - h, z - h)) /
                     (4 * h * h);
  expect(std::abs(j.v - f(r, z)) < 1e-14, "jet value");
  expect(std::abs(j.r - fr) < 1e-7, "jet d/dr");
  expect(std::abs(j.z - fz) < 1e-7, "jet d/dz");
  expect(std::abs(j.rr - frr) < 1e-5, "jet d2/dr2");
  expect(std::abs(j.zz - fzz) < 1e-5, "jet d2/dz2");
  expect(std::abs(j.rz - frz) < 1e-5, "jet d2/drdz");
}

void test_pcg() {
  // 1차원 Dirichlet 라플라시안 (SPD), 해 x_i = sin(i)
  const int n = 200;
  CsrMatrix a;
  a.n = n;
  a.row_start.push_back(0);
  for (int i = 0; i < n; ++i) {
    if (i > 0) { a.col.push_back(i - 1); a.val.push_back(-1.0); }
    a.col.push_back(i);
    a.val.push_back(2.0);
    if (i < n - 1) { a.col.push_back(i + 1); a.val.push_back(-1.0); }
    a.row_start.push_back(static_cast<int>(a.col.size()));
  }
  std::vector<double> x_true(n), b, x;
  for (int i = 0; i < n; ++i) x_true[i] = std::sin(static_cast<double>(i));
  a.multiply(x_true, b);
  IncompleteCholesky ic;
  ic.factor(a);
  const SolveResult res = pcg(a, ic, b, x, 1e-13, 10000);
  double err = 0.0;
  for (int i = 0; i < n; ++i) err = std::max(err, std::abs(x[i] - x_true[i]));
  expect(res.converged, "pcg converges");
  expect(err < 1e-8, "pcg solution accuracy");
  // 삼중대각에서 IC(0) 는 정확한 Cholesky → 한 번에 수렴
  expect(res.iterations <= 2, "IC(0) is exact for tridiagonal matrices");
}

// 임의의 매끈한 속도장을 투영한 뒤 이산 발산이 반올림 수준인지
void test_projection(Geometry geom, BoundaryKind rk, BoundaryKind zk, const std::string& label) {
  Problem pb;
  const double r0 = 0.0;
  pb.grid = Grid(geom, 24, 20, r0, 1.0, 0.0, 1.3);
  pb.bc.r_min.kind = (geom == Geometry::Axisymmetric) ? BoundaryKind::Axis : rk;
  pb.bc.r_max.kind = rk;
  pb.bc.z_min.kind = pb.bc.z_max.kind = zk;
  if (rk == BoundaryKind::Wall) {
    pb.bc.r_max.velocity = [](double, double z, double) {
      return std::array<double, 2>{0.0, std::sin(3 * z)};
    };
  }
  pb.material = std::make_shared<SmoothMaterial>(
      [](double r, double z, double) { return 2.0 + std::sin(5 * r) * std::cos(4 * z); },
      [](double, double, double) { return 0.01; }, false);
  FlowSolver s(pb);
  s.initialize(0.0, [](double r, double z) { return std::cos(2 * r + z); },
               [](double r, double z) { return std::sin(r - 3 * z); });
  expect(s.max_divergence() < 1e-9, "projection makes velocity divergence-free: " + label);
  expect(!s.stats().solver_failed, "projection solver converged: " + label);
  s.step(1e-3);
  expect(s.max_divergence() < 1e-9, "divergence-free after one step: " + label);
}

// 회귀 시험: 밀도비가 큰 두 층에서 압력 해법이 무너지지 않아야 한다
// (정칙화 위치를 첫 셀에 두었을 때 밀도비 1000 에서 발산했던 문제).
void test_large_density_ratio() {
  for (double ratio : {1e3, 1e6}) {
    Problem pb;
    pb.grid = Grid(Geometry::Planar, 8, 8, 0.0, 1.0, 0.0, 1.0);
    pb.bc.r_min.kind = pb.bc.r_max.kind = BoundaryKind::Wall;
    pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Wall;
    pb.material = std::make_shared<LayeredMaterial>(0.5, ratio, 1e-3, 1.0, 1.8e-5);
    pb.gravity = {0.0, -9.81};
    FlowSolver s(pb);
    s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
    s.step(1e-3);
    double umax = 0.0;
    for (double v : s.w().data()) umax = std::max(umax, std::abs(v));
    expect(!s.stats().solver_failed, "pressure solver converges at density ratio " +
                                         std::to_string(ratio));
    expect(umax < 1e-12, "rest state stays at rest at density ratio " + std::to_string(ratio));
  }
}

void test_boundary_values() {
  Problem pb;
  pb.grid = Grid(Geometry::Axisymmetric, 8, 6, 0.0, 1.0, 0.0, 1.0);
  pb.bc.r_min.kind = BoundaryKind::Axis;
  pb.bc.r_max.kind = BoundaryKind::Wall;
  pb.bc.r_max.velocity = [](double, double, double) { return std::array<double, 2>{0.0, 0.7}; };
  pb.bc.z_min.kind = pb.bc.z_max.kind = BoundaryKind::Wall;
  pb.material = std::make_shared<SmoothMaterial>([](double, double, double) { return 1.0; },
                                                 [](double, double, double) { return 1.0; },
                                                 false);
  FlowSolver s(pb);
  s.initialize(0.0, [](double, double) { return 0.3; }, [](double, double) { return 0.2; });
  bool axis_zero = true, wall_zero = true, bottom_zero = true;
  for (int j = 0; j < 6; ++j) {
    axis_zero = axis_zero && s.u()(0, j) == 0.0;
    wall_zero = wall_zero && s.u()(8, j) == 0.0;
  }
  for (int i = 0; i < 8; ++i) bottom_zero = bottom_zero && s.w()(i, 0) == 0.0 && s.w()(i, 6) == 0.0;
  expect(axis_zero, "axis radial velocity is zero");
  expect(wall_zero, "wall normal velocity imposed");
  expect(bottom_zero, "bottom/top normal velocity imposed");
  expect(!s.solves_u(0) && !s.solves_u(8) && s.solves_u(1), "solved r-face range");
}

void test_invalid_setups() {
  auto mat = std::make_shared<SmoothMaterial>([](double, double, double) { return 1.0; },
                                              [](double, double, double) { return 1.0; }, false);
  auto throws = [&](Problem pb, const std::string& label) {
    pb.material = mat;
    bool threw = false;
    try {
      FlowSolver s(pb);
    } catch (const std::invalid_argument&) {
      threw = true;
    }
    expect(threw, "rejects: " + label);
  };
  Problem a;
  a.grid = Grid(Geometry::Axisymmetric, 4, 4, 0.0, 1.0, 0.0, 1.0);
  a.bc.r_min.kind = BoundaryKind::Wall;
  throws(a, "axisymmetric r=0 without axis boundary");
  Problem b;
  b.grid = Grid(Geometry::Planar, 4, 4, 0.0, 1.0, 0.0, 1.0);
  b.bc.r_min.kind = BoundaryKind::Periodic;
  throws(b, "unpaired periodic boundary");
  Problem c;
  c.grid = Grid(Geometry::Planar, 4, 4, 0.0, 1.0, 0.0, 1.0);
  c.bc.r_min.kind = BoundaryKind::Axis;
  throws(c, "axis in planar geometry");
}

}  // namespace

int main() {
  test_jet();
  test_pcg();
  test_projection(Geometry::Planar, BoundaryKind::Wall, BoundaryKind::Wall, "planar walls");
  test_projection(Geometry::Planar, BoundaryKind::Periodic, BoundaryKind::Periodic,
                  "planar periodic");
  test_projection(Geometry::Axisymmetric, BoundaryKind::Wall, BoundaryKind::Wall,
                  "axisymmetric walls");
  test_projection(Geometry::Axisymmetric, BoundaryKind::Wall, BoundaryKind::Periodic,
                  "axisymmetric periodic z");
  test_large_density_ratio();
  test_boundary_values();
  test_invalid_setups();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
