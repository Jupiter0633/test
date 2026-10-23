// P1 단위 시험: Jet 도함수, 켤레기울기 해법, 보정 해법, 투영 후 발산, 경계 처리,
// 다중 스레드 결과가 1 스레드와 비트 단위로 같은지.

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <string>
#include <vector>

#include "common/parallel.h"
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

  // 도달 불가능한 허용치: 최대 반복까지 헛돌지 않고 정체로 멈춰야 한다
  std::vector<double> y;
  const SolveResult stalled = pcg(a, ic, b, y, 1e-30, 100000);
  expect(!stalled.converged, "unreachable tolerance is reported as not converged");
  expect(stalled.iterations < 2000, "stagnation stops the iteration early");
}

// 회귀 시험: 고운 격자에서 배정밀도 pcg 는 반올림 한계(≈ ε·|A||x|/|b|)에 막혀
// 1e-13 에 닿지 못하고 최대 반복까지 돈다 (P1 1-A(a) N≥64 의 출력 압력, 2026-10-22).
// 보정 해법은 허용치에 닿아야 하고, 닿을 수 없는 허용치에서는 빨리 실패해야 한다.
void test_refined_solve() {
  // 2차원 Dirichlet 라플라시안 128×128, 매끈한 해 → |A||x| / |b| ≈ N² 수준
  const int m = 128, n = m * m;
  CsrMatrix a;
  a.n = n;
  a.row_start.push_back(0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      const int k = i * m + j;
      if (i > 0) { a.col.push_back(k - m); a.val.push_back(-1.0); }
      if (j > 0) { a.col.push_back(k - 1); a.val.push_back(-1.0); }
      a.col.push_back(k);
      a.val.push_back(4.0);
      if (j < m - 1) { a.col.push_back(k + 1); a.val.push_back(-1.0); }
      if (i < m - 1) { a.col.push_back(k + m); a.val.push_back(-1.0); }
      a.row_start.push_back(static_cast<int>(a.col.size()));
    }
  }
  const double pi = 3.14159265358979323846;
  std::vector<double> x_true(n), b;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      x_true[i * m + j] = std::sin(pi * (i + 1) / (m + 1)) * std::sin(pi * (j + 1) / (m + 1)) + 1.0;
    }
  }
  a.multiply(x_true, b);
  IncompleteCholesky ic;
  ic.factor(a);

  std::vector<double> hi, lo;
  const SolveResult res = refined_solve(a, ic, b, hi, lo, 1e-13, 20000);
  expect(res.converged && res.relative_residual <= 1e-13, "refined solve reaches 1e-13 on 128x128");
  expect(res.iterations < 2000, "refined solve needs few iterations on 128x128");

  // 닿을 수 없는 허용치: 최대 반복까지 가지 않고 정체로 끝나야 한다
  hi.clear();
  lo.clear();
  const SolveResult stuck = refined_solve(a, ic, b, hi, lo, 1e-40, 1000000);
  expect(!stuck.converged && stuck.stagnated, "unreachable tolerance ends as stagnation");
  expect(stuck.iterations < 20000, "stagnation is detected early");

  // 따뜻한 시작: 이미 푼 해에서 다시 풀면 바로 끝난다
  const SolveResult again = refined_solve(a, ic, b, hi, lo, 1e-13, 20000);
  expect(again.converged && again.iterations == 0, "warm start from a converged solution");
}

// 격자 nr×nz (셀 c = i·nz + j) 의 변수 계수 압력형 행렬. periodic 이면 두 방향 모두 감는다.
// skew 이면 (i, j)–(i−1, j+1) 연결을 더해 파이프라인 순서로 나눌 수 없는 행렬을 만든다.
CsrMatrix grid_matrix(int nr, int nz, bool periodic, bool skew) {
  const int n = nr * nz;
  std::vector<std::vector<std::pair<int, double>>> rows(static_cast<std::size_t>(n));
  auto add = [&](int a, int b, double v) {
    for (auto& e : rows[a]) {
      if (e.first == b) {
        e.second += v;
        return;
      }
    }
    rows[a].push_back({b, v});
  };
  auto couple = [&](int a, int b, double coef) {
    add(a, a, coef);
    add(b, b, coef);
    add(a, b, -coef);
    add(b, a, -coef);
  };
  auto cell = [&](int i, int j) { return ((i + nr) % nr) * nz + (j + nz) % nz; };
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const double k = 1.0 + 0.9 * std::sin(0.37 * i + 0.11 * j * j);  // 밀도비 흉내
      if (i > 0 || periodic) couple(cell(i - 1, j), cell(i, j), k);
      if (j > 0 || periodic) couple(cell(i, j - 1), cell(i, j), 1.0 / (k + 0.05));
      if (skew && i > 0 && j + 1 < nz) couple(cell(i - 1, j + 1), cell(i, j), 0.1);
    }
  }
  add(n - 1, n - 1, rows[n - 1][0].second);  // 정칙화 (해석기와 같은 방식: 마지막 대각 두 배)
  CsrMatrix a;
  a.n = n;
  a.row_start.push_back(0);
  for (auto& r : rows) {
    std::sort(r.begin(), r.end());
    for (const auto& e : r) {
      a.col.push_back(e.first);
      a.val.push_back(e.second);
    }
    a.row_start.push_back(static_cast<int>(a.col.size()));
  }
  return a;
}

bool same_bits(const std::vector<double>& x, const std::vector<double>& y) {
  return x.size() == y.size() && std::memcmp(x.data(), y.data(), x.size() * sizeof(double)) == 0;
}

// 다중 스레드: 전처리·pcg·보정 해법이 스레드 수와 상관없이 1 스레드와 비트 단위로 같아야 한다.
void test_threads_bitwise() {
  struct Case {
    int nr, nz;
    bool periodic, skew, pipelined;
    const char* label;
  };
  const Case cases[] = {{70, 96, false, false, true, "walls"},
                        {64, 80, true, false, true, "periodic"},
                        {64, 80, false, true, false, "skew (pipeline refused)"}};
  for (const Case& c : cases) {
    const CsrMatrix a = grid_matrix(c.nr, c.nz, c.periodic, c.skew);
    std::vector<double> b(static_cast<std::size_t>(a.n));
    for (int k = 0; k < a.n; ++k) b[k] = std::sin(0.013 * k) + 0.3 * std::cos(0.7 * k);
    double mean = 0.0;
    for (double v : b) mean += v;
    for (double& v : b) v -= mean / a.n;

    hcw::parallel::set_threads(1);
    IncompleteCholesky ic1;
    ic1.factor(a, c.nz);
    std::vector<double> z1, x1, hi1, lo1;
    ic1.apply(b, z1);
    const SolveResult p1 = pcg(a, ic1, b, x1, 1e-10, 5000);
    const SolveResult r1 = refined_solve(a, ic1, b, hi1, lo1, 1e-13, 20000);
    expect(p1.converged && r1.converged, std::string("threads reference converges: ") + c.label);

    for (int t : {2, 3, 4, 7}) {
      hcw::parallel::set_threads(t);
      IncompleteCholesky ic;
      ic.factor(a, c.nz);
      const std::string tag = std::string(c.label) + ", threads " + std::to_string(t);
      expect((ic.pipeline_segments() > 0) == c.pipelined, "pipeline chosen as expected: " + tag);
      std::vector<double> z, x, hi, lo;
      ic.apply(b, z);
      expect(same_bits(z, z1), "IC(0) apply bitwise identical: " + tag);
      double side = 0.0;
      ic.apply(b, z, [&] { side = b[0] + b[1]; });
      expect(same_bits(z, z1) && side == b[0] + b[1], "IC(0) apply with side task: " + tag);
      const SolveResult p = pcg(a, ic, b, x, 1e-10, 5000);
      expect(same_bits(x, x1) && p.iterations == p1.iterations, "pcg bitwise identical: " + tag);
      const SolveResult r = refined_solve(a, ic, b, hi, lo, 1e-13, 20000);
      expect(same_bits(hi, hi1) && same_bits(lo, lo1) && r.iterations == r1.iterations &&
                 r.relative_residual == r1.relative_residual,
             "refined solve bitwise identical: " + tag);
      // 1 스레드로 인수분해한 것을 여러 스레드에서 써도 (직렬 풀이) 같다
      ic1.apply(b, z);
      expect(same_bits(z, z1), "serial IC(0) under threads: " + tag);
    }
  }

  // parallel_for 의 예외: 1 스레드에서 처음 나는 예외(가장 앞 인덱스)가 그대로 전달된다
  hcw::parallel::set_threads(4);
  std::string what;
  try {
    hcw::parallel::parallel_for(0, 10000, [](int lo, int hi) {
      for (int i = lo; i < hi; ++i) {
        if (i == 3000 || i == 9000) throw std::runtime_error("index " + std::to_string(i));
      }
    });
  } catch (const std::exception& e) {
    what = e.what();
  }
  expect(what == "index 3000", "parallel_for rethrows the first exception in index order");
  // 병렬 구간 안에서 다시 나누지 않는다 (작업 스레드 위에서 직렬로 돈다)
  int nested_ok = 1;
  hcw::parallel::for_each_thread([&](int) {
    if (hcw::parallel::available() != 1) nested_ok = 0;
  });
  expect(nested_ok == 1, "no nested parallel regions");
  hcw::parallel::set_threads(1);
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
  test_refined_solve();
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
  test_threads_bitwise();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
