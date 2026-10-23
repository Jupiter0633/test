// P5 단위 시험: Q8 형상함수, 스카이라인 해법, 요소 행렬, 하중, 정확해 식, 진동수 방정식, HHT 이론값.
// 검증 판정이 아니라 부품이 맞게 동작하는지 확인한다 (판정은 hcw_verify P5).

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "structure/assembly.h"
#include "structure/eigen.h"
#include "structure/loads.h"
#include "verify/p5_cases.h"
#include "verify/p5_exact.h"

using namespace hcw;
using namespace hcw::structure;

namespace {

int g_failures = 0, g_checks = 0;

void expect(bool ok, const std::string& what) {
  ++g_checks;
  if (!ok) {
    ++g_failures;
    std::cerr << "FAIL: " << what << "\n";
  }
}

void test_shape_functions() {
  double N[8], dx[8], dy[8];
  const double pts[3] = {-0.7, 0.1, 0.55};
  for (double xi : pts) {
    for (double eta : pts) {
      shape_functions(xi, eta, N, dx, dy);
      double s = 0, sx = 0, sy = 0;
      for (int a = 0; a < 8; ++a) {
        s += N[a];
        sx += dx[a];
        sy += dy[a];
      }
      expect(std::abs(s - 1) < 1e-14 && std::abs(sx) < 1e-14 && std::abs(sy) < 1e-14, "Q8 partition of unity");
      // 미분을 차분과 비교
      double Np[8], Nm[8], t1[8], t2[8];
      const double h = 1e-6;
      shape_functions(xi + h, eta, Np, t1, t2);
      shape_functions(xi - h, eta, Nm, t1, t2);
      for (int a = 0; a < 8; ++a) expect(std::abs((Np[a] - Nm[a]) / (2 * h) - dx[a]) < 1e-8, "dN/dxi");
      shape_functions(xi, eta + h, Np, t1, t2);
      shape_functions(xi, eta - h, Nm, t1, t2);
      for (int a = 0; a < 8; ++a) expect(std::abs((Np[a] - Nm[a]) / (2 * h) - dy[a]) < 1e-8, "dN/deta");
    }
  }
  // 절점에서 크로네커 델타
  const double nx[8] = {-1, 1, 1, -1, 0, 1, 0, -1}, ny[8] = {-1, -1, 1, 1, -1, 0, 1, 0};
  for (int b = 0; b < 8; ++b) {
    shape_functions(nx[b], ny[b], N, dx, dy);
    for (int a = 0; a < 8; ++a) expect(std::abs(N[a] - (a == b ? 1.0 : 0.0)) < 1e-15, "Q8 Kronecker delta");
  }
}

void test_skyline() {
  // 작은 대칭 양정치 띠 행렬
  const int n = 7;
  std::vector<int> top = {0, 0, 1, 0, 2, 4, 3};
  SkylineMatrix A(top);
  std::vector<std::vector<double>> dense(n, std::vector<double>(n, 0.0));
  for (int j = 0; j < n; ++j) {
    for (int i = top[j]; i <= j; ++i) {
      const double v = i == j ? 10.0 + j : 1.0 / (1.0 + i + 2 * j);
      A.add(i, j, v);
      dense[i][j] = dense[j][i] = v;
    }
  }
  std::vector<double> x0 = {1, -2, 3, 0.5, -1, 2, 0.25}, b(n, 0.0);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) b[i] += dense[i][j] * x0[j];
  }
  A.factor();
  A.solve(b);
  double err = 0;
  for (int i = 0; i < n; ++i) err = std::max(err, std::abs(b[i] - x0[i]));
  expect(err < 1e-13, "skyline LDL^T solve");
}

Mesh block(double r0, double r1, double z0, double z1, int nr, int nz) {
  MeshBuilder mb;
  mb.add_block(r0, r1, z0, z1, nr, nz);
  return mb.build();
}

void test_element_and_loads() {
  const Material mat;
  Model model(block(0.2, 1.0, 0.0, 0.5, 3, 2), mat);
  model.finalize();
  // 강체 z 이동은 힘이 없다
  Vec f;
  model.apply_K(rigid_z(model), f);
  double mx = 0;
  for (double v : f) mx = std::max(mx, std::abs(v));
  expect(mx < 1e-12, "rigid z translation gives zero elastic force");
  // 전체 질량 = ρ π (r1² − r0²) H
  expect(std::abs(total_mass(model) - kPi * (1.0 - 0.04) * 0.5) < 1e-12, "total mass");
  // 압력 합력 (r = 1 면) = p 2π r H, 방향 −r
  Vec F(model.num_dofs(), 0.0);
  add_pressure(model.mesh(), boundary_edges_where(model.mesh(), [](double r, double) { return r == 1.0; }), 2.0, F);
  double fr = 0;
  for (int n = 0; n < model.num_dofs() / 2; ++n) fr += F[2 * n];
  expect(std::abs(fr + 2.0 * kTwoPi * 1.0 * 0.5) < 1e-12, "pressure resultant on r = 1");
  // 경계 변 수: 둘레 2·(3+2) 개
  expect(boundary_edges(model.mesh()).size() == 10, "boundary edge count");
  // L 자 컵: 두 블록이 맞닿는 절점이 합쳐짐 (요소 4 개, 절점 수 확인)
  MeshBuilder mb;
  mb.add_block(0, 1, 0, 0.1, 2, 1);
  mb.add_block(0.5, 1, 0.1, 1, 1, 2);
  const Mesh cup = mb.build();
  expect(cup.nodes.size() == 13 + 13 - 3, "cup mesh merges shared nodes (" + std::to_string(cup.nodes.size()) + ")");
}

void test_exact_formulas() {
  // 제조해의 div σ 식을 응력 식의 4차 정확 중심 차분과 비교
  using namespace hcw::verify::p5exact;
  const double h = 1e-3;
  for (double r : {0.13, 0.5, 0.91}) {
    for (double z : {0.07, 0.4, 0.88}) {
      double sp[4], sm[4], sp2[4], sm2[4], s0[4], div[2];
      auto d = [&](int comp, bool dr) {
        if (dr) {
          mms_stress(r + h, z, 1.0, sp); mms_stress(r - h, z, 1.0, sm);
          mms_stress(r + 2 * h, z, 1.0, sp2); mms_stress(r - 2 * h, z, 1.0, sm2);
        } else {
          mms_stress(r, z + h, 1.0, sp); mms_stress(r, z - h, 1.0, sm);
          mms_stress(r, z + 2 * h, 1.0, sp2); mms_stress(r, z - 2 * h, 1.0, sm2);
        }
        return (8 * (sp[comp] - sm[comp]) - (sp2[comp] - sm2[comp])) / (12 * h);
      };
      mms_stress(r, z, 1.0, s0);
      mms_div_stress(r, z, 1.0, div);
      const double fr = d(0, true) + d(3, false) + (s0[0] - s0[1]) / r;
      const double fz = d(3, true) + d(2, false) + s0[3] / r;
      expect(std::abs(fr - div[0]) < 1e-8 && std::abs(fz - div[1]) < 1e-8, "MMS div sigma formula");
      // 응력 식을 변위 차분에서 다시 구해 비교 (E = 1, ν = 0.3)
      double up[2], um[2], u0[2];
      Material mat;
      mms_u(r + h, z, 1.0, up); mms_u(r - h, z, 1.0, um);
      const double err = (up[0] - um[0]) / (2 * h), dzr_r = (up[1] - um[1]) / (2 * h);
      mms_u(r, z + h, 1.0, up); mms_u(r, z - h, 1.0, um);
      const double ezz = (up[1] - um[1]) / (2 * h), dur_z = (up[0] - um[0]) / (2 * h);
      mms_u(r, z, 1.0, u0);
      const auto s = stress_from_strain(mat, err, u0[0] / r, ezz, dur_z + dzr_r);
      for (int c = 0; c < 4; ++c) {
        expect(std::abs(s[c] - s0[c]) < 1e-5 * (1 + std::abs(s0[c])),
               "MMS stress formula r=" + std::to_string(r) + " z=" + std::to_string(z) + " c=" + std::to_string(c) +
                   " fd=" + std::to_string(s[c]) + " exact=" + std::to_string(s0[c]));
      }
    }
  }
  // 원판 해: 윗면 σ_zz = −1, 아랫면 0, 중앙 처짐 = p a⁴/(64D)
  for (double hh : {0.1, 0.02}) {
    double s[4], u[2];
    plate_stress(0.4, hh / 2, hh, s);
    expect(std::abs(s[2] + 1) < 1e-12 && std::abs(s[3]) < 1e-12, "plate top traction");
    plate_stress(0.4, -hh / 2, hh, s);
    expect(std::abs(s[2]) < 1e-12 && std::abs(s[3]) < 1e-12, "plate bottom traction");
    plate_u(0, 0, hh, u);
    const double D = hh * hh * hh / (12 * (1 - 0.09));
    expect(std::abs(-u[1] * 64 * D - 1) < 1e-12, "plate centre deflection equals Kirchhoff value");
  }
}

void test_reference_values() {
  // tools/p5_reference_values.py 의 값과 비교
  const auto w1 = verify::hollow_cylinder_radial_omegas(0.5, 1.0, 0.3, 2);
  expect(w1.size() == 2 && std::abs(w1[0] - 1.4694732417) < 1e-9 && std::abs(w1[1] - 7.4639423187) < 1e-9,
         "hollow cylinder a=0.5 radial frequencies");
  const auto w2 = verify::hollow_cylinder_radial_omegas(0.9, 1.0, 0.3, 1);
  expect(w2.size() == 1 && std::abs(w2[0] - 1.1048413974) < 1e-9, "hollow cylinder a=0.9 radial frequency");
  for (double Om : {0.1, 1.0, 10.0}) {
    double rho, theta;
    verify::hht_dominant_root(Om, 0.0, 0.0, rho, theta);
    expect(std::abs(rho - 1) < 1e-13 && std::abs(theta - 2 * std::atan(Om / 2)) < 1e-12, "Newmark closed form");
  }
  for (double al : {-0.1, -0.3}) {
    double rho, theta;
    verify::hht_dominant_root(1e6, al, 0.0, rho, theta);
    expect(std::abs(rho - (1 + al) / (1 - al)) < 1e-6, "HHT spectral radius at infinity");
  }
}

void test_eigen_small() {
  // 막대 ν = 0, 고정–자유: 가장 낮은 진동수가 π/2 에 가까움 (Nz = 16)
  Material mat;
  mat.nu = 0.0;
  Model model(block(0, 0.1, 0, 1, 1, 16), mat);
  model.fix_axis();
  model.fix_where([](double, double z) { return z == 0.0; }, 1);
  model.finalize();
  const EigenResult er = lowest_modes(model, 2, 0.0);
  expect(er.converged, "subspace iteration converges");
  expect(std::abs(std::sqrt(er.omega2[0]) / (kPi / 2) - 1) < 1e-6, "rod first axial frequency");
}

}  // namespace

int main() {
  test_shape_functions();
  test_skyline();
  test_element_and_loads();
  test_exact_formulas();
  test_reference_values();
  test_eigen_small();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
