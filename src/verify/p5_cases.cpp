#include "verify/p5_cases.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <sstream>
#include <stdexcept>

#include "structure/assembly.h"
#include "structure/eigen.h"
#include "structure/loads.h"
#include "structure/time_integration.h"
#include "verify/p5_exact.h"
#include "verify/report_util.h"

namespace hcw::verify {

using namespace hcw::structure;

namespace {

const char* const kCriteriaCommit = "9bf931b";

// 정확 재현 상한 (03 §1.5 공통 정의)
constexpr double kExactDisp = 1e-10;
constexpr double kExactStress = 1e-9;
// 선형 풀이 잔차 상한
constexpr double kSolverResidual = 1e-10;

using Field2 = std::function<std::array<double, 2>(double, double)>;
using Field4 = std::function<std::array<double, 4>(double, double)>;

TestReport new_report(const std::string& id, const std::string& title) {
  TestReport rep;
  rep.id = id;
  rep.title = title;
  rep.phase = "P5";
  rep.criteria_section = "1.5";
  rep.criteria_commit = kCriteriaCommit;
  return rep;
}

Criterion make(const std::string& text, bool passed, const std::string& detail) {
  Criterion c;
  c.text = text;
  c.passed = passed;
  c.detail = detail;
  return c;
}

// 가장 고운 격자의 값 ≤ limit
Criterion finest_at_most(const std::vector<GridRun>& runs, const std::string& key, double limit) {
  const double v = runs.back().metrics.at(key);
  return make(key + ": finest grid <= " + fmt(limit), v <= limit, fmt(v));
}

// 모든 실행의 값 ≤ limit
Criterion all_at_most(const std::vector<GridRun>& runs, const std::string& key, double limit) {
  double worst = 0.0;
  bool ok = true;
  for (const auto& r : runs) {
    const auto it = r.metrics.find(key);
    if (it == r.metrics.end()) continue;
    worst = std::max(worst, it->second);
    ok = ok && it->second <= limit;
  }
  return make(key + ": every run <= " + fmt(limit), ok, "max " + fmt(worst));
}

void check_residual(GridRun& run, double res) {
  run.metrics["max_solver_residual"] = std::max(run.metrics["max_solver_residual"], res);
  if (!(res <= kSolverResidual) && run.failure.empty()) run.failure = "FAILED_SOLVER (residual " + fmt(res) + ")";
}

void check_finite(GridRun& run) {
  for (const auto& [k, v] : run.metrics) {
    if (!std::isfinite(v) && run.failure.empty()) run.failure = "FAILED_NONFINITE (" + k + ")";
  }
}

Vec nodal(const Mesh& m, const Field2& f) {
  Vec u(2 * m.nodes.size());
  for (std::size_t n = 0; n < m.nodes.size(); ++n) {
    const auto v = f(m.nodes[n].r, m.nodes[n].z);
    u[2 * n] = v[0];
    u[2 * n + 1] = v[1];
  }
  return u;
}

// 응력장 σ 의 표면력 σ·n
structure::TractionFn traction_of(const Field4& s) {
  return [s](double r, double z, double nr, double nz) {
    const auto v = s(r, z);
    return std::array<double, 2>{v[0] * nr + v[3] * nz, v[3] * nr + v[2] * nz};
  };
}

struct Norms {
  double disp = 0.0, stress = 0.0;
};

// 변위·응력 L2 상대오차 (3×3 가우스, 2πr 가중)
Norms l2_errors(const Model& model, const Vec& u, const Field2& ue, const Field4& se) {
  const Mesh& mesh = model.mesh();
  const auto& gp = gauss_points(3);
  const auto& gw = gauss_weights(3);
  double du = 0, nu = 0, ds = 0, ns = 0;
  for (int k = 0; k < static_cast<int>(mesh.elems.size()); ++k) {
    const auto& e = mesh.elems[k];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const PointEval p = evaluate(mesh, k, gp[i], gp[j]);
        const double w = gw[i] * gw[j] * p.detJ * kTwoPi * p.r;
        double uh[2] = {0, 0};
        for (int a = 0; a < 8; ++a) {
          uh[0] += p.N[a] * u[2 * e[a]];
          uh[1] += p.N[a] * u[2 * e[a] + 1];
        }
        const auto ux = ue(p.r, p.z);
        for (int c = 0; c < 2; ++c) {
          du += w * (uh[c] - ux[c]) * (uh[c] - ux[c]);
          nu += w * ux[c] * ux[c];
        }
        if (se) {
          const auto sh = element_stress(mesh, k, model.material(), u, gp[i], gp[j]);
          const auto sx = se(p.r, p.z);
          for (int c = 0; c < 4; ++c) {
            ds += w * (sh[c] - sx[c]) * (sh[c] - sx[c]);
            ns += w * sx[c] * sx[c];
          }
        }
      }
    }
  }
  Norms n;
  n.disp = std::sqrt(du / nu);
  n.stress = se ? std::sqrt(ds / ns) : 0.0;
  return n;
}

// 절점 변위 최대 오차 / max|u|. shift_z 만큼 z 성분을 빼고 비교 (관성 해제의 강체 이동 제거)
double max_nodal_error(const Mesh& mesh, const Vec& u, const Field2& ue, double shift_z = 0.0) {
  double err = 0.0, ref = 0.0;
  for (std::size_t n = 0; n < mesh.nodes.size(); ++n) {
    const auto x = ue(mesh.nodes[n].r, mesh.nodes[n].z);
    err = std::max(err, std::abs(u[2 * n] - x[0]));
    err = std::max(err, std::abs(u[2 * n + 1] - shift_z - x[1]));
    ref = std::max({ref, std::abs(x[0]), std::abs(x[1])});
  }
  return err / ref;
}

// 3×3 가우스점 응력 최대 오차 / ref (ref ≤ 0 이면 정확해 응력의 최댓값)
double max_gauss_stress_error(const Model& model, const Vec& u, const Field4& se, double ref = 0.0) {
  const Mesh& mesh = model.mesh();
  const auto& gp = gauss_points(3);
  double err = 0.0, mx = 0.0;
  for (int k = 0; k < static_cast<int>(mesh.elems.size()); ++k) {
    for (double xi : gp) {
      for (double eta : gp) {
        const PointEval p = evaluate(mesh, k, xi, eta);
        const auto sh = element_stress(mesh, k, model.material(), u, xi, eta);
        const auto sx = se(p.r, p.z);
        for (int c = 0; c < 4; ++c) {
          err = std::max(err, std::abs(sh[c] - sx[c]));
          mx = std::max(mx, std::abs(sx[c]));
        }
      }
    }
  }
  return err / (ref > 0.0 ? ref : mx);
}

// 실행 하나를 예외로부터 지킨다
template <class F>
GridRun guarded(const std::string& label, int n, F&& body) {
  GridRun run;
  run.label = label;
  run.n = n;
  const auto t0 = std::chrono::steady_clock::now();
  try {
    body(run);
  } catch (const std::exception& ex) {
    run.failure = std::string("FAILED_EXCEPTION (") + ex.what() + ")";
  }
  run.seconds = seconds_since(t0);
  check_finite(run);
  return run;
}

// 실행 묶음의 (★5) 기준: 단조 감소 + 관측 차수
void star5(TestReport& rep, const std::vector<GridRun>& runs, const std::string& key, double min_order,
           const std::string& prefix) {
  Criterion m = monotone(runs, key), o = order_at_least(runs, key, min_order);
  m.text = prefix + m.text;
  o.text = prefix + o.text;
  rep.criteria.push_back(m);
  rep.criteria.push_back(o);
}

bool any_failed(const std::vector<GridRun>& runs) {
  return std::any_of(runs.begin(), runs.end(), [](const GridRun& r) { return !r.failure.empty(); });
}

// 4-A(b)·4-G 제조해
std::array<double, 2> mms_u(double r, double z) {
  std::array<double, 2> o{};
  p5exact::mms_u(r, z, 1e-3, o.data());
  return o;
}
std::array<double, 4> mms_s(double r, double z) {
  std::array<double, 4> o{};
  p5exact::mms_stress(r, z, 1e-3, o.data());
  return o;
}
std::array<double, 2> mms_minus_div(double r, double z) {
  std::array<double, 2> o{};
  p5exact::mms_div_stress(r, z, 1e-3, o.data());
  return {-o[0], -o[1]};
}

}  // namespace

// ------------------------------------------------------------------ 4-A
TestReport run_4a(bool quick) {
  TestReport rep = new_report("4-A", "static exactness: (a) patch test, (b) manufactured solution");
  const Material mat;  // E = 1, ν = 0.3

  // (a) 찌그러진 4×4 격자의 패치 시험
  const double c1 = 1e-3, c2 = -2e-3;
  const auto s_patch = stress_from_strain(mat, c1, c1, c2, 0.0);
  const Field2 ue = [&](double r, double z) { return std::array<double, 2>{c1 * r, c2 * z}; };
  const Field4 se = [&](double, double) { return s_patch; };
  GridRun pa = guarded("(a) patch 4x4 distorted", 4, [&](GridRun& run) {
    MeshBuilder mb;
    mb.add_block(0, 1, 0, 1, 4, 4);
    Mesh mesh = mb.build();
    for (auto& nd : mesh.nodes) {
      const double x = nd.r * 4, y = nd.z * 4;
      const long i = std::lround(x), j = std::lround(y);
      if (std::abs(x - i) < 1e-9 && std::abs(y - j) < 1e-9 && i >= 1 && i <= 3 && j >= 1 && j <= 3) {
        const double s1 = static_cast<double>((i + 2 * j) % 3 - 1), s2 = static_cast<double>((2 * i + j) % 3 - 1);
        nd.r += 0.1 * 0.25 * s1;
        nd.z += 0.1 * 0.25 * s2;
      }
    }
    recenter_midside_nodes(mesh);
    Model model(mesh, mat);
    model.fix_axis();
    model.fix_where([](double, double z) { return z == 0.0; }, 1);
    model.finalize();
    Vec F(model.num_dofs(), 0.0);
    add_traction(model.mesh(), boundary_edges_where(model.mesh(), [](double r, double z) { return r == 1.0 || z == 1.0; }),
                 traction_of(se), F);
    double res = 0;
    const Vec u = solve_static(model, F, nodal(model.mesh(), ue), &res);
    check_residual(run, res);
    run.metrics["disp_max_error"] = max_nodal_error(model.mesh(), u, ue);
    run.metrics["stress_max_error"] = max_gauss_stress_error(model, u, se);
  });
  rep.runs.push_back(pa);

  // (b) 정적 제조해
  std::vector<GridRun> runs;
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{4, 8, 16, 32, 64};
  for (int n : ns) {
    runs.push_back(guarded("(b) MMS N=" + std::to_string(n), n, [&](GridRun& run) {
      MeshBuilder mb;
      mb.add_block(0, 1, 0, 1, n, n);
      Model model(mb.build(), mat);
      model.fix_axis();
      model.fix_where([](double r, double z) { return r == 1.0 || z == 0.0; }, 0);
      model.fix_where([](double r, double z) { return r == 1.0 || z == 0.0; }, 1);
      model.finalize();
      Vec F(model.num_dofs(), 0.0);
      add_body_force(model.mesh(), mms_minus_div, F);
      add_traction(model.mesh(), boundary_edges_where(model.mesh(), [](double, double z) { return z == 1.0; }),
                   traction_of(mms_s), F);
      double res = 0;
      const Vec u = solve_static(model, F, nodal(model.mesh(), mms_u), &res);
      check_residual(run, res);
      const Norms e = l2_errors(model, u, mms_u, mms_s);
      run.metrics["disp_L2"] = e.disp;
      run.metrics["stress_L2"] = e.stress;
    }));
  }
  rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());

  if (pa.failure.empty()) {
    rep.criteria.push_back(make("(a) patch: nodal displacement error / max|u| <= 1e-10",
                                pa.metrics["disp_max_error"] <= kExactDisp, fmt(pa.metrics["disp_max_error"])));
    rep.criteria.push_back(make("(a) patch: Gauss-point stress error / max|sigma| <= 1e-9",
                                pa.metrics["stress_max_error"] <= kExactStress, fmt(pa.metrics["stress_max_error"])));
  }
  if (!any_failed(runs)) {
    star5(rep, runs, "disp_L2", 2.7, "(b) ");
    star5(rep, runs, "stress_L2", 1.7, "(b) ");
  }
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-B
TestReport run_4b(bool quick) {
  TestReport rep = new_report("4-B", "thick cylinder under internal pressure (Lame, plane strain)");
  const Material mat;
  const double a = 0.5, b = 1.0, H = 0.25, p = 1.0;
  const double E = mat.E, nu = mat.nu, k = p * a * a / (b * b - a * a);
  const Field2 ue = [&](double r, double) {
    return std::array<double, 2>{(1 + nu) * k / E * ((1 - 2 * nu) * r + b * b / r), 0.0};
  };
  const Field4 se = [&](double r, double) {
    return std::array<double, 4>{k * (1 - b * b / (r * r)), k * (1 + b * b / (r * r)), 2 * nu * k, 0.0};
  };
  std::vector<GridRun> runs;
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{4, 8, 16, 32};
  for (int n : ns) {
    runs.push_back(guarded("Nr=" + std::to_string(n), n, [&](GridRun& run) {
      MeshBuilder mb;
      mb.add_block(a, b, 0, H, n, n / 2);
      Model model(mb.build(), mat);
      model.fix_where([&](double, double z) { return z == 0.0 || z == H; }, 1);
      model.finalize();
      Vec F(model.num_dofs(), 0.0);
      add_pressure(model.mesh(), boundary_edges_where(model.mesh(), [&](double r, double) { return r == a; }), p, F);
      double fr = 0.0;
      for (int nd = 0; nd < model.num_dofs() / 2; ++nd) fr += F[2 * nd];
      run.metrics["load_sum_error"] = std::abs(fr - p * kTwoPi * a * H) / (p * kTwoPi * a * H);
      double res = 0;
      const Vec u = solve_static(model, F, Vec(model.num_dofs(), 0.0), &res);
      check_residual(run, res);
      const Norms e = l2_errors(model, u, ue, se);
      run.metrics["disp_L2"] = e.disp;
      run.metrics["stress_L2"] = e.stress;
      const auto rec = recovered_nodal_stress(model.mesh(), mat, u);
      double hoop = 0.0, ura = 0.0;
      for (std::size_t nd = 0; nd < model.mesh().nodes.size(); ++nd) {
        if (model.mesh().nodes[nd].r != a) continue;
        hoop = std::max(hoop, std::abs(rec[nd][1] - se(a, 0)[1]) / std::abs(se(a, 0)[1]));
        ura = std::max(ura, std::abs(u[2 * nd] - ue(a, 0)[0]) / ue(a, 0)[0]);
      }
      run.metrics["hoop_inner_recovered"] = hoop;
      run.metrics["ur_inner_error"] = ura;
    }));
  }
  rep.runs = runs;
  if (!any_failed(runs)) {
    star5(rep, runs, "disp_L2", 2.7, "");
    star5(rep, runs, "hoop_inner_recovered", 1.7, "");
    rep.criteria.push_back(finest_at_most(runs, "hoop_inner_recovered", 0.005));
    rep.criteria.push_back(all_at_most(runs, "load_sum_error", 1e-12));
    rep.notes.push_back("report only: stress_L2, ur_inner_error");
  }
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-C
TestReport run_4c(bool quick) {
  TestReport rep = new_report("4-C", "circular plate bending, exact polynomial elasticity solution");
  const Material mat;
  struct Case {
    std::string name;
    double h;
    std::vector<std::pair<int, int>> meshes;
  };
  std::vector<Case> cases = {{"(a) h=0.1", 0.1, {{4, 2}, {8, 4}, {16, 8}, {32, 16}}},
                             {"(b) h=0.02", 0.02, {{10, 1}, {20, 2}, {40, 4}, {80, 8}}}};
  for (auto& c : cases) {
    if (quick) c.meshes.resize(2);
    const double h = c.h;
    const Field2 ue = [h](double r, double z) {
      std::array<double, 2> o{};
      p5exact::plate_u(r, z, h, o.data());
      return o;
    };
    const Field4 se = [h](double r, double z) {
      std::array<double, 4> o{};
      p5exact::plate_stress(r, z, h, o.data());
      return o;
    };
    std::vector<GridRun> runs;
    for (const auto& [nr, nz] : c.meshes) {
      runs.push_back(guarded(c.name + " Nr=" + std::to_string(nr) + " Nz=" + std::to_string(nz), nr, [&](GridRun& run) {
        MeshBuilder mb;
        mb.add_block(0, 1, -h / 2, h / 2, nr, nz);
        Model model(mb.build(), mat);
        model.fix_axis();
        model.fix_where([](double r, double) { return r == 1.0; }, 0);
        model.fix_where([](double r, double) { return r == 1.0; }, 1);
        model.finalize();
        Vec F(model.num_dofs(), 0.0);
        add_pressure(model.mesh(), boundary_edges_where(model.mesh(), [h](double, double z) { return z == h / 2; }), 1.0,
                     F);
        double res = 0;
        const Vec u = solve_static(model, F, nodal(model.mesh(), ue), &res);
        check_residual(run, res);
        const Norms e = l2_errors(model, u, ue, se);
        run.metrics["disp_L2"] = e.disp;
        run.metrics["stress_L2"] = e.stress;
        const auto rec = recovered_nodal_stress(model.mesh(), mat, u);
        for (std::size_t nd = 0; nd < model.mesh().nodes.size(); ++nd) {
          const auto& q = model.mesh().nodes[nd];
          if (q.r == 0.0 && q.z == 0.0) {
            run.metrics["center_deflection_error"] = std::abs(u[2 * nd + 1] - ue(0, 0)[1]) / std::abs(ue(0, 0)[1]);
          }
          if (q.r == 0.0 && q.z == -h / 2) {
            run.metrics["bottom_center_srr_recovered"] = rec[nd][0];
          }
        }
        run.metrics["bottom_center_srr_exact"] = se(0.0, -h / 2)[0];
      }));
    }
    rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());
    if (!any_failed(runs)) {
      star5(rep, runs, "disp_L2", 2.7, c.name + " ");
      star5(rep, runs, "stress_L2", 1.7, c.name + " ");
      Criterion f = finest_at_most(runs, "center_deflection_error", 1e-3);
      f.text = c.name + " " + f.text;
      rep.criteria.push_back(f);
    }
  }
  rep.notes.push_back("report only: bottom-centre recovered sigma_rr (exact 47.75 for h=0.1); "
                      "Kirchhoff p a^4/(64D) equals the exact centre deflection by construction");
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-E
std::vector<double> hollow_cylinder_radial_omegas(double a, double b, double nu, int count) {
  const double lam = nu / ((1 + nu) * (1 - 2 * nu)), mu = 1.0 / (2 * (1 + nu));
  const double c11 = lam + 2 * mu, cp = std::sqrt(c11);
  auto srr = [&](bool second, double k, double r) {
    const double x = k * r;
    const double f = second ? std::cyl_neumann(1.0, x) : std::cyl_bessel_j(1.0, x);
    const double f0 = second ? std::cyl_neumann(0.0, x) : std::cyl_bessel_j(0.0, x);
    const double fp = f0 - f / x;  // Z₁′(x) = Z₀(x) − Z₁(x)/x
    return c11 * k * fp + lam * f / r;
  };
  auto det = [&](double k) { return srr(false, k, a) * srr(true, k, b) - srr(true, k, a) * srr(false, k, b); };
  std::vector<double> out;
  double k0 = 0.05, d0 = det(k0);
  for (double k = k0 + 1e-3; k < 200.0 && static_cast<int>(out.size()) < count; k += 1e-3) {
    const double d1 = det(k);
    if ((d0 < 0) != (d1 < 0)) {
      double lo = k - 1e-3, hi = k, flo = d0;
      for (int it = 0; it < 200 && hi - lo > 1e-15 * hi; ++it) {
        const double mid = 0.5 * (lo + hi), fm = det(mid);
        if ((fm < 0) == (flo < 0)) {
          lo = mid;
          flo = fm;
        } else {
          hi = mid;
        }
      }
      out.push_back(0.5 * (lo + hi) * cp);
    }
    d0 = d1;
  }
  return out;
}

namespace {

// u 의 성분 comp 만 남긴 M-노름 (모드 φᵀMφ = 1 기준)
double component_fraction(const Model& model, const Vec& phi, int comp) {
  Vec x(phi.size(), 0.0), y;
  for (std::size_t i = comp; i < phi.size(); i += 2) x[i] = phi[i];
  model.apply_M(x, y);
  return std::sqrt(std::max(0.0, model.dot(x, y)));
}

}  // namespace

TestReport run_4e(bool quick) {
  TestReport rep = new_report("4-E", "natural frequencies: rod axial, hollow cylinder radial, (c) clamped plate report");
  struct Case {
    std::string name;
    double nu;
    std::function<Mesh(int)> mesh;
    std::function<void(Model&)> supports;
    std::vector<int> ns;
    std::vector<double> exact;
    int other_comp;  // 모드 식별: 이 성분이 반올림 수준이어야 함
  };
  const double kH1 = 0.125, kH2 = 0.05;
  std::vector<Case> cases;
  cases.push_back({"(a) rod", 0.0,
                   [](int n) {
                     MeshBuilder mb;
                     mb.add_block(0, 0.1, 0, 1, 1, n);
                     return mb.build();
                   },
                   [](Model& m) {
                     m.fix_axis();
                     m.fix_where([](double, double z) { return z == 0.0; }, 1);
                   },
                   {4, 8, 16, 32},
                   {kPi / 2, 3 * kPi / 2, 5 * kPi / 2},
                   0});
  cases.push_back({"(b1) a=0.5", 0.3,
                   [kH1](int n) {
                     MeshBuilder mb;
                     mb.add_block(0.5, 1, 0, kH1, n, n / 4);
                     return mb.build();
                   },
                   [kH1](Model& m) { m.fix_where([kH1](double, double z) { return z == 0.0 || z == kH1; }, 1); },
                   {4, 8, 16, 32},
                   hollow_cylinder_radial_omegas(0.5, 1.0, 0.3, 2),
                   1});
  cases.push_back({"(b2) a=0.9", 0.3,
                   [kH2](int n) {
                     MeshBuilder mb;
                     mb.add_block(0.9, 1, 0, kH2, n, n / 2);
                     return mb.build();
                   },
                   [kH2](Model& m) { m.fix_where([kH2](double, double z) { return z == 0.0 || z == kH2; }, 1); },
                   {2, 4, 8, 16},
                   hollow_cylinder_radial_omegas(0.9, 1.0, 0.3, 1),
                   1});
  for (auto& c : cases) {
    if (quick) c.ns.resize(2);
    std::vector<GridRun> runs;
    const int nm = static_cast<int>(c.exact.size());
    for (int n : c.ns) {
      runs.push_back(guarded(c.name + " N=" + std::to_string(n), n, [&](GridRun& run) {
        Material mat;
        mat.nu = c.nu;
        Model model(c.mesh(n), mat);
        c.supports(model);
        model.finalize();
        const EigenResult er = lowest_modes(model, nm, 0.0);
        check_residual(run, er.max_rel_residual);
        run.metrics["eigen_iterations"] = er.iterations;
        if (!er.converged && run.failure.empty()) run.failure = "FAILED_SOLVER (eigenvalues did not converge)";
        double other = 0.0;
        for (int m = 0; m < nm; ++m) {
          const double w = std::sqrt(er.omega2[m]);
          run.metrics["omega" + std::to_string(m + 1)] = w;
          run.metrics["omega" + std::to_string(m + 1) + "_error"] = std::abs(w - c.exact[m]) / c.exact[m];
          // 식별: 역반복으로 다듬은 모드의 다른 성분 (부분공간 반복의 덜 수렴한 벡터 성분을 없앰).
          // 판정하는 진동수는 다듬기 전 고유값 해법의 값이다.
          Vec phi = er.modes[m];
          refine_mode(model, er.omega2[m], phi);
          other = std::max(other, component_fraction(model, phi, c.other_comp));
        }
        run.metrics["mode_other_component"] = other;
      }));
    }
    rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());
    if (!any_failed(runs)) {
      std::ostringstream ex;
      ex.precision(12);
      for (double w : c.exact) ex << w << " ";
      rep.notes.push_back(c.name + " exact omega: " + ex.str());
      rep.criteria.push_back(make(c.name + " mode identification: other component (M-norm) <= 1e-8 on every grid",
                                  all_at_most(runs, "mode_other_component", 1e-8).passed,
                                  all_at_most(runs, "mode_other_component", 1e-8).detail));
      for (int m = 0; m < nm; ++m) {
        const std::string key = "omega" + std::to_string(m + 1) + "_error";
        star5(rep, runs, key, 3.5, c.name + " ");
        Criterion f = finest_at_most(runs, key, 1e-3);
        f.text = c.name + " " + f.text;
        rep.criteria.push_back(f);
      }
    }
  }
  conclude(rep, quick);
  return rep;
}

// 4-E(c) 보고만: 얇은 고정 원판과 Kirchhoff 값. 종료 조건에서 "보고서가 있음" 으로 따로 다루므로 별도 보고서로 낸다.
TestReport run_4e_c(bool quick) {
  TestReport rep = new_report("4-E(c)", "clamped thin plate frequencies vs Kirchhoff (report only)");
  // (c) 보고만: 얇은 고정 원판과 Kirchhoff 값
  const double lam2[3] = {10.2158262299, 39.7711482365, 89.1041439740};
  const std::vector<std::pair<int, int>> pm = quick ? std::vector<std::pair<int, int>>{{20, 2}}
                                                    : std::vector<std::pair<int, int>>{{40, 4}, {80, 8}};
  for (const auto& [nr, nz] : pm) {
    rep.runs.push_back(guarded("(c) clamped plate h=0.02 Nr=" + std::to_string(nr) + " (report only)", nr,
                               [&](GridRun& run) {
                                 const double h = 0.02;
                                 const Material mat;
                                 MeshBuilder mb;
                                 mb.add_block(0, 1, -h / 2, h / 2, nr, nz);
                                 Model model(mb.build(), mat);
                                 model.fix_axis();
                                 model.fix_where([](double r, double) { return r == 1.0; }, 0);
                                 model.fix_where([](double r, double) { return r == 1.0; }, 1);
                                 model.finalize();
                                 const EigenResult er = lowest_modes(model, 3, 0.0);
                                 check_residual(run, er.max_rel_residual);
                                 const double D = mat.E * h * h * h / (12 * (1 - mat.nu * mat.nu));
                                 for (int m = 0; m < 3; ++m) {
                                   const double wk = lam2[m] * std::sqrt(D / (mat.rho * h));
                                   run.metrics["omega" + std::to_string(m + 1) + "_over_kirchhoff"] =
                                       std::sqrt(er.omega2[m]) / wk;
                                 }
                               }));
  }
  rep.notes.push_back("report only: solid model vs Kirchhoff plate theory (thickness and clamp modelling differ)");
  rep.verdict = quick ? "QUICK_CHECK_ONLY" : "REPORT_ONLY";
  for (const auto& r : rep.runs) {
    if (!r.failure.empty()) {
      rep.verdict += " (" + r.failure.substr(0, r.failure.find(' ')) + ")";
      break;
    }
  }
  return rep;
}

// ------------------------------------------------------------------ 4-F
namespace {

Mesh cup_mesh() {
  MeshBuilder mb;
  mb.add_block(0, 1, 0, 0.1, 40, 4);
  mb.add_block(0.9, 1, 0.1, 1, 4, 36);
  return mb.build();
}

}  // namespace

TestReport run_4f(bool quick) {
  TestReport rep = new_report("4-F", "energy conservation of Newmark average acceleration (free-free cup)");
  const Material mat;
  const int steps = quick ? 200 : 2000;
  struct Variant {
    std::string name;
    double alpha;
    bool rayleigh;
  };
  const std::vector<Variant> variants = {
      {"Newmark undamped", 0.0, false}, {"HHT alpha=-0.1", -0.1, false}, {"Newmark + Rayleigh", 0.0, true}};
  for (const auto& var : variants) {
    rep.runs.push_back(guarded(var.name, steps, [&](GridRun& run) {
      Model model(cup_mesh(), mat);
      model.fix_axis();
      model.finalize();
      // 가장 낮은 탄성 모드. z 강체 이동 (ω = 0) 은 미리 알고 있으므로 빼고 구한다.
      const EigenResult er = lowest_modes(model, 1, -1.0, 2000, {rigid_z(model)});
      check_residual(run, er.max_rel_residual);
      if (!er.converged && run.failure.empty()) run.failure = "FAILED_SOLVER (eigenvalues did not converge)";
      const double w1 = std::sqrt(er.omega2[0]);
      run.metrics["omega1"] = w1;
      const double dt = 0.05 / w1;
      Rayleigh ray;
      if (var.rayleigh) ray.a1 = 2 * 0.02 / w1;  // 강성 비례, ζ₁ = 0.02
      HhtIntegrator ig(model, dt, var.alpha, ray);
      const Vec u0 = nodal(model.mesh(), mms_u);
      const Vec v0 = nodal(model.mesh(), [](double r, double z) {
        return std::array<double, 2>{1e-3 * r * std::cos(2 * kPi * z), 1e-3 * std::sin(kPi * r * r) * z};
      });
      const Vec zero(model.num_dofs(), 0.0);
      ig.initialize(u0, v0, zero, zero);
      const double E0 = ig.energy();
      double prev = E0, max_change = 0.0;
      int increases = 0;
      for (int s = 0; s < steps; ++s) {
        ig.step(zero, zero, zero, zero);
        const double En = ig.energy();
        max_change = std::max(max_change, std::abs(En - E0) / E0);
        if (En > prev * (1 + 1e-12)) ++increases;
        prev = En;
      }
      check_residual(run, ig.max_rel_residual());
      run.metrics["energy_max_relative_change"] = max_change;
      run.metrics["energy_final_over_initial"] = prev / E0;
      run.metrics["energy_increasing_steps"] = increases;
      run.metrics["dt_times_omega1"] = dt * w1;
    }));
  }
  if (!any_failed(rep.runs)) {
    const auto& r0 = rep.runs[0].metrics;
    rep.criteria.push_back(make("Newmark undamped: max |E_n - E_0| / E_0 <= 1e-8",
                                r0.at("energy_max_relative_change") <= 1e-8, fmt(r0.at("energy_max_relative_change"))));
    for (int k = 1; k <= 2; ++k) {
      const double inc = rep.runs[k].metrics.at("energy_increasing_steps");
      rep.criteria.push_back(make(rep.runs[k].label + ": energy never increases (E_{n+1} <= E_n (1 + 1e-12))", inc == 0,
                                  fmt(inc) + " increasing steps of " + std::to_string(steps)));
    }
  }
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-G
TestReport run_4g(bool quick) {
  TestReport rep = new_report("4-G", "dynamic manufactured solution (space and time refined together)");
  const Material mat;
  const double w = 2 * kPi, T = 1.0;
  auto g = [w](double t) { return std::sin(w * t) + std::cos(w * t); };
  auto gd = [w](double t) { return w * (std::cos(w * t) - std::sin(w * t)); };
  auto gdd = [w](double t) { return -w * w * (std::sin(w * t) + std::cos(w * t)); };
  for (double alpha : {0.0, -0.1}) {
    const std::string tag = alpha == 0.0 ? "Newmark" : "HHT alpha=-0.1";
    std::vector<GridRun> runs;
    const std::vector<int> ns = quick ? std::vector<int>{8, 16} : std::vector<int>{8, 16, 32, 64};
    for (int n : ns) {
      runs.push_back(guarded(tag + " N=" + std::to_string(n), n, [&](GridRun& run) {
        MeshBuilder mb;
        mb.add_block(0, 1, 0, 1, n, n);
        Model model(mb.build(), mat);
        model.fix_axis();
        model.fix_where([](double r, double z) { return r == 1.0 || z == 0.0; }, 0);
        model.fix_where([](double r, double z) { return r == 1.0 || z == 0.0; }, 1);
        model.finalize();
        // F(t) = Fs g(t) + Fi g''(t): Fs = 체적력 −div σ(U) + 윗면 표면력, Fi = ρU 의 일관 하중
        Vec Fs(model.num_dofs(), 0.0), Fi(model.num_dofs(), 0.0);
        add_body_force(model.mesh(), mms_minus_div, Fs);
        add_traction(model.mesh(), boundary_edges_where(model.mesh(), [](double, double z) { return z == 1.0; }),
                     traction_of(mms_s), Fs);
        add_body_force(model.mesh(), [&](double r, double z) {
          const auto u = mms_u(r, z);
          return std::array<double, 2>{mat.rho * u[0], mat.rho * u[1]};
        }, Fi);
        const Vec U = nodal(model.mesh(), mms_u);
        auto scaled = [&](double s) {
          Vec x(U.size());
          for (std::size_t i = 0; i < x.size(); ++i) x[i] = s * U[i];
          return x;
        };
        auto force = [&](double t) {
          Vec f(U.size());
          for (std::size_t i = 0; i < f.size(); ++i) f[i] = Fs[i] * g(t) + Fi[i] * gdd(t);
          return f;
        };
        const double dt = 0.05 * 8.0 / n;
        const int steps = static_cast<int>(std::lround(T / dt));
        HhtIntegrator ig(model, dt, alpha);
        ig.initialize(scaled(g(0)), scaled(gd(0)), force(0), scaled(gdd(0)));
        for (int s = 1; s <= steps; ++s) {
          const double t = s * dt;
          ig.step(force(t), scaled(g(t)), scaled(gd(t)), scaled(gdd(t)));
        }
        check_residual(run, ig.max_rel_residual());
        const double gT = g(steps * dt), gdT = gd(steps * dt);
        const Norms e = l2_errors(model, ig.u(), [&](double r, double z) {
          const auto u = mms_u(r, z);
          return std::array<double, 2>{u[0] * gT, u[1] * gT};
        }, nullptr);
        const Norms ev = l2_errors(model, ig.v(), [&](double r, double z) {
          const auto u = mms_u(r, z);
          return std::array<double, 2>{u[0] * gdT, u[1] * gdT};
        }, nullptr);
        run.metrics["disp_L2_at_T"] = e.disp;
        run.metrics["velocity_L2_at_T"] = ev.disp;
        run.metrics["steps"] = steps;
      }));
    }
    rep.runs.insert(rep.runs.end(), runs.begin(), runs.end());
    if (!any_failed(runs)) star5(rep, runs, "disp_L2_at_T", 1.7, tag + " ");
  }
  rep.notes.push_back("report only: velocity_L2_at_T");
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-H
void hht_dominant_root(double Omega, double alpha, double xi, double& rho, double& theta, double* spurious) {
  // 1 자유도 (ω = 1): 미지수 (u, v, a)_{n+1} = A⁻¹ B (u, v, a)_n  (tools/p5_reference_values.py 의 hht_matrix 와 같음)
  const double dt = Omega, ga = (1 - 2 * alpha) / 2, be = (1 - alpha) * (1 - alpha) / 4, K = 1.0, C = 2 * xi;
  double A[3][3] = {{1, 0, -be * dt * dt}, {0, 1, -ga * dt}, {(1 + alpha) * K, (1 + alpha) * C, 1}};
  double B[3][3] = {{1, dt, (0.5 - be) * dt * dt}, {0, 1, (1 - ga) * dt}, {alpha * K, alpha * C, 0}};
  // 가우스 소거로 X = A⁻¹ B
  for (int c = 0; c < 3; ++c) {
    int piv = c;
    for (int r = c + 1; r < 3; ++r) {
      if (std::abs(A[r][c]) > std::abs(A[piv][c])) piv = r;
    }
    std::swap(A[c], A[piv]);
    std::swap(B[c], B[piv]);
    for (int r = 0; r < 3; ++r) {
      if (r == c) continue;
      const double f = A[r][c] / A[c][c];
      for (int k = 0; k < 3; ++k) {
        A[r][k] -= f * A[c][k];
        B[r][k] -= f * B[c][k];
      }
    }
  }
  double X[3][3];
  for (int r = 0; r < 3; ++r) {
    for (int k = 0; k < 3; ++k) X[r][k] = B[r][k] / A[r][r];
  }
  const double tr = X[0][0] + X[1][1] + X[2][2];
  const double m2 = X[0][0] * X[1][1] - X[0][1] * X[1][0] + X[0][0] * X[2][2] - X[0][2] * X[2][0] +
                    X[1][1] * X[2][2] - X[1][2] * X[2][1];
  const double det = X[0][0] * (X[1][1] * X[2][2] - X[1][2] * X[2][1]) -
                     X[0][1] * (X[1][0] * X[2][2] - X[1][2] * X[2][0]) +
                     X[0][2] * (X[1][0] * X[2][1] - X[1][1] * X[2][0]);
  auto f = [&](double l) { return ((l - tr) * l + m2) * l - det; };
  double lo = -2.0, hi = 2.0;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    if ((f(mid) < 0) == (f(lo) < 0)) lo = mid;
    else hi = mid;
  }
  double l3 = 0.5 * (lo + hi);
  for (int it = 0; it < 3; ++it) {  // Newton 다듬기
    const double d = (3 * l3 - 2 * tr) * l3 + m2;
    if (d != 0.0) l3 -= f(l3) / d;
  }
  // (λ − λ₃)(λ² + pλ + q)
  const double p = l3 - tr, q = m2 + p * l3;
  rho = std::sqrt(q);
  theta = std::acos(std::clamp(-p / (2 * rho), -1.0, 1.0));
  if (spurious) *spurious = l3;
}

TestReport run_4h(bool quick) {
  TestReport rep = new_report("4-H", "time integration vs exact discrete single-mode solution");
  const Material mat;
  MeshBuilder mb;
  mb.add_block(0.5, 1, 0, 0.125, 16, 4);
  Model model(mb.build(), mat);
  model.fix_where([](double, double z) { return z == 0.0 || z == 0.125; }, 1);
  model.finalize();

  double w = 0.0;
  Vec phi;
  std::string setup_failure;
  try {
    const EigenResult er = lowest_modes(model, 2, 0.0);
    if (!er.converged) setup_failure = "FAILED_SOLVER (eigenvalues did not converge)";
    phi = er.modes[0];
    w = std::sqrt(refine_mode(model, er.omega2[0], phi));
    rep.notes.push_back("mode 1 omega (refined Rayleigh quotient) = " + fmt(w));
  } catch (const std::exception& ex) {
    setup_failure = std::string("FAILED_EXCEPTION (") + ex.what() + ")";
  }
  Vec Mphi;
  if (setup_failure.empty()) model.apply_M(phi, Mphi);

  struct Case {
    std::string name;
    double alpha;
    Rayleigh ray;
    bool step_load;
  };
  const double zeta = 0.02;
  const int steps = quick ? 200 : 600;
  const int skip = 50;
  for (double Om : {0.1, 1.0}) {
    std::vector<Case> cases = {{"(a) Newmark", 0.0, {}, false},
                               {"(b) HHT alpha=-0.1", -0.1, {}, false},
                               {"(b) HHT alpha=-0.3", -0.3, {}, false},
                               {"(c) Rayleigh mass", 0.0, {2 * zeta * w, 0.0}, false},
                               {"(c) Rayleigh stiffness", 0.0, {0.0, 2 * zeta / w}, false},
                               {"(d) Newmark step load", 0.0, {}, true}};
    for (const auto& c : cases) {
      const std::string label = c.name + " Omega=" + fmt(Om);
      GridRun run = guarded(label, steps, [&](GridRun& r) {
        if (!setup_failure.empty()) {
          r.failure = setup_failure;
          return;
        }
        const double dt = Om / w;
        HhtIntegrator ig(model, dt, c.alpha, c.ray);
        const Vec zero(model.num_dofs(), 0.0);
        Vec F = zero;
        const double qs = 1.0;  // 계단 하중 크기: F = ω² q_s M φ
        if (c.step_load) {
          for (std::size_t i = 0; i < F.size(); ++i) F[i] = w * w * qs * Mphi[i];
          ig.initialize(zero, zero, F, zero);
        } else {
          ig.initialize(phi, zero, zero, zero);
        }
        std::vector<double> q{model.dot(Mphi, ig.u())};
        double leak = 0.0, unorm = 0.0;
        Vec Mu;
        for (int s = 0; s < steps; ++s) {
          ig.step(F, zero, zero, zero);
          const double qn = model.dot(Mphi, ig.u());
          q.push_back(qn);
          Vec d = ig.u();
          for (std::size_t i = 0; i < d.size(); ++i) d[i] -= qn * phi[i];
          model.apply_M(d, Mu);
          leak = std::max(leak, std::sqrt(std::max(0.0, model.dot(d, Mu))));
          model.apply_M(ig.u(), Mu);
          unorm = std::max(unorm, std::sqrt(model.dot(ig.u(), Mu)));
        }
        check_residual(r, ig.max_rel_residual());
        r.metrics["leakage"] = leak / unorm;
        const double xi = c.ray.a0 / (2 * w) + c.ray.a1 * w / 2;
        double rho = 0, theta = 0, l3 = 0;
        hht_dominant_root(Om, c.alpha, xi, rho, theta, &l3);
        r.metrics["theory_rho"] = rho;
        r.metrics["theory_theta"] = theta;
        r.metrics["spurious_root"] = l3;
        if (c.step_load) {
          double err = 0.0, peak = 0.0;
          for (int n = 0; n <= steps; ++n) {
            err = std::max(err, std::abs(q[n] - qs * (1 - std::cos(n * theta))) / qs);
            peak = std::max(peak, q[n] / qs);
          }
          r.metrics["step_response_error"] = err;
          r.metrics["dynamic_amplification"] = peak;
        } else {
          // e_n = q_n² − q_{n−1} q_{n+1} = C² ρ^{2n} sin²θ
          double sxy = 0, sxx = 0;
          for (int n = skip; n + 2 <= steps; ++n) {
            const double en = q[n] * q[n] - q[n - 1] * q[n + 1];
            const double en1 = q[n + 1] * q[n + 1] - q[n] * q[n + 2];
            sxy += en * en1;
            sxx += en * en;
          }
          const double rm = std::sqrt(sxy / sxx);
          double num = 0, den = 0;
          for (int n = skip; n + 1 <= steps; ++n) {
            num += q[n] * (q[n + 1] + rm * rm * q[n - 1]);
            den += q[n] * q[n];
          }
          const double tm = std::acos(std::clamp(num / (2 * rm * den), -1.0, 1.0));
          r.metrics["measured_rho"] = rm;
          r.metrics["measured_theta"] = tm;
          r.metrics["rho_error"] = std::abs(rm - rho);
          r.metrics["theta_rel_error"] = std::abs(tm - theta) / theta;
        }
      });
      rep.runs.push_back(run);
    }
  }
  if (!any_failed(rep.runs)) {
    for (const auto& r : rep.runs) {
      if (r.metrics.count("step_response_error")) {
        rep.criteria.push_back(make(r.label + ": max |q_n - q_s(1 - cos n theta)| / q_s <= 1e-8",
                                    r.metrics.at("step_response_error") <= 1e-8, fmt(r.metrics.at("step_response_error"))));
      } else {
        rep.criteria.push_back(make(r.label + ": |rho_meas - rho| <= 1e-8", r.metrics.at("rho_error") <= 1e-8,
                                    fmt(r.metrics.at("rho_error"))));
        rep.criteria.push_back(make(r.label + ": |theta_meas - theta| / theta <= 1e-8",
                                    r.metrics.at("theta_rel_error") <= 1e-8, fmt(r.metrics.at("theta_rel_error"))));
      }
      rep.criteria.push_back(make(r.label + ": leakage to other modes <= 1e-8", r.metrics.at("leakage") <= 1e-8,
                                  fmt(r.metrics.at("leakage"))));
    }
  }
  conclude(rep, quick);
  return rep;
}

// ------------------------------------------------------------------ 4-I
TestReport run_4i(bool quick) {
  TestReport rep = new_report("4-I", "inertia loads: (a) prescribed motion, (b) inertia relief");
  (void)quick;  // 정확 재현 시험이라 축소판도 같은 격자
  const Material mat;
  const double R = 0.2, L = 1.0, A = 1.0, q = 1.0;
  const double E = mat.E, nu = mat.nu, rho = mat.rho;
  const double V = kPi * R * R * L;
  auto rod = [&]() {
    MeshBuilder mb;
    mb.add_block(0, R, 0, L, 2, 10);
    return mb.build();
  };
  // (a) 지정 운동: σ_zz = ρ A z
  const Field2 ua = [&](double r, double z) {
    return std::array<double, 2>{-nu * rho * A * z * r / E, rho * A * (z * z - L * L) / (2 * E) + nu * rho * A * r * r / (2 * E)};
  };
  const Field4 sa = [&](double, double z) { return std::array<double, 4>{0, 0, rho * A * z, 0}; };
  for (bool exact_top : {true, false}) {
    rep.runs.push_back(guarded(exact_top ? "(a1) top: exact displacement" : "(a2) top: clamped u = 0", 10, [&](GridRun& run) {
      Model model(rod(), mat);
      model.fix_axis();
      model.fix_where([&](double, double z) { return z == L; }, 0);
      model.fix_where([&](double, double z) { return z == L; }, 1);
      model.finalize();
      Vec F(model.num_dofs(), 0.0);
      add_inertia_z(model, A, F);
      double res = 0;
      const Vec up = exact_top ? nodal(model.mesh(), ua) : Vec(model.num_dofs(), 0.0);
      const Vec u = solve_static(model, F, up, &res);
      check_residual(run, res);
      const Vec Rz = reactions(model, u, F);
      double sum = 0.0;
      for (std::size_t n = 0; n < model.mesh().nodes.size(); ++n) {
        if (model.mesh().nodes[n].z == L) sum += Rz[2 * n + 1];
      }
      run.metrics["reaction_sum_error"] = std::abs(sum - rho * A * V) / (rho * A * V);
      if (exact_top) {
        run.metrics["disp_max_error"] = max_nodal_error(model.mesh(), u, ua);
        run.metrics["stress_max_error"] = max_gauss_stress_error(model, u, sa, rho * A * L);
      }
    }));
  }
  // (b) 관성 해제: σ_zz = −q (1 − z/L)
  const Field2 ub = [&](double r, double z) {
    return std::array<double, 2>{nu * q * r * (L - z) / (E * L), q * (z * z - 2 * L * z + nu * r * r) / (2 * E * L)};
  };
  const Field4 sb = [&](double, double z) { return std::array<double, 4>{0, 0, -q * (1 - z / L), 0}; };
  for (bool cup : {false, true}) {
    rep.runs.push_back(guarded(cup ? "(b2) cup, bottom pressure" : "(b) rod, end pressure", 10, [&](GridRun& run) {
      Model model(cup ? cup_mesh() : rod(), mat);
      model.fix_axis();
      // 강체 이동을 막는 한 절점: 막대는 윗면 축 (0, L), 컵은 측벽 윗 끝 (1, 1)
      const double fr = cup ? 1.0 : 0.0, fz = cup ? 1.0 : L;
      int fixed_node = -1;
      for (int n = 0; n < static_cast<int>(model.mesh().nodes.size()); ++n) {
        if (model.mesh().nodes[n].r == fr && model.mesh().nodes[n].z == fz) fixed_node = n;
      }
      if (fixed_node < 0) throw std::runtime_error("support node not found");
      model.fix(fixed_node, 1);
      model.finalize();
      Vec F(model.num_dofs(), 0.0);
      add_pressure(model.mesh(), boundary_edges_where(model.mesh(), [](double, double z) { return z == 0.0; }), q, F);
      double Fz = 0.0;
      for (int n = 0; n < model.num_dofs() / 2; ++n) Fz += F[2 * n + 1];
      const double Arel = add_inertia_relief_z(model, F);
      run.metrics["relief_acceleration"] = Arel;
      double res = 0;
      const Vec u = solve_static(model, F, Vec(model.num_dofs(), 0.0), &res);
      check_residual(run, res);
      const Vec Rz = reactions(model, u, F);
      run.metrics["support_reaction_over_load"] = std::abs(Rz[2 * fixed_node + 1]) / Fz;
      if (!cup) {
        // 질량 가중 평균 z 이동을 맞춘 뒤 비교
        Vec d = u;
        const Vec ex = nodal(model.mesh(), ub);
        for (std::size_t i = 0; i < d.size(); ++i) d[i] -= ex[i];
        Vec Md;
        model.apply_M(rigid_z(model), Md);
        const double c = model.dot(Md, d) / total_mass(model);
        run.metrics["disp_max_error"] = max_nodal_error(model.mesh(), u, ub, c);
        run.metrics["stress_max_error"] = max_gauss_stress_error(model, u, sb, q);
        run.metrics["relief_acceleration_error"] = std::abs(Arel - q / (rho * L)) / (q / (rho * L));
      }
    }));
  }
  if (!any_failed(rep.runs)) {
    const auto& a1 = rep.runs[0].metrics;
    rep.criteria.push_back(make("(a1) nodal displacement error / max|u| <= 1e-10", a1.at("disp_max_error") <= kExactDisp,
                                fmt(a1.at("disp_max_error"))));
    rep.criteria.push_back(make("(a1) Gauss-point stress error / (rho A L) <= 1e-9",
                                a1.at("stress_max_error") <= kExactStress, fmt(a1.at("stress_max_error"))));
    for (int k = 0; k < 2; ++k) {
      const double v = rep.runs[k].metrics.at("reaction_sum_error");
      rep.criteria.push_back(make(rep.runs[k].label + ": support reaction sum = rho A pi R^2 L within 1e-12", v <= 1e-12, fmt(v)));
    }
    const auto& b = rep.runs[2].metrics;
    rep.criteria.push_back(make("(b) support reaction <= 1e-12 * q pi R^2", b.at("support_reaction_over_load") <= 1e-12,
                                fmt(b.at("support_reaction_over_load"))));
    rep.criteria.push_back(make("(b) displacement (mean z shift removed) error / max|u| <= 1e-10",
                                b.at("disp_max_error") <= kExactDisp, fmt(b.at("disp_max_error"))));
    rep.criteria.push_back(make("(b) Gauss-point stress error / q <= 1e-9", b.at("stress_max_error") <= kExactStress,
                                fmt(b.at("stress_max_error"))));
    const double b2 = rep.runs[3].metrics.at("support_reaction_over_load");
    rep.criteria.push_back(make("(b2) cup: support reaction <= 1e-12 * total load", b2 <= 1e-12, fmt(b2)));
  }
  conclude(rep, quick);
  return rep;
}

}  // namespace hcw::verify
