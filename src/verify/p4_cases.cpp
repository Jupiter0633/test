#include "verify/p4_cases.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "body/cut_cell.h"
#include "body/force.h"
#include "body/shape.h"
#include "interface/two_phase.h"
#include "interface/vof_material.h"
#include "verify/report_util.h"

namespace hcw::verify {

using flow::Array2D;
using flow::Geometry;
using flow::Grid;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kRhoW = 1000.0, kMuW = 1.0e-3, kRhoA = 1.2, kMuA = 1.8e-5, kG = 9.81;
constexpr double kPatm = 1.0e5, kCw = 1480.0;
constexpr double kA = 0.05, kV = 2.0;  // 기준 물체 반지름 (반폭), 입수 속도 (§1.4 공통 조건)

TestReport p4_report(const std::string& id, const std::string& title) {
  TestReport rep;
  rep.id = id;
  rep.title = title;
  rep.phase = "P4";
  rep.criteria_section = "1.4";
  rep.criteria_commit = "b8ebdf3";  // §1.4 + 2차 시간 적분·압력 하한 (60d72a8 기준 그대로)
  return rep;
}

enum class BodyKind { Disk, ClosedCylinder, Wedge };

// 한 번의 입수 계산 조건. 길이는 a 단위.
struct EntrySetup {
  BodyKind kind = BodyKind::Disk;
  double beta_deg = 0.0;  // 쐐기 데드라이즈
  int n = 32;             // a/Δ
  double h0 = 0.1;        // 처음 간격 (바닥면 또는 꼭짓점 – 수면)
  double r_max = 3.0, z_min = -3.0, z_max = 1.5;
  double depth_end = 0.2;  // 계산 끝: 바닥면(꼭짓점) 이 처음 수면 아래 이 깊이
  double dt_factor = 1.0;
  double tol_factor = 1.0;
  double kappa = 1.4;  // 0 이면 비압축 공기
  bool cv = false;     // 3-F 제어면 운동량 수지
  double cv_rs = 2.0, cv_zb = -2.0, cv_zt = 1.25;
  // 개발용 그림 자료: 비어 있지 않으면 snapshot_interval [s] 마다 셀 값을 이 폴더에 CSV 로 쓴다
  std::string snapshot_dir;
  double snapshot_interval = 0.0;
};

struct EntryOut {
  GridRun run;
  std::vector<double> t, force, p_bottom, pocket_p, wetted;
  // 개발용 진단 (history.csv): 포켓 부피, 유효 셀 최저 절대압력, 하한 발동 셀 누적, 축 위 공기층 두께 (바닥면–수면)
  std::vector<double> pocket_vol, p_min_abs, floor_cells, gap_axis;
  // 3-F: 반 단계 시각의 표면 힘과 제어면 힘
  std::vector<double> fs_half, fcv_half;
  double impulse = 0.0, impulse_cv = 0.0;
};

std::shared_ptr<body::Shape> make_shape(BodyKind k, double beta_deg) {
  switch (k) {
    case BodyKind::Disk:
      return std::make_shared<body::Cylinder>(kA, 0.2 * kA);
    case BodyKind::ClosedCylinder:
      return std::make_shared<body::Cylinder>(kA, 2.0 * kA);
    case BodyKind::Wedge:
      return std::make_shared<body::Wedge>(kA, beta_deg * kPi / 180.0);
  }
  throw std::logic_error("unknown body kind");
}

// 같은 물·공기 분포의 이산 정수압 (물체 없음, 비압축으로 풀어 평균 0). P3 2-A 와 같은 방법.
Array2D hydrostatic(flow::Problem pb, const Array2D& f) {
  pb.body = nullptr;
  pb.air_kappa = 0.0;
  auto mat = std::make_shared<interface::VofMaterial>(kRhoW, 0.0, kRhoA, 0.0);
  interface::TwoPhaseSolver s(pb, mat, f);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; });
  Array2D p;
  s.flow().compute_pressure(p);
  return p;
}

double max_speed(const flow::FlowSolver& s) {
  double m = 0.0;
  for (double v : s.u().data()) m = std::max(m, std::abs(v));
  for (double v : s.w().data()) m = std::max(m, std::abs(v));
  return m;
}

int face_index(double x, double x0, double d) {
  const double q = (x - x0) / d;
  const long k = std::lround(q);
  if (std::abs(q - static_cast<double>(k)) > 1e-9) throw std::invalid_argument("control surface is not on a grid face");
  return static_cast<int>(k);
}

// 바깥 대기(영역 윗면) 와 이어지지 않은 공기 덩어리 중, 물체 아래(r < a, z < 바닥면) 에 닿는 것 = 갇힌 포켓.
// 공기 셀: 유체 비율 φ > 0.5 이고 유체 안 물 비율 f/φ < 0.5. 4-이웃 연결.
struct Pocket {
  bool found = false;
  double volume = 0.0, pressure = 0.0;
};

// 포켓 압력은 압력이 정의된 셀 (active) 만 평균한다.
Pocket find_pocket(const Grid& g, const Array2D& f, const Array2D& phi, const Array2D& p,
                   const std::vector<char>& active, double z_bottom) {
  const int nr = g.nr, nz = g.nz;
  auto idx = [&](int i, int j) { return static_cast<std::size_t>(i) * nz + j; };
  std::vector<char> air(static_cast<std::size_t>(nr) * nz, 0);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const double ph = phi(i, j);
      air[idx(i, j)] = ph > 0.5 && f(i, j) / ph < 0.5;
    }
  }
  std::vector<int> label(air.size(), -1);
  std::vector<std::size_t> stack;
  auto flood = [&](int i0, int j0, int lab) {
    stack.clear();
    stack.push_back(idx(i0, j0));
    label[idx(i0, j0)] = lab;
    while (!stack.empty()) {
      const std::size_t k = stack.back();
      stack.pop_back();
      const int i = static_cast<int>(k / nz), j = static_cast<int>(k % nz);
      const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
      for (int q = 0; q < 4; ++q) {
        const int ii = i + di[q], jj = j + dj[q];
        if (ii < 0 || jj < 0 || ii >= nr || jj >= nz) continue;
        const std::size_t kk = idx(ii, jj);
        if (!air[kk] || label[kk] >= 0) continue;
        label[kk] = lab;
        stack.push_back(kk);
      }
    }
  };
  for (int i = 0; i < nr; ++i) {
    if (air[idx(i, nz - 1)] && label[idx(i, nz - 1)] < 0) flood(i, nz - 1, 0);  // 대기
  }
  int next = 1;
  std::vector<char> is_pocket(1, 0);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      if (!air[idx(i, j)] || label[idx(i, j)] >= 0) continue;
      flood(i, j, next);
      is_pocket.push_back(0);
      ++next;
    }
  }
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const int lab = label[idx(i, j)];
      if (lab > 0 && g.rc(i) < kA && g.zc(j) < z_bottom) is_pocket[static_cast<std::size_t>(lab)] = 1;
    }
  }
  Pocket pk;
  double sp = 0.0, sv = 0.0;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nz; ++j) {
      const int lab = label[idx(i, j)];
      if (lab <= 0 || !is_pocket[static_cast<std::size_t>(lab)]) continue;
      const double v = std::max(0.0, phi(i, j) - f(i, j)) * g.cell_volume(i);
      pk.volume += v;
      if (active[idx(i, j)]) {
        sp += p(i, j) * v;
        sv += v;
      }
    }
  }
  pk.found = pk.volume > 0.0 && sv > 0.0;
  pk.pressure = pk.found ? sp / sv : 0.0;
  return pk;
}

// 3-F: 제어면 [0, rs] × [zb, zt] 의 유체 운동량 P (z 성분) 와 표면·체적 항 S
//   S = ∮ (−p n_z + τ_zn − ρ w (u·n)) dA + ∫ ρ g_z dV,   물체 힘 F_cv = S − dP/dt.
// 셀 질량 = (f ρ_w + (φ − f) ρ_a) V (고체 제외). z 면 운동량은 이웃 두 셀 질량의 절반씩.
struct CvTerms {
  double momentum = 0.0, surface = 0.0;
};

CvTerms cv_terms(const Grid& g, const flow::FlowSolver& fs, const Array2D& f, const Array2D& phi, double rs,
                 double zb, double zt) {
  const int is = face_index(rs, g.r0, g.dr), jb = face_index(zb, g.z0, g.dz), jt = face_index(zt, g.z0, g.dz);
  const Array2D& u = fs.u();
  const Array2D& w = fs.w();
  const Array2D& p = fs.pressure();
  auto mass = [&](int i, int j) {
    return (f(i, j) * kRhoW + std::max(0.0, phi(i, j) - f(i, j)) * kRhoA) * g.cell_volume(i);
  };
  auto fe = [&](int i, int j) { return phi(i, j) > 1e-12 ? std::clamp(f(i, j) / phi(i, j), 0.0, 1.0) : 0.0; };
  auto rho = [&](int i, int j) { return fe(i, j) * kRhoW + (1.0 - fe(i, j)) * kRhoA; };
  auto mu = [&](int i, int j) { return fe(i, j) * kMuW + (1.0 - fe(i, j)) * kMuA; };
  CvTerms c;
  double grav = 0.0;
  for (int i = 0; i < is; ++i) {
    for (int j = jb; j <= jt; ++j) {
      double m = 0.0;
      if (j - 1 >= jb) m += 0.5 * mass(i, j - 1);
      if (j < jt) m += 0.5 * mass(i, j);
      c.momentum += m * w(i, j);
    }
    for (int j = jb; j < jt; ++j) grav -= mass(i, j) * kG;
  }
  double s = 0.0;
  for (int i = 0; i < is; ++i) {
    for (int sgn : {+1, -1}) {
      const int j = sgn > 0 ? jt : jb;
      const double pf = 0.5 * (p(i, j - 1) + p(i, j));
      const double rf = 0.5 * (rho(i, j - 1) + rho(i, j)), mf = 0.5 * (mu(i, j - 1) + mu(i, j));
      const double tzz = 2.0 * mf * (w(i, j + 1) - w(i, j - 1)) / (2.0 * g.dz);
      s += sgn * (-pf + tzz - rf * w(i, j) * w(i, j)) * g.metric(g.rc(i)) * g.dr;
    }
  }
  for (int j = jb; j < jt; ++j) {
    const double rf = 0.5 * (rho(is - 1, j) + rho(is, j)), mf = 0.5 * (mu(is - 1, j) + mu(is, j));
    const double dudz = (u(is, j + 1) - u(is, j - 1)) / (2.0 * g.dz);
    const double dwdr = 0.5 * ((w(is, j) - w(is - 1, j)) + (w(is, j + 1) - w(is - 1, j + 1))) / g.dr;
    const double wc = 0.25 * (w(is - 1, j) + w(is, j) + w(is - 1, j + 1) + w(is, j + 1));
    s += (mf * (dudz + dwdr) - rf * wc * u(is, j)) * g.metric(g.rf(is)) * g.dz;
  }
  const double scale = g.geometry == Geometry::Axisymmetric ? 2.0 * kPi : 1.0;
  c.momentum *= scale;
  c.surface = (s + grav) * scale;
  return c;
}

// 국소 극대의 개수 (양 이웃보다 큰 표본)
int count_local_maxima(const std::vector<double>& x) {
  int c = 0;
  for (std::size_t k = 1; k + 1 < x.size(); ++k) {
    if (x[k] > x[k - 1] && x[k] > x[k + 1]) ++c;
  }
  return c;
}

EntryOut run_entry_impl(const EntrySetup& es) {
  EntryOut out;
  GridRun& run = out.run;
  run.n = es.n;
  const auto t0 = std::chrono::steady_clock::now();
  const bool planar = es.kind == BodyKind::Wedge;
  const auto shape = make_shape(es.kind, es.beta_deg);
  const double hh = shape->half_height();
  const double zc0 = es.h0 * kA + hh;
  const double depth_end = es.kind == BodyKind::Wedge ? es.depth_end * kA * std::tan(es.beta_deg * kPi / 180.0)
                                                      : es.depth_end * kA;
  const double t_end = (es.h0 * kA + depth_end) / kV;

  flow::Problem pb;
  const int nr = static_cast<int>(std::lround(es.r_max * es.n));
  const int nz = static_cast<int>(std::lround((es.z_max - es.z_min) * es.n));
  pb.grid = Grid(planar ? Geometry::Planar : Geometry::Axisymmetric, nr, nz, 0.0, es.r_max * kA, es.z_min * kA,
                 es.z_max * kA);
  pb.bc.r_min.kind = flow::BoundaryKind::Axis;  // 축대칭 축, 평면에서는 대칭면
  pb.gravity = {0.0, -kG};
  const Grid g = pb.grid;

  // 처음 물: z < 0 (셀마다 정확한 비율), 고체 부분은 0
  Array2D f0(g.nr, g.nz);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) f0(i, j) = std::clamp((0.0 - g.zf(j)) / g.dz, 0.0, 1.0);
  }
  const Array2D p_hyd = hydrostatic(pb, f0);  // 물체는 공기 중 (고체를 공기로 보고 정수압)
  const Array2D phi0 = body::fluid_fraction(g, *shape, zc0);
  for (int i = 0; i < g.nr; ++i) {
    for (int j = 0; j < g.nz; ++j) f0(i, j) = std::min(f0(i, j), phi0(i, j));
  }

  auto bd = std::make_shared<flow::Body>();
  bd->shape = shape;
  bd->motion.z = [zc0](double t) { return zc0 - kV * t; };
  bd->motion.w = [](double) { return -kV; };
  bd->motion.a = [](double) { return 0.0; };
  bd->interpolation_order = 1;
  pb.body = bd;
  pb.air_kappa = es.kappa;

  flow::SolverOptions opt;
  opt.convection = flow::Convection::Consistent;
  opt.time_order = 2;            // 03 §1.4 2차 시간 적분 (b8ebdf3)
  opt.air_pressure_floor = 1e3;  // 03 §1.4 공기 압력 하한
  opt.projection_tolerance *= es.tol_factor;
  opt.pressure_tolerance *= es.tol_factor;
  auto mat = std::make_shared<interface::VofMaterial>(kRhoW, kMuW, kRhoA, kMuA);
  interface::TwoPhaseSolver s(pb, mat, f0, true, opt);
  s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; }, false);
  Array2D p0 = p_hyd;
  if (es.kappa > 0.0) {
    for (double& v : p0.data()) v += kPatm;
  }
  s.flow().set_pressure(p0);
  const double water0 = s.vof().water_volume();

  const double h = std::min(g.dr, g.dz);
  const double rho_c_v = kRhoW * kCw * kV;
  double peak_force = 0.0, peak_bottom = -1e300, peak_s1 = -1e300, peak_s2 = -1e300, peak_point = -1e300;
  double pocket_max_volume = 0.0, side_peak = -1e300;
  bool pocket_seen = false;
  long retries = 0;
  double next_snapshot = 0.0;
  int snapshot_count = 0;
  double f_prev = 0.0, t_prev = 0.0;
  CvTerms cv_prev;
  bool have_prev = false;

  auto gauge_ref = [&]() {
    // 압축성: 대기압. 비압축: 영역 오른쪽 위 셀 (먼 공기) 의 압력
    return es.kappa > 0.0 ? kPatm : s.flow().pressure()(g.nr - 1, g.nz - 1);
  };
  auto measure = [&]() {
    const double t = s.time();
    const double zc = zc0 - kV * t;
    // 압력이 정의된 셀 = 이번 단계 압력 미지수 (P2 와 같이 pressure_geometry; 물체 옆에서 모든 면을
    // 운동량이 정한 셀은 주변 유효 셀의 최소제곱 1차식으로 외삽된다)
    const body::CutGeometry& geom = s.flow().pressure_geometry();
    const Array2D& p = s.flow().pressure();
    if (geom.pieces.empty()) {  // 첫 단계 전: 물체는 공기 중, 대기압 → 모든 값 0
      out.t.push_back(t);
      out.force.push_back(0.0);
      out.p_bottom.push_back(0.0);
      out.pocket_p.push_back(0.0);
      out.wetted.push_back(0.0);
      out.pocket_vol.push_back(0.0);
      out.p_min_abs.push_back(0.0);
      out.floor_cells.push_back(0.0);
      out.gap_axis.push_back(0.0);
      if (es.cv) cv_prev = cv_terms(g, s.flow(), s.vof().f(), s.fluid_fraction(t), es.cv_rs * kA, es.cv_zb * kA, es.cv_zt * kA);
      f_prev = 0.0;
      t_prev = t;
      have_prev = true;
      return;
    }
    const double ref = gauge_ref();
    // 힘은 대기압 기준 (게이지) 압력으로 적분한다. 영역 밖으로 나간 물체 부분 (3-I 원통 윗면) 에는
    // 대기압이 걸린다고 보는 것과 같다. 물체 전체가 영역 안이면 결과는 같다 (닫힌 표면의 상수 압력 합력 0).
    Array2D p_gauge = p;
    for (double& v : p_gauge.data()) v -= ref;
    const double force = body::axial_pressure_force(g, geom, p_gauge) +
                         body::axial_viscous_force(g, geom, *shape, zc, s.flow().w(), kMuW, -kV);
    auto bottom = [&](double rmax) {
      return [=](const body::BoundaryPiece& pc) { return pc.nz_area > 0.0 && pc.z < zc && pc.r <= rmax; };
    };
    const double pb_all = body::average_surface_pressure(g, geom, p, bottom(kA * (1.0 + 1e-9))) - ref;
    const double ps1 = body::average_surface_pressure(g, geom, p, bottom(0.1 * kA)) - ref;
    const double ps2 = body::average_surface_pressure(g, geom, p, bottom(0.2 * kA)) - ref;
    const double ppt = body::max_surface_pressure(g, geom, p) - ref;
    const Array2D phi = s.fluid_fraction(t);
    const Pocket pk = find_pocket(g, s.vof().f(), phi, p, geom.active, zc - hh);
    if (!es.snapshot_dir.empty() && t >= next_snapshot - 1e-12) {
      // 셀마다: 중심 r, z, 유체 비율 φ, 유체 안 물 비율, 게이지 압력, 셀 중심 속도 (u, w)
      char name[64];
      std::snprintf(name, sizeof(name), "snap_%03d.csv", snapshot_count++);
      std::ofstream o(std::filesystem::path(es.snapshot_dir) / name);
      o << "# t=" << t << " zc=" << zc << " half_height=" << hh << " a=" << kA << " planar=" << (planar ? 1 : 0)
        << " nr=" << g.nr << " nz=" << g.nz << " dr=" << g.dr << " dz=" << g.dz << " z0=" << g.z0 << "\n";
      o << "r,z,phi,water,p_gauge,u,w,active\n";
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) {
          const double ph = phi(i, j);
          const double fw = ph > 1e-12 ? std::clamp(s.vof().f()(i, j) / ph, 0.0, 1.0) : 0.0;
          const bool act = geom.active[static_cast<std::size_t>(i) * g.nz + j] != 0;
          const double uc = 0.5 * (s.flow().u()(i, j) + s.flow().u()(i + 1, j));
          const double wc = 0.5 * (s.flow().w()(i, j) + s.flow().w()(i, j + 1));
          o << g.rc(i) << "," << g.zc(j) << "," << ph << "," << fw << "," << (act ? p_gauge(i, j) : 0.0) << "," << uc
            << "," << wc << "," << (act ? 1 : 0) << "\n";
        }
      }
      next_snapshot += es.snapshot_interval;
    }
    // 젖은 반폭 (3-A): 바닥면 조각 중, 그 셀 또는 바로 아래 셀의 유체 안 물 비율이 0.5 를 넘는 가장 바깥 r
    double wet = 0.0;
    auto wet_cell = [&](int i, int j) {
      if (j < 0) return false;
      const double ph = phi(i, j);
      return ph > 1e-12 && s.vof().f()(i, j) / ph > 0.5;
    };
    for (const auto& pc : geom.pieces) {
      if (!(pc.nz_area > 0.0)) continue;
      if (wet_cell(pc.i, pc.j) || wet_cell(pc.i, pc.j - 1)) wet = std::max(wet, pc.r);
    }
    if (es.kind == BodyKind::ClosedCylinder) {
      // 옆벽: 법선이 수평인 조각 (n_z 면적 0) 의 압력 평균 — 3-I 보고
      // 옆벽 r ≈ a 의 조각 (법선이 수평, n_z 면적 0) 의 압력 평균 — 3-I 보고
      const double side = body::average_surface_pressure_by_length(g, geom, p, [&](const body::BoundaryPiece& pc) {
        return pc.nz_area == 0.0 && std::abs(pc.r - kA) < 0.5 * g.dr;
      }) - ref;
      side_peak = std::max(side_peak, side);
    }
    out.t.push_back(t);
    out.force.push_back(force);
    out.p_bottom.push_back(pb_all);
    out.pocket_p.push_back(pk.found ? pk.pressure - ref : 0.0);
    out.wetted.push_back(wet);
    {
      double pmin = 1e300;
      for (std::size_t k = 0; k < geom.active.size(); ++k) {
        if (geom.active[k]) pmin = std::min(pmin, p.data()[k]);
      }
      // 축 옆 첫 열에서 바닥면 바로 아래부터 내려가며 물 비율이 처음 0.5 를 넘는 셀의 윗면 = 수면
      const double zb = zc - hh;
      double gap = 0.0;
      for (int j = g.nz - 1; j >= 0; --j) {
        if (g.zc(j) >= zb) continue;
        const double ph = phi(0, j);
        if (ph > 0.5 && s.vof().f()(0, j) / ph > 0.5) {
          gap = zb - (g.zf(j) + s.vof().f()(0, j) / ph * g.dz);
          break;
        }
      }
      out.pocket_vol.push_back(pk.found ? pk.volume : 0.0);
      out.p_min_abs.push_back(pmin + (es.kappa > 0.0 ? 0.0 : kPatm - ref));
      out.floor_cells.push_back(static_cast<double>(s.flow().stats().pressure_floor_cells));
      out.gap_axis.push_back(gap);
    }
    if (pk.found) {
      pocket_seen = true;
      pocket_max_volume = std::max(pocket_max_volume, pk.volume);
    }
    peak_force = std::max(peak_force, std::abs(force));
    peak_bottom = std::max(peak_bottom, pb_all);
    peak_s1 = std::max(peak_s1, ps1);
    peak_s2 = std::max(peak_s2, ps2);
    peak_point = std::max(peak_point, ppt);
    if (have_prev) {
      const double dt = t - t_prev;
      out.impulse += 0.5 * (f_prev + force) * dt;
      if (es.cv) {
        const CvTerms c = cv_terms(g, s.flow(), s.vof().f(), phi, es.cv_rs * kA, es.cv_zb * kA, es.cv_zt * kA);
        const double fcv = 0.5 * (cv_prev.surface + c.surface) - (c.momentum - cv_prev.momentum) / dt;
        out.fs_half.push_back(0.5 * (f_prev + force));
        out.fcv_half.push_back(fcv);
        out.impulse_cv += fcv * dt;
        cv_prev = c;
      }
    } else if (es.cv) {
      cv_prev = cv_terms(g, s.flow(), s.vof().f(), phi, es.cv_rs * kA, es.cv_zb * kA, es.cv_zt * kA);
    }
    f_prev = force;
    t_prev = t;
    have_prev = true;
  };

  try {
    measure();
    while (s.time() < t_end - 1e-15) {
      const double vm = max_speed(s.flow());
      double dt = 0.1 * h / kV;
      if (vm > 0.0) dt = std::min(dt, 0.25 * h / vm);
      dt = std::min(dt * es.dt_factor, t_end - s.time());
      // 한 단계 안에서 속도가 급히 커져 VOF 이동 한계 (CFL 0.5) 를 넘거나 공기 압력이 0 이하가 되면
      // 그 단계를 되돌리고 Δt 를 절반으로 다시 계산한다 (최대 8 번). 03 §1.4 시간 간격 보충 (2026-10-23).
      for (int halvings = 0;; ++halvings) {
        const interface::TwoPhaseSolver backup = s;
        const interface::VofMaterial mat_backup = *mat;
        try {
          s.step(dt);
          break;
        } catch (const std::exception& e) {
          const std::string msg = e.what();
          const bool retry = msg.find("CFL") != std::string::npos ||
                             msg.find("non-positive absolute pressure") != std::string::npos;
          if (!retry || halvings >= 8) throw;
          s = backup;
          *mat = mat_backup;
          dt *= 0.5;
          ++retries;
        }
      }
      measure();
      if (s.flow().stats().solver_failed || !std::isfinite(out.force.back())) break;
    }
    run.metrics["end_time_over_t_end"] = s.time() / t_end;
    run.metrics["impulse"] = out.impulse;
    run.metrics["peak_force"] = peak_force;
    run.metrics["peak_bottom_avg_pressure"] = peak_bottom;
    run.metrics["peak_sensor_0.1a_pressure"] = peak_s1;
    run.metrics["peak_sensor_0.2a_pressure"] = peak_s2;
    run.metrics["peak_point_pressure"] = peak_point;
    run.metrics["peak_point_over_rho_c_V"] = peak_point / rho_c_v;
    run.metrics["pocket_found"] = pocket_seen ? 1.0 : 0.0;
    run.metrics["pocket_max_volume_over_a3"] =
        pocket_max_volume * (planar ? 1.0 : 2.0 * kPi) / (planar ? kA * kA : kA * kA * kA);
    run.metrics["pocket_pressure_local_maxima"] = static_cast<double>(count_local_maxima(out.pocket_p));
    run.metrics["clipped_water_over_initial"] = s.clipped_water() / water0;
    run.metrics["steps"] = static_cast<double>(s.flow().stats().steps);
    run.metrics["step_retries"] = static_cast<double>(retries);
    run.metrics["pressure_floor_cells"] = static_cast<double>(s.flow().stats().pressure_floor_cells);
    if (es.kind == BodyKind::ClosedCylinder) run.metrics["peak_side_wall_pressure"] = side_peak;
    run.seconds = seconds_since(t0);
    finish_run(run, s.flow());
    if (s.time() < t_end - 1e-12 && run.failure.empty()) run.failure = "FAILED_UNSTABLE";
  } catch (const std::exception& e) {
    run.failure = std::string("FAILED_EXCEPTION ") + e.what();
  }
  return out;
}

// 준비 단계 (정수압, 물체 설정) 의 예외도 실행 실패로 기록한다
EntryOut run_entry(const EntrySetup& es) {
  try {
    return run_entry_impl(es);
  } catch (const std::exception& e) {
    EntryOut o;
    o.run.n = es.n;
    o.run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    return o;
  }
}

// 3-D: 한 점 최대 압력과 ρ_w c_w V 비교 (모든 입수 실행의 보고서에 넣는다)
void note_3d(TestReport& rep, const GridRun& r) {
  if (!r.failure.empty() || !r.metrics.count("peak_point_over_rho_c_V")) return;
  const double x = r.metrics.at("peak_point_over_rho_c_V");
  rep.notes.push_back("3-D " + r.label + ": peak point pressure / (rho_w c_w V) = " + fmt(x) +
                      (x > 1.0 ? "  -> OUT_OF_SCOPE (water compressibility not modelled)" : "  (within)"));
}

// 시간 이력 표본 (보고용, 10 개 남짓)
std::string history(const EntryOut& o, const std::vector<double>& y, double scale) {
  std::string s;
  const std::size_t n = o.t.size();
  const std::size_t stride = std::max<std::size_t>(1, n / 10);
  for (std::size_t k = 0; k < n; k += stride) {
    s += (s.empty() ? "" : ", ") + fmt(o.t[k] * 1e3) + "ms:" + fmt(y[k] * scale);
  }
  return s;
}

std::string label_of(const EntrySetup& es) { return "a/dx=" + std::to_string(es.n); }

bool runs_ok(const std::vector<GridRun>& runs) {
  return std::all_of(runs.begin(), runs.end(), [](const GridRun& r) { return r.failure.empty(); });
}

}  // namespace

bool gci_fine(double phi1, double phi2, double phi3, double& gci, double& p) {
  const double e21 = phi2 - phi1, e32 = phi3 - phi2;
  if (e32 == 0.0) return false;
  const double ratio = e21 / e32;
  if (!(ratio > 0.0 && ratio < 1.0)) return false;
  p = std::clamp(std::log(e32 / e21) / std::log(2.0), 0.5, 3.0);
  gci = 1.25 * std::abs(e21 / phi1) / (std::pow(2.0, p) - 1.0);
  return true;
}

// ---------------------------------------------------------------- 개발용
TestReport entry_dev_study(int n, double dt_factor) {
  TestReport rep = p4_report("dev-dt", "Disk entry time-step study (development only, not a verdict)");
  EntrySetup es;
  es.n = n;
  es.dt_factor = dt_factor;
  es.cv = true;
  EntryOut o = run_entry(es);
  o.run.label = label_of(es) + " dt x " + fmt(dt_factor);
  if (o.run.failure.empty()) {
    double dmax = 0.0, fpeak = 0.0;
    for (std::size_t k = 0; k < o.fs_half.size(); ++k) {
      dmax = std::max(dmax, std::abs(o.fs_half[k] - o.fcv_half[k]));
      fpeak = std::max(fpeak, std::abs(o.fs_half[k]));
    }
    o.run.metrics["max_force_difference_over_peak"] = fpeak > 0.0 ? dmax / fpeak : 0.0;
    o.run.metrics["impulse_cv"] = o.impulse_cv;
    o.run.metrics["impulse_difference_relative"] =
        o.impulse != 0.0 ? std::abs(o.impulse - o.impulse_cv) / std::abs(o.impulse) : 0.0;
    rep.notes.push_back("bottom area-averaged pressure (kPa): " + history(o, o.p_bottom, 1e-3));
    rep.notes.push_back("pocket mean pressure (kPa): " + history(o, o.pocket_p, 1e-3));
    rep.notes.push_back("force F_s (N): " + history(o, o.force, 1.0));
    std::string hs;
    const std::size_t stride = std::max<std::size_t>(1, o.fs_half.size() / 12);
    for (std::size_t k = 0; k < o.fs_half.size(); k += stride) {
      hs += (hs.empty() ? "" : ", ") + fmt(o.t[k + 1] * 1e3) + "ms:" + fmt(o.fs_half[k]) + "/" + fmt(o.fcv_half[k]);
    }
    rep.notes.push_back("F_s / F_cv (N) at half steps: " + hs);
  }
  rep.runs.push_back(o.run);
  rep.verdict = o.run.failure.empty() ? "DEV_ONLY" : "DEV_ONLY (" + o.run.failure.substr(0, o.run.failure.find(' ')) + ")";
  return rep;
}

TestReport entry_dev_snapshots(int n, double dt_factor, const std::string& dir, double interval) {
  TestReport rep = p4_report("dev-snap", "Disk entry snapshots for figures (development only, not a verdict)");
  std::filesystem::create_directories(dir);
  EntrySetup es;
  es.n = n;
  es.dt_factor = dt_factor;
  es.snapshot_dir = dir;
  es.snapshot_interval = interval;
  EntryOut o = run_entry(es);
  o.run.label = label_of(es) + " dt x " + fmt(dt_factor);
  std::ofstream h(std::filesystem::path(dir) / "history.csv");
  h << "t,force,p_bottom,p_pocket,pocket_vol,p_min_abs,floor_cells,gap_axis\n";
  for (std::size_t k = 0; k < o.t.size(); ++k) {
    h << o.t[k] << "," << o.force[k] << "," << o.p_bottom[k] << "," << o.pocket_p[k] << "," << o.pocket_vol[k] << ","
      << o.p_min_abs[k] << "," << o.floor_cells[k] << "," << o.gap_axis[k] << "\n";
  }
  rep.runs.push_back(o.run);
  rep.notes.push_back("snapshots and history.csv written to " + dir);
  rep.verdict = o.run.failure.empty() ? "DEV_ONLY" : "DEV_ONLY (" + o.run.failure.substr(0, o.run.failure.find(' ')) + ")";
  return rep;
}

// ---------------------------------------------------------------- 3-A
TestReport run_3a(bool quick) {
  TestReport rep = p4_report("3-A", "Planar wedge entry at constant speed (report only)");
  for (double beta : {10.0, 20.0, 30.0}) {
    EntrySetup es;
    es.kind = BodyKind::Wedge;
    es.beta_deg = beta;
    es.n = quick ? 16 : 32;  // a/Δ = 8 은 쐐기 높이가 1–3 칸이라 빠른 점검으로도 너무 거칠다
    es.depth_end = 0.5;  // 꼭짓점이 쐐기 높이의 절반 깊이까지
    EntryOut o = run_entry(es);
    o.run.label = "deadrise " + fmt(beta) + " deg, " + label_of(es);
    if (o.run.failure.empty()) {
      std::vector<double> cf(o.t.size()), cp(o.p_bottom.size());
      for (std::size_t k = 0; k < o.t.size(); ++k) {
        cf[k] = o.wetted[k] > 0.0 ? o.force[k] / (kRhoW * kV * kV * o.wetted[k]) : 0.0;
      }
      o.run.metrics["peak_pressure_coefficient"] = o.run.metrics.at("peak_point_pressure") / (0.5 * kRhoW * kV * kV);
      rep.notes.push_back(o.run.label + " force coefficient F/(rho V^2 c) (c: wetted half-width): " +
                          history(o, cf, 1.0));
      rep.notes.push_back(o.run.label + " wetted half-width / a: " + history(o, o.wetted, 1.0 / kA));
      note_3d(rep, o.run);
    }
    rep.runs.push_back(o.run);
  }
  rep.notes.push_back("report only (03 section 1.4): Wagner [R-02] and Zhao & Faltinsen [R-03] not verified against "
                      "the original, so no reference values are shown");
  rep.notes.push_back("planar half wedge (symmetry plane at x = 0), half-width a = 0.05 m, V = 2 m/s, apex gap 0.1a, "
                      "end when the apex is half the wedge height below the initial surface; force per unit depth");
  rep.verdict = quick ? "QUICK_CHECK_ONLY" : "REPORT_ONLY";
  for (const auto& r : rep.runs) {
    if (!r.failure.empty()) rep.verdict += " (" + r.failure.substr(0, r.failure.find(' ')) + ")";
  }
  return rep;
}

// ---------------------------------------------------------------- 3-B
TestReport run_3b(bool quick) {
  TestReport rep = p4_report("3-B", "Flat disk entry with air cushioning (qualitative)");
  EntrySetup es;
  es.n = quick ? 8 : 32;
  EntryOut o = run_entry(es);
  o.run.label = label_of(es);
  if (o.run.failure.empty()) {
    rep.notes.push_back("bottom area-averaged pressure (kPa): " + history(o, o.p_bottom, 1e-3));
    rep.notes.push_back("pocket mean pressure (kPa, 0 = no pocket): " + history(o, o.pocket_p, 1e-3));
    rep.notes.push_back("force (N): " + history(o, o.force, 1.0));
    note_3d(rep, o.run);
  }
  rep.runs.push_back(o.run);
  if (runs_ok(rep.runs)) {
    const auto& m = o.run.metrics;
    Criterion c1;
    c1.text = "an air pocket is trapped under the bottom at some time";
    c1.passed = m.at("pocket_found") > 0.5;
    c1.detail = "max pocket volume / a^3 = " + fmt(m.at("pocket_max_volume_over_a3"));
    rep.criteria.push_back(c1);
    Criterion c2;
    c2.text = "pocket mean pressure has at least two local maxima (oscillates)";
    c2.passed = m.at("pocket_pressure_local_maxima") >= 2.0;
    c2.detail = fmt(m.at("pocket_pressure_local_maxima")) + " local maxima";
    rep.criteria.push_back(c2);
    Criterion c3;
    c3.text = "peak point pressure < rho_w c_w V";
    c3.passed = m.at("peak_point_over_rho_c_V") < 1.0;
    c3.detail = "ratio " + fmt(m.at("peak_point_over_rho_c_V"));
    rep.criteria.push_back(c3);
  }
  rep.notes.push_back("disk radius a = 0.05 m, thickness 0.2a, V = 2 m/s, gap 0.1a, isentropic air kappa = 1.4; "
                      "qualitative basis [R-05, R-06] (originals not verified)");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 3-C
TestReport run_3c(bool quick) {
  TestReport rep = p4_report("3-C", "Numerical independence of disk entry (grid GCI judged, others reported)");
  const std::vector<int> ns = quick ? std::vector<int>{4, 8, 16} : std::vector<int>{8, 16, 32, 64};
  std::vector<GridRun> grid_runs;
  for (int n : ns) {
    EntrySetup es;
    es.n = n;
    EntryOut o = run_entry(es);
    o.run.label = label_of(es);
    note_3d(rep, o.run);
    grid_runs.push_back(o.run);
  }
  rep.runs = grid_runs;
  if (runs_ok(grid_runs)) {
    const std::size_t k = grid_runs.size();
    for (const char* key : {"impulse", "peak_bottom_avg_pressure", "peak_force", "peak_point_pressure"}) {
      const double p1 = grid_runs[k - 1].metrics.at(key), p2 = grid_runs[k - 2].metrics.at(key),
                   p3 = grid_runs[k - 3].metrics.at(key);
      double gci = 0.0, p = 0.0;
      const bool mono = gci_fine(p1, p2, p3, gci, p);
      const std::string key_s = key;
      const bool judged = key_s == "impulse" || key_s == "peak_bottom_avg_pressure";
      const double limit = key_s == "impulse" ? 0.05 : 0.10;
      const std::string detail = mono ? "GCI_fine = " + fmt(gci) + ", p = " + fmt(p)
                                      : "not monotone: " + fmt(p3) + ", " + fmt(p2) + ", " + fmt(p1);
      if (judged) {
        Criterion c;
        c.text = key_s + ": monotone on the three finest grids and GCI_fine <= " + fmt(limit);
        c.passed = mono && gci <= limit;
        c.detail = detail;
        rep.criteria.push_back(c);
      } else {
        rep.notes.push_back(key_s + " (report only): " + detail);
      }
    }
  }
  // 보고만: a/Δ = 32 (quick 8) 기준 대비 변화율
  const int nb = quick ? 8 : 32;
  const GridRun* base = nullptr;
  for (const auto& r : grid_runs) {
    if (r.n == nb && r.failure.empty()) base = &r;
  }
  struct Variant {
    std::string name;
    std::function<void(EntrySetup&)> set;
  };
  const std::vector<Variant> variants = {
      {"dt x 0.5", [](EntrySetup& e) { e.dt_factor = 0.5; }},
      {"gap 0.05a", [](EntrySetup& e) { e.h0 = 0.05; }},
      {"gap 0.2a", [](EntrySetup& e) { e.h0 = 0.2; }},
      {"domain r<=4a, z>=-4a", [](EntrySetup& e) { e.r_max = 4.0; e.z_min = -4.0; }},
      {"Poisson tolerance x 10", [](EntrySetup& e) { e.tol_factor = 10.0; }},
  };
  for (const auto& v : variants) {
    EntrySetup es;
    es.n = nb;
    v.set(es);
    EntryOut o = run_entry(es);
    if (!o.run.failure.empty()) {
      rep.notes.push_back("variant " + v.name + ": " + o.run.failure);
      continue;
    }
    note_3d(rep, o.run);
    if (!base) continue;
    std::string s;
    for (const char* key : {"impulse", "peak_bottom_avg_pressure", "peak_force"}) {
      const double b = base->metrics.at(key), x = o.run.metrics.at(key);
      s += std::string(s.empty() ? "" : ", ") + key + " " + fmt((x - b) / b);
    }
    rep.notes.push_back("variant " + v.name + " (relative change vs a/dx=" + std::to_string(nb) + "): " + s);
  }
  rep.notes.push_back("GCI per Celik et al. [R-25]: safety factor 1.25, ratio 2, observed order clipped to [0.5, 3]; "
                      "peak force and peak point pressure GCI are report only");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 3-E
TestReport run_3e(bool quick) {
  TestReport rep = p4_report("3-E", "Air model sensitivity for disk entry (report only)");
  for (double kappa : {1.4, 1.0, 0.0}) {
    EntrySetup es;
    es.n = quick ? 8 : 32;
    es.kappa = kappa;
    EntryOut o = run_entry(es);
    o.run.label = kappa > 0.0 ? "isentropic kappa=" + fmt(kappa) : std::string("incompressible air");
    if (kappa == 1.0) o.run.label = "isothermal kappa=1";
    note_3d(rep, o.run);
    rep.runs.push_back(o.run);
  }
  rep.notes.push_back("3-B conditions; incompressible air pressures are relative to the top-right (far air) cell");
  rep.verdict = quick ? "QUICK_CHECK_ONLY" : "REPORT_ONLY";
  for (const auto& r : rep.runs) {
    if (!r.failure.empty()) rep.verdict += " (" + r.failure.substr(0, r.failure.find(' ')) + ")";
  }
  return rep;
}

// ---------------------------------------------------------------- 3-F
TestReport run_3f(bool quick) {
  TestReport rep = p4_report("3-F", "Momentum balance: surface force vs control-volume force");
  rep.criteria_commit = "eecd38c";  // 제어면 윗면 1.25a (변경 기록 03 §1.4)
  const std::vector<int> ns = quick ? std::vector<int>{8, 16} : std::vector<int>{32, 64};
  for (int n : ns) {
    EntrySetup es;
    es.n = n;
    es.cv = true;
    EntryOut o = run_entry(es);
    o.run.label = label_of(es);
    if (o.run.failure.empty()) {
      double dmax = 0.0, fpeak = 0.0;
      for (std::size_t k = 0; k < o.fs_half.size(); ++k) {
        dmax = std::max(dmax, std::abs(o.fs_half[k] - o.fcv_half[k]));
        fpeak = std::max(fpeak, std::abs(o.fs_half[k]));
      }
      o.run.metrics["max_force_difference_over_peak"] = fpeak > 0.0 ? dmax / fpeak : 0.0;
      o.run.metrics["impulse_cv"] = o.impulse_cv;
      std::string hs;
      const std::size_t stride = std::max<std::size_t>(1, o.fs_half.size() / 12);
      for (std::size_t k = 0; k < o.fs_half.size(); k += stride) {
        hs += (hs.empty() ? "" : ", ") + fmt(o.t[k + 1] * 1e3) + "ms:" + fmt(o.fs_half[k]) + "/" + fmt(o.fcv_half[k]);
      }
      rep.notes.push_back(o.run.label + " F_s / F_cv (N) at half steps: " + hs);
      o.run.metrics["impulse_difference_relative"] =
          o.impulse != 0.0 ? std::abs(o.impulse - o.impulse_cv) / std::abs(o.impulse) : 0.0;
    }
    rep.runs.push_back(o.run);
  }
  if (runs_ok(rep.runs)) {
    Criterion a;
    a.text = "max |F_s - F_cv| <= 1 % of peak |F_s| on every grid";
    a.passed = true;
    Criterion b;
    b.text = "|I_s - I_cv| <= 0.5 % of |I_s| on every grid";
    b.passed = true;
    for (const auto& r : rep.runs) {
      const double x = r.metrics.at("max_force_difference_over_peak"), y = r.metrics.at("impulse_difference_relative");
      a.passed = a.passed && x <= 0.01;
      b.passed = b.passed && y <= 0.005;
      a.detail += (a.detail.empty() ? "" : ", ") + fmt(x);
      b.detail += (b.detail.empty() ? "" : ", ") + fmt(y);
    }
    rep.criteria.push_back(a);
    rep.criteria.push_back(b);
  }
  rep.notes.push_back("control surface r = 2a, z = -2a .. 1.25a; F_cv = (stress - momentum flux + gravity on the "
                      "fluid) - dP/dt at half steps; cell mass from the water/air volumes outside the body");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 3-G
TestReport run_3g(bool quick) {
  TestReport rep = p4_report("3-G", "Static buoyancy of a submerged closed cylinder");
  const std::vector<int> ns = quick ? std::vector<int>{4, 8} : std::vector<int>{16, 32};
  const int steps = quick ? 20 : 200;
  const double buoy = kRhoW * kG * kPi * kA * kA * 2.0 * kA;
  for (int n : ns) {
    GridRun run;
    run.n = n;
    run.label = "a/dx=" + std::to_string(n);
    const auto t0 = std::chrono::steady_clock::now();
    try {
      flow::Problem pb;
      pb.grid = Grid(Geometry::Axisymmetric, 3 * n, static_cast<int>(std::lround(6.5 * n)), 0.0, 3.0 * kA, -5.0 * kA,
                     1.5 * kA);
      pb.bc.r_min.kind = flow::BoundaryKind::Axis;
      pb.gravity = {0.0, -kG};
      const Grid g = pb.grid;
      const auto shape = std::make_shared<body::Cylinder>(kA, 2.0 * kA);
      const double zc = -2.0 * kA;  // 윗면 −a
      Array2D f(g.nr, g.nz);
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) f(i, j) = std::clamp((0.0 - g.zf(j)) / g.dz, 0.0, 1.0);
      }
      const Array2D p_hyd = hydrostatic(pb, f);  // 고체를 물로 보고 정수압
      const Array2D phi = body::fluid_fraction(g, *shape, zc);
      for (int i = 0; i < g.nr; ++i) {
        for (int j = 0; j < g.nz; ++j) f(i, j) = std::min(f(i, j), phi(i, j));
      }
      auto bd = std::make_shared<flow::Body>();
      bd->shape = shape;
      bd->motion.z = [zc](double) { return zc; };
      bd->interpolation_order = 1;
      pb.body = bd;
      pb.air_kappa = 1.4;
      flow::SolverOptions opt;
      opt.convection = flow::Convection::Consistent;
  opt.time_order = 2;            // 03 §1.4 2차 시간 적분 (b8ebdf3)
  opt.air_pressure_floor = 1e3;  // 03 §1.4 공기 압력 하한
      auto mat = std::make_shared<interface::VofMaterial>(kRhoW, kMuW, kRhoA, kMuA);
      interface::TwoPhaseSolver s(pb, mat, f, true, opt);
      s.initialize(0.0, [](double, double) { return 0.0; }, [](double, double) { return 0.0; }, false);
      Array2D p0 = p_hyd;
      for (double& v : p0.data()) v += kPatm;
      s.flow().set_pressure(p0);
      const double dt = 0.1 * g.dz / std::sqrt(kG * kA);
      double vmax = 0.0;
      for (int k = 0; k < steps; ++k) {
        s.step(dt);
        vmax = std::max(vmax, max_speed(s.flow()));
        if (s.flow().stats().solver_failed) break;
      }
      const double force = body::axial_pressure_force(g, s.flow().pressure_geometry(), s.flow().pressure()) +
                           body::axial_viscous_force(g, s.flow().geometry(), *shape, zc, s.flow().w(), kMuW, 0.0);
      run.metrics["force"] = force;
      run.metrics["buoyancy_exact"] = buoy;
      run.metrics["force_relative_error"] = std::abs(force - buoy) / buoy;
      run.metrics["max_speed_over_run"] = vmax;
      run.seconds = seconds_since(t0);
      finish_run(run, s.flow());
    } catch (const std::exception& e) {
      run.failure = std::string("FAILED_EXCEPTION ") + e.what();
    }
    rep.runs.push_back(run);
  }
  if (runs_ok(rep.runs)) {
    Criterion c;
    c.text = "force_relative_error <= 0.01 on every grid";
    c.passed = true;
    for (const auto& r : rep.runs) {
      c.passed = c.passed && r.metrics.at("force_relative_error") <= 0.01;
      c.detail += (c.detail.empty() ? "" : ", ") + fmt(r.metrics.at("force_relative_error"));
    }
    rep.criteria.push_back(c);
  }
  rep.notes.push_back("closed cylinder radius a, height 2a, top at depth a; buoyancy = rho_w g pi a^2 (2a) = " +
                      fmt(buoy) + " N; still water from the discrete hydrostatic pressure + 1e5 Pa, dt = 0.1 dx / "
                      "sqrt(g a), " + std::to_string(steps) + " steps; max speed is report only");
  conclude(rep, quick);
  return rep;
}

// ---------------------------------------------------------------- 3-I
TestReport run_3i(bool quick) {
  TestReport rep = p4_report("3-I", "Side-wall effect: disk vs closed cylinder (report only)");
  std::vector<EntryOut> outs;
  for (BodyKind k : {BodyKind::Disk, BodyKind::ClosedCylinder}) {
    EntrySetup es;
    es.kind = k;
    es.n = quick ? 8 : 32;
    EntryOut o = run_entry(es);
    o.run.label = k == BodyKind::Disk ? "disk (thickness 0.2a)" : "closed cylinder (height 2a)";
    note_3d(rep, o.run);
    rep.runs.push_back(o.run);
    outs.push_back(o);
  }
  if (runs_ok(rep.runs)) {
    for (const char* key : {"impulse", "peak_bottom_avg_pressure", "peak_force"}) {
      const double d = rep.runs[0].metrics.at(key), c = rep.runs[1].metrics.at(key);
      rep.notes.push_back(std::string(key) + ": cylinder vs disk relative difference " + fmt((c - d) / d));
    }
  }
  rep.notes.push_back("3-B conditions; the closed cylinder is computed as a solid cylinder (rigid body, D-12)");
  rep.verdict = quick ? "QUICK_CHECK_ONLY" : "REPORT_ONLY";
  for (const auto& r : rep.runs) {
    if (!r.failure.empty()) rep.verdict += " (" + r.failure.substr(0, r.failure.find(' ')) + ")";
  }
  return rep;
}

}  // namespace hcw::verify
