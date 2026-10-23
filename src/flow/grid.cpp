#include "flow/grid.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hcw::flow {

Grid::Grid(Geometry g, int nr_, int nz_, double r0_, double r1_, double z0_, double z1_)
    : geometry(g), nr(nr_), nz(nz_), r0(r0_), r1(r1_), z0(z0_), z1(z1_) {
  if (nr < 2 || nz < 2) throw std::invalid_argument("grid needs at least 2 cells per direction");
  if (!(r1 > r0) || !(z1 > z0)) throw std::invalid_argument("grid extent must be positive");
  if (g == Geometry::Axisymmetric && r0 < 0.0) {
    throw std::invalid_argument("axisymmetric grid must have r0 >= 0");
  }
  hr_ = (r1 - r0) / nr;
  hz_ = (z1 - z0) / nz;
  dr_min_ = hr_;
  dz_min_ = hz_;
}

Grid Grid::stretched(Geometry g, const std::vector<double>& rf, const std::vector<double>& zf) {
  if (rf.size() < 3 || zf.size() < 3) throw std::invalid_argument("grid needs at least 2 cells per direction");
  Grid out(g, static_cast<int>(rf.size()) - 1, static_cast<int>(zf.size()) - 1, rf.front(), rf.back(), zf.front(),
           zf.back());
  for (std::size_t k = 1; k < rf.size(); ++k) {
    if (!(rf[k] > rf[k - 1])) throw std::invalid_argument("r faces must increase");
  }
  for (std::size_t k = 1; k < zf.size(); ++k) {
    if (!(zf[k] > zf[k - 1])) throw std::invalid_argument("z faces must increase");
  }
  out.uniform_ = false;
  out.rfa_ = rf;
  out.zfa_ = zf;
  out.dr_min_ = 1e300;
  out.dz_min_ = 1e300;
  for (int i = 0; i < out.nr; ++i) out.dr_min_ = std::min(out.dr_min_, out.drc(i));
  for (int j = 0; j < out.nz; ++j) out.dz_min_ = std::min(out.dz_min_, out.dzc(j));
  return out;
}

namespace {

int find_cell(const std::vector<double>& f, double x) {
  const auto it = std::upper_bound(f.begin(), f.end(), x);
  return static_cast<int>(it - f.begin()) - 1;
}

int find_face(const std::vector<double>& f, double x) {
  const int k = std::clamp(find_cell(f, x), 0, static_cast<int>(f.size()) - 2);
  const double h = f[static_cast<std::size_t>(k) + 1] - f[static_cast<std::size_t>(k)];
  if (std::abs(x - f[static_cast<std::size_t>(k)]) <= 1e-9 * h) return k;
  if (std::abs(x - f[static_cast<std::size_t>(k) + 1]) <= 1e-9 * h) return k + 1;
  return -1;
}

int uniform_face(double x, double x0, double h) {
  const double q = (x - x0) / h;
  const long k = std::lround(q);
  return std::abs(q - static_cast<double>(k)) > 1e-9 ? -1 : static_cast<int>(k);
}

}  // namespace

int Grid::cell_r(double r) const {
  return uniform_ ? static_cast<int>(std::floor((r - r0) / hr_)) : find_cell(rfa_, r);
}
int Grid::cell_z(double z) const {
  return uniform_ ? static_cast<int>(std::floor((z - z0) / hz_)) : find_cell(zfa_, z);
}
int Grid::face_r(double r) const { return uniform_ ? uniform_face(r, r0, hr_) : find_face(rfa_, r); }
int Grid::face_z(double z) const { return uniform_ ? uniform_face(z, z0, hz_) : find_face(zfa_, z); }

namespace {

// core 바깥 길이 len 을 채우는 칸 크기 (core 쪽부터). k 번째 칸 = min(h_max, h_core q^k).
// 칸 수 N 은 q = ratio 로 len 에 닿는 가장 작은 수, q 는 [1/ratio, ratio] 에서 이분법으로 합이 len 이 되게 정한다
// (이웃 크기 비 ≤ ratio 유지). len 이 너무 짧아 q = 1/ratio 로도 넘치면 그때만 비례로 줄인다.
std::vector<double> outer_sizes(double len, double h_core, double ratio, double h_max) {
  std::vector<double> h;
  if (len <= 0.0) return h;
  auto fill = [&](double q, int n) {
    std::vector<double> v(static_cast<std::size_t>(n));
    double size = h_core, sum = 0.0;
    for (int k = 0; k < n; ++k) {
      size = std::min(h_max, size * q);
      v[static_cast<std::size_t>(k)] = size;
      sum += size;
    }
    return std::make_pair(v, sum);
  };
  int n = 0;
  for (double size = h_core, sum = 0.0; sum < len; ++n) {
    size = std::min(h_max, size * ratio);
    sum += size;
  }
  double lo_q = 1.0 / ratio, hi_q = ratio;
  if (fill(lo_q, n).second > len) {
    hi_q = lo_q;
  } else {
    for (int it = 0; it < 200; ++it) {
      const double mid = 0.5 * (lo_q + hi_q);
      (fill(mid, n).second > len ? hi_q : lo_q) = mid;
    }
  }
  auto [v, sum] = fill(hi_q, n);
  for (double& x : v) x *= len / sum;
  return v;
}

// 단조 3 차 Hermite 보간 (Fritsch–Carlson). x 는 균일 매개변수 s_k = k/N 의 값.
double pchip(const std::vector<double>& y, double s) {
  const int n = static_cast<int>(y.size()) - 1;
  const double t = s * n;
  int k = std::clamp(static_cast<int>(std::floor(t)), 0, n - 1);
  const double u = t - k;
  auto slope = [&](int m) { return y[static_cast<std::size_t>(m) + 1] - y[static_cast<std::size_t>(m)]; };
  auto deriv = [&](int m) {
    if (m == 0) return slope(0);
    if (m == n) return slope(n - 1);
    const double a = slope(m - 1), b = slope(m);
    if (a * b <= 0.0) return 0.0;
    return 2.0 * a * b / (a + b);  // 조화 평균 (같은 간격)
  };
  const double y0 = y[static_cast<std::size_t>(k)], y1 = y[static_cast<std::size_t>(k) + 1];
  const double d0 = deriv(k), d1 = deriv(k + 1);
  const double h00 = (1 + 2 * u) * (1 - u) * (1 - u), h10 = u * (1 - u) * (1 - u);
  const double h01 = u * u * (3 - 2 * u), h11 = u * u * (u - 1);
  return h00 * y0 + h10 * d0 + h01 * y1 + h11 * d1;
}

}  // namespace

std::vector<double> stretched_faces(double lo, double hi, double c0, double c1, double h_core, double ratio,
                                    double h_max, int refine, const std::vector<double>& pins) {
  if (!(lo <= c0 && c0 < c1 && c1 <= hi) || !(h_core > 0.0) || !(ratio >= 1.0) || !(h_max >= h_core)) {
    throw std::invalid_argument("stretched_faces: bad arguments");
  }
  // core 칸 크기는 정확히 h_core: 칸 수를 올림하고 core 를 양쪽으로 고르게 넓힌다 (영역 안으로 밀어 넣음)
  int n_core = std::max(1, static_cast<int>(std::ceil((c1 - c0) / h_core - 1e-9)));
  const int n_fit = static_cast<int>(std::floor((hi - lo) / h_core + 1e-9));
  if (n_core >= n_fit) {  // 영역 전체가 core
    n_core = std::max(1, static_cast<int>(std::lround((hi - lo) / h_core)));
    c0 = lo;
    c1 = hi;
  } else {
    const double grow = 0.5 * (n_core * h_core - (c1 - c0));
    c0 -= grow;
    c1 += grow;
    if (c0 < lo) c1 += lo - c0, c0 = lo;
    if (c1 > hi) c0 -= c1 - hi, c1 = hi;
  }
  const double hc = (c1 - c0) / n_core;
  // 바깥: core 끝에서 영역 끝까지, 고정점 (pins) 마다 구간을 나눠 각 구간을 정확히 채운다.
  // 다음 구간은 앞 구간 마지막 칸 크기에서 이어서 커진다.
  auto side = [&](double from, double to) {
    std::vector<double> stops;
    for (double p : pins) {
      if ((p - from) * (to - p) > 0.0) stops.push_back(p);
    }
    std::sort(stops.begin(), stops.end(), [&](double x, double y) { return std::abs(x - from) < std::abs(y - from); });
    stops.push_back(to);
    std::vector<double> sizes;
    double start = from, h_prev = hc;
    for (double stop : stops) {
      const auto seg = outer_sizes(std::abs(stop - start), h_prev, ratio, h_max);
      sizes.insert(sizes.end(), seg.begin(), seg.end());
      if (!seg.empty()) h_prev = seg.back();
      start = stop;
    }
    return std::make_pair(sizes, stops);
  };
  const auto [below, below_stops] = side(c0, lo);
  const auto [above, above_stops] = side(c1, hi);
  std::vector<double> base;
  double x = lo;
  base.push_back(x);
  for (auto it = below.rbegin(); it != below.rend(); ++it) {
    x += *it;
    base.push_back(x);
  }
  base.back() = c0;
  for (int k = 1; k <= n_core; ++k) base.push_back(c0 + (c1 - c0) * k / n_core);
  for (double h : above) {
    x = base.back() + h;
    base.push_back(x);
  }
  base.back() = hi;
  // 고정점은 반올림 오차 없이 정확한 값으로 둔다
  for (const auto* stops : {&below_stops, &above_stops}) {
    for (double p : *stops) {
      auto it = std::min_element(base.begin(), base.end(), [p](double u, double v) { return std::abs(u - p) < std::abs(v - p); });
      *it = p;
    }
  }
  if (refine <= 0) return base;
  const int n = static_cast<int>(base.size()) - 1;
  const int m = n << refine;
  std::vector<double> out(static_cast<std::size_t>(m) + 1);
  for (int k = 0; k <= m; ++k) out[static_cast<std::size_t>(k)] = pchip(base, static_cast<double>(k) / m);
  out.front() = lo;
  out.back() = hi;
  return out;
}

}  // namespace hcw::flow
