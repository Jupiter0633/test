#pragma once

// 2차원 균일 엇갈림(MAC) 격자.
// 첫 번째 좌표는 r (축대칭) 또는 x (평면), 두 번째 좌표는 z (또는 y).
// 압력은 셀 중심, u(첫 번째 성분)는 r 면, w(두 번째 성분)는 z 면에 놓인다.

#include <cstddef>
#include <vector>

namespace hcw::flow {

enum class Geometry { Planar, Axisymmetric };

struct Grid {
  Geometry geometry = Geometry::Planar;
  int nr = 0, nz = 0;
  double r0 = 0, r1 = 1, z0 = 0, z1 = 1;
  double dr = 1, dz = 1;

  Grid() = default;
  Grid(Geometry g, int nr_, int nz_, double r0_, double r1_, double z0_, double z1_);

  double rc(int i) const { return r0 + (i + 0.5) * dr; }  // 셀 중심
  double rf(int i) const { return r0 + i * dr; }          // r 면
  double zc(int j) const { return z0 + (j + 0.5) * dz; }
  double zf(int j) const { return z0 + j * dz; }

  // 축대칭이면 r, 평면이면 1. 체적·면적 가중치 (2π 는 모든 항에 공통이라 생략).
  double metric(double r) const { return geometry == Geometry::Axisymmetric ? r : 1.0; }
  double cell_volume(int i) const { return metric(rc(i)) * dr * dz; }
};

// 행 우선 2차원 배열. (i, j) → i * nj + j
class Array2D {
 public:
  Array2D() = default;
  Array2D(int ni, int nj, double value = 0.0)
      : ni_(ni), nj_(nj), a_(static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj), value) {}

  int ni() const { return ni_; }
  int nj() const { return nj_; }
  double& operator()(int i, int j) { return a_[index(i, j)]; }
  double operator()(int i, int j) const { return a_[index(i, j)]; }
  std::vector<double>& data() { return a_; }
  const std::vector<double>& data() const { return a_; }
  void fill(double v) { a_.assign(a_.size(), v); }

 private:
  std::size_t index(int i, int j) const {
    return static_cast<std::size_t>(i) * static_cast<std::size_t>(nj_) + static_cast<std::size_t>(j);
  }
  int ni_ = 0, nj_ = 0;
  std::vector<double> a_;
};

}  // namespace hcw::flow
