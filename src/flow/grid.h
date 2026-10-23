#pragma once

// 2차원 엇갈림(MAC) 격자. 균일 격자와 늘어나는 격자 (칸 경계 위치를 따로 주는 텐서곱 격자, P4 03 §1.4) 를 지원한다.
// 첫 번째 좌표는 r (축대칭) 또는 x (평면), 두 번째 좌표는 z (또는 y).
// 압력은 셀 중심, u(첫 번째 성분)는 r 면, w(두 번째 성분)는 z 면에 놓인다.
// 셀 중심은 양쪽 면의 가운데이다.
// 균일 격자는 예전과 똑같은 식 (r0 + i·Δ 등) 으로 계산해 기존 결과가 비트 단위로 같다.

#include <cstddef>
#include <vector>

namespace hcw::flow {

enum class Geometry { Planar, Axisymmetric };

struct Grid {
  Geometry geometry = Geometry::Planar;
  int nr = 0, nz = 0;
  double r0 = 0, r1 = 1, z0 = 0, z1 = 1;

  Grid() = default;
  // 균일 격자
  Grid(Geometry g, int nr_, int nz_, double r0_, double r1_, double z0_, double z1_);
  // 늘어나는 격자: 칸 경계 위치 (단조 증가, 크기 nr+1, nz+1)
  static Grid stretched(Geometry g, const std::vector<double>& r_faces, const std::vector<double>& z_faces);

  bool uniform() const { return uniform_; }

  double rf(int i) const { return uniform_ ? r0 + i * hr_ : rfa_[static_cast<std::size_t>(i)]; }  // r 면
  double zf(int j) const { return uniform_ ? z0 + j * hz_ : zfa_[static_cast<std::size_t>(j)]; }
  double rc(int i) const { return uniform_ ? r0 + (i + 0.5) * hr_ : 0.5 * (rf(i) + rf(i + 1)); }  // 셀 중심
  double zc(int j) const { return uniform_ ? z0 + (j + 0.5) * hz_ : 0.5 * (zf(j) + zf(j + 1)); }

  // 셀 폭
  double drc(int i) const { return uniform_ ? hr_ : rf(i + 1) - rf(i); }
  double dzc(int j) const { return uniform_ ? hz_ : zf(j + 1) - zf(j); }
  // 면을 사이에 둔 두 셀 중심 사이 거리 (r 면 i, z 면 j). 영역 경계 면 (0, n) 은 붙은 셀의 폭
  // (벽·축의 유령 셀이 그 셀의 거울상이라 중심 거리가 셀 폭과 같다).
  double drn(int i) const {
    if (uniform_) return hr_;
    if (i <= 0) return drc(0);
    if (i >= nr) return drc(nr - 1);
    return rc(i) - rc(i - 1);
  }
  double dzn(int j) const {
    if (uniform_) return hz_;
    if (j <= 0) return dzc(0);
    if (j >= nz) return dzc(nz - 1);
    return zc(j) - zc(j - 1);
  }
  // 좌표가 들어 있는 셀 번호 (영역 밖이면 −1 또는 n 쪽 값이 나올 수 있음 — 호출하는 쪽에서 자름)
  int cell_r(double r) const;
  int cell_z(double z) const;
  // 좌표가 면 위에 있으면 그 면 번호, 아니면 −1 (상대 허용 1e−9 × 이웃 칸 크기)
  int face_r(double r) const;
  int face_z(double z) const;
  double dr_min() const { return dr_min_; }
  double dz_min() const { return dz_min_; }

  // 축대칭이면 r, 평면이면 1. 체적·면적 가중치 (2π 는 모든 항에 공통이라 생략).
  double metric(double r) const { return geometry == Geometry::Axisymmetric ? r : 1.0; }
  double cell_volume(int i, int j) const { return metric(rc(i)) * drc(i) * dzc(j); }

 private:
  bool uniform_ = true;
  double hr_ = 1, hz_ = 1;
  double dr_min_ = 1, dz_min_ = 1;
  std::vector<double> rfa_, zfa_;
};

// 늘어나는 격자 한 방향의 칸 경계 (03 §1.4).
// 가장 거친 격자: [lo, hi] 에서 core [c0, c1] 은 크기 h_core 로 같고, 바깥은 이웃 크기 비 ratio (≤ 1.05) 로
// 커지며 크기 상한 h_max. 영역 끝에 맞도록 바깥 칸들을 같은 비율로 늘이거나 줄인다.
// refine 번 세분: 가장 거친 격자의 경계를 매개변수 s = k/N 에 대해 단조 3 차 보간 (PCHIP) 한 사상에서
// s = k/(N·2^refine) 로 다시 잡는다. 균일 core 안의 경계는 정확히 등간격으로 나뉜다.
std::vector<double> stretched_faces(double lo, double hi, double c0, double c1, double h_core, double ratio,
                                    double h_max, int refine);

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
