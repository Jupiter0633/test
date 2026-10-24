#include "flow/poisson.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "common/parallel.h"

namespace hcw::flow {
namespace {

// 내적은 원래 순서 그대로 한 스레드가 더한다 (나누면 반올림이 달라져 결과가 바뀜).
double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// 파이프라인 구간 하나의 최소 행 수와 파이프라인을 쓰는 최소 미지수 수.
// 이보다 작으면 동기화 비용이 계산보다 커서 직렬로 푼다 (결과는 어느 쪽이든 같음).
constexpr int kMinSegmentRows = 8;
constexpr int kMinPipelineRows = 4096;
// 파이프라인 삼각 풀이에서 한 스레드가 섞어 푸는 열 수
constexpr int kInterleave = 4;

}  // namespace

void CsrMatrix::multiply(const std::vector<double>& x, std::vector<double>& y) const {
  y.resize(static_cast<std::size_t>(n));
  parallel::parallel_for(0, n, [&](int lo, int hi) {
    for (int i = lo; i < hi; ++i) {
      double s = 0.0;
      for (int k = row_start[i]; k < row_start[i + 1]; ++k) s += val[k] * x[col[k]];
      y[i] = s;
    }
  });
}

void IncompleteCholesky::Flags::reset(std::size_t size) {
  n = size;
  a = size ? std::make_unique<std::atomic<std::uint32_t>[]>(size * kStride) : nullptr;
  for (std::size_t k = 0; k < size * kStride; ++k) a[k].store(0, std::memory_order_relaxed);
}

void IncompleteCholesky::factor(const CsrMatrix& a, int inner) {
  n_ = a.n;
  row_start_.assign(static_cast<std::size_t>(n_) + 1, 0);
  col_.clear();
  val_.clear();
  diag_.assign(static_cast<std::size_t>(n_), 0.0);

  // 하삼각 패턴 복사 (대각 제외)
  for (int i = 0; i < n_; ++i) {
    for (int k = a.row_start[i]; k < a.row_start[i + 1]; ++k) {
      if (a.col[k] < i) {
        col_.push_back(a.col[k]);
        val_.push_back(a.val[k]);
      } else if (a.col[k] == i) {
        diag_[i] = a.val[k];
      }
    }
    row_start_[i + 1] = static_cast<int>(col_.size());
  }

  // L_ij = (A_ij - Σ_{m<j} L_im L_jm) / L_jj,  L_ii = sqrt(A_ii - Σ L_im²)  (패턴 제한)
  for (int i = 0; i < n_; ++i) {
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) {
      const int j = col_[p];
      double s = val_[p];
      int q = row_start_[j];
      for (int pp = row_start_[i]; pp < p; ++pp) {
        while (q < row_start_[j + 1] && col_[q] < col_[pp]) ++q;
        if (q < row_start_[j + 1] && col_[q] == col_[pp]) s -= val_[pp] * val_[q];
      }
      val_[p] = s / diag_[j];
    }
    double d = diag_[i];
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) d -= val_[p] * val_[p];
    if (!(d > 0.0)) {
      throw std::runtime_error("incomplete Cholesky breakdown at row " + std::to_string(i));
    }
    diag_[i] = std::sqrt(d);
  }
  work_.assign(static_cast<std::size_t>(n_), 0.0);
  plan_pipeline(inner);
}

void IncompleteCholesky::plan_pipeline(int inner) {
  seg_ = 0;
  lt_start_.clear();
  lt_row_.clear();
  lt_val_.clear();
  const int t = parallel::threads();
  if (t < 2 || inner <= 0 || n_ < kMinPipelineRows || n_ % inner != 0) return;
  // 스레드가 3 개 이상이면 하나는 남겨 apply(r, z, side) 의 side (pcg 의 내적) 를 맡긴다
  const int segs = std::min(t >= 3 ? t - 1 : t, inner / kMinSegmentRows);
  if (segs < 2) return;
  inner_ = inner;
  ncol_ = n_ / inner;
  seg_begin_.assign(static_cast<std::size_t>(segs) + 1, 0);
  for (int s = 0; s <= segs; ++s) {
    seg_begin_[s] = static_cast<int>(static_cast<long long>(inner) * s / segs);
  }
  std::vector<int> seg_of(static_cast<std::size_t>(inner));
  for (int s = 0; s < segs; ++s) {
    for (int j = seg_begin_[s]; j < seg_begin_[s + 1]; ++j) seg_of[j] = s;
  }
  // 앞 풀이: 행 i 는 L 의 열 q (< i) 를 읽는다 → q 의 (열, 구간) 이 둘 다 i 의 것 이하여야 한다.
  // 같은 구간의 앞 열이면 열 안 위치도 i 의 것 이하여야 한다 (열 섞기, apply_segment 참조).
  // 뒤 풀이: 행 q 는 L 의 행 i (> q) 를 읽는다 → 같은 조건을 반대로 본 것이라 따로 검사할 필요 없다.
  for (int i = 0; i < n_; ++i) {
    const int ci = i / inner, si = seg_of[i % inner];
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) {
      const int q = col_[p];
      const int cq = q / inner, sq = seg_of[q % inner];
      const bool ok = cq <= ci && sq <= si && (cq == ci || sq < si || q % inner <= i % inner);
      if (!ok) return;  // 이 순서로는 못 나눔 → 직렬
    }
  }
  // Lᵀ (열 q 마다 행 i 내림차순): 직렬 뒤 풀이의 흩뿌리기 z[q] −= L_iq z[i] (i 내림차순) 와
  // 같은 순서로 모아 빼기 위해
  lt_start_.assign(static_cast<std::size_t>(n_) + 1, 0);
  for (std::size_t p = 0; p < col_.size(); ++p) ++lt_start_[col_[p] + 1];
  for (int q = 0; q < n_; ++q) lt_start_[q + 1] += lt_start_[q];
  lt_row_.assign(col_.size(), 0);
  lt_val_.assign(col_.size(), 0.0);
  std::vector<int> fill(lt_start_.begin(), lt_start_.end() - 1);
  for (int i = n_ - 1; i >= 0; --i) {
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) {
      const int q = col_[p];
      lt_row_[fill[q]] = i;
      lt_val_[fill[q]] = val_[p];
      ++fill[q];
    }
  }
  fwd_done_.reset(static_cast<std::size_t>(ncol_) * segs);
  bwd_done_.reset(static_cast<std::size_t>(ncol_) * segs);
  epoch_ = 0;
  seg_ = segs;
}

void IncompleteCholesky::apply(const std::vector<double>& r, std::vector<double>& z) const {
  run(r, z, nullptr);
}

void IncompleteCholesky::apply(const std::vector<double>& r, std::vector<double>& z,
                               const std::function<void()>& side) const {
  run(r, z, &side);
}

void IncompleteCholesky::run(const std::vector<double>& r, std::vector<double>& z,
                             const std::function<void()>* side) const {
  z.resize(static_cast<std::size_t>(n_));
  const int t = parallel::available();
  if (seg_ == 0 || t < seg_) {
    if (side) (*side)();
    apply_serial(r, z);
    return;
  }
  if (side && t == seg_) (*side)();  // 남는 스레드가 없다
  if (++epoch_ == 0) ++epoch_;       // 0 은 "아직 안 함" 표시라 건너뛴다
  const std::uint32_t e = epoch_;
  parallel::for_each_thread([&](int tid) {
    if (tid < seg_) {
      apply_segment(tid, e, r, z);
    } else if (tid == seg_ && side) {
      (*side)();
    }
  });
}

void IncompleteCholesky::apply_segment(int s, std::uint32_t epoch, const std::vector<double>& r,
                                       std::vector<double>& z) const {
  auto wait = [&](const Flags& f, std::size_t k) {
    while (f.at(k).load(std::memory_order_acquire) != epoch) parallel::detail::pause();
  };
  const int j0 = seg_begin_[s], len = seg_begin_[s + 1] - j0;
  auto flag = [&](int c) { return static_cast<std::size_t>(c) * seg_ + s; };
  // 한 스레드 안에서도 열 kInterleave 개를 한 행씩 어긋나게 섞어 푼다. 한 열 안의 행들은 앞 행의
  // 나눗셈 결과를 기다려야 해서 CPU 가 놀지만, 이웃 열의 행은 따로 계산할 수 있기 때문이다.
  // 단계 k 에서 묶음의 m 번째 열은 열 안 위치 k − m 의 행을 푼다 (그 행이 읽는 앞 열의 같은 위치는
  // 이미 단계 k − 1 에서 끝났다; plan_pipeline 이 이 조건을 검사함).
  // L y = r : 열 0 → 끝. 열 묶음은 구간 s−1 이 같은 묶음을 끝낸 뒤 푼다.
  for (int c0 = 0; c0 < ncol_; c0 += kInterleave) {
    const int g = std::min(kInterleave, ncol_ - c0);
    if (s > 0) wait(fwd_done_, flag(c0 + g - 1) - 1);
    for (int k = 0; k < len + g - 1; ++k) {
      for (int m = 0; m < g; ++m) {
        const int jj = k - m;
        if (jj < 0 || jj >= len) continue;
        const int i = (c0 + m) * inner_ + j0 + jj;
        double v = r[i];
        for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) v -= val_[p] * work_[col_[p]];
        work_[i] = v / diag_[i];
      }
    }
    for (int m = 0; m < g; ++m) fwd_done_.at(flag(c0 + m)).store(epoch, std::memory_order_release);
  }
  // Lᵀ z = y : 열 끝 → 0, 열 안에서도 위 → 아래. 열 묶음은 구간 s+1 이 같은 묶음을 끝낸 뒤 푼다.
  for (int c0 = ncol_ - 1; c0 >= 0; c0 -= kInterleave) {
    const int g = std::min(kInterleave, c0 + 1);
    if (s + 1 < seg_) wait(bwd_done_, flag(c0 - g + 1) + 1);
    for (int k = 0; k < len + g - 1; ++k) {
      for (int m = 0; m < g; ++m) {
        const int jj = len - 1 - (k - m);
        if (jj < 0 || jj >= len) continue;
        const int q = (c0 - m) * inner_ + j0 + jj;
        double v = work_[q];
        for (int p = lt_start_[q]; p < lt_start_[q + 1]; ++p) v -= lt_val_[p] * z[lt_row_[p]];
        z[q] = v / diag_[q];
      }
    }
    for (int m = 0; m < g; ++m) bwd_done_.at(flag(c0 - m)).store(epoch, std::memory_order_release);
  }
}

void IncompleteCholesky::apply_serial(const std::vector<double>& r, std::vector<double>& z) const {
  // L y = r
  for (int i = 0; i < n_; ++i) {
    double s = r[i];
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) s -= val_[p] * work_[col_[p]];
    work_[i] = s / diag_[i];
  }
  // Lᵀ z = y
  for (int i = 0; i < n_; ++i) z[i] = work_[i];
  for (int i = n_ - 1; i >= 0; --i) {
    z[i] /= diag_[i];
    for (int p = row_start_[i]; p < row_start_[i + 1]; ++p) z[col_[p]] -= val_[p] * z[i];
  }
}

SolveResult pcg(const CsrMatrix& a, const IncompleteCholesky& m, const std::vector<double>& b,
                std::vector<double>& x, double tolerance, int max_iterations) {
  SolveResult result;
  const std::size_t n = static_cast<std::size_t>(a.n);
  x.resize(n, 0.0);
  const double b_norm = std::sqrt(dot(b, b));
  if (b_norm == 0.0) {
    x.assign(n, 0.0);
    result.converged = true;
    return result;
  }

  std::vector<double> r(n), z(n), p(n), q(n);
  a.multiply(x, q);
  for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - q[i];
  double r_norm = std::sqrt(dot(r, r));
  // 초기 추정이 0 보다 나쁘면 버린다. 그대로 쓰면 ||b|| 기준 허용치가
  // 초기 잔차 대비 기계 정밀도 아래로 내려가 수렴할 수 없게 된다.
  if (r_norm > b_norm) {
    x.assign(n, 0.0);
    r = b;
    r_norm = b_norm;
  }
  if (r_norm <= tolerance * b_norm) {
    result.converged = true;
    result.relative_residual = r_norm / b_norm;
    return result;
  }
  m.apply(r, z);
  p = z;
  double rz = dot(r, z);
  // 정체 감지: 최근 kStall 회 동안 최저 잔차가 개선되지 않으면 반올림 한계에
  // 도달한 것으로 보고 멈춘다 (허용치 미달 → converged = false).
  constexpr int kStall = 500;
  double best = r_norm;
  int best_it = 0;

  for (int it = 1; it <= max_iterations; ++it) {
    a.multiply(p, q);
    const double alpha = rz / dot(p, q);
    parallel::parallel_for(0, a.n, [&](int lo, int hi) {
      for (int i = lo; i < hi; ++i) {
        x[i] += alpha * p[i];
        r[i] -= alpha * q[i];
      }
    });
    // ||r||² (한 스레드가 원래 순서로 더함) 와 전처리 z = M⁻¹ r 은 서로 상관없으므로 동시에 한다.
    // 아래에서 r 을 참 잔차로 바꾸면 z 를 다시 구한다. 수렴·정체로 끝나면 z 는 버린다.
    // 어느 경우든 쓰이는 z 는 직렬 코드와 같은 r 로 구한 것이라 결과가 같다.
    double rr = 0.0;
    m.apply(r, z, [&] { rr = dot(r, r); });
    r_norm = std::sqrt(rr);
    result.iterations = it;
    bool z_stale = false;
    if (r_norm <= tolerance * b_norm) {
      // 누적 반올림을 피하려고 참 잔차로 확인한다.
      z_stale = true;
      a.multiply(x, q);
      parallel::parallel_for(0, a.n, [&](int lo, int hi) {
        for (int i = lo; i < hi; ++i) r[i] = b[i] - q[i];
      });
      r_norm = std::sqrt(dot(r, r));
      if (r_norm <= tolerance * b_norm) {
        result.converged = true;
        break;
      }
    }
    if (r_norm < 0.999 * best) {
      best = r_norm;
      best_it = it;
    } else if (it - best_it >= kStall) {
      break;
    }
    if (z_stale) m.apply(r, z);
    const double rz_new = dot(r, z);
    const double beta = rz_new / rz;
    rz = rz_new;
    parallel::parallel_for(0, a.n, [&](int lo, int hi) {
      for (int i = lo; i < hi; ++i) p[i] = z[i] + beta * p[i];
    });
  }
  result.relative_residual = r_norm / b_norm;
  return result;
}

namespace {

// s + e = a + b 를 정확히 (Knuth TwoSum)
inline void two_sum(double a, double b, double& s, double& e) {
  s = a + b;
  const double bb = s - a;
  e = (a - (s - bb)) + (b - bb);
}

// r = b − A(hi + lo). 행마다 배정밀도 두 배 수준으로 누적한다 (Ogita–Rump–Oishi Dot2 방식).
double accurate_residual(const CsrMatrix& a, const std::vector<double>& b,
                         const std::vector<double>& hi, const std::vector<double>& lo,
                         std::vector<double>& r) {
  r.resize(static_cast<std::size_t>(a.n));
  parallel::parallel_for(0, a.n, [&](int row_lo, int row_hi) {
    for (int i = row_lo; i < row_hi; ++i) {
      double s = b[i], c = 0.0;
      for (int k = a.row_start[i]; k < a.row_start[i + 1]; ++k) {
        const double av = a.val[k];
        const double xh = hi[a.col[k]];
        const double p = -av * xh;
        const double pe = std::fma(-av, xh, -p);  // 곱의 반올림 오차
        double e;
        two_sum(s, p, s, e);
        c += e + pe - av * lo[a.col[k]];
      }
      r[i] = s + c;
    }
  });
  double sum_sq = 0.0;  // 합은 원래 순서로 한 스레드가
  for (int i = 0; i < a.n; ++i) sum_sq += r[i] * r[i];
  return std::sqrt(sum_sq);
}

}  // namespace

SolveResult refined_solve(const CsrMatrix& a, const IncompleteCholesky& m,
                          const std::vector<double>& b, std::vector<double>& hi,
                          std::vector<double>& lo, double tolerance, int max_iterations) {
  // 내부 pcg 가 한 번에 줄이는 비율. 배정밀도 반올림 한계보다 충분히 커야 한다.
  constexpr double kInnerTolerance = 1e-6;
  // 바깥 반복에서 잔차가 이 비율보다 덜 줄면 정체로 본다.
  constexpr double kMinReduction = 0.5;
  constexpr int kMaxStalls = 3;

  SolveResult result;
  const std::size_t n = static_cast<std::size_t>(a.n);
  hi.resize(n, 0.0);
  lo.resize(n, 0.0);
  const double b_norm = std::sqrt(dot(b, b));
  if (b_norm == 0.0) {
    hi.assign(n, 0.0);
    lo.assign(n, 0.0);
    result.converged = true;
    return result;
  }

  std::vector<double> r, d;
  double r_norm = accurate_residual(a, b, hi, lo, r);
  // 초기 추정이 0 보다 나쁘면 버린다 (pcg 와 같은 이유).
  if (r_norm > b_norm) {
    hi.assign(n, 0.0);
    lo.assign(n, 0.0);
    r = b;
    r_norm = b_norm;
  }

  int stalls = 0;
  while (r_norm > tolerance * b_norm) {
    const int budget = max_iterations - result.iterations;
    if (budget <= 0) break;
    d.assign(n, 0.0);
    const SolveResult inner = pcg(a, m, r, d, kInnerTolerance, budget);
    result.iterations += inner.iterations;
    result.refinements += 1;
    // (hi, lo) += d
    parallel::parallel_for(0, a.n, [&](int row_lo, int row_hi) {
      for (int i = row_lo; i < row_hi; ++i) {
        double s, e;
        two_sum(hi[i], d[i], s, e);
        e += lo[i];
        hi[i] = s + e;
        lo[i] = e - (hi[i] - s);
      }
    });
    const double previous = r_norm;
    r_norm = accurate_residual(a, b, hi, lo, r);
    if (r_norm > kMinReduction * previous) {
      if (++stalls >= kMaxStalls) {
        result.stagnated = true;
        break;
      }
    } else {
      stalls = 0;
    }
  }
  result.relative_residual = r_norm / b_norm;
  result.converged = r_norm <= tolerance * b_norm;
  return result;
}

}  // namespace hcw::flow
