#include "flow/poisson.h"

#include <cmath>
#include <stdexcept>
#include <string>

namespace hcw::flow {
namespace {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

}  // namespace

void CsrMatrix::multiply(const std::vector<double>& x, std::vector<double>& y) const {
  y.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double s = 0.0;
    for (int k = row_start[i]; k < row_start[i + 1]; ++k) s += val[k] * x[col[k]];
    y[i] = s;
  }
}

void IncompleteCholesky::factor(const CsrMatrix& a) {
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
}

void IncompleteCholesky::apply(const std::vector<double>& r, std::vector<double>& z) const {
  z.resize(static_cast<std::size_t>(n_));
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

  for (int it = 1; it <= max_iterations; ++it) {
    a.multiply(p, q);
    const double alpha = rz / dot(p, q);
    for (std::size_t i = 0; i < n; ++i) {
      x[i] += alpha * p[i];
      r[i] -= alpha * q[i];
    }
    r_norm = std::sqrt(dot(r, r));
    result.iterations = it;
    if (r_norm <= tolerance * b_norm) {
      // 누적 반올림을 피하려고 참 잔차로 확인한다.
      a.multiply(x, q);
      for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - q[i];
      r_norm = std::sqrt(dot(r, r));
      if (r_norm <= tolerance * b_norm) {
        result.converged = true;
        break;
      }
    }
    m.apply(r, z);
    const double rz_new = dot(r, z);
    const double beta = rz_new / rz;
    rz = rz_new;
    for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
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
  double sum_sq = 0.0;
  for (int i = 0; i < a.n; ++i) {
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
    sum_sq += r[i] * r[i];
  }
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
    for (std::size_t i = 0; i < n; ++i) {
      double s, e;
      two_sum(hi[i], d[i], s, e);
      e += lo[i];
      hi[i] = s + e;
      lo[i] = e - (hi[i] - s);
    }
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
