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
  // 정체 감지: 최근 kStall 회 동안 최저 잔차가 개선되지 않으면 반올림 한계에
  // 도달한 것으로 보고 멈춘다 (허용치 미달 → converged = false).
  constexpr int kStall = 500;
  double best = r_norm;
  int best_it = 0;

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
    if (r_norm < 0.999 * best) {
      best = r_norm;
      best_it = it;
    } else if (it - best_it >= kStall) {
      break;
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

}  // namespace hcw::flow
