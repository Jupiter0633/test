#include "flow/linear.h"

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

double norm(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

}  // namespace

void Ilu0::factor(const CsrMatrix& a) {
  n_ = a.n;
  row_start_ = a.row_start;
  col_ = a.col;
  val_ = a.val;
  diag_pos_.assign(static_cast<std::size_t>(n_), -1);
  for (int i = 0; i < n_; ++i) {
    for (int k = row_start_[i]; k < row_start_[i + 1]; ++k) {
      if (col_[k] == i) diag_pos_[i] = k;
    }
    if (diag_pos_[i] < 0) throw std::runtime_error("ILU(0): missing diagonal at row " + std::to_string(i));
  }
  std::vector<int> where(static_cast<std::size_t>(n_), -1);
  for (int i = 0; i < n_; ++i) {
    for (int k = row_start_[i]; k < row_start_[i + 1]; ++k) where[col_[k]] = k;
    for (int k = row_start_[i]; k < row_start_[i + 1] && col_[k] < i; ++k) {
      const int j = col_[k];
      const double pivot = val_[diag_pos_[j]];
      if (pivot == 0.0) throw std::runtime_error("ILU(0): zero pivot at row " + std::to_string(j));
      val_[k] /= pivot;
      for (int kk = diag_pos_[j] + 1; kk < row_start_[j + 1]; ++kk) {
        const int w = where[col_[kk]];
        if (w >= 0) val_[w] -= val_[k] * val_[kk];
      }
    }
    for (int k = row_start_[i]; k < row_start_[i + 1]; ++k) where[col_[k]] = -1;
    if (val_[diag_pos_[i]] == 0.0) throw std::runtime_error("ILU(0): zero pivot at row " + std::to_string(i));
  }
}

void Ilu0::apply(const std::vector<double>& r, std::vector<double>& z) const {
  z.resize(static_cast<std::size_t>(n_));
  for (int i = 0; i < n_; ++i) {  // L y = r (단위 대각)
    double s = r[i];
    for (int k = row_start_[i]; k < diag_pos_[i]; ++k) s -= val_[k] * z[col_[k]];
    z[i] = s;
  }
  for (int i = n_ - 1; i >= 0; --i) {  // U z = y
    double s = z[i];
    for (int k = diag_pos_[i] + 1; k < row_start_[i + 1]; ++k) s -= val_[k] * z[col_[k]];
    z[i] = s / val_[diag_pos_[i]];
  }
}

SolveResult bicgstab(const CsrMatrix& a, const Ilu0& m, const std::vector<double>& b,
                     std::vector<double>& x, double tolerance, int max_iterations) {
  constexpr int kStall = 200;
  SolveResult result;
  const std::size_t n = static_cast<std::size_t>(a.n);
  x.resize(n, 0.0);
  const double b_norm = norm(b);
  if (b_norm == 0.0) {
    x.assign(n, 0.0);
    result.converged = true;
    return result;
  }
  std::vector<double> r(n), r0, p(n, 0.0), v(n, 0.0), s(n), t(n), ph(n), sh(n), ax(n);
  a.multiply(x, ax);
  for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - ax[i];
  if (norm(r) > b_norm) {  // 초기 추정이 0 보다 나쁘면 버린다
    x.assign(n, 0.0);
    r = b;
  }
  double r_norm = norm(r);
  r0 = r;
  double rho = 1.0, alpha = 1.0, omega = 1.0, best = r_norm;
  int best_it = 0;
  for (int it = 1; it <= max_iterations && r_norm > tolerance * b_norm; ++it) {
    result.iterations = it;
    const double rho_new = dot(r0, r);
    if (rho_new == 0.0) break;
    const double beta = (rho_new / rho) * (alpha / omega);
    rho = rho_new;
    for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + beta * (p[i] - omega * v[i]);
    m.apply(p, ph);
    a.multiply(ph, v);
    const double r0v = dot(r0, v);
    if (r0v == 0.0) break;
    alpha = rho / r0v;
    for (std::size_t i = 0; i < n; ++i) s[i] = r[i] - alpha * v[i];
    m.apply(s, sh);
    a.multiply(sh, t);
    const double tt = dot(t, t);
    omega = tt > 0.0 ? dot(t, s) / tt : 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      x[i] += alpha * ph[i] + omega * sh[i];
      r[i] = s[i] - omega * t[i];
    }
    r_norm = norm(r);
    if (r_norm <= tolerance * b_norm || it % 50 == 0) {  // 참 잔차로 확인·교체
      a.multiply(x, ax);
      for (std::size_t i = 0; i < n; ++i) r[i] = b[i] - ax[i];
      r_norm = norm(r);
    }
    if (r_norm < 0.999 * best) {
      best = r_norm;
      best_it = it;
    } else if (it - best_it >= kStall) {
      break;
    }
    if (omega == 0.0) break;
  }
  result.relative_residual = r_norm / b_norm;
  result.converged = r_norm <= tolerance * b_norm;
  return result;
}

}  // namespace hcw::flow
