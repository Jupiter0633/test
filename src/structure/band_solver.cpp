#include "structure/band_solver.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace hcw::structure {

SkylineMatrix::SkylineMatrix(const std::vector<int>& top) : top_(top), start_(top.size() + 1, 0) {
  for (std::size_t j = 0; j < top.size(); ++j) {
    if (top[j] < 0 || top[j] > static_cast<int>(j)) throw std::invalid_argument("SkylineMatrix: bad profile");
    start_[j + 1] = start_[j] + (j - static_cast<std::size_t>(top[j]) + 1);
  }
  values_.assign(start_.back(), 0.0);
}

void SkylineMatrix::add(int i, int j, double v) {
  if (i > j) std::swap(i, j);
  if (i < top_[j]) throw std::logic_error("SkylineMatrix::add outside the profile");
  values_[at(i, j)] += v;
}

double SkylineMatrix::get(int i, int j) const {
  if (i > j) std::swap(i, j);
  if (i < top_[j]) return 0.0;
  return values_[at(i, j)];
}

void SkylineMatrix::factor() {
  const int n = size();
  for (int j = 0; j < n; ++j) {
    const int tj = top_[j];
    double* cj = values_.data() + start_[j];  // cj[i − tj] = (i, j) 자리
    // g_ij = a_ij − Σ_{k<i} l_ki g_kj
    for (int i = tj + 1; i < j; ++i) {
      const int ti = top_[i];
      const double* ci = values_.data() + start_[i];
      const int k0 = std::max(ti, tj);
      double s = 0.0;
      for (int k = k0; k < i; ++k) s += ci[k - ti] * cj[k - tj];
      cj[i - tj] -= s;
    }
    // l_ij = g_ij / d_i, d_j = a_jj − Σ g_ij l_ij
    const double ajj = cj[j - tj];
    double d = ajj;
    for (int i = tj; i < j; ++i) {
      const double g = cj[i - tj];
      const double l = g / values_[at(i, i)];
      d -= g * l;
      cj[i - tj] = l;
    }
    if (!(std::abs(d) > 1e-14 * std::abs(ajj))) throw std::runtime_error("SkylineMatrix: singular pivot");
    cj[j - tj] = d;
  }
  factored_ = true;
}

void SkylineMatrix::solve(std::vector<double>& b) const {
  if (!factored_) throw std::logic_error("SkylineMatrix::solve before factor");
  const int n = size();
  // L y = b
  for (int j = 0; j < n; ++j) {
    const int tj = top_[j];
    const double* cj = values_.data() + start_[j];
    double s = 0.0;
    for (int k = tj; k < j; ++k) s += cj[k - tj] * b[k];
    b[j] -= s;
  }
  for (int j = 0; j < n; ++j) b[j] /= values_[at(j, j)];
  // Lᵀ x = y
  for (int j = n - 1; j >= 0; --j) {
    const int tj = top_[j];
    const double* cj = values_.data() + start_[j];
    const double xj = b[j];
    for (int k = tj; k < j; ++k) b[k] -= cj[k - tj] * xj;
  }
}

}  // namespace hcw::structure
