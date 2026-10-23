#include "structure/eigen.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace hcw::structure {

void jacobi_eigen(int n, std::vector<double>& a, std::vector<double>& v) {
  v.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int i = 0; i < n; ++i) v[i * n + i] = 1.0;
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0, tot = 0.0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        tot += a[i * n + j] * a[i * n + j];
        if (i != j) off += a[i * n + j] * a[i * n + j];
      }
    }
    if (off <= 1e-30 * tot) return;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = a[p * n + q];
        if (apq == 0.0) continue;
        const double theta = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; ++k) {  // 열 p, q
          const double akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {  // 행 p, q
          const double apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = v[k * n + p], vkq = v[k * n + q];
          v[k * n + p] = c * vkp - s * vkq;
          v[k * n + q] = s * vkp + c * vkq;
        }
      }
    }
  }
}

namespace {

// 작은 일반 고유값 문제 Kr q = λ Mr q (Mr 양의 정부호). 고유값 오름차순, q 는 Mr-정규 직교.
void small_generalized(int n, const std::vector<double>& Kr, const std::vector<double>& Mr, std::vector<double>& lam,
                       std::vector<double>& Q) {
  // Mr = L Lᵀ (Cholesky)
  std::vector<double> L(static_cast<std::size_t>(n) * n, 0.0);
  for (int j = 0; j < n; ++j) {
    double d = Mr[j * n + j];
    for (int k = 0; k < j; ++k) d -= L[j * n + k] * L[j * n + k];
    if (!(d > 0.0)) throw std::runtime_error("subspace mass matrix is not positive definite");
    L[j * n + j] = std::sqrt(d);
    for (int i = j + 1; i < n; ++i) {
      double s = Mr[i * n + j];
      for (int k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
      L[i * n + j] = s / L[j * n + j];
    }
  }
  // C = L⁻¹ Kr L⁻ᵀ
  auto lsolve_cols = [&](std::vector<double>& X) {  // X ← L⁻¹ X (열마다)
    for (int c = 0; c < n; ++c) {
      for (int i = 0; i < n; ++i) {
        double s = X[i * n + c];
        for (int k = 0; k < i; ++k) s -= L[i * n + k] * X[k * n + c];
        X[i * n + c] = s / L[i * n + i];
      }
    }
  };
  std::vector<double> C = Kr;
  lsolve_cols(C);  // L⁻¹ Kr
  std::vector<double> Ct(C.size());
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) Ct[i * n + j] = C[j * n + i];
  }
  lsolve_cols(Ct);  // L⁻¹ (L⁻¹ Kr)ᵀ = L⁻¹ Kr L⁻ᵀ (Kr 대칭)
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      const double m = 0.5 * (Ct[i * n + j] + Ct[j * n + i]);
      Ct[i * n + j] = Ct[j * n + i] = m;
    }
  }
  std::vector<double> V;
  jacobi_eigen(n, Ct, V);
  std::vector<int> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) { return Ct[a * n + a] < Ct[b * n + b]; });
  lam.resize(n);
  Q.assign(static_cast<std::size_t>(n) * n, 0.0);
  for (int c = 0; c < n; ++c) {
    lam[c] = Ct[idx[c] * n + idx[c]];
    // q = L⁻ᵀ v
    for (int i = n - 1; i >= 0; --i) {
      double s = V[i * n + idx[c]];
      for (int k = i + 1; k < n; ++k) s -= L[k * n + i] * Q[k * n + c];
      Q[i * n + c] = s / L[i * n + i];
    }
  }
}

}  // namespace

EigenResult lowest_modes(const Model& model, int count, double shift, int max_iterations,
                         const std::vector<Vec>& known) {
  EigenResult res;
  const int nf = model.num_free();
  const int q = std::min(nf, std::max(2 * count, count + 8));
  SkylineMatrix A = model.assemble(1.0, -shift);
  A.factor();

  // 처음 벡터: 정해진 의사 난수 (재현 가능)
  std::vector<Vec> X(q, Vec(nf));
  unsigned long long state = 88172645463325252ULL;
  for (auto& x : X) {
    for (double& v : x) {
      state ^= state << 13;
      state ^= state >> 7;
      state ^= state << 17;
      v = static_cast<double>(state % 2000001ULL) / 1000000.0 - 1.0;
    }
  }
  std::vector<double> lam_old(count, 0.0);
  Vec full(model.num_dofs(), 0.0), y;
  // 미리 아는 모드 (자유 자유도 부분) 와 그 M 곱
  std::vector<Vec> kf, kMf;
  std::vector<double> kn;
  for (const auto& k : known) {
    model.apply_M(k, y);
    kf.push_back(model.restrict_free(k));
    kMf.push_back(model.restrict_free(y));
    kn.push_back(model.dot(kf.back(), kMf.back()));
  }
  auto deflate = [&](Vec& x) {
    for (std::size_t m = 0; m < kf.size(); ++m) {
      const double c = model.dot(kMf[m], x) / kn[m];
      for (std::size_t i = 0; i < x.size(); ++i) x[i] -= c * kf[m][i];
    }
  };
  for (auto& x : X) deflate(x);
  std::vector<Vec> Xb(q), KX(q), MX(q);
  for (int it = 1; it <= max_iterations; ++it) {
    for (int c = 0; c < q; ++c) {
      model.scatter_free(X[c], full);
      model.apply_M(full, y);
      double rr = 0.0;
      Xb[c] = solve_checked(model, A, 1.0, -shift, model.restrict_free(y), &rr);
      res.max_rel_residual = std::max(res.max_rel_residual, rr);
      deflate(Xb[c]);
      model.scatter_free(Xb[c], full);
      model.apply_K(full, y);
      KX[c] = model.restrict_free(y);
      model.apply_M(full, y);
      MX[c] = model.restrict_free(y);
    }
    std::vector<double> Kr(static_cast<std::size_t>(q) * q), Mr(Kr.size());
    for (int i = 0; i < q; ++i) {
      for (int j = 0; j < q; ++j) {
        Kr[i * q + j] = model.dot(Xb[i], KX[j]);
        Mr[i * q + j] = model.dot(Xb[i], MX[j]);
      }
    }
    for (int i = 0; i < q; ++i) {
      for (int j = i + 1; j < q; ++j) {
        Kr[i * q + j] = Kr[j * q + i] = 0.5 * (Kr[i * q + j] + Kr[j * q + i]);
        Mr[i * q + j] = Mr[j * q + i] = 0.5 * (Mr[i * q + j] + Mr[j * q + i]);
      }
    }
    std::vector<double> lam, Q;
    small_generalized(q, Kr, Mr, lam, Q);
    for (int c = 0; c < q; ++c) {
      std::fill(X[c].begin(), X[c].end(), 0.0);
      for (int k = 0; k < q; ++k) {
        const double w = Q[k * q + c];
        for (int i = 0; i < nf; ++i) X[c][i] += w * Xb[k][i];
      }
    }
    // 수렴: 구하려는 고유값 각각의 상대 변화 ≤ 1e−12. 0 에 가까운 (강체) 값은 가장 큰 구하려는 값에 대해.
    const double scale = std::abs(lam[count - 1]);
    bool done = it > 1;
    for (int c = 0; c < count; ++c) {
      const double ref = std::max(std::abs(lam[c]), 1e-6 * scale);
      if (std::abs(lam[c] - lam_old[c]) > 1e-12 * ref) done = false;
      lam_old[c] = lam[c];
    }
    res.iterations = it;
    if (done) {
      res.converged = true;
      break;
    }
  }
  res.omega2 = lam_old;
  res.modes.resize(count);
  for (int c = 0; c < count; ++c) {
    res.modes[c].assign(model.num_dofs(), 0.0);
    model.scatter_free(X[c], res.modes[c]);
  }
  return res;
}

double refine_mode(const Model& model, double omega2, Vec& mode, int iterations) {
  const double sigma = omega2 * (1.0 - 1e-4);
  SkylineMatrix A = model.assemble(1.0, -sigma);
  A.factor();
  Vec y, x = model.restrict_free(mode);
  Vec full(model.num_dofs(), 0.0);
  for (int it = 0; it < iterations; ++it) {
    model.scatter_free(x, full);
    model.apply_M(full, y);
    Vec rhs = model.restrict_free(y);
    A.solve(rhs);
    x = rhs;
    model.scatter_free(x, full);
    model.apply_M(full, y);
    const double nrm = std::sqrt(model.dot(full, y));
    for (double& v : x) v /= nrm;
  }
  mode.assign(model.num_dofs(), 0.0);
  model.scatter_free(x, mode);
  model.apply_K(mode, y);
  const double k = model.dot(mode, y);
  model.apply_M(mode, y);
  return k / model.dot(mode, y);
}

}  // namespace hcw::structure
