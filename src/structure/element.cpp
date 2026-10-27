#include "structure/element.h"

#include <cmath>
#include <stdexcept>

namespace hcw::structure {

const std::vector<double>& gauss_points(int n) {
  static const std::vector<double> p2{-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
  static const std::vector<double> p3{-std::sqrt(0.6), 0.0, std::sqrt(0.6)};
  if (n == 2) return p2;
  if (n == 3) return p3;
  throw std::invalid_argument("gauss_points: n must be 2 or 3");
}

const std::vector<double>& gauss_weights(int n) {
  static const std::vector<double> w2{1.0, 1.0};
  static const std::vector<double> w3{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  if (n == 2) return w2;
  if (n == 3) return w3;
  throw std::invalid_argument("gauss_weights: n must be 2 or 3");
}

void shape_functions(double xi, double eta, double* N, double* dNxi, double* dNeta) {
  static const double cx[4] = {-1, 1, 1, -1}, cy[4] = {-1, -1, 1, 1};
  for (int a = 0; a < 4; ++a) {
    const double x = cx[a], y = cy[a];
    N[a] = 0.25 * (1 + xi * x) * (1 + eta * y) * (xi * x + eta * y - 1);
    dNxi[a] = 0.25 * x * (1 + eta * y) * (2 * xi * x + eta * y);
    dNeta[a] = 0.25 * y * (1 + xi * x) * (xi * x + 2 * eta * y);
  }
  // 중간 절점 4 (0,−1), 6 (0,1): ½(1−ξ²)(1+ηy); 5 (1,0), 7 (−1,0): ½(1+ξx)(1−η²)
  for (int a : {4, 6}) {
    const double y = a == 4 ? -1.0 : 1.0;
    N[a] = 0.5 * (1 - xi * xi) * (1 + eta * y);
    dNxi[a] = -xi * (1 + eta * y);
    dNeta[a] = 0.5 * y * (1 - xi * xi);
  }
  for (int a : {5, 7}) {
    const double x = a == 5 ? 1.0 : -1.0;
    N[a] = 0.5 * (1 + xi * x) * (1 - eta * eta);
    dNxi[a] = 0.5 * x * (1 - eta * eta);
    dNeta[a] = -eta * (1 + xi * x);
  }
}

PointEval evaluate(const Mesh& mesh, int elem, double xi, double eta) {
  PointEval p;
  double dxi[8], deta[8];
  shape_functions(xi, eta, p.N.data(), dxi, deta);
  const auto& e = mesh.elems[elem];
  double j11 = 0, j12 = 0, j21 = 0, j22 = 0;  // j11 = ∂r/∂ξ, j12 = ∂z/∂ξ, j21 = ∂r/∂η, j22 = ∂z/∂η
  for (int a = 0; a < 8; ++a) {
    const Node& n = mesh.nodes[e[a]];
    p.r += p.N[a] * n.r;
    p.z += p.N[a] * n.z;
    j11 += dxi[a] * n.r;
    j12 += dxi[a] * n.z;
    j21 += deta[a] * n.r;
    j22 += deta[a] * n.z;
  }
  p.detJ = j11 * j22 - j12 * j21;
  if (!(p.detJ > 0.0)) throw std::runtime_error("structure element has non-positive Jacobian");
  for (int a = 0; a < 8; ++a) {
    p.dNr[a] = (j22 * dxi[a] - j12 * deta[a]) / p.detJ;
    p.dNz[a] = (-j21 * dxi[a] + j11 * deta[a]) / p.detJ;
  }
  return p;
}

std::array<double, 4> stress_from_strain(const Material& mat, double err, double ett, double ezz, double grz) {
  const double lam = mat.lambda(), mu = mat.mu(), tr = err + ett + ezz;
  return {lam * tr + 2 * mu * err, lam * tr + 2 * mu * ett, lam * tr + 2 * mu * ezz, mu * grz};
}

void element_matrices(const Mesh& mesh, int elem, const Material& mat, double* K, double* M) {
  for (int i = 0; i < 256; ++i) K[i] = M[i] = 0.0;
  const double lam = mat.lambda(), mu = mat.mu();
  const double D[4][4] = {{lam + 2 * mu, lam, lam, 0},
                          {lam, lam + 2 * mu, lam, 0},
                          {lam, lam, lam + 2 * mu, 0},
                          {0, 0, 0, mu}};
  const auto& gp = gauss_points(3);
  const auto& gw = gauss_weights(3);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      const PointEval p = evaluate(mesh, elem, gp[i], gp[j]);
      const double w = gw[i] * gw[j] * p.detJ * kTwoPi * p.r;
      // B (4×16): 행 = (ε_rr, ε_θθ, ε_zz, γ_rz)
      double B[4][16] = {};
      for (int a = 0; a < 8; ++a) {
        B[0][2 * a] = p.dNr[a];
        B[1][2 * a] = p.N[a] / p.r;
        B[2][2 * a + 1] = p.dNz[a];
        B[3][2 * a] = p.dNz[a];
        B[3][2 * a + 1] = p.dNr[a];
      }
      double DB[4][16];
      for (int s = 0; s < 4; ++s) {
        for (int c = 0; c < 16; ++c) {
          double v = 0;
          for (int t = 0; t < 4; ++t) v += D[s][t] * B[t][c];
          DB[s][c] = v;
        }
      }
      for (int r = 0; r < 16; ++r) {
        for (int c = 0; c < 16; ++c) {
          double v = 0;
          for (int s = 0; s < 4; ++s) v += B[s][r] * DB[s][c];
          K[16 * r + c] += w * v;
        }
      }
      for (int a = 0; a < 8; ++a) {
        for (int b = 0; b < 8; ++b) {
          const double m = w * mat.rho * p.N[a] * p.N[b];
          M[16 * (2 * a) + 2 * b] += m;
          M[16 * (2 * a + 1) + 2 * b + 1] += m;
        }
      }
    }
  }
}

std::array<double, 4> element_stress(const Mesh& mesh, int elem, const Material& mat, const std::vector<double>& u,
                                     double xi, double eta) {
  const PointEval p = evaluate(mesh, elem, xi, eta);
  const auto& e = mesh.elems[elem];
  double err = 0, ett = 0, ezz = 0, grz = 0;
  for (int a = 0; a < 8; ++a) {
    const double ur = u[2 * e[a]], uz = u[2 * e[a] + 1];
    err += p.dNr[a] * ur;
    ett += p.N[a] * ur;
    ezz += p.dNz[a] * uz;
    grz += p.dNz[a] * ur + p.dNr[a] * uz;
  }
  return stress_from_strain(mat, err, ett / p.r, ezz, grz);
}

std::vector<std::array<double, 4>> recovered_nodal_stress(const Mesh& mesh, const Material& mat,
                                                          const std::vector<double>& u) {
  const std::size_t nn = mesh.nodes.size();
  std::vector<std::array<double, 4>> sum(nn, {0, 0, 0, 0});
  std::vector<int> cnt(nn, 0);
  static const double cx[4] = {-1, 1, 1, -1}, cy[4] = {-1, -1, 1, 1};
  const double g = 1.0 / std::sqrt(3.0);
  for (int k = 0; k < static_cast<int>(mesh.elems.size()); ++k) {
    std::array<double, 4> sg[4];
    for (int q = 0; q < 4; ++q) sg[q] = element_stress(mesh, k, mat, u, cx[q] * g, cy[q] * g);
    for (int c = 0; c < 4; ++c) {
      // 가우스점 좌표를 ±1 로 늘린 쌍일차 보간을 꼭짓점 (s, t) = (±√3, ±√3) 에서 계산
      const double s = cx[c] * std::sqrt(3.0), t = cy[c] * std::sqrt(3.0);
      std::array<double, 4> v{0, 0, 0, 0};
      for (int q = 0; q < 4; ++q) {
        const double w = 0.25 * (1 + s * cx[q]) * (1 + t * cy[q]);
        for (int m = 0; m < 4; ++m) v[m] += w * sg[q][m];
      }
      const int n = mesh.elems[k][c];
      for (int m = 0; m < 4; ++m) sum[n][m] += v[m];
      ++cnt[n];
    }
  }
  for (std::size_t n = 0; n < nn; ++n) {
    if (cnt[n] > 0) {
      for (auto& v : sum[n]) v /= cnt[n];
    }
  }
  for (const auto& e : mesh.elems) {
    for (int s = 0; s < 4; ++s) {
      const auto l = edge_local_nodes(s);
      for (int m = 0; m < 4; ++m) sum[e[l[1]]][m] = 0.5 * (sum[e[l[0]]][m] + sum[e[l[2]]][m]);
    }
  }
  return sum;
}

}  // namespace hcw::structure
