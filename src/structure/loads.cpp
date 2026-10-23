#include "structure/loads.h"

#include <cmath>

#include "structure/element.h"

namespace hcw::structure {

void add_traction(const Mesh& mesh, const std::vector<Edge>& edges, const TractionFn& traction, Vec& F) {
  const auto& gp = gauss_points(3);
  const auto& gw = gauss_weights(3);
  for (const auto& ed : edges) {
    const auto l = edge_local_nodes(ed.side);
    const auto& e = mesh.elems[ed.elem];
    for (int g = 0; g < 3; ++g) {
      const double s = gp[g];
      // 변 위 2차 형상함수 (시작, 중간, 끝)
      const double N[3] = {0.5 * s * (s - 1), 1 - s * s, 0.5 * s * (s + 1)};
      const double dN[3] = {s - 0.5, -2 * s, s + 0.5};
      double r = 0, z = 0, dr = 0, dz = 0;
      for (int a = 0; a < 3; ++a) {
        const Node& n = mesh.nodes[e[l[a]]];
        r += N[a] * n.r;
        z += N[a] * n.z;
        dr += dN[a] * n.r;
        dz += dN[a] * n.z;
      }
      const double len = std::sqrt(dr * dr + dz * dz);
      // 반시계로 도는 변의 바깥 법선 = 접선을 시계 방향으로 90° 돌린 것
      const double nr = dz / len, nz = -dr / len;
      const auto t = traction(r, z, nr, nz);
      const double w = gw[g] * len * kTwoPi * r;
      for (int a = 0; a < 3; ++a) {
        const int node = e[l[a]];
        F[2 * node] += w * N[a] * t[0];
        F[2 * node + 1] += w * N[a] * t[1];
      }
    }
  }
}

void add_pressure(const Mesh& mesh, const std::vector<Edge>& edges, double p, Vec& F) {
  add_traction(mesh, edges,
               [p](double, double, double nr, double nz) { return std::array<double, 2>{-p * nr, -p * nz}; }, F);
}

void add_body_force(const Mesh& mesh, const BodyFn& body, Vec& F) {
  const auto& gp = gauss_points(3);
  const auto& gw = gauss_weights(3);
  for (int k = 0; k < static_cast<int>(mesh.elems.size()); ++k) {
    const auto& e = mesh.elems[k];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        const PointEval p = evaluate(mesh, k, gp[i], gp[j]);
        const auto b = body(p.r, p.z);
        const double w = gw[i] * gw[j] * p.detJ * kTwoPi * p.r;
        for (int a = 0; a < 8; ++a) {
          F[2 * e[a]] += w * p.N[a] * b[0];
          F[2 * e[a] + 1] += w * p.N[a] * b[1];
        }
      }
    }
  }
}

Vec rigid_z(const Model& model) {
  Vec ez(model.num_dofs(), 0.0);
  for (int n = 0; n < model.num_dofs() / 2; ++n) ez[2 * n + 1] = 1.0;
  return ez;
}

double total_mass(const Model& model) {
  const Vec ez = rigid_z(model);
  Vec Mez;
  model.apply_M(ez, Mez);
  return model.dot(ez, Mez);
}

void add_inertia_z(const Model& model, double A, Vec& F) {
  const Vec ez = rigid_z(model);
  Vec Mez;
  model.apply_M(ez, Mez);
  for (std::size_t i = 0; i < F.size(); ++i) F[i] -= A * Mez[i];
}

double add_inertia_relief_z(const Model& model, Vec& F) {
  double Fz = 0.0;
  for (int n = 0; n < model.num_dofs() / 2; ++n) Fz += F[2 * n + 1];
  const double A = Fz / total_mass(model);
  add_inertia_z(model, A, F);
  return A;
}

}  // namespace hcw::structure
