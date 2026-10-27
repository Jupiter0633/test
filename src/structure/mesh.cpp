#include "structure/mesh.h"

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace hcw::structure {

std::array<int, 3> edge_local_nodes(int side) {
  static const int t[4][3] = {{0, 4, 1}, {1, 5, 2}, {2, 6, 3}, {3, 7, 0}};
  return {t[side][0], t[side][1], t[side][2]};
}

void MeshBuilder::add_block(double r0, double r1, double z0, double z1, int nr, int nz) {
  if (nr < 1 || nz < 1 || !(r1 > r0) || !(z1 > z0) || r0 < 0.0) {
    throw std::invalid_argument("MeshBuilder::add_block: bad block");
  }
  blocks_.push_back({r0, r1, z0, z1, nr, nz});
}

Mesh MeshBuilder::build() const {
  Mesh mesh;
  // 같은 위치의 절점을 합치는 열쇠: 좌표를 1e−9 단위로 반올림
  std::map<std::pair<long long, long long>, int> index;
  auto node_at = [&](double r, double z) {
    const auto key = std::make_pair(std::llround(r * 1e9), std::llround(z * 1e9));
    const auto it = index.find(key);
    if (it != index.end()) return it->second;
    const int id = static_cast<int>(mesh.nodes.size());
    mesh.nodes.push_back({r, z});
    index.emplace(key, id);
    return id;
  };
  for (const auto& b : blocks_) {
    // 절점 격자 (2nr+1) × (2nz+1), 위치는 반 요소 간격
    auto lattice = [&](int i, int j) {
      const double r = i == 2 * b.nr ? b.r1 : b.r0 + (b.r1 - b.r0) * i / (2.0 * b.nr);
      const double z = j == 2 * b.nz ? b.z1 : b.z0 + (b.z1 - b.z0) * j / (2.0 * b.nz);
      return node_at(r, z);
    };
    for (int j = 0; j < b.nz; ++j) {
      for (int i = 0; i < b.nr; ++i) {
        const int i0 = 2 * i, j0 = 2 * j;
        mesh.elems.push_back({lattice(i0, j0), lattice(i0 + 2, j0), lattice(i0 + 2, j0 + 2), lattice(i0, j0 + 2),
                              lattice(i0 + 1, j0), lattice(i0 + 2, j0 + 1), lattice(i0 + 1, j0 + 2),
                              lattice(i0, j0 + 1)});
      }
    }
  }
  return mesh;
}

std::vector<Edge> boundary_edges(const Mesh& mesh) {
  // 변은 중간 절점으로 구별된다 (중간 절점은 많아야 두 요소가 공유).
  std::map<int, int> count;
  for (const auto& e : mesh.elems) {
    for (int s = 0; s < 4; ++s) ++count[e[edge_local_nodes(s)[1]]];
  }
  std::vector<Edge> out;
  for (int k = 0; k < static_cast<int>(mesh.elems.size()); ++k) {
    for (int s = 0; s < 4; ++s) {
      if (count[mesh.elems[k][edge_local_nodes(s)[1]]] == 1) out.push_back({k, s});
    }
  }
  return out;
}

std::vector<Edge> boundary_edges_where(const Mesh& mesh, const std::function<bool(double, double)>& pred) {
  std::vector<Edge> out;
  for (const auto& ed : boundary_edges(mesh)) {
    bool all = true;
    for (int a : edge_local_nodes(ed.side)) {
      const auto& n = mesh.nodes[mesh.elems[ed.elem][a]];
      all = all && pred(n.r, n.z);
    }
    if (all) out.push_back(ed);
  }
  return out;
}

void recenter_midside_nodes(Mesh& mesh) {
  for (const auto& e : mesh.elems) {
    for (int s = 0; s < 4; ++s) {
      const auto l = edge_local_nodes(s);
      const Node& a = mesh.nodes[e[l[0]]];
      const Node& b = mesh.nodes[e[l[2]]];
      mesh.nodes[e[l[1]]] = {0.5 * (a.r + b.r), 0.5 * (a.z + b.z)};
    }
  }
}

}  // namespace hcw::structure
