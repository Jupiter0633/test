#pragma once

// P5 구조 격자: r–z 단면의 8 절점 2차 사각형 (Q8) 요소.
// 기준: docs/design/03_verification_plan.md §1.5 (기준 고정 커밋 9bf931b), 요소 결정 D-28.
//
// 요소 안 절점 순서 (ξ, η):
//   꼭짓점 0 (−1,−1), 1 (1,−1), 2 (1,1), 3 (−1,1)  — r 를 가로, z 를 세로로 볼 때 반시계 방향
//   중간 4 (0,−1), 5 (1,0), 6 (0,1), 7 (−1,0)
// 변 k 의 세 절점: 변 0 = (0,4,1), 1 = (1,5,2), 2 = (2,6,3), 3 = (3,7,0). 반시계로 돈다.

#include <array>
#include <functional>
#include <vector>

namespace hcw::structure {

struct Node {
  double r = 0.0, z = 0.0;
};

struct Mesh {
  std::vector<Node> nodes;
  std::vector<std::array<int, 8>> elems;
};

// 요소의 한 변 (경계 하중, 경계 조건 지정용)
struct Edge {
  int elem = 0;
  int side = 0;  // 0..3
};

// 변 side 의 세 절점 (요소 안 번호: 시작 꼭짓점, 중간, 끝 꼭짓점)
std::array<int, 3> edge_local_nodes(int side);

// 직사각형 블록 여러 개를 붙여 격자를 만든다. 같은 위치의 절점은 하나로 합친다.
class MeshBuilder {
 public:
  // [r0, r1] × [z0, z1] 을 nr × nz 요소로 균일 분할
  void add_block(double r0, double r1, double z0, double z1, int nr, int nz);
  Mesh build() const;

 private:
  struct Block {
    double r0, r1, z0, z1;
    int nr, nz;
  };
  std::vector<Block> blocks_;
};

// 한 요소에만 속한 변 (바깥 경계)
std::vector<Edge> boundary_edges(const Mesh& mesh);

// 바깥 경계 변 가운데 세 절점이 모두 pred(r, z) 를 만족하는 변
std::vector<Edge> boundary_edges_where(const Mesh& mesh, const std::function<bool(double, double)>& pred);

// 중간 절점을 양 끝 꼭짓점의 가운데로 다시 놓는다 (꼭짓점을 옮긴 뒤 변을 직선으로 유지).
void recenter_midside_nodes(Mesh& mesh);

}  // namespace hcw::structure
