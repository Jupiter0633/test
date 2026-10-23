#pragma once

// 구조 모델: 격자 + 재료 + 지정 자유도. 요소 행렬을 저장해 두고
// (1) 전역 행렬 곱 K x, M x 를 요소별로 계산하고 (모든 자유도, 반력 계산에도 씀),
// (2) 자유 자유도끼리의 c_K K + c_M M 을 스카이라인 행렬로 조립한다.
// 자유도 번호는 절점 그래프의 역 Cuthill–McKee 순서로 정해 띠 폭을 줄인다.

#include <functional>
#include <vector>

#include "structure/band_solver.h"
#include "structure/element.h"
#include "structure/mesh.h"

namespace hcw::structure {

using Vec = std::vector<double>;

class Model {
 public:
  Model(Mesh mesh, Material mat);

  const Mesh& mesh() const { return mesh_; }
  const Material& material() const { return mat_; }
  int num_dofs() const { return 2 * static_cast<int>(mesh_.nodes.size()); }

  // 자유도 지정 (값은 풀 때 전역 벡터로 따로 준다). finalize() 전에만.
  void fix(int node, int comp);
  // pred(r, z) 인 모든 절점의 성분 comp (0 = r, 1 = z) 를 지정
  void fix_where(const std::function<bool(double, double)>& pred, int comp);
  // 축 (r = 0) 절점의 u_r 지정 (축대칭 조건)
  void fix_axis();
  bool is_fixed(int dof) const { return fixed_[dof]; }
  const std::vector<int>& fixed_dofs() const { return fixed_list_; }

  // 요소 행렬 계산, 자유도 번호, 스카이라인 모양 결정
  void finalize();
  int num_free() const { return num_free_; }
  int equation(int dof) const { return eq_[dof]; }  // 지정 자유도는 −1

  // 전역 (모든 자유도) 곱: y = K x, y = M x
  void apply_K(const Vec& x, Vec& y) const;
  void apply_M(const Vec& x, Vec& y) const;
  // y = cK K x + cM M x
  void apply(double cK, double cM, const Vec& x, Vec& y) const;

  // 자유 자유도끼리의 cK K + cM M (분해 전)
  SkylineMatrix assemble(double cK, double cM) const;

  // 전역 ↔ 자유 벡터
  Vec restrict_free(const Vec& full) const;
  void scatter_free(const Vec& free, Vec& full) const;  // 자유 자리만 덮어씀

  double dot(const Vec& a, const Vec& b) const;

 private:
  Mesh mesh_;
  Material mat_;
  std::vector<bool> fixed_;
  std::vector<int> fixed_list_;
  std::vector<int> eq_;
  std::vector<int> dof_of_eq_;
  std::vector<double> ke_, me_;  // 요소마다 256 개
  std::vector<int> top_;
  int num_free_ = 0;
  bool finalized_ = false;
};

// 정적 문제: 자유 자유도에서 K u = F, 지정 자유도 값은 u_prescribed (전역 벡터의 지정 자리).
// 상대 잔차 ||K u − F||_free / ||F − K u_p||_free 를 rel_residual 에 쓴다.
// 1e−10 을 넘으면 호출한 쪽이 FAILED_SOLVER 로 처리한다.
Vec solve_static(const Model& model, const Vec& F, const Vec& u_prescribed, double* rel_residual);

// 분해된 행렬과 그 행렬의 cK, cM 으로 자유 자유도 해를 구하고 잔차를 확인한다.
// rhs, 해 모두 자유 자유도 벡터.
Vec solve_checked(const Model& model, const SkylineMatrix& factored, double cK, double cM, const Vec& rhs,
                  double* rel_residual);

// 반력 R = K u − F (전역 벡터, 자유 자리 값은 평형 잔차)
Vec reactions(const Model& model, const Vec& u, const Vec& F);

}  // namespace hcw::structure
