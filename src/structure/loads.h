#pragma once

// 하중 → 절점력 (일관 하중). 모두 2πr 를 곱한 전체 원둘레 값이며 전역 벡터 F 에 더한다.
// - 경계 변의 표면력 t(r, z): ∫ N t 2πr ds (변마다 3 점 가우스)
// - 압력 p: 표면력 −p n (n 은 바깥 법선)
// - 체적력 b(r, z): ∫ N b 2πr dA (3×3 가우스)
// - 지정 운동 관성력: 체적력 −ρA (A 는 z 방향 강체 가속도)
// - 관성 해제: 합력과 질량에서 A 를 구해 −ρA 를 더함 (D-11: 두 방식 모두 제공)

#include <array>
#include <functional>
#include <vector>

#include "structure/assembly.h"
#include "structure/mesh.h"

namespace hcw::structure {

// traction(r, z, n_r, n_z) → (t_r, t_z)
using TractionFn = std::function<std::array<double, 2>(double, double, double, double)>;
using BodyFn = std::function<std::array<double, 2>(double, double)>;

void add_traction(const Mesh& mesh, const std::vector<Edge>& edges, const TractionFn& traction, Vec& F);
void add_pressure(const Mesh& mesh, const std::vector<Edge>& edges, double p, Vec& F);
void add_body_force(const Mesh& mesh, const BodyFn& body, Vec& F);

// z 방향 강체 이동 벡터 e_z (모든 절점 u_z = 1)
Vec rigid_z(const Model& model);

// 전체 질량 m = e_zᵀ M e_z
double total_mass(const Model& model);

// 지정 운동: 몸체 좌표계에서 z 가속도 A 의 관성력 −A M e_z 를 F 에 더한다 (일관 질량)
void add_inertia_z(const Model& model, double A, Vec& F);

// 관성 해제: F 의 z 합력 / 질량 = A 를 구하고 관성력을 더한다. A 를 돌려준다.
double add_inertia_relief_z(const Model& model, Vec& F);

}  // namespace hcw::structure
