#pragma once

// 축대칭 Q8 요소: 강성·일관 질량·응력. 선형 탄성, 소변형 (docs/design/02_design_document.md §3.9).
// 미지수 순서: 요소 절점 a 의 (u_r, u_z) = 요소 자유도 2a, 2a+1. 전역 자유도는 2·절점 + 성분.
// 모든 체적·면적분에 2πr 를 곱한다 (원둘레 한 바퀴 전체의 값).
// 응력 성분 순서: (σ_rr, σ_θθ, σ_zz, σ_rz).

#include <array>
#include <vector>

#include "structure/mesh.h"

namespace hcw::structure {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

struct Material {
  double E = 1.0;
  double nu = 0.3;
  double rho = 1.0;
  double lambda() const { return E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)); }
  double mu() const { return E / (2.0 * (1.0 + nu)); }
};

// 1 차원 가우스 점 (n = 2, 3)
const std::vector<double>& gauss_points(int n);
const std::vector<double>& gauss_weights(int n);

// 요소 안 한 점 (ξ, η) 의 형상함수와 그 r, z 미분
struct PointEval {
  double r = 0.0, z = 0.0, detJ = 0.0;
  std::array<double, 8> N{}, dNr{}, dNz{};
};

void shape_functions(double xi, double eta, double* N, double* dNxi, double* dNeta);
PointEval evaluate(const Mesh& mesh, int elem, double xi, double eta);

// 16×16 행 우선 강성·질량 (3×3 가우스, D-28)
void element_matrices(const Mesh& mesh, int elem, const Material& mat, double* K, double* M);

// 전역 변위 u (크기 2·절점) 에서 요소 점 (ξ, η) 의 응력
std::array<double, 4> element_stress(const Mesh& mesh, int elem, const Material& mat, const std::vector<double>& u,
                                     double xi, double eta);

// 변위로부터 계산한 응력 (해석해의 변형률을 넣을 때 쓰는 공용 식)
std::array<double, 4> stress_from_strain(const Material& mat, double err, double ett, double ezz, double grz);

// 절점 회복 응력 (03 §1.5 공통 정의): 요소 2×2 가우스점 응력을 꼭짓점까지 1차식으로 외삽하고
// 이웃 요소 값을 평균한다. 중간 절점은 양 끝 꼭짓점 (평균된) 값의 평균.
std::vector<std::array<double, 4>> recovered_nodal_stress(const Mesh& mesh, const Material& mat,
                                                          const std::vector<double>& u);

}  // namespace hcw::structure
