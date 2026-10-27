#pragma once

// 고유진동수: K φ = ω² M φ 의 가장 낮은 모드들 (지정 자유도는 0).
// 이동-역반복 부분공간 반복: X ← (K − σM)⁻¹ M X, 그리고 Rayleigh–Ritz.
// 멈춤 조건 (03 §1.5 공통 정의): 구하려는 모든 고유값의 상대 변화 ≤ 1e−12.

#include <string>
#include <vector>

#include "structure/assembly.h"

namespace hcw::structure {

struct EigenResult {
  std::vector<double> omega2;  // 오름차순 ω²
  std::vector<Vec> modes;      // 전역 벡터, φᵀMφ = 1
  int iterations = 0;
  bool converged = false;
  SolveCheck check;  // 선형 풀이 확인값의 최댓값
};

// shift σ 는 K − σM 이 양의 정부호가 되도록 (지지가 없으면 음수) 준다.
// known 은 미리 아는 모드 (예: 지지가 없을 때 z 강체 이동). 매 반복에서 이들과 M-직교가 되게 빼서
// 구하는 모드에서 제외한다.
EigenResult lowest_modes(const Model& model, int count, double shift, int max_iterations = 2000,
                         const std::vector<Vec>& known = {});

// 한 모드를 역반복으로 다듬는다 (이동 σ = λ(1 − 1e−4)). 4-H 의 "한 모드만 들뜬" 처음 조건용.
// 결과 φ 는 φᵀMφ = 1, 돌려주는 값은 Rayleigh 몫 φᵀKφ.
double refine_mode(const Model& model, double omega2, Vec& mode, int iterations = 6);

// 작은 대칭 행렬 고유값 (Jacobi). a 는 n×n 행 우선, 결과 고유벡터는 v 의 열.
void jacobi_eigen(int n, std::vector<double>& a, std::vector<double>& v);

}  // namespace hcw::structure
