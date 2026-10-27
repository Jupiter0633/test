#pragma once

// P5 검증 시험 (단계 4: 구조 동해석, 유체 없음). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.5 (기준 고정 커밋 9bf931b) 를 그대로 따른다.
// 4-A 패치·정적 제조해, 4-B Lamé 원통, 4-C 원판 굽힘, 4-E 고유진동수, 4-F 에너지 보존,
// 4-G 동적 제조해, 4-H 한 모드 이산 정확해, 4-I 관성력·관성 해제. 4-D 는 미실행 (D-26).

#include <array>

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_4a(bool quick);
TestReport run_4b(bool quick);
TestReport run_4c(bool quick);
TestReport run_4e(bool quick);
TestReport run_4e_c(bool quick);  // 보고만
TestReport run_4f(bool quick);
TestReport run_4g(bool quick);
TestReport run_4h(bool quick);
TestReport run_4i(bool quick);

// 두꺼운 원통 (안 a, 바깥 b) 의 평면 변형 반경 진동수 ω (E = ρ = 1, ν) 가운데 낮은 쪽부터 count 개.
// 진동수 방정식은 03 §1.5 4-E(b) (Bessel J₁, Y₁).
std::vector<double> hollow_cylinder_radial_omegas(double a, double b, double nu, int count);

// 1 자유도 HHT-α (α = 0 → Newmark) 의 증폭 행렬 지배 근 ρ e^{±iθ}. Ω = ωΔt, ξ = 감쇠비.
// 또 하나의 (실수) 근을 spurious 에 쓴다.
void hht_dominant_root(double Omega, double alpha, double xi, double& rho, double& theta, double* spurious = nullptr);

}  // namespace hcw::verify
