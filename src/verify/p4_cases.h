#pragma once

// P4 검증 시험 (단계 3: 지정 운동의 정수면 입수, 강체). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.4 (기준 고정 커밋 60d72a8, 3-F 제어면 변경 eecd38c) 를 따른다.
// 3-A 쐐기 (보고만), 3-B 원판 정성, 3-C 격자 수렴 (GCI), 3-D 상한 (각 보고서 안), 3-E 공기 모델 (보고만),
// 3-F 운동량 수지, 3-G 정지 부력, 3-I 옆벽 효과 (보고만). 3-H 는 미실행.

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_3a(bool quick);
TestReport run_3b(bool quick);
TestReport run_3c(bool quick);
TestReport run_3e(bool quick);
TestReport run_3f(bool quick);
TestReport run_3g(bool quick);
TestReport run_3i(bool quick);

// Celik 등의 GCI_fine (안전계수 1.25, 격자비 2). phi1: 가장 고운 격자, phi3: 가장 거친 격자.
// 차이 비 R = (phi2 − phi1)/(phi3 − phi2) 가 0 < R < 1 이 아니면 false.
// 관측 차수 p 는 [0.5, 3] 으로 자른다 (03 §1.4).
bool gci_fine(double phi1, double phi2, double phi3, double& gci, double& p);

}  // namespace hcw::verify
