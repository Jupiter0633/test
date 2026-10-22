#pragma once

// P3 검증 시험 (단계 2: 자유수면). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.3 (기준 고정 커밋 d74d3e8) 를 그대로 따른다.
// 2-A 정지 수면, 2-B 순수 수송, 2-C 정재파, 2-D 수주 붕괴 (보고만), 2-F 공기 스프링.

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_2b(bool quick);
TestReport run_2a(bool quick);
TestReport run_2c(bool quick);
TestReport run_2d(bool quick);  // 보고만
TestReport run_2f(bool quick);

// 축 위 중심 (0, zc), 반지름 R 인 구가 셀 [r0, r1] × [z0, z1] 에서 차지하는 r 가중 체적 (2π 제외)
double sphere_cell_volume(double r0, double r1, double z0, double z1, double zc, double radius);

}  // namespace hcw::verify
