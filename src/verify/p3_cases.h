#pragma once

// P3 검증 시험 (단계 2: 자유수면). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.3 (기준 고정 커밋 d74d3e8) 를 그대로 따른다.
// 2-A 정지 수면, 2-B 순수 수송, 2-C 정재파, 2-D 수주 붕괴 (안정성), 2-F 공기 스프링, 2-H 제한자 대류 정확도.
// 2-D·2-H·(lim) 재실행의 기준은 같은 절의 추가분 (기준 고정 8ffed34).

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_2b(bool quick);
TestReport run_2a(bool quick, bool limiter = false);
TestReport run_2c(bool quick, bool limiter = false);
TestReport run_2d(bool quick);  // 2026-10-23 판정 시험으로 변경 (제한자 대류)
TestReport run_2f(bool quick, bool limiter = false);
// 제한자 대류의 정확도 (1-C, 1-A(b) 재풀이). limiter = true 인 2-A·2-C·2-F 는 id 에 "(lim)" 이 붙는다.
TestReport run_2h(bool quick);

// 축 위 중심 (0, zc), 반지름 R 인 구가 셀 [r0, r1] × [z0, z1] 에서 차지하는 r 가중 체적 (2π 제외)
double sphere_cell_volume(double r0, double r1, double z0, double z1, double zc, double radius);

}  // namespace hcw::verify
