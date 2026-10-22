#pragma once

// P3 검증 시험 (단계 2: 자유수면). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.3 (기준 고정 커밋 d74d3e8) 를 그대로 따른다.
// 지금은 2-B (순수 수송) 만 있다. 2-A, 2-C, 2-D, 2-F 는 2상 유동 결합 후 추가한다.

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_2b(bool quick);

// 축 위 중심 (0, zc), 반지름 R 인 구가 셀 [r0, r1] × [z0, z1] 에서 차지하는 r 가중 체적 (2π 제외)
double sphere_cell_volume(double r0, double r1, double z0, double z1, double zc, double radius);

}  // namespace hcw::verify
