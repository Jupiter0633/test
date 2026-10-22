#pragma once

// P2 검증 시험 1-E, 1-F, 1-G, 1-H. 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.2 (기준 고정 커밋 93cf47d, 1-H 변경 d4c2574) 를 그대로 따른다.
// 1-E(b) (이동 중 힘 이력)는 보고만 한다.

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_1e(bool quick);
TestReport run_1f(bool quick);
TestReport run_1g(bool quick);
TestReport run_1h(bool quick);

}  // namespace hcw::verify
