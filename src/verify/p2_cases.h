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

// 개발용 (판정 아님): 1-G (정지 구 Stokes 항력) 와 1-E (가속하는 구, (b) 이동 구 포함) 를
// 절단 셀 스위치 (SolverOptions::cut_cell, docs/design/07_cut_cell_momentum_design.md) 를 켜거나 꺼서 계산한다.
TestReport dev_p2_cut_cell(bool quick, bool cut_cell);

}  // namespace hcw::verify
