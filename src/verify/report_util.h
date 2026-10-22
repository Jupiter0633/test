#pragma once

// 검증 시험 공용 도구: 시간 측정, 실행 마무리, 판정 기준(단조 감소, 관측 차수), 최종 판정.
// P1 (p1_cases.cpp) 에서 옮겨 P2 와 같이 쓴다. 판정 논리는 P1 판정 때와 같다.

#include <chrono>
#include <string>
#include <vector>

#include "flow/solver.h"
#include "verify/p1_cases.h"

namespace hcw::verify {

double seconds_since(std::chrono::steady_clock::time_point t0);

// 종료 시각까지 Δt 로 진행하고 마지막 단계만 남은 시간으로 맞춘다.
void advance_to(flow::FlowSolver& s, double t_end, double dt);

// 해법 통계를 옮기고 해법 실패·비유한 값이면 실패 사유를 적는다.
void finish_run(GridRun& run, const flow::FlowSolver& s);

std::string fmt(double v);

// 오차가 격자 세분화마다 줄어드는지
Criterion monotone(const std::vector<GridRun>& runs, const std::string& key);

// 가장 고운 두 격자의 관측 차수 ≥ min
Criterion order_at_least(const std::vector<GridRun>& runs, const std::string& key, double min);

// 실행 실패가 있으면 그 사유, 아니면 기준 전부 통과 시 VERIFIED, 아니면 FAILED_CRITERIA.
void conclude(TestReport& rep, bool quick);

}  // namespace hcw::verify
