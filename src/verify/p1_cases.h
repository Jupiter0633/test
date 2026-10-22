#pragma once

// P1 검증 시험 1-A ~ 1-D. 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.1 (기준 고정 커밋 35a77ac) 을 그대로 따른다.

#include <map>
#include <string>
#include <vector>

#include "flow/solver.h"

namespace hcw::verify {

struct GridRun {
  std::string label;                      // 예: "N=64", "planar N=32 aligned"
  int n = 0;
  std::map<std::string, double> metrics;  // 오차 등
  flow::SolverStats stats;
  double seconds = 0.0;
  std::string failure;                    // 실행 실패 사유 (비어 있으면 정상)
};

struct Criterion {
  std::string text;
  bool passed = false;
  std::string detail;
};

struct TestReport {
  std::string id;
  std::string title;
  // 기준 위치: docs/design/03_verification_plan.md 의 절과 기준 고정 커밋
  std::string phase = "P1";
  std::string criteria_section = "1.1";
  std::string criteria_commit = "35a77ac";
  std::vector<GridRun> runs;
  std::vector<Criterion> criteria;
  std::string verdict;  // VERIFIED / FAILED_<사유> / QUICK_CHECK_ONLY
  std::vector<std::string> notes;
};

// quick = true 이면 작은 격자로만 실행 (동작 확인용, 판정 아님).
TestReport run_1a(bool variable_density, bool quick);
TestReport run_1b(bool quick);
TestReport run_1c(bool quick);
TestReport run_1d(bool quick);

// 한 격자 실행 (조건은 위 시험과 같고 대류 방식만 고름). P3 2-H 가 새 대류 방식으로 다시 쓴다.
// run_1a_grid 는 예외를 그대로 던진다.
GridRun run_1a_grid(int n, bool variable_density, flow::Convection convection);
GridRun run_1c_grid(int n, flow::Convection convection);

// p_obs = ln(e_coarse / e_fine) / ln 2
double observed_order(double e_coarse, double e_fine);

}  // namespace hcw::verify
