#pragma once

// P4 검증 시험 (단계 3: 지정 운동의 정수면 입수, 강체). 조건과 통과 기준은
// docs/design/03_verification_plan.md §1.4 (기준 고정 커밋 60d72a8, 3-F 제어면 변경 eecd38c,
// 늘어나는 격자·Δt × ¼ d753d5a / acd6130) 를 따른다.
// 3-A 쐐기 (보고만), 3-B 원판 정성, 3-C 격자 수렴 (GCI), 3-D 상한 (각 보고서 안), 3-E 공기 모델 (보고만),
// 3-F 운동량 수지, 3-G 정지 부력, 3-I 옆벽 효과 (보고만). 3-H 는 미실행.

#include "verify/p1_cases.h"

namespace hcw::verify {

TestReport run_3a(bool quick);
TestReport run_3b(bool quick);
TestReport run_3c(bool quick);
TestReport run_3e(bool quick);
TestReport run_3f(bool quick);
// plan 이 있으면 늘어나는 격자 (S-5)
TestReport run_3g(bool quick, const GridPlan* plan = nullptr);
// 늘어나는 격자 시험 S-1 – S-5 (03 §1.4). id: "S-1" … "S-5"
TestReport run_s(const std::string& id, bool quick);
TestReport run_3i(bool quick);

// 개발용 (판정 아님, 결과 저장 안 함): 3-B 조건 원판 입수를 a/Δ = n, Δt 규칙 × dt_factor 로 계산한다.
// 3-F 의 제어면 힘도 함께 구해 충격량·첨두값·두 힘의 차이와 시간 이력을 돌려준다 (verdict = DEV_ONLY).
// 시간 간격 수렴 실험 (P4-analysis-2026-10-23.md) 용.
// stretched = true 이면 P4 새 조건의 늘어나는 격자 (n = a/Δ_core). dt_factor 는 그대로 곱한다.
// water_floor = false 이면 물 압력 하한을 끈다 (3-F 원인 비교용)
TestReport entry_dev_study(int n, double dt_factor, bool stretched = false, bool water_floor = true);

// 개발용 그림 자료: 위와 같은 원판 입수를 계산하며 interval [s] 마다 셀 값 (snap_NNN.csv) 과
// 끝에 시간 이력 (history.csv: t, 힘, 바닥 평균 압력, 포켓 압력) 을 dir 에 쓴다. 그림은 tools/plot_entry.py.
// water_floor = false 이면 물 압력 하한을 끈다 (개발 비교용)
TestReport entry_dev_snapshots(int n, double dt_factor, const std::string& dir, double interval,
                               bool water_floor = true, bool stretched = false);

// Celik 등의 GCI_fine (안전계수 1.25, 격자비 2). phi1: 가장 고운 격자, phi3: 가장 거친 격자.
// 차이 비 R = (phi2 − phi1)/(phi3 − phi2) 가 0 < R < 1 이 아니면 false.
// 관측 차수 p 는 [0.5, 3] 으로 자른다 (03 §1.4).
bool gci_fine(double phi1, double phi2, double phi3, double& gci, double& p);

}  // namespace hcw::verify
