# 인수인계 문서 (새 세션은 이 파일부터 읽으세요)

작성: 2026-10-21, 클라우드 세션에서 노트북 로컬 세션으로 옮기기 위해 작성.
이 대화의 내용은 새 세션에 전달되지 않습니다. **필요한 모든 결정과 규칙은 이 저장소 문서에 있습니다.**

---

## 1. 프로젝트 한 줄 요약

논문용 정량 분석을 위한 독립 입수충격 해석기 `hcw`. 축이 수직인 **바닥이 막힌 중공 원통**이 평평한 바닥판부터 입수할 때 압력 분포, 충격력, 충격량, 구조응력 시간 이력을 계산한다. 1차 버전은 2차원 축대칭 정수면 입수이고, 최종 목표는 3차원 규칙파 중 입수다.

## 2. 사용자에 대해 (진행 방식)

- **개발 경험이 거의 없습니다.** 쉬운 한국어로 설명하고, 용어는 처음 나올 때 풀어 쓰고, 명령은 그대로 복사해 실행할 수 있게 주세요.
- 결정이 필요한 것은 선택지·장단점·추천을 표로 제시하고 승인을 받은 뒤 진행합니다.
- 설치·설정 변경은 사용자가 직접 합니다. 필요하면 단계별로 안내합니다.
- 프로그램 출력 메시지는 Windows 콘솔 한글 깨짐을 피하려고 **영어**로 합니다. 문서와 주석은 한국어입니다.

## 3. 반드시 지킬 규칙 (사용자가 처음에 정한 원칙)

1. **코드가 실행되는 것과 물리적 검증 통과를 구분합니다.** 상태: `RUN_COMPLETED` / `RUN_HEALTHY` / `VERIFIED` / `VALIDATED` / `FAILED_<사유>` (`docs/design/02_design_document.md` §5).
2. **통과 기준은 실행 전에 커밋으로 고정합니다.** 결과를 본 뒤 기준을 완화하지 않습니다. 불통과는 `FAILED_*` 로 남기고 지우지 않습니다. 기준 변경이 꼭 필요하면 사유·날짜·변경 전 판정을 함께 남깁니다.
3. **실제로 확인한 문헌만 근거로 씁니다.** 확인 수준은 `docs/design/06_references.md`. 원문을 못 본 수치(표, 그래프 값)는 통과 기준으로 쓰지 않습니다.
4. **개발용 임시 입력과 논문 조건을 분리합니다.** `cases/dev/` (`DEV_` 접두어) / `cases/verification/` / `cases/paper/` (SHA-256 잠금). `hcw` 프로그램이 이를 강제합니다.
5. **논문 결과는 워크스테이션에서만** (`config/machines.ini` 에 `workstation` 으로 등록된 기기 + 수정 없는 커밋 빌드). 노트북은 `development` 로 등록되어 있습니다.
6. 각 단계는 "구현 완료"와 "검증 통과"를 별도로 확인하고, 판정 기록을 `validation/records/` 에 남깁니다 (템플릿: `TEMPLATE.md`).

## 4. 환경

| 기기 | 역할 | 상태 |
|---|---|---|
| 노트북 `DESKTOP-4JJOTVL` | **개발 전용** | Windows 11 Pro, i7-1165G7 (4코어/8스레드), 16 GB, GPU 없음(CUDA 불가). Visual Studio Community 2026 (MSVC 19.51, CMake 4.3.1), Git for Windows, GitHub Desktop, Python 3.12 |
| 워크스테이션 | 논문 계산 | **Windows** (사용자 확정). CPU·RAM·GPU 미확인, 현재 접근 불가 |

- 빌드는 반드시 **"x64 Native Tools Command Prompt for VS 2026"** 에서 (Developer PowerShell 은 32비트라 빌드 설정이 일부러 오류를 냄).
- 저장소 위치(노트북): `%USERPROFILE%\Documents\GitHub\test`
- 작업 브랜치: `claude/hollow-cylinder-water-entry-design-c9vq6i` (main 에 아직 병합하지 않음, PR 없음).

```bat
cmake --preset windows-x64-release
cmake --build --preset windows-x64-release
ctest --preset windows-x64-release
out\build\windows-x64-release\hcw.exe info
out\build\windows-x64-release\hcw_verify.exe all
```

## 5. 확정된 결정 (`docs/design/05_open_decisions.md` §A)

| 항목 | 결정 |
|---|---|
| 수치해석 방법 | 고정 격자 유한체적 + 기하학적 VOF + 날카로운 가상경계법(IBM) (SPH 미채택) |
| 공기 압축성 | 압축성 공기 + 비압축 물 (압력식에 압축 항을 넣되 0 으로 두면 비압축). 물 압축성은 ρ_w c_w V 상한 감시 |
| 이동 경계 | IBM + 물체 고정 좌표계 교차검증 모드 (축대칭 지정운동 한정) |
| 언어 | C++17 + CMake (+ 이후 OpenMP, 3D 에서 MS-MPI), 후처리는 Python |
| 구조 모델 | 자체 축대칭 고체 유한요소. 상용 FE 라이선스 없음 → 교차검증 수단 미정 (D-26) |
| 개발 도구 | Visual Studio Community 2026 |
| 검증 1-C | Ghia 리드구동 공동 표를 원문으로 확인하지 못해 **Taylor–Green 와류로 교체** (승인됨). Ghia 비교는 원문 확보 시까지 보류 (사용자는 원문 없음) |

## 6. 진행 현황

| 단계 | 상태 | 근거 |
|---|---|---|
| 설계 | 완료 | `docs/design/` 01–06 |
| P0 (빌드 틀, 입력 검사, 잠금, 실행 기록) | **완료** — Linux·Windows 노트북 모두 시험 18/18 통과 | `validation/records/P0-exit-2026-10-21.md` |
| P1 (단상·변수밀도 유동, 압력 해법) | **구현 완료 (커밋 f330e84), 판정 대기** | 기준 고정 커밋 `35a77ac`, `03_verification_plan.md` §1.1 |
| P2 이후 | 미착수 | `04_development_phases.md` |

### P1 세부 상태
- 클라우드(Linux)에서 `hcw_verify all` 전체 실행 → **2시간 한도로 중단, 판정 결과 없음.** 원인: N=256 에서 고정 기준의 압력 해법 허용치(1e−12/1e−13)가 배정밀도 반올림 한계 아래라 도달 불가. 기록: `validation/records/P1-attempt-2026-10-21-linux.md`.
- **P1 판정 전에 사용자 결정 필요: D-27** (`05_open_decisions.md`). 결정 없이 기준(허용치)을 바꾸지 말 것.
- 정체 감지 수정 후 N=256 은 단계당 약 6.5 초(클라우드 기준)로 진행되지만 `FAILED_SOLVER` 로 표시됨.
- 축소판(동작 확인용, 판정 아님)에서 1-A(b) 압력 L∞ 차수가 거친 격자에서 1.46 (기준 1.5). 전체 격자 결과로 판정하며, 불통과면 그대로 `FAILED_CRITERIA` 로 기록.
- 축소판 확인 중 찾아 고친 버그 2개 (커밋 f330e84 에 포함, 회귀 시험 있음):
  1. 압력 정칙화 셀을 첫 셀에 두면 밀도비 1000 에서 IC(0) 마지막 피벗이 무너짐 → 마지막 셀로 이동.
  2. 이전 해를 초기 추정으로 쓸 때 0 보다 나쁜 추정이면 수렴 불가 → 0 에서 재시작.

### P1 구현상 알려진 한계 (다음 단계에서 다룰 것)
- 압력 해법은 IC(0) 전처리 켤레기울기법 — 256² 정도까지만 실용적. 큰 격자에는 다중격자로 교체 필요.
- 균일 격자만 지원. 국소 세분화는 **P4 로 이동** (코드 작성 전 결정, `03` §1.1 범위 메모).
- 점성항 명시 처리(공기상에서 시간 간격 제한이 커질 수 있음), 대류 2차 중심차분(자유수면용 제한자 방식은 P3 에서 별도 제조해 시험).

## 7. 다음 할 일 (순서대로)

1. GitHub Desktop: **Fetch origin → Pull origin**.
2. 노트북에서 빌드 + `ctest` (빠른 시험 23개 통과 확인).
3. **D-27 결정을 먼저 받은 뒤**, 노트북에서 `hcw_verify.exe all` 실행 (수십 분 이상). `results\verification\P1\...\report.txt` 를 근거로 판정 기록 `validation/records/P1-exit-<날짜>.md` 작성. 기준은 절대 바꾸지 말 것.
4. P1 판정에 따라:
   - 모두 `VERIFIED` → `04_development_phases.md` 진행 현황 갱신, P2 로.
   - `FAILED_*` → 원인 분석과 수정. 수정 후 다시 판정하되 **이전 불통과 기록은 남김**.
5. **P2 (가상경계법, 이동 물체)**: 코드 작성 전에 1-E ~ 1-H 의 기준안을 쉬운 설명과 함께 사용자에게 제시하고 승인 후 커밋으로 고정.
   - 주의: 1-E 의 구 부가질량 (1/2)ρV 와 1-F 의 자유수면 위 원판 "절반" 값은 **원전 미확인** (`06_references.md` R-21). 무한유체 원판 (8/3)ρa³ 만 2차 출처로 확인됨. 기준에 쓰기 전에 출처를 확인하거나, 확인 불가하면 사용자에게 알리고 대체안을 제시.
   - P2 진입 전 결정 필요: D-05 (측벽 두께 해상 방식).

## 8. 주요 문서 위치

| 문서 | 내용 |
|---|---|
| `docs/design/02_design_document.md` | 방법 비교, 항목별 처리·한계, 실행/검증 구분, 입력 분리 |
| `docs/design/03_verification_plan.md` | 단계별 검증 계획. §1.1 = P1 고정 기준 |
| `docs/design/04_development_phases.md` | 단계표 P0–P9 와 진행 현황 |
| `docs/design/05_open_decisions.md` | 확정·미결정 항목 |
| `docs/design/06_references.md` | 문헌과 확인 수준 |
| `docs/P0_build_guide.md`, `docs/P1_verification_guide.md` | 노트북 빌드·검증 안내 |
| `validation/records/` | 판정 기록 |

## 9. 새 세션 첫 메시지 예시

> docs/HANDOFF.md 를 읽고 이어서 진행해. 먼저 노트북에서 빌드하고 P1 검증을 돌려줘.
