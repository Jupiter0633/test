# 판정 기록: <단계>-<시험 ID> (<YYYY-MM-DD>)

> `docs/design/03_verification_plan.md` §0.4 형식. 통과·불통과 모두 같은 형식으로 남깁니다.

| 항목 | 값 |
|---|---|
| 시험 ID | 예: 1-A |
| 판정 기준 커밋 | 기준이 실행 전에 커밋된 해시 |
| 코드 커밋 / 트리 상태 | `run_meta.txt` 의 git_commit / git_tree_state |
| 입력 파일 / SHA-256 | |
| 기기 / 역할 | |
| 컴파일러 / 구성 | |

## 기준 (실행 전 고정)

- 지표:
- 통과 기준과 근거 유형 (T/U/P/C):

## 결과

- 측정값:
- 판정: `VERIFIED` / `VALIDATED` / `FAILED_<사유>` / 실행 성공만(`RUN_COMPLETED`)

## 비고

- 기준을 변경했다면: 변경 사유, 날짜, 변경 전 기준에 의한 판정 (삭제하지 않음)
