# hcw — 바닥이 막힌 중공 원통의 입수충격 해석기

논문용 정량 분석을 위한 독립 해석기입니다. 축이 수직인 중공 원통이 평평한 바닥판부터 입수할 때의 압력 분포, 충격력, 충격량, 구조응력 시간 이력을 계산하는 것이 목표입니다.

## 현재 상태

| 단계 | 내용 | 상태 |
|---|---|---|
| 설계 | 방법 선택, 검증 계획 | 완료 (`docs/design/`) |
| **P0** | 빌드 틀, 입력 검사, 논문 조건 잠금, 실행 기록 | 완료. Linux·Windows(노트북) 시험 18/18 통과 |
| P1 이후 | 유체·구조 계산 | 미착수 |

**아직 계산 기능은 없습니다.** `hcw prepare` 는 입력과 출력 규칙만 확인하고 `PREPARED` 상태를 기록합니다.

## 문서

- 설계 문서: [`docs/design/00_README.md`](docs/design/00_README.md)
- 노트북 빌드·시험 안내: [`docs/P0_build_guide.md`](docs/P0_build_guide.md)

## 폴더 구조

```
src/        해석기 소스 (C++17)
tests/      자동 시험
cases/      입력 파일: dev/ (개발용, DEV_), verification/ (검증용), paper/ (논문 조건, 잠금)
config/     기기 등록부 (machines.ini)
results/    계산 결과 (git 제외)
validation/ 검증 판정 기록
docs/       설계 문서와 안내
```
