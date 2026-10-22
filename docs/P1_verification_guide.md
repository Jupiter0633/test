# P1 검증 실행 안내 (노트북, Windows)

P1 은 물 한 가지(또는 주어진 밀도 분포)의 유동을 푸는 단계입니다. 자유수면과 원통은 아직 없습니다.
검증 시험 4종(1-A ~ 1-D)의 조건과 통과 기준은 `docs/design/03_verification_plan.md` §1.1 에 **코드 작성 전에 고정**되어 있습니다 (커밋 35a77ac).

## 1. 최신 코드 받기와 빌드

GitHub Desktop 에서 **Fetch origin → Pull origin** 후, "x64 Native Tools Command Prompt for VS 2026" 에서:

```bat
cd /d "%USERPROFILE%\Documents\GitHub\test"
cmake --preset windows-x64-release
cmake --build --preset windows-x64-release
ctest --preset windows-x64-release
```

`ctest` 는 빠른 시험만 돌립니다 (단위 시험 + 검증 시험 축소판). 축소판은 **동작 확인용이며 판정이 아닙니다.**

## 2. 전체 검증 실행 (판정)

```bat
out\build\windows-x64-release\hcw_verify.exe all
```

- 시험별로 실행할 수도 있습니다: `1A`, `1A-a`, `1A-b`, `1B`, `1C`, `1D`
- 가장 고운 격자(256×256) 때문에 **수십 분 이상** 걸릴 수 있습니다. 창을 닫지 마세요.
- 결과는 화면과 `results\verification\P1\<시험>\<시각>\report.txt` 에 저장됩니다.

## 3. 결과 읽는 법

| 표시 | 뜻 |
|---|---|
| `[PASS]` / `[FAIL]` | 고정된 기준 하나하나의 통과 여부 |
| `p_obs` | 관측 수렴 차수. 격자를 절반으로 줄일 때 오차가 줄어드는 정도 (2 에 가까우면 설계대로) |
| `verdict : VERIFIED` | 모든 기준 통과 |
| `verdict : FAILED_CRITERIA` | 기준 하나 이상 불통과 — 결과를 그대로 기록하고 원인을 분석합니다 |
| `verdict : FAILED_SOLVER` | 압력 방정식이 허용치까지 수렴하지 못함 |
| `QUICK_CHECK_ONLY` | 축소판 실행. 판정 아님 |

`report.txt` 의 첫 부분에 코드 커밋, 컴파일러, 기기 이름이 기록됩니다. 판정 기록(`validation/records/`)은 이 파일을 근거로 작성합니다.

## 4. 정식 판정 전에는 빌드 폴더를 지우고 새로 빌드 (2026-10-22 추가)

```bat
rmdir /s /q out\build\windows-x64-release
cmake --preset windows-x64-release
cmake --build --preset windows-x64-release
```

- 이유: 한국어 Visual Studio 의 컴파일러는 헤더 의존 정보(`/showIncludes`)를 한국어로 출력합니다.
  - 빌드 폴더를 처음 만들 때 CMake 가 기록한 문구의 인코딩이 실제 출력과 어긋나면, Ninja 가 헤더 의존성을 기록하지 못합니다.
  - 그러면 헤더만 바뀐 경우 그 헤더를 쓰는 파일이 다시 컴파일되지 않습니다. 이전 코드와 새 코드가 섞인 프로그램이 만들어질 수 있습니다.
- 2026-10-22 노트북에서 실제로 확인했습니다. 빌드 폴더를 새로 만든 뒤에는 의존성이 정상 기록되는 것도 확인했습니다.
- 정식 판정은 반드시 새로 만든 빌드 폴더에서, 수정 없는 커밋으로 실행합니다.
