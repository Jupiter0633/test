# P0 빌드·시험 안내 (노트북, Windows)

P0 는 계산을 하지 않습니다. 해석기의 바탕 틀이 노트북에서 빌드되고 규칙대로 동작하는지 확인하는 단계입니다.

- 빌드: 소스 코드 → 실행 파일 `hcw.exe`
- 시험: 자동 시험 18개 (`ctest`)

준비물: Visual Studio Community 2026 ("C++를 사용한 데스크톱 개발"), GitHub Desktop 에서 최신 내용 받기 (**Fetch origin → Pull origin**).

---

## 방법 1. 명령 창 (권장 — 결과를 그대로 복사해 보내기 쉬움)

1. 시작 메뉴에서 **"x64 Native Tools Command Prompt for VS 2026"** 을 엽니다.
   - 반드시 **x64** 가 붙은 창이어야 합니다. "Developer PowerShell" 은 32비트(x86)로 열려서, 빌드 설정이 일부러 오류를 냅니다.
2. 저장소 폴더로 이동합니다. (GitHub Desktop 의 **Repository → Show in Explorer** 로 연 폴더 경로)
   ```bat
   cd /d "%USERPROFILE%\Documents\GitHub\test"
   ```
3. 빌드 설정 → 빌드 → 시험을 차례로 실행합니다.
   ```bat
   cmake --preset windows-x64-release
   cmake --build --preset windows-x64-release
   ctest --preset windows-x64-release
   ```
4. 마지막 줄이 `100% tests passed, 0 tests failed out of 18` 이면 성공입니다.
5. 프로그램 정보를 확인합니다.
   ```bat
   out\build\windows-x64-release\hcw.exe info
   ```

## 방법 2. Visual Studio 화면

1. Visual Studio 실행 → **로컬 폴더 열기** → 저장소 폴더(`test`) 선택.
2. 위쪽 구성 목록에서 **Windows x64 Release** 를 고릅니다.
3. 메뉴 **빌드 → 모두 빌드**.
4. 메뉴 **테스트 → CTest 실행** (또는 테스트 탐색기에서 모두 실행).

메뉴 이름은 Visual Studio 버전에 따라 조금 다를 수 있습니다. 다르면 화면을 보내 주세요.

---

## 노트북 등록 (한 번만)

1. `hcw.exe info` 출력의 `machine name` 값을 확인합니다 (예: `DESKTOP-AB12CD`).
2. `config\machines.ini` 의 `[machines]` 아래에 **소문자로** 한 줄 추가합니다.
   ```ini
   desktop-ab12cd = development
   ```
3. 노트북은 반드시 `development` 입니다. `workstation` 으로 등록하면 논문 조건 실행 차단 장치가 무력화됩니다.

등록하지 않아도 개발용 실행에는 지장이 없습니다 (메타데이터에 `role = unregistered` 로 기록됨).

---

## 동작 확인 예시

```bat
out\build\windows-x64-release\hcw.exe check cases\dev\DEV_example.case
out\build\windows-x64-release\hcw.exe prepare cases\dev\DEV_example.case
out\build\windows-x64-release\hcw.exe prepare cases\dev\DEV_example.case --output-dir results\paper\x
```

| 명령 | 기대 결과 |
|---|---|
| `check` | `ok: DEV_example (purpose = development, ...)` |
| `prepare` | `results\dev\DEV_example\<시각>\` 폴더 생성, `status: PREPARED (no calculation was performed)` |
| `prepare ... --output-dir results\paper\x` | `refused: nothing was written` (종료 코드 3) — 개발용 입력은 논문 결과 폴더에 쓸 수 없음 |

## 문제가 생기면

| 증상 | 원인·조치 |
|---|---|
| `hcw must be built as a 64-bit (x64) program` | x86 창에서 실행함 → "x64 Native Tools" 창 사용. 이미 만든 `out` 폴더를 지우고 다시 |
| `'ninja' ... not found` 또는 CMake 생성기 오류 | Visual Studio 설치에서 "Windows용 C++ CMake 도구" 가 빠졌을 수 있음 |
| 시험 일부 실패 | 출력 전체를 복사해 보내 주세요. 실패를 무시하고 넘어가지 않습니다 |
| 경로에 한글이 있어 오류 | 오류 메시지와 폴더 경로를 보내 주세요 |

출력 메시지는 Windows 콘솔 한글 깨짐을 피하려고 영어로 나옵니다.
