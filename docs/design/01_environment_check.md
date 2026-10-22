# 01. 실행 환경 확인 (읽기 전용)

확인일: 2026-10-21
방법: 읽기 전용 명령만 사용 (`/etc/os-release`, `uname`, `lscpu`, `free`, `nvidia-smi`, `lspci`, 각 도구의 `--version`, `ldconfig -p`). 설치나 설정 변경은 하지 않았습니다.

## 0. 중요 — 이 결과는 "워크스테이션"이 아닙니다

이번 작업은 **클라우드의 일회성 컨테이너(KVM 가상머신)** 에서 실행되었습니다.
아래 수치는 이 컨테이너의 사양이며, **논문 계산에 쓸 워크스테이션의 사양이 아닙니다.**
워크스테이션 사양은 §3의 명령을 사용자가 직접 실행해 기록해야 하며, 그 전까지 워크스테이션 사양은 **미확인** 상태입니다.

## 1. 이번 세션 컨테이너에서 확인한 사양

| 항목 | 확인 결과 |
|---|---|
| OS | Ubuntu 24.04.5 LTS, Linux 6.18.44 (x86_64), KVM 하이퍼바이저 |
| CPU | Intel Xeon Processor @ 2.80 GHz, 4 vCPU (소켓 1 × 코어 4, 코어당 스레드 1), AVX2/AVX-512 지원 |
| RAM | 15 GiB, 스왑 없음 |
| GPU | 없음 (`nvidia-smi` 없음, `lspci`에 VGA/3D 장치 없음) |
| 디스크 | 세션 할당량 기준 약 31 GB 여유 |
| 컴파일러 | gcc/g++ 13.3.0, clang 있음, **gfortran 없음** |
| 빌드 도구 | CMake 3.28.3, make, ninja, git 2.43.0 |
| 병렬 | libgomp(OpenMP) 있음, **MPI 없음**, **CUDA(nvcc) 없음** |
| 수치 라이브러리 | libblas/liblapack/FFTW3 런타임 라이브러리만 확인. Eigen 헤더, PETSc, hypre, HDF5 헤더 없음 |
| Python | 3.13.16, NumPy 2.5.3, **SciPy 없음, matplotlib 없음** |

해석: 이 컨테이너는 문서 작성과 소규모 2차원 단위 시험 정도에만 적합합니다. 논문 계산(특히 3차원)은 워크스테이션에서 수행해야 합니다.

## 2. 설계에 영향을 주는 워크스테이션 정보 (미확인 → 사용자 확인 필요)

| 필요한 정보 | 설계에 미치는 영향 |
|---|---|
| OS (Linux / Windows / WSL2) | 빌드 체계(CMake 공통), MPI 구현 선택 |
| 물리 코어 수, 소켓 수 | OpenMP만으로 충분한지, MPI 영역분할이 필요한지 |
| RAM 용량 | 3차원 격자 상한. 대략 셀당 1–2 kB를 가정하면 64 GB에서 3천만–6천만 셀 수준 (가정치, 구현 후 실측 필요) |
| GPU 모델·메모리·CUDA 버전 | GPU 이식을 1차 범위에 넣을지 여부 (현재 설계는 GPU를 요구하지 않음) |
| 설치된 컴파일러·MPI·PETSc/hypre | 압력 Poisson 해법을 외부 라이브러리로 할지 자체 다중격자로 할지 |
| 상용 FE 소프트웨어(ABAQUS/ANSYS 등) 사용 가능 여부 | 구조 해석기의 교차검증 수단 |

## 3. 워크스테이션에서 실행할 읽기 전용 확인 명령

Linux:
```bash
cat /etc/os-release; uname -a
lscpu; free -h; df -h
nvidia-smi; lspci | grep -i -E 'vga|3d'
gcc --version; g++ --version; gfortran --version; cmake --version
mpicc --version; mpirun --version; nvcc --version
python3 -c "import numpy, scipy, matplotlib; print(numpy.__version__, scipy.__version__, matplotlib.__version__)"
```

Windows (PowerShell):
```powershell
Get-ComputerInfo -Property OsName,OsVersion,CsProcessors,CsTotalPhysicalMemory
Get-CimInstance Win32_Processor | Select Name,NumberOfCores,NumberOfLogicalProcessors
Get-CimInstance Win32_VideoController | Select Name,AdapterRAM,DriverVersion
nvidia-smi
```

결과는 이 파일 §4에 붙여 넣어 기록합니다.

## 4. 사용자 기기 확인 결과

### 4.1 노트북 (사용자가 PowerShell 로 직접 실행, 2026-10-21)

**역할(사용자 확정, 2026-10-21): 프로그램 개발 전용. 논문 계산용 워크스테이션과는 별개 기기입니다.**

| 항목 | 확인 결과 |
|---|---|
| OS | Microsoft Windows 11 Pro, 버전 10.0.26200 |
| CPU | 11th Gen Intel Core i7-1165G7 @ 2.80 GHz, 물리 코어 4 / 논리 프로세서 8 (저전력 노트북용 CPU) |
| RAM | 16,873,873,408 바이트 (약 15.7 GiB) |
| GPU | Intel Iris Xe Graphics (내장 그래픽, 시스템 메모리 공유). AdapterRAM 1 GiB 로 표시되나 이 값은 내장 그래픽에서 실제 사용 가능량과 다를 수 있음 |
| NVIDIA/CUDA | 없음 (`nvidia-smi` 명령 없음) |
| WSL | 설치되어 있지 않음 (`wsl --status` 실행 시 설치 안내가 떴으나 취소되어 아무것도 설치되지 않음) |
| C/C++ 컴파일러 | 없음 (`cl`, `gcc`, `g++` 모두 찾을 수 없음) |
| CMake | 없음 |
| MPI | 없음 (`mpiexec` 없음) |
| Python | 3.12.10 (`C:\Users\User\AppData\Local\Programs\Python\Python312\python.exe`). `WindowsApps\python.exe` 는 Microsoft Store 실행 별칭 |
| NumPy/SciPy/matplotlib | 미확인 |

설계에 미치는 영향:
- 2차원 축대칭 개발과 단계 1–2 검증 시험에는 충분할 것으로 예상합니다 (구현 후 실측 필요).
- **3차원 규칙파 해석(P9)의 논문 계산에는 부적합합니다.** 코어 4개, 메모리 약 16 GB, CUDA GPU 없음. 3차원 단계에는 별도 워크스테이션이 필요합니다.
- CUDA GPU 가 없으므로 GPU 이식은 이 기기에서 시험할 수 없습니다 (D-16).
- Windows 이므로 빌드 체계는 CMake 공통으로 두고, 개발 도구는 WSL2(Ubuntu) 또는 Visual Studio(MSVC) 중 선택이 필요합니다 (D-07 에 추가).

- **현재 C++ 코드를 빌드할 수 없습니다.** P0(빌드 체계) 착수 전에 C++ 개발 도구 설치가 필요합니다 (D-24). 설치는 사용자 결정 사항이며 이 문서 작성 과정에서는 설치하지 않았습니다.

주의: `wsl --status` 는 WSL 이 없으면 설치 안내 프롬프트를 띄웁니다. 읽기 전용 확인에는 `where.exe wsl` 이나 `Get-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux` (관리자 권한 필요) 를 쓰는 편이 안전합니다.

### 4.1.1 노트북 개발 도구 설치 결과 (사용자 확인, 2026-10-21)

| 항목 | 결과 |
|---|---|
| Visual Studio | Community 2026, v18.11.1 ("C++를 사용한 데스크톱 개발" 워크로드: MSVC 빌드 도구, C++ CMake 도구, Windows 11 SDK 10.0.26100, JIT 디버거, AddressSanitizer, vcpkg) |
| C++ 컴파일러 | MSVC `cl` 19.51.36260 — 기본 Developer PowerShell 은 **x86** 대상으로 열림 |
| CMake | 4.3.1-msvc1 (Visual Studio 포함본) |
| GitHub Desktop | 설치·저장소 연결 완료 (설계 브랜치 추적) |
| Git for Windows | 설치 완료 (명령줄 PATH 등록). 빌드 커밋 기록 정상 확인 |
| 기기 이름 | DESKTOP-4JJOTVL, 하드웨어 스레드 8 — `config/machines.ini` 에 development 로 등록 |

주의: 해석기는 **64비트(x64)로 빌드**해야 합니다 (32비트는 메모리 약 4 GB 한계). CMake 의 Visual Studio 생성기를 `-A x64` 로 쓰거나, 명령줄 빌드는 "x64 Native Tools" 환경에서 수행합니다. 빌드 설정(P0)에서 x64 를 강제하고 32비트 빌드는 오류로 처리합니다.

### 4.2 논문 계산용 워크스테이션

- 노트북과 별개 기기 (사용자 확정).
- OS: **Windows** (사용자 확정, 2026-10-21).
- 영향: 개발(노트북)과 계산(워크스테이션)이 모두 Windows 이므로 MSVC + CMake 체계로 통일 가능 (D-24). MPI 가 필요해지면 Windows 에서는 MS-MPI 가 후보.
- 확인일: **2026-10-22** (사용자가 §3 Windows 명령을 워크스테이션에서 실행, 출력 원문을 대화로 전달)
- 결과:

| 항목 | 값 |
|---|---|
| 컴퓨터 이름 | `DESKTOP-G0U5CPK` (`hcw info` 의 machine name 으로 재확인 예정) |
| OS | Microsoft Windows 11 Pro for Workstations, 10.0.26200 |
| CPU | Intel Xeon w5-2455X, 물리 코어 12, 논리 프로세서 24 |
| 메모리 | 67 907 026 944 바이트 (약 63 GiB) |
| GPU | NVIDIA RTX 4000 Ada Generation, 메모리 20 475 MiB (`nvidia-smi`; `Win32_VideoController` 의 AdapterRAM 4 GiB 는 32비트 필드 한계로 잘린 값), 드라이버 596.71, CUDA 13.2 (드라이버 지원 버전) |

- 해석 (판정 아님, 계획용):
  - **코어 12개**: 지금 해석기는 단일 스레드라 노트북과 속도 차이가 크지 않습니다. 대신 정식 검증 여러 개를 동시에 돌릴 수 있고, OpenMP 병렬화(D-07)를 넣으면 2차원·3차원 계산에서 큰 이득이 예상됩니다 (구현 후 실측).
  - **메모리 약 63 GiB**: §2 의 가정치(셀당 1–2 kB)로 3차원 약 3천만–6천만 셀 규모. 구현 후 실측으로 확정합니다.
  - **GPU 20 GiB (CUDA 지원)**: GPU 이식은 현재 설계 범위 밖입니다 (D-16 에서 P9 전에 결정).

### 4.3 개발 기기와 계산 기기가 다르기 때문에 생기는 요구사항

1. **이식성:** 워크스테이션 OS 가 미정이므로 코드는 표준 C++ 와 CMake 만 사용하고, OS 전용 API·경로 표기·컴파일러 확장을 쓰지 않습니다. 외부 라이브러리는 워크스테이션에 설치 가능한 것으로 제한합니다.
2. **노트북에서 도는 크기:** 모든 검증 시험(단계 1–6)은 노트북에서 수 분–수십 분 안에 끝나는 축소 격자 버전을 함께 둡니다. 큰 격자 수렴 시험은 워크스테이션에서 수행합니다.
3. **재검증:** 노트북에서 받은 `VERIFIED` 판정은 개발 확인용입니다. 논문 계산(P8) 전에 **워크스테이션 빌드로 단계 1–6 검증 시험을 다시 실행**해 같은 판정이 나오는지 확인하고, 두 기기의 결과 차이를 기록합니다.
4. **실행 기기 기록:** 모든 결과 메타데이터에 기기 이름, OS, 컴파일러와 버전, 최적화 옵션, 스레드 수를 남깁니다. `results/paper/` 출력은 워크스테이션으로 등록된 기기에서만 허용합니다.
