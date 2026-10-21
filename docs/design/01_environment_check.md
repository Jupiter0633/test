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

## 4. 워크스테이션 확인 결과 (미기입)

- 확인일: 미기입
- 결과: 미기입
