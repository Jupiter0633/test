cmake_minimum_required(VERSION 3.21)

project(hcw
  VERSION 0.0.1
  DESCRIPTION "Water-entry impact solver for a bottom-closed hollow cylinder"
  LANGUAGES CXX)

# 해석기는 4 GB 이상의 메모리를 다뤄야 하므로 64비트 빌드만 허용한다.
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  message(FATAL_ERROR
    "hcw must be built as a 64-bit (x64) program. "
    "Use the 'windows-x64-*' preset or the x64 Native Tools environment.")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
  add_compile_options(/W4 /permissive- /utf-8 /EHsc)
  add_compile_definitions(_CRT_SECURE_NO_WARNINGS NOMINMAX)
else()
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# 빌드할 때마다 git 커밋과 변경 여부를 기록한다 (결과 메타데이터용).
set(HCW_BUILD_INFO_HEADER "${CMAKE_BINARY_DIR}/generated/hcw_build_info.h")
add_custom_target(hcw_build_info
  COMMAND ${CMAKE_COMMAND}
    -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
    -DOUTPUT=${HCW_BUILD_INFO_HEADER}
    -P ${CMAKE_SOURCE_DIR}/cmake/BuildInfo.cmake
  BYPRODUCTS ${HCW_BUILD_INFO_HEADER}
  COMMENT "Recording git commit for run metadata")

# 다중 스레드 실행 (선택, 기본 1). 결과는 스레드 수와 상관없이 같다 (src/common/parallel.h).
add_library(hcw_parallel STATIC src/common/parallel.cpp)
target_include_directories(hcw_parallel PUBLIC ${CMAKE_SOURCE_DIR}/src)

add_library(hcw_core STATIC
  src/common/sha256.cpp
  src/common/text.cpp
  src/io/case_file.cpp
  src/io/case_policy.cpp
  src/io/run_metadata.cpp)
target_include_directories(hcw_core
  PUBLIC ${CMAKE_SOURCE_DIR}/src
  PRIVATE ${CMAKE_BINARY_DIR}/generated)
target_compile_definitions(hcw_core PRIVATE
  HCW_COMPILER_ID="${CMAKE_CXX_COMPILER_ID}"
  HCW_COMPILER_VERSION="${CMAKE_CXX_COMPILER_VERSION}"
  HCW_BUILD_CONFIG="$<CONFIG>"
  HCW_VERSION="${PROJECT_VERSION}")
add_dependencies(hcw_core hcw_build_info)
target_link_libraries(hcw_core PUBLIC hcw_parallel)
# 순서 의존만으로는 Ninja 가 헤더 내용이 바뀐 빌드에서 run_metadata 를 다시 컴파일하지 않을 수 있다
# (2026-10-22 확인: 작업 트리가 바뀐 뒤 첫 빌드가 "clean" 으로 기록됨). 명시적 의존을 건다.
set_source_files_properties(src/io/run_metadata.cpp PROPERTIES OBJECT_DEPENDS ${HCW_BUILD_INFO_HEADER})

add_executable(hcw src/app/main.cpp)
target_link_libraries(hcw PRIVATE hcw_core)

# P1: 단상·변수밀도 유동 해석기, P2: 물체 (형상, 컷셀, 표면 힘), P3: 기하 VOF
add_library(hcw_flow STATIC
  src/flow/grid.cpp
  src/flow/material.cpp
  src/flow/poisson.cpp
  src/flow/linear.cpp
  src/flow/solver.cpp
  src/body/shape.cpp
  src/body/cut_cell.cpp
  src/body/force.cpp
  src/interface/plic.cpp
  src/interface/vof.cpp
  src/interface/vof_material.cpp
  src/interface/two_phase.cpp)
target_include_directories(hcw_flow PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(hcw_flow PUBLIC hcw_parallel)

# 검증 시험 (docs/design/03_verification_plan.md §1.1 P1, §1.2 P2)
add_library(hcw_verify_cases STATIC
  src/verify/report_util.cpp
  src/verify/p1_cases.cpp
  src/verify/potential.cpp
  src/verify/p2_cases.cpp
  src/verify/p3_cases.cpp
  src/verify/p4_cases.cpp)
target_link_libraries(hcw_verify_cases PUBLIC hcw_flow)

add_executable(hcw_verify src/app/verify_main.cpp)
target_link_libraries(hcw_verify PRIVATE hcw_verify_cases hcw_core)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
