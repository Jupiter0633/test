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

add_executable(hcw src/app/main.cpp)
target_link_libraries(hcw PRIVATE hcw_core)

# P1: 단상·변수밀도 유동 해석기
add_library(hcw_flow STATIC
  src/flow/grid.cpp
  src/flow/material.cpp
  src/flow/poisson.cpp
  src/flow/solver.cpp)
target_include_directories(hcw_flow PUBLIC ${CMAKE_SOURCE_DIR}/src)

# P1 검증 시험 (docs/design/03_verification_plan.md §1.1)
add_library(hcw_verify_p1 STATIC src/verify/p1_cases.cpp)
target_link_libraries(hcw_verify_p1 PUBLIC hcw_flow)

add_executable(hcw_verify src/app/verify_main.cpp)
target_link_libraries(hcw_verify PRIVATE hcw_verify_p1 hcw_core)

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
