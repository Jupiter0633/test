# hcw 를 실행하고 종료 코드(및 선택적으로 파일 존재/부재)를 확인한다.
# 입력: HCW, ARGS('|' 구분), EXPECT, WORKDIR, ABSENT(선택), PRESENT_GLOB(선택)

if(ARGS STREQUAL "")
  set(arg_list "")
else()
  string(REPLACE "|" ";" arg_list "${ARGS}")
endif()

execute_process(
  COMMAND "${HCW}" ${arg_list}
  WORKING_DIRECTORY "${WORKDIR}"
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)

message("exit code: ${rc}\nstdout:\n${out}\nstderr:\n${err}")

if(NOT rc STREQUAL "${EXPECT}")
  message(FATAL_ERROR "expected exit code ${EXPECT}, got ${rc}")
endif()

if(NOT ABSENT STREQUAL "" AND EXISTS "${ABSENT}")
  message(FATAL_ERROR "path must not exist after this command: ${ABSENT}")
endif()

if(NOT PRESENT_GLOB STREQUAL "")
  file(GLOB found "${PRESENT_GLOB}")
  if(found STREQUAL "")
    message(FATAL_ERROR "expected a file matching ${PRESENT_GLOB}")
  endif()
endif()
