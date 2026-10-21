add_executable(hcw_unit_tests test_p0.cpp)
target_link_libraries(hcw_unit_tests PRIVATE hcw_core)
add_test(NAME unit_p0 COMMAND hcw_unit_tests)

add_executable(hcw_unit_tests_p1 test_p1.cpp)
target_link_libraries(hcw_unit_tests_p1 PRIVATE hcw_flow)
add_test(NAME unit_p1 COMMAND hcw_unit_tests_p1)

# P1 검증 시험의 축소판 (동작 확인용, 판정 아님). 전체 판정은 hcw_verify all 로 따로 실행.
foreach(t 1A 1B 1C 1D)
  add_test(NAME verify_quick_${t} COMMAND hcw_verify ${t} --quick)
  set_tests_properties(verify_quick_${t} PROPERTIES LABELS quick)
endforeach()

# 명령줄 시험: 매 실행마다 시험용 프로젝트 루트를 새로 복사한다.
set(HCW_TEST_ROOT "${CMAKE_CURRENT_BINARY_DIR}/cli_root")
add_test(NAME cli_setup
  COMMAND ${CMAKE_COMMAND}
    -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/fixtures/cli_root
    -DDST=${HCW_TEST_ROOT}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/copy_fixture.cmake)
set_tests_properties(cli_setup PROPERTIES FIXTURES_SETUP cli_root)

# hcw_cli_test(<이름> <기대 종료코드> <인자...> [ABSENT <경로>] [PRESENT_GLOB <패턴>])
function(hcw_cli_test name expect)
  cmake_parse_arguments(T "" "ABSENT;PRESENT_GLOB;DEPENDS" "" ${ARGN})
  list(JOIN T_UNPARSED_ARGUMENTS "|" joined)
  add_test(NAME cli_${name}
    COMMAND ${CMAKE_COMMAND}
      -DHCW=$<TARGET_FILE:hcw>
      -DARGS=${joined}
      -DEXPECT=${expect}
      -DWORKDIR=${HCW_TEST_ROOT}
      -DABSENT=${T_ABSENT}
      -DPRESENT_GLOB=${T_PRESENT_GLOB}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_expect.cmake)
  set_tests_properties(cli_${name} PROPERTIES FIXTURES_REQUIRED cli_root)
  if(T_DEPENDS)
    set_tests_properties(cli_${name} PROPERTIES DEPENDS cli_${T_DEPENDS})
  endif()
endfunction()

set(R ${HCW_TEST_ROOT})

hcw_cli_test(no_arguments 1)
hcw_cli_test(info 0 info)
hcw_cli_test(check_dev 0 check ${R}/cases/dev/DEV_example.case)
hcw_cli_test(check_dev_without_prefix 2 check ${R}/cases/dev/NoPrefix.case)
hcw_cli_test(check_dev_in_paper_folder 2 check ${R}/cases/paper/DEV_misplaced.case)
hcw_cli_test(check_missing_file 5 check ${R}/cases/dev/DEV_nope.case)
hcw_cli_test(prepare_dev 0 prepare ${R}/cases/dev/DEV_example.case
  PRESENT_GLOB "${R}/results/dev/DEV_example/*/run_meta.txt")

# P0 종료 조건: 개발·검증용 입력은 results/paper 에 쓸 수 없어야 한다.
hcw_cli_test(prepare_dev_into_paper_refused 3
  prepare ${R}/cases/dev/DEV_example.case --output-dir ${R}/results/paper/DEV_example
  ABSENT ${R}/results/paper)
hcw_cli_test(prepare_dev_into_paper_dotdot_refused 3
  prepare ${R}/cases/dev/DEV_example.case --output-dir ${R}/results/dev/../paper/x
  ABSENT ${R}/results/paper)
hcw_cli_test(prepare_verification_into_paper_refused 3
  prepare ${R}/cases/verification/V_example.case --output-dir ${R}/results/paper/V_example
  ABSENT ${R}/results/paper)

hcw_cli_test(lock_dev_refused 3 lock ${R}/cases/dev/DEV_example.case)
hcw_cli_test(lock_paper 0 lock ${R}/cases/paper/P_example.case)
hcw_cli_test(lock_paper_again_refused 3 lock ${R}/cases/paper/P_example.case
  DEPENDS lock_paper)
hcw_cli_test(check_paper_locked 0 check ${R}/cases/paper/P_example.case DEPENDS lock_paper)
hcw_cli_test(check_paper_tampered 3 check ${R}/cases/paper/P_tampered.case)
# 시험 기기는 등록부에 없으므로 논문 조건 실행은 거부되어야 한다.
hcw_cli_test(prepare_paper_unregistered_machine_refused 3
  prepare ${R}/cases/paper/P_example.case
  ABSENT ${R}/results/paper DEPENDS lock_paper)
