#pragma once

// 모든 실행 결과 폴더에 남기는 메타데이터 (run_meta.txt).
// 실행 상태 정의: docs/design/02_design_document.md §5

#include <filesystem>
#include <string>
#include <vector>

#include "io/case_policy.h"

namespace hcw {

// PREPARED: 입력·정책 검사만 통과하고 계산은 하지 않음 (P0 의 유일한 상태).
// 이후 단계: RUN_COMPLETED, RUN_HEALTHY, VERIFIED, VALIDATED, FAILED_<사유>.
enum class RunStatus { Prepared };

std::string run_status_name(RunStatus s);

struct BuildInfo {
  std::string version;
  std::string git_commit;
  std::string git_tree_state;
  std::string compiler;
  std::string build_config;
  std::string platform;
};

BuildInfo current_build_info();

struct RunMetadata {
  RunStatus status = RunStatus::Prepared;
  std::string created_utc;
  std::string machine;
  MachineRole machine_role = MachineRole::Unregistered;
  unsigned hardware_threads = 0;
  BuildInfo build;
  LoadedCase input;
  LockState lock = LockState::Missing;
  std::vector<std::string> notes;
};

// run_meta.txt 를 쓰고 입력 파일 사본(input.case)을 함께 남긴다.
bool write_run_folder(const std::filesystem::path& dir, const RunMetadata& meta,
                      std::string& error);

}  // namespace hcw
