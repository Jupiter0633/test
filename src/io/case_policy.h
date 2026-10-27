#pragma once

// 입력 목적(개발/검증/논문) 구분, 논문 조건 잠금, 출력 위치 허용 규칙.
// 근거: docs/design/02_design_document.md §5–6, 01_environment_check.md §4.3

#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include "io/case_file.h"

namespace hcw {

namespace fs = std::filesystem;

enum class Purpose { Development, Verification, Paper };

std::optional<Purpose> parse_purpose(const std::string& text);
std::string purpose_name(Purpose p);       // development / verification / paper
std::string purpose_directory(Purpose p);  // dev / verification / paper

struct Diagnostics {
  std::vector<std::string> errors;
  std::vector<std::string> warnings;
  bool ok() const { return errors.empty(); }
};

// 저장소 최상위 표식 파일. 이 파일이 있는 폴더를 프로젝트 루트로 본다.
inline constexpr const char* kRootMarker = ".hcw-root";

std::optional<fs::path> find_project_root(const fs::path& start);

// 읽어 들인 입력 파일과 그 해시.
struct LoadedCase {
  fs::path path;
  std::string sha256;  // 정규화된 텍스트 기준
  CaseFile data;
  Purpose purpose = Purpose::Development;
  std::string name;
};

// 입력 파일을 읽고 형식·위치·이름 규칙을 검사한다.
std::optional<LoadedCase> load_case(const fs::path& case_path, const fs::path& root,
                                    Diagnostics& diag);

// 논문 조건 잠금 파일: <입력파일>.lock
enum class LockState { Missing, Match, Mismatch, Invalid };

fs::path lock_path_for(const fs::path& case_path);
LockState check_lock(const LoadedCase& c, std::string& detail);
bool write_lock(const LoadedCase& c, std::string& error);

// 기기 등록부: <루트>/config/machines.ini 의 [machines] 섹션
//   <컴퓨터 이름> = development | workstation
enum class MachineRole { Unregistered, Development, Workstation };

std::string current_machine_name();
MachineRole lookup_machine(const fs::path& root, const std::string& machine,
                           Diagnostics& diag);
std::string machine_role_name(MachineRole r);

// child 가 parent 안(또는 같은 폴더)에 있는지. 존재하지 않는 경로도 처리.
bool is_within(const fs::path& child, const fs::path& parent);

fs::path default_output_dir(const fs::path& root, const LoadedCase& c,
                            const std::string& timestamp);

// 출력 허용 여부를 판단하는 순수 함수 (단위 시험 대상).
struct OutputRequest {
  Purpose purpose = Purpose::Development;
  fs::path root;
  fs::path output_dir;
  LockState lock = LockState::Missing;
  MachineRole machine = MachineRole::Unregistered;
  std::string git_commit;      // "unknown" 이면 출처 불명
  std::string git_tree_state;  // clean / modified / unknown
};

Diagnostics check_output_policy(const OutputRequest& req);

}  // namespace hcw
