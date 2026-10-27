#include "io/run_metadata.h"

#include <fstream>

#include "hcw_build_info.h"

namespace hcw {
namespace {

std::string lock_state_name(LockState s) {
  switch (s) {
    case LockState::Missing: return "missing";
    case LockState::Match: return "match";
    case LockState::Mismatch: return "mismatch";
    case LockState::Invalid: return "invalid";
  }
  return "unknown";
}

}  // namespace

std::string run_status_name(RunStatus s) {
  switch (s) {
    case RunStatus::Prepared: return "PREPARED";
  }
  return "UNKNOWN";
}

BuildInfo current_build_info() {
  BuildInfo b;
  b.version = HCW_VERSION;
  b.git_commit = HCW_GIT_COMMIT;
  b.git_tree_state = HCW_GIT_TREE_STATE;
  b.compiler = std::string(HCW_COMPILER_ID) + " " + HCW_COMPILER_VERSION;
  b.build_config = HCW_BUILD_CONFIG;
#if defined(_WIN32)
  b.platform = "windows";
#elif defined(__linux__)
  b.platform = "linux";
#elif defined(__APPLE__)
  b.platform = "macos";
#else
  b.platform = "other";
#endif
  b.platform += sizeof(void*) == 8 ? "-x64" : "-x86";
  return b;
}

bool write_run_folder(const std::filesystem::path& dir, const RunMetadata& meta,
                      std::string& error) {
  std::error_code ec;
  if (std::filesystem::exists(dir, ec)) {
    error = "output folder already exists: " + dir.string();
    return false;
  }
  std::filesystem::create_directories(dir, ec);
  if (ec) {
    error = "cannot create output folder " + dir.string() + ": " + ec.message();
    return false;
  }

  std::ofstream out(dir / "run_meta.txt", std::ios::binary);
  if (!out) {
    error = "cannot write run_meta.txt";
    return false;
  }
  out << "# hcw run metadata\n";
  if (meta.input.purpose == Purpose::Development) {
    out << "# WATERMARK: DEVELOPMENT INPUT - NOT PAPER CONDITIONS\n";
  }
  out << "[run]\n"
      << "status = " << run_status_name(meta.status) << "\n"
      << "created_utc = " << meta.created_utc << "\n";
  for (std::size_t i = 0; i < meta.notes.size(); ++i) {
    out << "note_" << i + 1 << " = " << meta.notes[i] << "\n";
  }
  out << "\n[input]\n"
      << "name = " << meta.input.name << "\n"
      << "purpose = " << purpose_name(meta.input.purpose) << "\n"
      << "file = " << meta.input.path.generic_string() << "\n"
      << "sha256 = " << meta.input.sha256 << "\n"
      << "lock = " << lock_state_name(meta.lock) << "\n"
      << "\n[build]\n"
      << "version = " << meta.build.version << "\n"
      << "git_commit = " << meta.build.git_commit << "\n"
      << "git_tree_state = " << meta.build.git_tree_state << "\n"
      << "compiler = " << meta.build.compiler << "\n"
      << "config = " << meta.build.build_config << "\n"
      << "platform = " << meta.build.platform << "\n"
      << "\n[machine]\n"
      << "name = " << meta.machine << "\n"
      << "role = " << machine_role_name(meta.machine_role) << "\n"
      << "hardware_threads = " << meta.hardware_threads << "\n"
      << "solver_threads = " << meta.solver_threads << "\n";
  if (!out) {
    error = "failed while writing run_meta.txt";
    return false;
  }
  out.close();

  std::filesystem::copy_file(meta.input.path, dir / "input.case", ec);
  if (ec) {
    error = "cannot copy input file: " + ec.message();
    return false;
  }
  return true;
}

}  // namespace hcw
