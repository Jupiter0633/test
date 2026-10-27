// hcw: 바닥이 막힌 중공 원통 입수충격 해석기 — P0 명령줄 진입점.
// P0 에서는 계산을 하지 않는다. 입력 검사, 논문 조건 잠금, 출력 허용 판단,
// 실행 메타데이터 기록만 수행한다.

#include <filesystem>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "common/parallel.h"
#include "common/text.h"
#include "io/case_policy.h"
#include "io/run_metadata.h"

namespace fs = std::filesystem;
using namespace hcw;

namespace {

enum ExitCode {
  kOk = 0,
  kUsage = 1,
  kInvalidInput = 2,
  kPolicyRefused = 3,
  kIoError = 5,
};

void print_usage() {
  std::cout <<
      "hcw - water-entry impact solver (P0: input checks and run bookkeeping only)\n"
      "\n"
      "usage:\n"
      "  hcw info                               show build and machine information\n"
      "  hcw check   <case-file>                validate a case file\n"
      "  hcw lock    <case-file>                lock a purpose = paper case (SHA-256)\n"
      "  hcw prepare <case-file> [--output-dir <dir>]\n"
      "                                         check output policy and create a run folder\n"
      "\n"
      "exit codes: 0 ok, 1 usage, 2 invalid input, 3 refused by policy, 5 file error\n";
}

void print_diagnostics(const Diagnostics& d) {
  for (const auto& w : d.warnings) std::cout << "warning: " << w << "\n";
  for (const auto& e : d.errors) std::cerr << "error: " << e << "\n";
}

struct Context {
  fs::path root;
  LoadedCase input;
};

// 입력 파일을 찾고 읽는다. 실패하면 종료 코드를 돌려준다.
int load(const std::string& case_arg, Context& ctx) {
  const fs::path case_path(case_arg);
  std::error_code ec;
  if (!fs::is_regular_file(case_path, ec)) {
    std::cerr << "error: case file not found: " << case_arg << "\n";
    return kIoError;
  }
  const auto root = find_project_root(case_path);
  if (!root) {
    std::cerr << "error: no project root ('" << kRootMarker
              << "' file) found above " << case_arg << "\n";
    return kInvalidInput;
  }
  Diagnostics diag;
  auto loaded = load_case(case_path, *root, diag);
  print_diagnostics(diag);
  if (!loaded) return kInvalidInput;
  ctx.root = *root;
  ctx.input = std::move(*loaded);
  return kOk;
}

int cmd_info() {
  const BuildInfo b = current_build_info();
  std::cout << "version          " << b.version << "\n"
            << "git commit       " << b.git_commit << " (" << b.git_tree_state << ")\n"
            << "compiler         " << b.compiler << " [" << b.build_config << "]\n"
            << "platform         " << b.platform << "\n"
            << "machine name     " << current_machine_name() << "\n"
            << "hardware threads " << std::thread::hardware_concurrency() << "\n";
  return kOk;
}

int cmd_check(const std::string& case_arg) {
  Context ctx;
  if (int rc = load(case_arg, ctx); rc != kOk) return rc;
  std::cout << "ok: " << ctx.input.name << " (purpose = " << purpose_name(ctx.input.purpose)
            << ", sha256 = " << ctx.input.sha256 << ")\n";
  if (ctx.input.purpose == Purpose::Paper) {
    std::string detail;
    const LockState lock = check_lock(ctx.input, detail);
    std::cout << "lock: " << detail << "\n";
    if (lock == LockState::Mismatch || lock == LockState::Invalid) return kPolicyRefused;
  }
  return kOk;
}

int cmd_lock(const std::string& case_arg) {
  Context ctx;
  if (int rc = load(case_arg, ctx); rc != kOk) return rc;
  std::string error;
  if (!write_lock(ctx.input, error)) {
    std::cerr << "error: " << error << "\n";
    return kPolicyRefused;
  }
  std::cout << "locked: " << lock_path_for(ctx.input.path).string() << "\n";
  return kOk;
}

int cmd_prepare(const std::string& case_arg, const std::string& output_arg) {
  Context ctx;
  if (int rc = load(case_arg, ctx); rc != kOk) return rc;

  RunMetadata meta;
  meta.created_utc = utc_now_iso();
  meta.machine = current_machine_name();
  meta.hardware_threads = std::thread::hardware_concurrency();
  {
    std::string thread_error;
    const int n = parallel::threads_from_environment(thread_error);
    if (n < 1) {
      std::cerr << "error: " << thread_error << "\n";
      return kInvalidInput;
    }
    meta.solver_threads = n;
  }
  meta.build = current_build_info();
  meta.input = ctx.input;

  Diagnostics registry_diag;
  meta.machine_role = lookup_machine(ctx.root, meta.machine, registry_diag);
  print_diagnostics(registry_diag);
  if (!registry_diag.ok()) return kInvalidInput;

  std::string lock_detail;
  meta.lock = check_lock(ctx.input, lock_detail);

  fs::path out_dir = output_arg.empty()
                         ? default_output_dir(ctx.root, ctx.input, utc_now_compact())
                         : fs::path(output_arg);

  OutputRequest req;
  req.purpose = ctx.input.purpose;
  req.root = ctx.root;
  req.output_dir = out_dir;
  req.lock = meta.lock;
  req.machine = meta.machine_role;
  req.git_commit = meta.build.git_commit;
  req.git_tree_state = meta.build.git_tree_state;
  const Diagnostics policy = check_output_policy(req);
  if (!policy.ok()) {
    print_diagnostics(policy);
    if (ctx.input.purpose == Purpose::Paper) std::cerr << "lock: " << lock_detail << "\n";
    std::cerr << "refused: nothing was written\n";
    return kPolicyRefused;
  }

  // 같은 초에 두 번 실행해도 기존 폴더를 덮어쓰지 않도록 번호를 붙인다.
  if (output_arg.empty()) {
    const fs::path base = out_dir;
    std::error_code ec;
    for (int n = 2; fs::exists(out_dir, ec); ++n) {
      out_dir = base;
      out_dir += "_" + std::to_string(n);
    }
  }

  meta.notes.push_back("P0 build: no solver yet; input and output policy checks only");
  std::string error;
  if (!write_run_folder(out_dir, meta, error)) {
    std::cerr << "error: " << error << "\n";
    return kIoError;
  }
  std::cout << "prepared: " << out_dir.string() << "\n"
            << "status: " << run_status_name(meta.status)
            << " (no calculation was performed)\n";
  return kOk;
}

}  // namespace

int main(int argc, char** argv) {
  const std::vector<std::string> args(argv + 1, argv + argc);
  if (args.empty() || args[0] == "help" || args[0] == "--help" || args[0] == "-h") {
    print_usage();
    return args.empty() ? kUsage : kOk;
  }
  const std::string& cmd = args[0];

  if (cmd == "info" && args.size() == 1) return cmd_info();
  if (cmd == "check" && args.size() == 2) return cmd_check(args[1]);
  if (cmd == "lock" && args.size() == 2) return cmd_lock(args[1]);
  if (cmd == "prepare") {
    if (args.size() == 2) return cmd_prepare(args[1], "");
    if (args.size() == 4 && args[2] == "--output-dir") return cmd_prepare(args[1], args[3]);
  }
  std::cerr << "error: invalid command line\n\n";
  print_usage();
  return kUsage;
}
