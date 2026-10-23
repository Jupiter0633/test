// hcw_verify: 검증 시험 실행기 (P1, P2, P3).
// 기준: docs/design/03_verification_plan.md §1.1 (P1, 고정 커밋 35a77ac), §1.2 (P2, 고정 커밋 93cf47d)
// 결과는 화면과 results/verification/<단계>/<시험>/<UTC 시각>/report.txt 에 남긴다.

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "common/parallel.h"
#include "common/text.h"
#include "io/case_policy.h"
#include "io/run_metadata.h"
#include "verify/p1_cases.h"
#include "verify/p2_cases.h"
#include "verify/p3_cases.h"
#include "verify/p4_cases.h"

namespace fs = std::filesystem;
using namespace hcw;
using namespace hcw::verify;

namespace {

std::string render(const TestReport& rep, const std::string& header) {
  std::ostringstream o;
  o.precision(6);
  o << header;
  o << "test     : " << rep.id << " - " << rep.title << "\n"
    << "criteria : docs/design/03_verification_plan.md section " << rep.criteria_section
    << " (locked in commit " << rep.criteria_commit << ")\n\n";

  // 지표 이름 모음
  std::vector<std::string> keys;
  for (const auto& run : rep.runs) {
    for (const auto& [k, v] : run.metrics) {
      if (std::find(keys.begin(), keys.end(), k) == keys.end()) keys.push_back(k);
    }
  }
  for (const auto& run : rep.runs) {
    o << "[" << run.label << "]";
    if (!run.failure.empty()) o << "  " << run.failure;
    o << "\n";
    for (const auto& k : keys) {
      const auto it = run.metrics.find(k);
      if (it != run.metrics.end()) o << "  " << k << " = " << it->second << "\n";
    }
    o << "  poisson: solves " << run.stats.poisson_solves << ", total iterations "
      << run.stats.poisson_iterations << ", refinements " << run.stats.poisson_refinements
      << ", max iterations " << run.stats.max_poisson_iterations << ", max relative residual "
      << run.stats.max_relative_residual << "\n"
      << "  steps " << run.stats.steps << ", wall time " << run.seconds << " s\n";
  }
  // 연속 격자 사이 관측 차수 (참고용)
  // 1-D 는 격자 세분화 시험이 아니므로 차수를 출력하지 않는다. P2 는 시험 안에서 따로 보고한다.
  if (rep.runs.size() >= 2 && rep.id != "1-D" && rep.phase == "P1") {
    o << "\nobserved orders between successive grids:\n";
    for (const auto& k : keys) {
      if (k == "steps" || k == "final_step_change") continue;
      o << "  " << k << ":";
      for (std::size_t i = 0; i + 1 < rep.runs.size(); ++i) {
        const auto a = rep.runs[i].metrics.find(k), b = rep.runs[i + 1].metrics.find(k);
        if (a == rep.runs[i].metrics.end() || b == rep.runs[i + 1].metrics.end()) continue;
        o << " " << observed_order(a->second, b->second);
      }
      o << "\n";
    }
  }
  o << "\ncriteria:\n";
  const bool any_failed = std::any_of(rep.runs.begin(), rep.runs.end(),
                                      [](const GridRun& r) { return !r.failure.empty(); });
  if (rep.criteria.empty()) {
    o << (any_failed ? "  (not evaluated because a run failed)\n" : "  (none: report only)\n");
  }
  for (const auto& c : rep.criteria) {
    o << "  [" << (c.passed ? "PASS" : "FAIL") << "] " << c.text << "  (" << c.detail << ")\n";
  }
  for (const auto& n : rep.notes) o << "note: " << n << "\n";
  o << "\nverdict  : " << rep.verdict << "\n";
  return o.str();
}

std::string header_text(bool quick) {
  const BuildInfo b = current_build_info();
  std::ostringstream o;
  o << "hcw_verify\n"
    << "created  : " << utc_now_iso() << "\n"
    << "build    : " << b.git_commit << " (" << b.git_tree_state << "), " << b.compiler << " ["
    << b.build_config << "], " << b.platform << "\n"
    << "machine  : " << current_machine_name() << ", hardware threads "
    << std::thread::hardware_concurrency() << "\n"
    << "threads  : " << parallel::describe() << "\n"
    << "mode     : " << (quick ? "QUICK (reduced grids, not a verification verdict)" : "FULL")
    << "\n";
  return o.str();
}

void save(const std::string& phase, const std::string& id, const std::string& text, bool quick) {
  const auto root = find_project_root(fs::current_path());
  if (!root) {
    std::cout << "(not saved: run from inside the repository to save a report)\n";
    return;
  }
  std::string folder = id;
  for (char& c : folder) {
    if (c == '(' || c == ')') c = '_';
  }
  const fs::path dir = *root / "results" / "verification" / phase / folder /
                       (utc_now_compact() + (quick ? "_quick" : ""));
  std::error_code ec;
  fs::create_directories(dir, ec);
  std::ofstream out(dir / "report.txt", std::ios::binary);
  out << text;
  std::cout << "saved: " << (dir / "report.txt").string() << "\n";
}

void usage() {
  std::cout << "usage: hcw_verify <test | P1 | P2 | P3 | P4> [--quick] [--threads N]\n"
               "  P1 tests: 1A, 1A-a, 1A-b, 1B, 1C, 1D   ('all' = P1, kept for compatibility)\n"
               "  P2 tests: 1E, 1F, 1G, 1H\n"
               "  P3 tests: 2A, 2B, 2C, 2D, 2F, 2H, 2A-B, 2C-B, 2F-B (-B: new convection scheme B)\n"
               "  P4 tests: 2T, 3A, 3E, 3I (report only), 3B, 3C, 3F, 3G, S1-S5 (stretched grids)   (3-D is reported inside the entry tests)\n"
               "  --quick   reduced grids for a fast functional check (no verdict)\n"
               "  --threads N   use N CPU threads (default 1, or the HCW_THREADS environment\n"
               "                variable). Results are bitwise identical for any N.\n"
               "exit code: 0 if every requested test is VERIFIED (or quick mode), 4 otherwise\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  bool quick = false;
  std::string thread_error;
  int threads = parallel::threads_from_environment(thread_error);
  for (auto it = args.begin(); it != args.end();) {
    if (*it == "--quick") {
      quick = true;
      it = args.erase(it);
    } else if (*it == "--threads") {
      it = args.erase(it);
      char* end = nullptr;
      const long n = it == args.end() ? 0 : std::strtol(it->c_str(), &end, 10);
      if (n < 1 || n > 256 || *end != '\0') {
        std::cerr << "error: --threads needs an integer from 1 to 256\n";
        return 1;
      }
      threads = static_cast<int>(n);
      thread_error.clear();
      it = args.erase(it);
    } else {
      ++it;
    }
  }
  if (threads < 1) {
    std::cerr << "error: " << thread_error << "\n";
    return 1;
  }
  parallel::set_threads(threads);
  if (args.size() != 1) {
    usage();
    return 1;
  }
  const std::string what = args[0];
  std::vector<std::string> tests;
  if (what == "all" || what == "P1") tests = {"1A-a", "1A-b", "1B", "1C", "1D"};
  else if (what == "P2") tests = {"1E", "1F", "1G", "1H"};
  else if (what == "P3") tests = {"2A", "2B", "2C", "2D", "2F", "2H", "2A-B", "2C-B", "2F-B"};
  else if (what == "P4") tests = {"S1", "S2", "S3", "S4", "S5", "2T", "3A", "3B", "3C", "3E", "3F", "3G", "3I"};
  else if (what == "1A") tests = {"1A-a", "1A-b"};
  else if (what == "1A-a" || what == "1A-b" || what == "1B" || what == "1C" || what == "1D" ||
           what == "1E" || what == "1F" || what == "1G" || what == "1H" || what == "2A" || what == "2B" || what == "2C" || what == "2D" || what == "2F" ||
           what == "2H" || what == "2T" || what == "2A-B" || what == "2C-B" || what == "2F-B" || what == "3A" ||
           what == "3B" || what == "3C" || what == "3E" || what == "3F" || what == "3G" || what == "3I" ||
           what == "S1" || what == "S2" || what == "S3" || what == "S4" || what == "S5")
    tests = {what};
  else {
    usage();
    return 1;
  }

  const std::string header = header_text(quick);
  bool all_ok = true;
  for (const auto& t : tests) {
    std::cout << "running " << t << (quick ? " (quick)" : "") << " ...\n" << std::flush;
    TestReport rep;
    if (t == "1A-a") rep = run_1a(false, quick);
    if (t == "1A-b") rep = run_1a(true, quick);
    if (t == "1B") rep = run_1b(quick);
    if (t == "1C") rep = run_1c(quick);
    if (t == "1D") rep = run_1d(quick);
    if (t == "1E") rep = run_1e(quick);
    if (t == "1F") rep = run_1f(quick);
    if (t == "1G") rep = run_1g(quick);
    if (t == "1H") rep = run_1h(quick);
    if (t == "2B") rep = run_2b(quick);
    if (t == "2A") rep = run_2a(quick);
    if (t == "2C") rep = run_2c(quick);
    if (t == "2D") rep = run_2d(quick);
    if (t == "2F") rep = run_2f(quick);
    if (t == "2H") rep = run_2h(quick);
    if (t == "2T") rep = run_2t(quick);
    if (t == "2A-B") rep = run_2a(quick, true);
    if (t == "2C-B") rep = run_2c(quick, true);
    if (t == "2F-B") rep = run_2f(quick, true);
    if (t == "3A") rep = run_3a(quick);
    if (t == "3B") rep = run_3b(quick);
    if (t == "3C") rep = run_3c(quick);
    if (t == "3E") rep = run_3e(quick);
    if (t == "3F") rep = run_3f(quick);
    if (t == "3G") rep = run_3g(quick);
    if (t == "3I") rep = run_3i(quick);
    if (t.size() == 2 && t[0] == 'S') rep = run_s(std::string("S-") + t[1], quick);
    const std::string text = render(rep, header);
    std::cout << text;
    if (!quick) save(rep.phase, rep.id, text, quick);
    std::cout << "\n";
    if (!quick && rep.verdict != "VERIFIED" && rep.verdict != "REPORT_ONLY") all_ok = false;
    if (quick && rep.verdict != "QUICK_CHECK_ONLY") all_ok = false;
  }
  return all_ok ? 0 : 4;
}
