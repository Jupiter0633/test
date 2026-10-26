// 개발용 실험 프로그램 (판정 아님, 결과 파일을 남기지 않는다).
// 정식 판정은 hcw_verify 로만 한다. 여기 결과는 화면 출력뿐이며 기록은 사람이 옮긴다.

#include <cstdlib>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

#include "common/parallel.h"
#include "common/text.h"
#include "io/case_policy.h"
#include "io/run_metadata.h"
#include "verify/p2_cases.h"
#include "verify/p4_cases.h"

using namespace hcw;

namespace {

void usage() {
  std::cout << "usage: hcw_dev entry-dt <a/dx> <dt factor> [--threads N]\n"
               "       hcw_dev entry-snap <a/dx> <dt factor> <out dir> <interval ms> [no-water-floor] [--threads N]\n"
               "  entry-dt  disk entry (P4 3-B conditions) with the time-step rule scaled by <dt factor>;\n"
               "            prints impulse, peaks and the 3-F surface vs control-volume force (development only)\n"
               "  entry-snap  same run, writing cell snapshots and history.csv for tools/plot_entry.py\n"
               "  --stretched  P4 stretched grid (a/dx is then a/dx_core); add it to either command\n"
               "  --history <file>  entry-dt only: write the half-step F_s / F_cv history as CSV\n"
               "  --cut-cell  entry-dt only: cut-cell moving boundary switch (development)\n"
               "       hcw_dev p2-cut <quick|full> <on|off>  1-G and 1-E with the cut-cell switch (development)\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  std::string thread_error;
  int threads = parallel::threads_from_environment(thread_error);
  bool stretched = false;  // P4 새 조건의 늘어나는 격자 (n = a/Δ_core)
  bool water_floor = true;  // --no-water-floor: 물 압력 하한 끔 (entry-dt 개발 비교용)
  bool cut_cell = false;  // --cut-cell: 절단 셀 이동 경계 스위치 (07 설계, 개발 비교용)
  std::string history_csv;  // --history <file>: entry-dt 의 반 단계 F_s / F_cv 전체 이력
  for (auto it = args.begin(); it != args.end();) {
    if (*it == "--no-water-floor") {
      water_floor = false;
      it = args.erase(it);
    } else if (*it == "--cut-cell") {
      cut_cell = true;
      it = args.erase(it);
    } else if (*it == "--history") {
      it = args.erase(it);
      if (it == args.end()) {
        std::cerr << "error: --history needs a file name\n";
        return 1;
      }
      history_csv = *it;
      it = args.erase(it);
    } else if (*it == "--stretched") {
      stretched = true;
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
      it = args.erase(it);
    } else {
      ++it;
    }
  }
  if (threads < 1) threads = 1;
  parallel::set_threads(threads);
  const bool snap = !args.empty() && args[0] == "entry-snap";
  const bool p2cut = args.size() == 3 && args[0] == "p2-cut" && (args[1] == "quick" || args[1] == "full") &&
                     (args[2] == "on" || args[2] == "off");
  if (!((args.size() == 3 && args[0] == "entry-dt") || ((args.size() == 5 || args.size() == 6) && snap) || p2cut)) {
    usage();
    return 1;
  }
  const int n = p2cut ? 4 : std::atoi(args[1].c_str());
  const double f = p2cut ? 1.0 : std::atof(args[2].c_str());
  if (n < 4 || !(f > 0.0) || (snap && !(std::atof(args[4].c_str()) > 0.0))) {
    usage();
    return 1;
  }
  const BuildInfo b = current_build_info();
  std::cout << "hcw_dev (development only, not a verification verdict)\n"
            << "created  : " << utc_now_iso() << "\n"
            << "build    : " << b.git_commit << " (" << b.git_tree_state << "), " << b.compiler << " ["
            << b.build_config << "], " << b.platform << "\n"
            << "machine  : " << current_machine_name() << ", hardware threads "
            << std::thread::hardware_concurrency() << "\n"
            << "threads  : " << parallel::describe() << "\n";
  const verify::TestReport rep =
      p2cut ? verify::dev_p2_cut_cell(args[1] == "quick", args[2] == "on")
      : snap ? verify::entry_dev_snapshots(n, f, args[3], 1e-3 * std::atof(args[4].c_str()),
                                         !(args.size() == 6 && args[5] == "no-water-floor"), stretched)
           : verify::entry_dev_study(n, f, stretched, water_floor, history_csv, cut_cell);
  for (const auto& r : rep.runs) {
    std::cout << "[" << r.label << "]" << (r.failure.empty() ? "" : "  " + r.failure) << "\n";
    for (const auto& [k, v] : r.metrics) std::cout << "  " << k << " = " << v << "\n";
    std::cout << "  steps " << r.stats.steps << ", wall time " << r.seconds << " s\n";
  }
  for (const auto& note : rep.notes) std::cout << "note: " << note << "\n";
  std::cout << "verdict  : " << rep.verdict << "\n";
  return rep.verdict == "DEV_ONLY" ? 0 : 4;
}
