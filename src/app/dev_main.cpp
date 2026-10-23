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
#include "verify/p4_cases.h"

using namespace hcw;

namespace {

void usage() {
  std::cout << "usage: hcw_dev entry-dt <a/dx> <dt factor> [--threads N]\n"
               "  entry-dt  disk entry (P4 3-B conditions) with the time-step rule scaled by <dt factor>;\n"
               "            prints impulse, peaks and the 3-F surface vs control-volume force (development only)\n";
}

}  // namespace

int main(int argc, char** argv) {
  std::vector<std::string> args(argv + 1, argv + argc);
  std::string thread_error;
  int threads = parallel::threads_from_environment(thread_error);
  for (auto it = args.begin(); it != args.end();) {
    if (*it == "--threads") {
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
  if (args.size() != 3 || args[0] != "entry-dt") {
    usage();
    return 1;
  }
  const int n = std::atoi(args[1].c_str());
  const double f = std::atof(args[2].c_str());
  if (n < 4 || !(f > 0.0)) {
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
  const verify::TestReport rep = verify::entry_dev_study(n, f);
  for (const auto& r : rep.runs) {
    std::cout << "[" << r.label << "]" << (r.failure.empty() ? "" : "  " + r.failure) << "\n";
    for (const auto& [k, v] : r.metrics) std::cout << "  " << k << " = " << v << "\n";
    std::cout << "  steps " << r.stats.steps << ", wall time " << r.seconds << " s\n";
  }
  for (const auto& note : rep.notes) std::cout << "note: " << note << "\n";
  std::cout << "verdict  : " << rep.verdict << "\n";
  return rep.verdict == "DEV_ONLY" ? 0 : 4;
}
