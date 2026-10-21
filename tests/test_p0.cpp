// P0 단위 시험: 해시, 텍스트 정규화, 입력 파일 해석, 입력 규칙, 잠금, 출력 허용 규칙.
// 외부 시험 프레임워크 없이 동작한다. 실패 항목 수를 종료 코드로 돌려준다.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "common/sha256.h"
#include "common/text.h"
#include "io/case_file.h"
#include "io/case_policy.h"

namespace fs = std::filesystem;
using namespace hcw;

namespace {

int g_failures = 0;
int g_checks = 0;

void expect(bool condition, const std::string& what) {
  ++g_checks;
  if (!condition) {
    ++g_failures;
    std::cerr << "FAIL: " << what << "\n";
  }
}

void write_text(const fs::path& p, const std::string& text) {
  fs::create_directories(p.parent_path());
  std::ofstream out(p, std::ios::binary);
  out << text;
}

// 시험마다 새 임시 프로젝트 루트를 만든다.
fs::path make_root(const std::string& tag) {
  const fs::path root = fs::temp_directory_path() / ("hcw_test_" + tag);
  fs::remove_all(root);
  fs::create_directories(root);
  write_text(root / kRootMarker, "");
  return root;
}

std::string case_text(const std::string& name, const std::string& purpose) {
  return "[case]\nname = " + name + "\npurpose = " + purpose + "\n";
}

void test_sha256() {
  // NIST FIPS 180-4 예제 (빈 문자열, "abc", 448비트 메시지, 'a' 100만 개)
  expect(sha256_hex("") ==
             "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
         "sha256 empty");
  expect(sha256_hex("abc") ==
             "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
         "sha256 abc");
  expect(sha256_hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
             "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1",
         "sha256 448-bit");
  expect(sha256_hex(std::string(1000000, 'a')) ==
             "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0",
         "sha256 million a");
  // 채움(padding) 경계 길이. 기댓값은 Python hashlib 으로 독립 계산.
  expect(sha256_hex(std::string(55, 'a')) ==
             "9f4390f8d30c2dd92ec9f095b65e2b9ae9b0a925a5258e241c9f1e910f734318",
         "sha256 55 bytes");
  expect(sha256_hex(std::string(56, 'a')) ==
             "b35439a4ac6f0948b6d6f9e3c6af0f5f590ce20f1bde7090ef7970686ec6738a",
         "sha256 56 bytes");
  expect(sha256_hex(std::string(63, 'a')) ==
             "7d3e74a05d7db15bce4ad9ec0658ea98e3f06eeecf16b4c6fff2da457ddc2f34",
         "sha256 63 bytes");
  expect(sha256_hex(std::string(64, 'a')) ==
             "ffe054fe7ae0cb6dc65c3af9b61d5209f439851db43d0ba5997337df154668eb",
         "sha256 64 bytes");
}

void test_normalize() {
  expect(normalize_text("a\r\nb\r\n") == "a\nb\n", "CRLF to LF");
  expect(normalize_text("a\rb") == "a\nb", "CR to LF");
  expect(normalize_text("\xEF\xBB\xBF[case]\n") == "[case]\n", "BOM removed");
  expect(normalize_text("x\n") == "x\n", "LF unchanged");
}

void test_parser() {
  std::vector<std::string> errors;
  const CaseFile ok = CaseFile::parse(
      "# comment\n[case]\nname = A_1  # trailing comment\npurpose = \"paper\"\n\n[geometry]\n"
      "outer_radius = 0.1\n",
      errors);
  expect(errors.empty(), "valid file parses without errors");
  expect(ok.get("case", "name") == std::optional<std::string>("A_1"), "value with comment");
  expect(ok.get("case", "purpose") == std::optional<std::string>("paper"), "quoted value");
  expect(ok.get("geometry", "outer_radius") == std::optional<std::string>("0.1"),
         "second section");
  expect(!ok.get("case", "missing").has_value(), "missing key");

  const std::vector<std::pair<std::string, std::string>> bad = {
      {"name = x\n", "key outside section"},
      {"[case]\nname\n", "missing '='"},
      {"[case]\nName = x\n", "uppercase key"},
      {"[case]\nname =\n", "empty value"},
      {"[case]\nname = a\nname = b\n", "duplicate key"},
      {"[case]\n[case]\n", "duplicate section"},
      {"[Case]\n", "uppercase section"},
      {"[case\n", "unterminated section"},
  };
  for (const auto& [text, label] : bad) {
    std::vector<std::string> errs;
    CaseFile::parse(text, errs);
    expect(!errs.empty(), "parser rejects: " + label);
  }
}

void test_load_case_rules() {
  const fs::path root = make_root("load");
  Diagnostics d;

  write_text(root / "cases/dev/DEV_a.case", case_text("DEV_a", "development"));
  expect(load_case(root / "cases/dev/DEV_a.case", root, d).has_value() && d.ok(),
         "valid development case");

  struct Bad {
    std::string rel, text, label;
  };
  const std::vector<Bad> bad = {
      {"cases/dev/a.case", case_text("a", "development"), "development without DEV_"},
      {"cases/paper/DEV_b.case", case_text("DEV_b", "paper"), "paper with DEV_ prefix"},
      {"cases/paper/DEV_c.case", case_text("DEV_c", "development"), "development in paper dir"},
      {"cases/dev/DEV_d.case", case_text("DEV_x", "development"), "name differs from file"},
      {"cases/dev/DEV_e.case", case_text("DEV_e", "draft"), "unknown purpose"},
      {"cases/dev/DEV_f.case", "[case]\nname = DEV_f\n", "missing purpose"},
      {"cases/dev/DEV_g.case", "[case]\npurpose = development\n", "missing name"},
      {"cases/dev/DEV h.case", case_text("DEV h", "development"), "space in name"},
      {"other/DEV_i.case", case_text("DEV_i", "development"), "outside cases folder"},
  };
  for (const auto& b : bad) {
    write_text(root / b.rel, b.text);
    Diagnostics diag;
    expect(!load_case(root / b.rel, root, diag).has_value() && !diag.ok(),
           "load_case rejects: " + b.label);
  }

  Diagnostics warn;
  write_text(root / "cases/verification/V_a.case",
             case_text("V_a", "verification") + "extra = 1\n[geometry]\nr = 1\n");
  expect(load_case(root / "cases/verification/V_a.case", root, warn).has_value() &&
             warn.warnings.size() == 2,
         "unused keys and sections give warnings, not errors");

  expect(find_project_root(root / "cases/dev/DEV_a.case") == root, "project root found");
  fs::remove_all(root);
}

void test_lock() {
  const fs::path root = make_root("lock");
  const fs::path file = root / "cases/paper/P_a.case";
  write_text(file, case_text("P_a", "paper"));

  Diagnostics d;
  auto c = load_case(file, root, d);
  expect(c.has_value(), "paper case loads");
  if (!c) return;

  std::string detail, error;
  expect(check_lock(*c, detail) == LockState::Missing, "lock missing before locking");
  expect(write_lock(*c, error), "lock written");
  expect(check_lock(*c, detail) == LockState::Match, "lock matches after locking");
  expect(!write_lock(*c, error), "re-locking is refused");

  // 줄바꿈만 CRLF 로 바뀐 경우(Windows git 변환)는 같은 입력으로 본다.
  write_text(file, "[case]\r\nname = P_a\r\npurpose = paper\r\n");
  auto crlf = load_case(file, root, d);
  expect(crlf && check_lock(*crlf, detail) == LockState::Match, "CRLF copy still matches");

  // 내용이 바뀌면 불일치.
  write_text(file, case_text("P_a", "paper") + "description = changed\n");
  auto changed = load_case(file, root, d);
  expect(changed && check_lock(*changed, detail) == LockState::Mismatch,
         "modified case is detected");

  write_text(lock_path_for(file), "[lock]\nsha256 = short\n");
  expect(changed && check_lock(*changed, detail) == LockState::Invalid, "malformed lock");

  const fs::path dev = root / "cases/dev/DEV_a.case";
  write_text(dev, case_text("DEV_a", "development"));
  auto dev_case = load_case(dev, root, d);
  expect(dev_case && !write_lock(*dev_case, error), "development cases cannot be locked");
  fs::remove_all(root);
}

void test_machine_registry() {
  const fs::path root = make_root("machines");
  Diagnostics d;
  expect(lookup_machine(root, "PC1", d) == MachineRole::Unregistered,
         "missing registry means unregistered");
  write_text(root / "config/machines.ini",
             "[machines]\nlaptop_1 = development\nws_1 = workstation\nbad_1 = server\n");
  Diagnostics d2;
  expect(lookup_machine(root, "WS_1", d2) == MachineRole::Workstation,
         "workstation lookup is case-insensitive");
  expect(lookup_machine(root, "laptop_1", d2) == MachineRole::Development, "development role");
  expect(lookup_machine(root, "other", d2) == MachineRole::Unregistered, "unknown machine");
  Diagnostics d3;
  lookup_machine(root, "bad_1", d3);
  expect(!d3.ok(), "invalid role is an error");
  fs::remove_all(root);
}

void test_is_within() {
  const fs::path base = fs::temp_directory_path() / "hcw_within";
  expect(is_within(base / "results/paper/x", base / "results/paper"), "inside");
  expect(is_within(base / "results/paper", base / "results/paper"), "same folder");
  expect(!is_within(base / "results/paper2/x", base / "results/paper"), "sibling prefix");
  expect(!is_within(base / "results/dev/x", base / "results/paper"), "other folder");
  expect(is_within(base / "results/dev/../paper/x", base / "results/paper"), "'..' resolved");
  expect(!is_within(base / "results", base / "results/paper"), "parent is not inside");
}

void test_output_policy() {
  const fs::path root = fs::temp_directory_path() / "hcw_policy";
  OutputRequest ok_paper;
  ok_paper.purpose = Purpose::Paper;
  ok_paper.root = root;
  ok_paper.output_dir = root / "results/paper/P_a/t";
  ok_paper.lock = LockState::Match;
  ok_paper.machine = MachineRole::Workstation;
  ok_paper.git_commit = "abc123";
  ok_paper.git_tree_state = "clean";
  expect(check_output_policy(ok_paper).ok(), "fully valid paper request allowed");

  auto refused = [&](OutputRequest r, const std::string& label) {
    expect(!check_output_policy(r).ok(), "refused: " + label);
  };
  OutputRequest r = ok_paper;
  r.lock = LockState::Missing;
  refused(r, "paper without lock");
  r = ok_paper;
  r.lock = LockState::Mismatch;
  refused(r, "paper with modified case");
  r = ok_paper;
  r.machine = MachineRole::Development;
  refused(r, "paper on development machine");
  r = ok_paper;
  r.machine = MachineRole::Unregistered;
  refused(r, "paper on unregistered machine");
  r = ok_paper;
  r.git_commit = "unknown";
  refused(r, "paper without commit");
  r = ok_paper;
  r.git_tree_state = "modified";
  refused(r, "paper from modified sources");
  r = ok_paper;
  r.output_dir = root / "results/dev/P_a";
  refused(r, "paper outside results/paper");

  for (Purpose p : {Purpose::Development, Purpose::Verification}) {
    OutputRequest d;
    d.purpose = p;
    d.root = root;
    d.git_commit = "unknown";
    d.git_tree_state = "modified";
    d.output_dir = root / "results/dev/x";
    expect(check_output_policy(d).ok(), purpose_name(p) + " may write outside paper results");
    d.output_dir = root / "results/paper/x";
    refused(d, purpose_name(p) + " into results/paper");
  }
}

}  // namespace

int main() {
  test_sha256();
  test_normalize();
  test_parser();
  test_load_case_rules();
  test_lock();
  test_machine_registry();
  test_is_within();
  test_output_policy();
  std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
  return g_failures == 0 ? 0 : 1;
}
