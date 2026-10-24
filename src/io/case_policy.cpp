#include "io/case_policy.h"

#include <cstdlib>
#include <fstream>
#include <set>

#include "common/sha256.h"
#include "common/text.h"

#ifndef _WIN32
#include <unistd.h>
#endif

namespace hcw {
namespace {

const char* kDevPrefix = "DEV_";
const char* kLockSection = "lock";

bool valid_case_name(const std::string& s) {
  if (s.empty()) return false;
  for (char c : s) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
                    (c >= '0' && c <= '9') || c == '_' || c == '-';
    if (!ok) return false;
  }
  return true;
}

fs::path normalized(const fs::path& p) {
  std::error_code ec;
  fs::path abs = fs::absolute(p, ec);
  if (ec) abs = p;
  fs::path canon = fs::weakly_canonical(abs, ec);
  return (ec ? abs : canon).lexically_normal();
}

bool same_path(const fs::path& a, const fs::path& b) {
#ifdef _WIN32
  return to_lower(normalized(a).generic_string()) == to_lower(normalized(b).generic_string());
#else
  return normalized(a) == normalized(b);
#endif
}

}  // namespace

std::optional<Purpose> parse_purpose(const std::string& text) {
  if (text == "development") return Purpose::Development;
  if (text == "verification") return Purpose::Verification;
  if (text == "paper") return Purpose::Paper;
  return std::nullopt;
}

std::string purpose_name(Purpose p) {
  switch (p) {
    case Purpose::Development: return "development";
    case Purpose::Verification: return "verification";
    case Purpose::Paper: return "paper";
  }
  return "unknown";
}

std::string purpose_directory(Purpose p) {
  switch (p) {
    case Purpose::Development: return "dev";
    case Purpose::Verification: return "verification";
    case Purpose::Paper: return "paper";
  }
  return "unknown";
}

std::optional<fs::path> find_project_root(const fs::path& start) {
  fs::path dir = normalized(start);
  std::error_code ec;
  if (!fs::is_directory(dir, ec)) dir = dir.parent_path();
  while (!dir.empty()) {
    if (fs::exists(dir / kRootMarker, ec)) return dir;
    const fs::path parent = dir.parent_path();
    if (parent == dir) break;
    dir = parent;
  }
  return std::nullopt;
}

std::optional<LoadedCase> load_case(const fs::path& case_path, const fs::path& root,
                                    Diagnostics& diag) {
  std::string raw;
  if (!read_file(case_path, raw)) {
    diag.errors.push_back("cannot read case file: " + case_path.string());
    return std::nullopt;
  }
  const std::string text = normalize_text(raw);

  LoadedCase c;
  c.path = normalized(case_path);
  c.sha256 = sha256_hex(text);
  c.data = CaseFile::parse(text, diag.errors);

  const auto name = c.data.get("case", "name");
  const auto purpose_text = c.data.get("case", "purpose");
  if (!c.data.has_section("case")) diag.errors.push_back("missing [case] section");
  if (!name) {
    diag.errors.push_back("[case] name is required");
  } else if (!valid_case_name(*name)) {
    diag.errors.push_back("[case] name '" + *name + "' may only contain A-Z a-z 0-9 _ -");
  } else {
    c.name = *name;
  }
  std::optional<Purpose> purpose;
  if (!purpose_text) {
    diag.errors.push_back("[case] purpose is required (development | verification | paper)");
  } else {
    purpose = parse_purpose(*purpose_text);
    if (!purpose) {
      diag.errors.push_back("[case] purpose '" + *purpose_text +
                            "' is invalid (development | verification | paper)");
    }
  }

  if (purpose && !c.name.empty()) {
    c.purpose = *purpose;
    if (c.path.stem().string() != c.name) {
      diag.errors.push_back("file name '" + c.path.filename().string() +
                            "' must match [case] name '" + c.name + "'");
    }
    const bool dev_name = starts_with(c.name, kDevPrefix);
    if (*purpose == Purpose::Development && !dev_name) {
      diag.errors.push_back("development cases must be named with the 'DEV_' prefix");
    }
    if (*purpose != Purpose::Development && dev_name) {
      diag.errors.push_back("'DEV_' prefix is reserved for development cases");
    }
    const fs::path expected_dir = root / "cases" / purpose_directory(*purpose);
    if (!same_path(c.path.parent_path(), expected_dir)) {
      diag.errors.push_back("purpose '" + purpose_name(*purpose) + "' cases must be stored in " +
                            expected_dir.string());
    }
  }

  // P0 는 [case] 외의 섹션을 아직 해석하지 않는다. 무시됨을 알린다.
  static const std::set<std::string> case_keys = {"name", "purpose", "description"};
  for (const auto& [section, entries] : c.data.sections()) {
    if (section != "case") {
      diag.warnings.push_back("section [" + section + "] is not used yet and will be ignored");
      continue;
    }
    for (const auto& [key, entry] : entries) {
      if (case_keys.count(key) == 0) {
        diag.warnings.push_back("line " + std::to_string(entry.line) + ": unknown key '" + key +
                                "' in [case] will be ignored");
      }
    }
  }

  if (!diag.ok()) return std::nullopt;
  return c;
}

fs::path lock_path_for(const fs::path& case_path) {
  fs::path p = case_path;
  p += ".lock";
  return p;
}

LockState check_lock(const LoadedCase& c, std::string& detail) {
  const fs::path lock = lock_path_for(c.path);
  std::error_code ec;
  if (!fs::exists(lock, ec)) {
    detail = "no lock file (" + lock.filename().string() + ")";
    return LockState::Missing;
  }
  std::string raw;
  if (!read_file(lock, raw)) {
    detail = "cannot read lock file";
    return LockState::Invalid;
  }
  std::vector<std::string> errors;
  const CaseFile lockfile = CaseFile::parse(normalize_text(raw), errors);
  const auto hash = lockfile.get(kLockSection, "sha256");
  if (!errors.empty() || !hash || hash->size() != 64) {
    detail = "lock file is malformed";
    return LockState::Invalid;
  }
  if (*hash != c.sha256) {
    detail = "case file changed after it was locked (locked " + hash->substr(0, 12) +
             "..., now " + c.sha256.substr(0, 12) + "...)";
    return LockState::Mismatch;
  }
  detail = "lock matches";
  return LockState::Match;
}

bool write_lock(const LoadedCase& c, std::string& error) {
  if (c.purpose != Purpose::Paper) {
    error = "only purpose = paper cases can be locked";
    return false;
  }
  const fs::path lock = lock_path_for(c.path);
  std::error_code ec;
  if (fs::exists(lock, ec)) {
    error = "lock file already exists: " + lock.string() +
            " (a locked case must not be re-locked; delete the lock deliberately and record why)";
    return false;
  }
  std::ofstream out(lock, std::ios::binary);
  if (!out) {
    error = "cannot write lock file: " + lock.string();
    return false;
  }
  out << "# hcw case lock. Do not edit. Records the SHA-256 of the normalized case file.\n"
      << "[lock]\n"
      << "case = " << c.path.filename().string() << "\n"
      << "sha256 = " << c.sha256 << "\n"
      << "locked_utc = " << utc_now_iso() << "\n"
      << "machine = " << current_machine_name() << "\n";
  if (!out) {
    error = "failed while writing lock file";
    return false;
  }
  return true;
}

std::string current_machine_name() {
#ifdef _WIN32
  if (const char* name = std::getenv("COMPUTERNAME")) return name;
  return "unknown";
#else
  char buf[256] = {};
  if (gethostname(buf, sizeof(buf) - 1) == 0 && buf[0] != '\0') return buf;
  return "unknown";
#endif
}

MachineRole lookup_machine(const fs::path& root, const std::string& machine,
                           Diagnostics& diag) {
  const fs::path registry = root / "config" / "machines.ini";
  std::string raw;
  if (!read_file(registry, raw)) {
    diag.warnings.push_back("machine registry not found: " + registry.string());
    return MachineRole::Unregistered;
  }
  std::vector<std::string> errors;
  const CaseFile file = CaseFile::parse(normalize_text(raw), errors);
  for (const auto& e : errors) diag.errors.push_back("machines.ini " + e);

  const auto section = file.sections().find("machines");
  if (section == file.sections().end()) return MachineRole::Unregistered;
  // 컴퓨터 이름은 대소문자를 구분하지 않는다 (Windows 는 대문자로 보고함).
  const std::string wanted = to_lower(machine);
  for (const auto& [key, entry] : section->second) {
    if (key != wanted) continue;
    if (entry.value == "workstation") return MachineRole::Workstation;
    if (entry.value == "development") return MachineRole::Development;
    diag.errors.push_back("machines.ini line " + std::to_string(entry.line) + ": role '" +
                          entry.value + "' is invalid (development | workstation)");
  }
  return MachineRole::Unregistered;
}

std::string machine_role_name(MachineRole r) {
  switch (r) {
    case MachineRole::Development: return "development";
    case MachineRole::Workstation: return "workstation";
    case MachineRole::Unregistered: return "unregistered";
  }
  return "unknown";
}

bool is_within(const fs::path& child, const fs::path& parent) {
  const fs::path c = normalized(child);
  const fs::path p = normalized(parent);
  auto ci = c.begin();
  for (auto pi = p.begin(); pi != p.end(); ++pi, ++ci) {
    if (pi->empty()) continue;  // 끝의 '/' 로 생기는 빈 요소
    if (ci == c.end()) return false;
#ifdef _WIN32
    if (to_lower(ci->string()) != to_lower(pi->string())) return false;
#else
    if (*ci != *pi) return false;
#endif
  }
  return true;
}

fs::path default_output_dir(const fs::path& root, const LoadedCase& c,
                            const std::string& timestamp) {
  return root / "results" / purpose_directory(c.purpose) / c.name / timestamp;
}

Diagnostics check_output_policy(const OutputRequest& req) {
  Diagnostics d;
  const fs::path paper_results = req.root / "results" / "paper";
  const bool into_paper = is_within(req.output_dir, paper_results);

  if (req.purpose != Purpose::Paper) {
    if (into_paper) {
      d.errors.push_back("purpose '" + purpose_name(req.purpose) +
                         "' may not write into " + paper_results.string());
    }
    return d;
  }

  if (!into_paper) {
    d.errors.push_back("paper results must be written under " + paper_results.string());
  }
  if (req.lock != LockState::Match) {
    d.errors.push_back("paper case is not locked or was modified after locking");
  }
  if (req.machine != MachineRole::Workstation) {
    d.errors.push_back("paper results may only be produced on a machine registered as "
                       "'workstation' in config/machines.ini");
  }
  if (req.git_commit == "unknown" || req.git_commit.empty()) {
    d.errors.push_back("build has no git commit recorded; paper results need a traceable build");
  }
  if (req.git_tree_state != "clean") {
    d.errors.push_back("build was made from modified sources (git tree '" + req.git_tree_state +
                       "'); paper results need a clean commit");
  }
  return d;
}

}  // namespace hcw
