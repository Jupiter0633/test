#include "io/case_file.h"

#include <sstream>

#include "common/text.h"

namespace hcw {
namespace {

bool valid_identifier(std::string_view s) {
  if (s.empty()) return false;
  for (char c : s) {
    const bool ok = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') || c == '_' || c == '-';
    if (!ok) return false;
  }
  return true;
}

std::string line_error(int line, const std::string& message) {
  return "line " + std::to_string(line) + ": " + message;
}

}  // namespace

CaseFile CaseFile::parse(std::string_view text, std::vector<std::string>& errors) {
  CaseFile file;
  std::istringstream in{std::string(text)};
  std::string raw;
  std::string current;
  int line = 0;

  while (std::getline(in, raw)) {
    ++line;
    const std::size_t hash = raw.find('#');
    const std::string content = trim(hash == std::string::npos ? raw : raw.substr(0, hash));
    if (content.empty()) continue;

    if (content.front() == '[') {
      if (content.back() != ']') {
        errors.push_back(line_error(line, "section header must end with ']'"));
        continue;
      }
      const std::string name = trim(content.substr(1, content.size() - 2));
      if (!valid_identifier(name)) {
        errors.push_back(line_error(line, "invalid section name '" + name +
                                              "' (use lowercase letters, digits, '_', '-')"));
        continue;
      }
      if (file.sections_.count(name) != 0) {
        errors.push_back(line_error(line, "duplicate section [" + name + "]"));
      }
      current = name;
      file.sections_[current];
      continue;
    }

    const std::size_t eq = content.find('=');
    if (eq == std::string::npos) {
      errors.push_back(line_error(line, "expected 'key = value'"));
      continue;
    }
    if (current.empty()) {
      errors.push_back(line_error(line, "key outside of any [section]"));
      continue;
    }
    const std::string key = trim(content.substr(0, eq));
    std::string value = trim(content.substr(eq + 1));
    if (!valid_identifier(key)) {
      errors.push_back(line_error(line, "invalid key '" + key +
                                            "' (use lowercase letters, digits, '_', '-')"));
      continue;
    }
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
      value = value.substr(1, value.size() - 2);
    }
    if (value.empty()) {
      errors.push_back(line_error(line, "empty value for '" + key + "'"));
      continue;
    }
    Section& section = file.sections_[current];
    if (section.count(key) != 0) {
      errors.push_back(line_error(line, "duplicate key '" + key + "' in [" + current + "]"));
      continue;
    }
    section[key] = CaseEntry{value, line};
  }
  return file;
}

bool CaseFile::has_section(const std::string& section) const {
  return sections_.count(section) != 0;
}

std::optional<std::string> CaseFile::get(const std::string& section,
                                         const std::string& key) const {
  const auto s = sections_.find(section);
  if (s == sections_.end()) return std::nullopt;
  const auto k = s->second.find(key);
  if (k == s->second.end()) return std::nullopt;
  return k->second.value;
}

}  // namespace hcw
