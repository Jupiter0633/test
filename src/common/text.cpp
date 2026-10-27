#include "common/text.h"

#include <cctype>
#include <ctime>
#include <fstream>
#include <iterator>

namespace hcw {

std::string trim(std::string_view s) {
  std::size_t b = 0, e = s.size();
  while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
  while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
  return std::string(s.substr(b, e - b));
}

std::string to_lower(std::string_view s) {
  std::string out(s);
  for (char& c : out) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return out;
}

bool starts_with(std::string_view s, std::string_view prefix) {
  return s.size() >= prefix.size() && s.substr(0, prefix.size()) == prefix;
}

std::string normalize_text(std::string_view raw) {
  if (starts_with(raw, "\xEF\xBB\xBF")) raw.remove_prefix(3);
  std::string out;
  out.reserve(raw.size());
  for (std::size_t i = 0; i < raw.size(); ++i) {
    if (raw[i] == '\r') {
      out.push_back('\n');
      if (i + 1 < raw.size() && raw[i + 1] == '\n') ++i;
    } else {
      out.push_back(raw[i]);
    }
  }
  return out;
}

bool read_file(const std::filesystem::path& path, std::string& out) {
  std::ifstream in(path, std::ios::binary);
  if (!in) return false;
  out.assign(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
  return !in.bad();
}

namespace {

std::tm utc_now_tm() {
  const std::time_t now = std::time(nullptr);
  std::tm tm{};
#ifdef _WIN32
  gmtime_s(&tm, &now);
#else
  gmtime_r(&now, &tm);
#endif
  return tm;
}

std::string format_now(const char* fmt) {
  const std::tm tm = utc_now_tm();
  char buf[32];
  std::strftime(buf, sizeof(buf), fmt, &tm);
  return buf;
}

}  // namespace

std::string utc_now_iso() { return format_now("%Y-%m-%dT%H:%M:%SZ"); }
std::string utc_now_compact() { return format_now("%Y%m%dT%H%M%SZ"); }

}  // namespace hcw
