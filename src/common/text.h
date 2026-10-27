#pragma once

#include <filesystem>
#include <string>
#include <string_view>

namespace hcw {

std::string trim(std::string_view s);
std::string to_lower(std::string_view s);
bool starts_with(std::string_view s, std::string_view prefix);

// UTF-8 BOM 제거 + CRLF/CR → LF. Windows 메모장·git 자동 줄바꿈 변환과
// 무관하게 같은 파일이 같은 해시를 갖도록 해시 전에 적용한다.
std::string normalize_text(std::string_view raw);

// 파일 전체를 바이트 그대로 읽는다. 실패 시 false.
bool read_file(const std::filesystem::path& path, std::string& out);

// 현재 UTC 시각. iso: 2026-10-21T12:34:56Z, compact: 20261021T123456Z
std::string utc_now_iso();
std::string utc_now_compact();

}  // namespace hcw
