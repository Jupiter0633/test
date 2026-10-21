#pragma once

#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace hcw {

// 입력 파일 형식 (INI 유사):
//   # 주석 (줄 어디서든 '#' 이후는 무시)
//   [section]
//   key = value
// 키와 섹션 이름은 소문자·숫자·'_'·'-' 만 허용 (Windows 컴퓨터 이름의 '-' 때문). 섹션 밖의 키, 중복 키, 빈 값은 오류.
struct CaseEntry {
  std::string value;
  int line = 0;
};

class CaseFile {
 public:
  using Section = std::map<std::string, CaseEntry>;

  // 정규화된 텍스트(normalize_text 적용 후)를 해석한다. 오류는 errors 에 누적.
  static CaseFile parse(std::string_view text, std::vector<std::string>& errors);

  bool has_section(const std::string& section) const;
  std::optional<std::string> get(const std::string& section, const std::string& key) const;
  const std::map<std::string, Section>& sections() const { return sections_; }

 private:
  std::map<std::string, Section> sections_;
};

}  // namespace hcw
