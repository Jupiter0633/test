#pragma once

#include <string>
#include <string_view>

namespace hcw {

// FIPS 180-4 SHA-256. 입력 파일 잠금(해시)에 사용한다.
std::string sha256_hex(std::string_view data);

}  // namespace hcw
