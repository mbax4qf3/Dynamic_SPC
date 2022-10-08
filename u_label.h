#ifndef SPC_U_LABEL_H_
#define SPC_U_LABEL_H_

#include <cstdint>
#include <string>
#include <vector>

#include "macros.h"

namespace spc {
using uint32_t = std::uint32_t;
using uint64_t = std::uint64_t;

uint32_t constexpr kNumVBits = 25; //23
uint32_t constexpr kNumDBits = 10; //10
uint32_t constexpr kNumCBits = 29; //31
uint32_t constexpr kUBC = (static_cast<uint32_t>(1) << kNumCBits) - 1;

struct LabelEntry final {
  uint64_t v_d_c;
};

inline LabelEntry LEMerge(const uint32_t v,
                          const uint32_t d,
                          const uint32_t c) {
  return {(((static_cast<uint64_t>(v) << kNumDBits) | (d)) << kNumCBits) | (c)};
}

inline uint32_t LEExtractV(const LabelEntry& le) {
  return static_cast<uint32_t>(le.v_d_c >> (kNumDBits + kNumCBits));
}

inline uint32_t LEExtractD(const LabelEntry& le) {
  uint32_t mask = (static_cast<uint32_t>(1) << kNumDBits) - 1;
  return static_cast<uint32_t>((le.v_d_c >> kNumCBits) & mask);
}

inline uint32_t LEExtractC(const LabelEntry& le) {
  uint32_t mask = (static_cast<uint32_t>(1) << kNumCBits) - 1;
  return static_cast<uint32_t>(le.v_d_c & mask);
}

inline void NormalV(const uint32_t v) {
  if (v >= (static_cast<uint32_t>(1) << kNumVBits)) {
    std::string msg = "too many vertices: " + std::to_string(v);
    ASSERT_INFO(false, msg.c_str());
  }
}

inline void NormalD(const uint32_t d) {
  if (d >= (static_cast<uint32_t>(1) << kNumDBits)) {
    std::string msg = "large distance: " + std::to_string(d);
    ASSERT_INFO(false, msg.c_str());
  }
}

inline void NormalC(const uint32_t c) {
  if (c >= (static_cast<uint32_t>(1) << kNumCBits)) {
    std::string msg = "large count: " + std::to_string(c);
    ASSERT_INFO(false, msg.c_str());
  }
}

using Graph = std::vector<std::vector<uint32_t>>;
using Label = std::vector<std::vector<LabelEntry>>;
} // namespace spc

#endif
