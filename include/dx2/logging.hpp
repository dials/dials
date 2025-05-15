#pragma once

#include <cstdlib>
#include <fmt/color.h>
#include <fmt/core.h>
#include <string_view>

namespace dx2_log {

// Returns true if DX2_DEBUG is set
inline bool is_debug_enabled() { return std::getenv("DX2_DEBUG") != nullptr; }

inline void info(std::string_view msg) {
  fmt::print(fg(fmt::color::green) | fmt::emphasis::bold, "[INFO] ");
  fmt::print("{}\n", msg);
}

inline void warning(std::string_view msg) {
  fmt::print(stderr, fg(fmt::color::yellow) | fmt::emphasis::bold,
             "[WARNING] ");
  fmt::print(stderr, "{}\n", msg);
}

inline void error(std::string_view msg) {
  fmt::print(stderr, fg(fmt::color::red) | fmt::emphasis::bold, "[ERROR] ");
  fmt::print(stderr, "{}\n", msg);
}

inline void debug(std::string_view msg) {
  if (is_debug_enabled()) {
    fmt::print(fg(fmt::color::blue) | fmt::emphasis::bold, "[DEBUG] ");
    fmt::print("{}\n", msg);
  }
}

// Optionally, overloads that support format-style args:
template <typename... Args>
inline void debug(fmt::format_string<Args...> fmt_str, Args &&...args) {
  if (is_debug_enabled()) {
    fmt::print(fg(fmt::color::blue) | fmt::emphasis::bold, "[DEBUG] ");
    fmt::print(fmt_str, std::forward<Args>(args)...);
    fmt::print("\n");
  }
}

} // namespace dx2_log