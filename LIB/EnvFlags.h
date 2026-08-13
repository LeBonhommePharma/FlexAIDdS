// EnvFlags.h — one parser for FLEXAIDDS_* boolean environment switches
//
// The codebase had grown three incompatible conventions for the same idea:
//
//   std::atoi(e) != 0                      (BindingMode.cpp, cluster.cpp)
//   e && e[0] != '\0'                      (BindingMode.cpp election weight)
//   std::getenv(name) != nullptr           (FLEXAIDDS_NO_TENCOM, debug flags)
//
// They disagree on the values a scientist is most likely to type: under atoi,
// FLEXAIDDS_ELECT_LEGACY_ACF=true parses as 0 and silently selects the OPPOSITE
// arm, while under the presence-only convention FLEXAIDDS_NO_TENCOM=0 still
// turns the feature off. An A/B control that silently means its own negation is
// worse than no control at all.
//
// env_bool() accepts 1/true/yes/on and 0/false/no/off, case-insensitively, and
// returns the supplied default when the variable is unset, empty, or
// unparseable — so an unrecognised value can never flip an arm by accident.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstdlib>
#include <cstring>
#include <cctype>

namespace flexaids {

/// Parse a FLEXAIDDS_* boolean switch. Unset/empty/unparseable → `fallback`.
inline bool env_bool(const char* name, bool fallback = false) noexcept
{
    const char* raw = std::getenv(name);
    if (raw == nullptr) return fallback;

    // Trim surrounding whitespace; shell exports pick it up more often than
    // you would like (FLEXAIDDS_X=" 1").
    while (*raw != '\0' && std::isspace(static_cast<unsigned char>(*raw))) ++raw;
    if (*raw == '\0') return fallback;

    char buf[16];
    std::size_t n = 0;
    while (raw[n] != '\0' && n < sizeof(buf) - 1 &&
           !std::isspace(static_cast<unsigned char>(raw[n]))) {
        buf[n] = static_cast<char>(std::tolower(static_cast<unsigned char>(raw[n])));
        ++n;
    }
    buf[n] = '\0';

    if (std::strcmp(buf, "1") == 0 || std::strcmp(buf, "true") == 0 ||
        std::strcmp(buf, "yes") == 0 || std::strcmp(buf, "on") == 0)
        return true;
    if (std::strcmp(buf, "0") == 0 || std::strcmp(buf, "false") == 0 ||
        std::strcmp(buf, "no") == 0 || std::strcmp(buf, "off") == 0)
        return false;

    return fallback;
}

/// FLEXAIDDS_PARALLEL_REPRODUCE — DEFAULT OFF (METHODOLOGY §1).
/// Explicit 1/true/on opts in; unset/empty/0 stays serial inline eval.
inline bool parallel_reproduce_enabled() noexcept
{
    return env_bool("FLEXAIDDS_PARALLEL_REPRODUCE", false);
}

}  // namespace flexaids
