// PdbCoords.h — Shared strict PDB coordinate decoder (RMSD + PoseBust)
//
// One finite, fail-closed parser for FlexAID compact-negative XYZ spans and
// normal fixed-column PDB coordinates. Used by DatasetRunner RMSD and PoseBust
// loaders so both validators consume identical geometry.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <string>

namespace flexaids {
namespace pdb_coords {

/// Parse three finite floats from PDB columns 31–54 (0-based 30, len 24).
/// Handles FlexAID compact negatives such as " -0.635 -80.275-146.614".
/// Returns false if fewer/more than three numbers or any non-finite value.
inline bool parse_xyz_span(const std::string& line, std::array<float, 3>& xyz) {
    if (line.size() < 54) return false;
    const std::string span = line.substr(30, 24);
    const char* p = span.c_str();
    const char* end = p + span.size();
    float vals[3] = {0.f, 0.f, 0.f};
    int n = 0;
    while (p < end && n < 3) {
        while (p < end && std::isspace(static_cast<unsigned char>(*p))) ++p;
        if (p >= end) break;
        char* next = nullptr;
        errno = 0;
        const float v = std::strtof(p, &next);
        if (next == p || errno == ERANGE || !std::isfinite(v)) return false;
        vals[n++] = v;
        p = next;
    }
    if (n != 3) return false;
    // Trailing non-space junk in the span is a parse failure (fail-closed).
    while (p < end) {
        if (!std::isspace(static_cast<unsigned char>(*p))) return false;
        ++p;
    }
    xyz = {vals[0], vals[1], vals[2]};
    return true;
}

/// Fixed 8-char field must consume the whole trimmed field as one finite float.
inline bool parse_fixed_float(const std::string& field, float& out) {
    std::string s = field;
    // trim
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())))
        s.erase(s.begin());
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())))
        s.pop_back();
    if (s.empty()) return false;
    char* end = nullptr;
    errno = 0;
    out = std::strtof(s.c_str(), &end);
    if (end == s.c_str() || errno == ERANGE || !std::isfinite(out)) return false;
    while (*end) {
        if (!std::isspace(static_cast<unsigned char>(*end))) return false;
        ++end;
    }
    return true;
}

}  // namespace pdb_coords
}  // namespace flexaids
