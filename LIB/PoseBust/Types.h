// Types.h — Locked contract for native PoseBust (C++26)
//
// Clean-room PoseBusters-compatible pose validation for FlexAIDdS.
// Apache-2.0. No posebusters/RDKit source.
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace flexaids::posebust {

// ─── 3D vector ───────────────────────────────────────────────────────────────

struct Vec3 {
    float x = 0.f;
    float y = 0.f;
    float z = 0.f;
};

inline float dot(const Vec3& a, const Vec3& b) noexcept {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 cross(const Vec3& a, const Vec3& b) noexcept {
    return Vec3{a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x};
}
inline float norm(const Vec3& v) noexcept { return std::sqrt(dot(v, v)); }
inline float dist2(const Vec3& a, const Vec3& b) noexcept {
    const float dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}
inline float dist(const Vec3& a, const Vec3& b) noexcept { return std::sqrt(dist2(a, b)); }
inline Vec3 operator-(const Vec3& a, const Vec3& b) noexcept {
    return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 operator+(const Vec3& a, const Vec3& b) noexcept {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 operator*(const Vec3& a, float s) noexcept {
    return Vec3{a.x * s, a.y * s, a.z * s};
}

// ─── Molecule graph ──────────────────────────────────────────────────────────

struct Atom {
    int         id         = 0;   // 1-based external index when from SDF
    std::string element;          // "C", "Cl", …
    float       x          = 0.f;
    float       y          = 0.f;
    float       z          = 0.f;
    int         atomic_num = 0;   // Z
    bool        is_h       = false;

    // Convenience aliases used by check modules
    [[nodiscard]] Vec3 pos() const noexcept { return Vec3{x, y, z}; }
    void set_pos(const Vec3& p) noexcept { x = p.x; y = p.y; z = p.z; }
    [[nodiscard]] int Z() const noexcept { return atomic_num; }
};

// Bond order: 1 single, 2 double, 3 triple, 4 aromatic (MDL).
struct Bond {
    int a     = 0;  // 0-based
    int b     = 0;
    int order = 1;
};

struct Molecule {
    std::string                   name;
    std::vector<Atom>             atoms;
    std::vector<Bond>             bonds;
    std::vector<std::vector<int>> adj;

    void build_adjacency() {
        adj.assign(atoms.size(), {});
        for (const Bond& bond : bonds) {
            if (bond.a < 0 || bond.b < 0) continue;
            if (static_cast<std::size_t>(bond.a) >= atoms.size() ||
                static_cast<std::size_t>(bond.b) >= atoms.size())
                continue;
            adj[static_cast<std::size_t>(bond.a)].push_back(bond.b);
            adj[static_cast<std::size_t>(bond.b)].push_back(bond.a);
        }
    }

    [[nodiscard]] int n_heavy() const {
        int n = 0;
        for (const Atom& a : atoms)
            if (!a.is_h) ++n;
        return n;
    }

    [[nodiscard]] bool empty() const { return atoms.empty(); }
};

// ─── Validation report (locked contract) ─────────────────────────────────────

struct CheckItem {
    std::string key;     // stable machine id: "no_internal_clash"
    std::string label;   // human: "Internal steric clash"
    bool        passed = false;
    std::string detail;
    // Optional continuous diagnostics (NaN if unused)
    float metric    = std::numeric_limits<float>::quiet_NaN();
    float threshold = std::numeric_limits<float>::quiet_NaN();
    int   n_checked = 0;
    int   n_failed  = 0;
};

struct PoseBustReport {
    std::vector<CheckItem> checks;
    bool                   ran     = false;
    std::string            backend;  // "native" | "skipped" | "error"
    std::string            error;    // non-empty if ran failed to execute

    [[nodiscard]] bool all_passed() const {
        if (!ran || !error.empty()) return false;
        for (const CheckItem& c : checks)
            if (!c.passed) return false;
        return !checks.empty();  // empty after run ⇒ not a pass
    }
    [[nodiscard]] bool success_pb() const { return all_passed(); }

    [[nodiscard]] int n_pass() const {
        int n = 0;
        for (const CheckItem& c : checks)
            if (c.passed) ++n;
        return n;
    }
    [[nodiscard]] int n_fail() const {
        int n = 0;
        for (const CheckItem& c : checks)
            if (!c.passed) ++n;
        return n;
    }
    [[nodiscard]] int n_checks() const { return static_cast<int>(checks.size()); }

    [[nodiscard]] std::string failed_keys_csv() const {
        std::string out;
        for (const CheckItem& c : checks) {
            if (c.passed) continue;
            if (!out.empty()) out += ';';
            out += c.key;
        }
        return out;
    }

    // Continuous summaries (filled by Engine when available)
    float min_lig_prot_dist = std::numeric_limits<float>::quiet_NaN();
    float volume_overlap    = std::numeric_limits<float>::quiet_NaN();
};

// Suite selection (dock = protein-conditioned; redock adds identity; mol = ligand-only)
enum class Suite { Dock, Redock, Mol };

// Backend selection for DatasetRunner
enum class Backend { Native, Off };

}  // namespace flexaids::posebust
