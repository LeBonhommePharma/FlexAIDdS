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
#include <string_view>
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
    // Aromaticity carried SEPARATELY from `order` so a ring can be written with
    // its Kekule orders (1/2 — sanitizable by RDKit) without losing the aromatic
    // marking that FlexAIDdS typing needs (C.ar/N.ar VCT types).
    //
    // Set from either encoding on load: legacy MDL order 4, or the
    // FLEXAIDDS_AROMATIC_BONDS SD tag written alongside Kekule orders. write_sdf
    // re-emits the tag, so a kekulized ligand survives a load/write round trip
    // with both its sanitizability AND its aromaticity intact.
    bool aromatic = false;
};

struct Molecule {
    std::string                   name;
    std::vector<Atom>             atoms;
    std::vector<Bond>             bonds;
    std::vector<std::vector<int>> adj;

    /// True when the bond block was FABRICATED rather than read.
    ///
    /// Three loader paths invent connectivity because their input carries none:
    /// infer_bonds() (covalent-radius proximity, order ALWAYS 1) and the PDB
    /// CONECT path (connectivity only, "order unknown"). An aromatic ring then
    /// arrives as a saturated all-single cage -- the HUP/1GPK failure mode the
    /// ligand extractor's own comment names.
    ///
    /// Six check sites read Bond::order (ChecksChemistry 161-164, 303, 649;
    /// ChecksGeometry 128-138, 322, 682), covering bond lengths, bond angles,
    /// aromatic-ring flatness, double-bond flatness and double-bond
    /// stereochemistry. Run against invented orders they return a verdict about
    /// chemistry the file never contained -- a false pass or a false failure,
    /// with nothing recording that the orders were guessed.
    ///
    /// So those checks are marked `skipped` instead: NOT ASSESSED, never a
    /// failure. That is the same rule the admission contract applies to a target
    /// PoseBusters cannot score, and `skipped` is already excluded from both
    /// n_passed() and n_failed(), so a not-assessed check can neither inflate
    /// nor deflate a rate. Connectivity- and distance-only checks are unaffected
    /// and still run.
    ///
    /// Cleared when real topology arrives: a parsed SDF bond block, or a
    /// successful assign_topology_from_reference().
    bool topology_inferred = false;

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
    /// Not computed / not applicable (missing InChI binary, no cofactor
    /// entities in an apo crop, …). Ignored by all_passed() / n_fail() /
    /// failed_keys_csv() so a skipped key cannot inflate a native pass or a
    /// native fail. JSON still emits the row.
    bool        skipped = false;
    // Optional continuous diagnostics (NaN if unused)
    float metric    = std::numeric_limits<float>::quiet_NaN();
    float threshold = std::numeric_limits<float>::quiet_NaN();
    int   n_checked = 0;
    int   n_failed  = 0;
};

// NativePoseQC diagnostic subset (extraction + protein clash/volume).
// NOT the authoritative PoseBusters gate — that is upstream `bust` via BustCli.
// Soft chemistry/geometry heuristics remain in `checks` for diagnostics only.
// Diagnostic subset keys use upstream PoseBusters names (True = pass).
inline constexpr const char* kNativeQcDiagnosticKeys[] = {
    "mol_pred_loaded",
    "mol_cond_loaded",
    "all_atoms_connected",
    "internal_steric_clash",
    "minimum_distance_to_protein",
    "protein-ligand_maximum_distance",
    "volume_overlap_with_protein",
};
/// Checks whose verdict is a function of Bond::order, and which therefore CANNOT
/// be assessed when Molecule::topology_inferred is true (orders fabricated).
///
/// Each entry traced to the site that reads the order:
///   bond_lengths                    ChecksGeometry:440  <- ideal length per order
///   bond_angles                     ChecksGeometry:500  <- ideal_angle_at_center 128-138
///   aromatic_ring_flatness          ring_has_mdl_aromatic_bond 316-327 (order == 4)
///   double_bond_flatness            ChecksGeometry:682  (order < 2 -> skip)
///   non-aromatic_ring_non-flatness  aromatic determination, same order test
///   double_bond_stereochemistry     ChecksChemistry:649 (order < 2 || == 4)
///   no_radicals                     ChecksChemistry:547 <- bond_order_sum 160-172
///   inchi_convertible               ChecksChemistry:445 <- gated on sanity.ok,
///                                   and an InChI derived from invented orders
///                                   describes a molecule the file never held
///
/// DELIBERATELY ABSENT, because they read connectivity or coordinates only and
/// stay assessable: all_atoms_connected, internal_steric_clash,
/// minimum_distance_to_protein, protein-ligand_maximum_distance,
/// volume_overlap_with_protein, molecular_formula, molecular_bonds, mol_*_loaded.
///
/// `sanitization` is ABSENT BY A DIFFERENT ROUTE and the reason is recorded here
/// so the omission is not silent. It reads Bond::order at ChecksChemistry:303
/// (bad_order), but a fabricated order is always 1 -- a VALID order -- so the
/// test cannot fail on inferred topology; it would merely report the bond block
/// as validated when none was read. That component is therefore neutralised at
/// the source line instead of skipping the whole check, which keeps
/// native_sanity's other components (non-finite coordinates, unknown elements,
/// bad atom indices) assessable. Skipping the key outright would have
/// suppressed real defects.
///
/// This list is verified to FIRE by mutation test, which is not the same as
/// verified COMPLETE: a check added later that reads Bond::order must be added
/// here too, or it will score fabricated chemistry silently.
inline constexpr const char* kOrderDependentKeys[] = {
    "bond_lengths",
    "bond_angles",
    "aromatic_ring_flatness",
    "double_bond_flatness",
    "non-aromatic_ring_non-flatness",
    "double_bond_stereochemistry",
    "no_radicals",
    "inchi_convertible",
};

struct PoseBustReport {
    std::vector<CheckItem> checks;
    bool                   ran     = false;
    std::string            backend;  // "native_pose_qc" | "bust_cli" | "skipped" | "error"
    std::string            error;    // hard execution failure only (not soft warnings)
    std::string            warning;  // soft diagnostics (e.g. topology assign miss)
    int                    n_ligand_atoms = 0;
    int                    n_protein_atoms_cropped = 0;

    [[nodiscard]] const CheckItem* find_check(std::string_view key) const {
        for (const CheckItem& c : checks)
            if (c.key == key) return &c;
        return nullptr;
    }

    [[nodiscard]] bool all_passed() const {
        if (!ran || !error.empty()) return false;
        int n_scored = 0;
        for (const CheckItem& c : checks) {
            if (c.skipped) continue;
            ++n_scored;
            if (!c.passed) return false;
        }
        return n_scored > 0;
    }

    /// Full NativePoseQC suite (diagnostic / parity target). Not claim gate.
    [[nodiscard]] bool success_pb_full() const { return all_passed(); }

    /// NativePoseQC diagnostic subset only (extract + clash/volume).
    /// Missing keys ⇒ fail closed. NEVER use as DatasetRunner.success_pb.
    [[nodiscard]] bool native_qc_diagnostic_pass() const {
        if (!ran || !error.empty() || checks.empty()) return false;
        for (const char* key : kNativeQcDiagnosticKeys) {
            const CheckItem* c = find_check(key);
            if (!c || !c->passed) return false;
        }
        return true;
    }

    /// @deprecated Alias of native_qc_diagnostic_pass — not the claim gate.
    [[nodiscard]] bool success_pb_campaign() const {
        return native_qc_diagnostic_pass();
    }

    /// @deprecated Do not map to DockingResult.success_pb (that is rmsd∧bust).
    [[nodiscard]] bool success_pb() const { return native_qc_diagnostic_pass(); }

    [[nodiscard]] int n_pass() const {
        int n = 0;
        for (const CheckItem& c : checks)
            if (!c.skipped && c.passed) ++n;
        return n;
    }
    [[nodiscard]] int n_fail() const {
        int n = 0;
        for (const CheckItem& c : checks)
            if (!c.skipped && !c.passed) ++n;
        return n;
    }
    [[nodiscard]] int n_skipped() const {
        int n = 0;
        for (const CheckItem& c : checks)
            if (c.skipped) ++n;
        return n;
    }
    [[nodiscard]] int n_checks() const { return static_cast<int>(checks.size()); }

    [[nodiscard]] std::string failed_keys_csv() const {
        std::string out;
        for (const CheckItem& c : checks) {
            if (c.skipped || c.passed) continue;
            if (!out.empty()) out += ';';
            out += c.key;
        }
        return out;
    }

    [[nodiscard]] std::string failed_campaign_keys_csv() const {
        std::string out;
        for (const char* key : kNativeQcDiagnosticKeys) {
            const CheckItem* c = find_check(key);
            if (c && c->passed) continue;
            if (!out.empty()) out += ';';
            out += key;
            if (!c) out += "(missing)";
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
// Native  = NativePoseQC clean-room diagnostic suite
// BustCli = official upstream PoseBusters CLI (default → pb_pass)
// Off     = skip
enum class Backend { BustCli, Native, Off };

}  // namespace flexaids::posebust
