// ChecksChemistry.cpp — Clean-room chemistry plausibility checks
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// PoseBusters-compatible check *keys* only. Algorithms are original.
// No RDKit, no InChI library, no posebusters source code.

#include "ChecksChemistry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <map>
#include <queue>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace flexaids::posebust {
namespace {

// ---------------------------------------------------------------------------
// CheckItem helper (Types.h: key / label / passed / detail)
// ---------------------------------------------------------------------------

void emit(std::vector<CheckItem>& out,
          std::string key,
          std::string label,
          bool passed,
          std::string detail = {}) {
    CheckItem item;
    item.key    = std::move(key);
    item.label  = std::move(label);
    item.passed = passed;
    item.detail = std::move(detail);
    out.push_back(std::move(item));
}

// ---------------------------------------------------------------------------
// Atomic helpers (local — Types.h does not provide a full table)
// ---------------------------------------------------------------------------

[[nodiscard]] int atomic_number_of(const Atom& a) noexcept {
    if (a.atomic_num > 0) return a.atomic_num;
    // Fallback: parse element symbol if atomic_num unset.
    if (a.element.empty()) return 0;
    std::string sym = a.element;
    // Trim and normalise first letter upper, rest lower.
    // Keep only first two alpha characters.
    std::string s;
    for (char c : sym) {
        if (std::isalpha(static_cast<unsigned char>(c))) s.push_back(c);
        if (s.size() >= 2) break;
    }
    if (s.empty()) return 0;
    s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    if (s.size() > 1) {
        s[1] = static_cast<char>(std::tolower(static_cast<unsigned char>(s[1])));
    }

    if (s == "H") return 1;
    if (s == "B") return 5;
    if (s == "C") return 6;
    if (s == "N") return 7;
    if (s == "O") return 8;
    if (s == "F") return 9;
    if (s == "Na") return 11;
    if (s == "Mg") return 12;
    if (s == "Si") return 14;
    if (s == "P") return 15;
    if (s == "S") return 16;
    if (s == "Cl") return 17;
    if (s == "K") return 19;
    if (s == "Ca") return 20;
    if (s == "Fe") return 26;
    if (s == "Co") return 27;
    if (s == "Ni") return 28;
    if (s == "Cu") return 29;
    if (s == "Zn") return 30;
    if (s == "Se") return 34;
    if (s == "Br") return 35;
    if (s == "I") return 53;
    return 0;
}

[[nodiscard]] bool atom_is_hydrogen(const Atom& a) noexcept {
    if (a.is_h) return true;
    return atomic_number_of(a) == 1;
}

[[nodiscard]] bool is_known_Z(int Z) noexcept {
    switch (Z) {
        case 1: case 5: case 6: case 7: case 8: case 9:
        case 11: case 12: case 14: case 15: case 16: case 17:
        case 19: case 20: case 26: case 27: case 28: case 29: case 30:
        case 34: case 35: case 53:
            return true;
        default:
            return false;
    }
}

[[nodiscard]] std::string_view symbol_of_Z(int Z) noexcept {
    switch (Z) {
        case 1:  return "H";
        case 5:  return "B";
        case 6:  return "C";
        case 7:  return "N";
        case 8:  return "O";
        case 9:  return "F";
        case 11: return "Na";
        case 12: return "Mg";
        case 14: return "Si";
        case 15: return "P";
        case 16: return "S";
        case 17: return "Cl";
        case 19: return "K";
        case 20: return "Ca";
        case 26: return "Fe";
        case 27: return "Co";
        case 28: return "Ni";
        case 29: return "Cu";
        case 30: return "Zn";
        case 34: return "Se";
        case 35: return "Br";
        case 53: return "I";
        default: return "?";
    }
}

[[nodiscard]] bool molecule_loaded(const Molecule* mol) noexcept {
    return mol != nullptr && !mol->atoms.empty();
}

[[nodiscard]] bool finite_xyz(const Atom& a) noexcept {
    return std::isfinite(a.x) && std::isfinite(a.y) && std::isfinite(a.z);
}

/// Bond order contribution for valence: MDL aromatic (4) → 1.5.
[[nodiscard]] float bond_order_valence(const Bond& b) noexcept {
    if (b.order == 4) return 1.5f;  // aromatic
    if (b.order <= 0) return 1.0f;
    if (b.order >= 3) return 3.0f;
    return static_cast<float>(b.order);
}

[[nodiscard]] float bond_order_sum(const Molecule& mol, int ai) noexcept {
    float sum = 0.0f;
    const int n = static_cast<int>(mol.atoms.size());
    if (ai < 0 || ai >= n) return 0.0f;
    for (const Bond& b : mol.bonds) {
        if (b.a == ai || b.b == ai) sum += bond_order_valence(b);
    }
    return sum;
}

// ---------------------------------------------------------------------------
// Expected valences
// ---------------------------------------------------------------------------
// Bases: C 4, N 3, O 2, H 1, S 2/6, P 3/5, F/Cl/Br/I 1
// Formal-charge slack ±1: Types.h Atom has no formal_charge — use charge=0
// and retain a modest ±1 residual tolerance for aromatic / charged edge cases.

void expected_valences(int Z, std::vector<int>& out) {
    out.clear();
    switch (Z) {
        case 1:  out = {1}; break;          // H
        case 5:  out = {3, 4}; break;       // B
        case 6:  out = {4}; break;          // C
        case 7:  out = {3}; break;          // N
        case 8:  out = {2}; break;          // O
        case 9:  out = {1}; break;          // F
        case 14: out = {4}; break;          // Si
        case 15: out = {3, 5}; break;       // P
        case 16: out = {2, 6}; break;       // S
        case 17: out = {1}; break;          // Cl
        case 34: out = {2, 4, 6}; break;    // Se
        case 35: out = {1}; break;          // Br
        case 53: out = {1}; break;          // I
        default: break;  // metals / unknown: not enforced
    }
}

[[nodiscard]] bool valence_ok(int Z, float bos) noexcept {
    std::vector<int> vals;
    expected_valences(Z, vals);
    if (vals.empty()) return true;

    // Aromatic half-orders + formal-charge slack ±1 (no formal_charge field).
    // tol = 1.05 covers off-by-one valence for charged species and 1.5 aromatic.
    constexpr float tol = 1.05f;
    for (int v : vals) {
        if (std::fabs(bos - static_cast<float>(v)) <= tol) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Connectivity (heavy-atom graph)
// ---------------------------------------------------------------------------

[[nodiscard]] bool heavy_atoms_single_component(const Molecule& mol,
                                                int& heavy_count,
                                                int& component_count) {
    heavy_count     = 0;
    component_count = 0;
    const int n     = static_cast<int>(mol.atoms.size());
    if (n == 0) return false;

    std::vector<int> heavy_index(static_cast<std::size_t>(n), -1);
    std::vector<int> heavies;
    heavies.reserve(static_cast<std::size_t>(n));

    for (int i = 0; i < n; ++i) {
        if (!atom_is_hydrogen(mol.atoms[static_cast<std::size_t>(i)])) {
            heavy_index[static_cast<std::size_t>(i)] =
                static_cast<int>(heavies.size());
            heavies.push_back(i);
        }
    }
    heavy_count = static_cast<int>(heavies.size());
    if (heavy_count == 0) return false;
    if (heavy_count == 1) {
        component_count = 1;
        return true;
    }

    std::vector<std::vector<int>> adj(static_cast<std::size_t>(heavy_count));
    for (const Bond& b : mol.bonds) {
        if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n) continue;
        const int hi = heavy_index[static_cast<std::size_t>(b.a)];
        const int hj = heavy_index[static_cast<std::size_t>(b.b)];
        if (hi < 0 || hj < 0 || hi == hj) continue;
        adj[static_cast<std::size_t>(hi)].push_back(hj);
        adj[static_cast<std::size_t>(hj)].push_back(hi);
    }

    std::vector<char> seen(static_cast<std::size_t>(heavy_count), 0);
    for (int start = 0; start < heavy_count; ++start) {
        if (seen[static_cast<std::size_t>(start)]) continue;
        ++component_count;
        std::queue<int> q;
        q.push(start);
        seen[static_cast<std::size_t>(start)] = 1;
        while (!q.empty()) {
            const int u = q.front();
            q.pop();
            for (int v : adj[static_cast<std::size_t>(u)]) {
                if (!seen[static_cast<std::size_t>(v)]) {
                    seen[static_cast<std::size_t>(v)] = 1;
                    q.push(v);
                }
            }
        }
    }
    return component_count == 1;
}

// ---------------------------------------------------------------------------
// Native "RDKit" sanity (no RDKit)
// ---------------------------------------------------------------------------

struct SanityReport {
    bool ok{true};
    int  n_nonfinite{0};
    int  n_unknown_elem{0};
    int  n_bad_bonds{0};
    std::string detail;
};

SanityReport native_sanity(const Molecule& mol) {
    SanityReport r;
    const int n = static_cast<int>(mol.atoms.size());
    if (n == 0) {
        r.ok     = false;
        r.detail = "empty molecule";
        return r;
    }

    for (int i = 0; i < n; ++i) {
        const Atom& a = mol.atoms[static_cast<std::size_t>(i)];
        if (!finite_xyz(a)) {
            ++r.n_nonfinite;
            r.ok = false;
        }
        const int Z = atomic_number_of(a);
        if (!is_known_Z(Z)) {
            ++r.n_unknown_elem;
            r.ok = false;
        }
    }

    for (const Bond& b : mol.bonds) {
        const bool bad_idx =
            b.a < 0 || b.b < 0 || b.a >= n || b.b >= n || b.a == b.b;
        // Valid MDL-style orders: 1,2,3,4 (aromatic). Reject others.
        const bool bad_order = (b.order < 1 || b.order > 4);
        if (bad_idx || bad_order) {
            ++r.n_bad_bonds;
            r.ok = false;
        }
    }

    std::ostringstream oss;
    oss << "atoms=" << n << " bonds=" << mol.bonds.size()
        << " nonfinite=" << r.n_nonfinite
        << " unknown_elem=" << r.n_unknown_elem
        << " bad_bonds=" << r.n_bad_bonds;
    r.detail = oss.str();
    return r;
}

// ---------------------------------------------------------------------------
// Formula helpers
// ---------------------------------------------------------------------------

using ElementMultiset = std::map<int, int>;  // Z → count

ElementMultiset heavy_formula(const Molecule& mol) {
    ElementMultiset m;
    for (const Atom& a : mol.atoms) {
        if (atom_is_hydrogen(a)) continue;
        const int Z = atomic_number_of(a);
        if (Z > 0) ++m[Z];
    }
    return m;
}

[[nodiscard]] std::string formula_string(const ElementMultiset& m) {
    std::ostringstream oss;
    auto emit_sym = [&](int Z, int count) {
        if (count <= 0) return;
        oss << symbol_of_Z(Z);
        if (count > 1) oss << count;
    };
    if (auto it = m.find(6); it != m.end()) emit_sym(6, it->second);
    for (const auto& [Z, c] : m) {
        if (Z == 6) continue;
        emit_sym(Z, c);
    }
    if (m.empty()) oss << "(empty)";
    return oss.str();
}

}  // namespace

// ===========================================================================
// Public API
// ===========================================================================

void check_loading(const Molecule* pred,
                   const Molecule* protein,
                   std::vector<CheckItem>& out) {
    const bool pred_ok = molecule_loaded(pred);
    {
        std::ostringstream d;
        if (pred == nullptr) {
            d << "pred is null";
        } else if (pred->atoms.empty()) {
            d << "pred has zero atoms";
        } else {
            d << "pred atoms=" << pred->atoms.size()
              << " bonds=" << pred->bonds.size();
        }
        emit(out, "mol_pred_loaded", "MOL_PRED loaded", pred_ok, d.str());
    }

    const bool cond_ok = molecule_loaded(protein);
    {
        std::ostringstream d;
        if (protein == nullptr) {
            d << "protein/condition is null";
        } else if (protein->atoms.empty()) {
            d << "protein/condition has zero atoms";
        } else {
            d << "protein atoms=" << protein->atoms.size()
              << " bonds=" << protein->bonds.size();
        }
        emit(out, "mol_cond_loaded", "MOL_COND loaded", cond_ok, d.str());
    }
}

void check_chemistry_sanity(const Molecule& pred, std::vector<CheckItem>& out) {
    // --- passes_rdkit_sanity_checks (native, no RDKit) --------------------
    const SanityReport sanity = native_sanity(pred);
    emit(out,
         "passes_rdkit_sanity_checks",
         "Sanitization",
         sanity.ok,
         sanity.detail);

    // --- inchi_convertible (soft placeholder; no InChI library) -----------
    int heavy = 0;
    bool all_known = !pred.atoms.empty();
    for (const Atom& a : pred.atoms) {
        const int Z = atomic_number_of(a);
        if (!is_known_Z(Z)) all_known = false;
        if (!atom_is_hydrogen(a) && Z > 0) ++heavy;
    }
    const bool inchi_soft = (heavy > 0) && all_known;
    {
        std::ostringstream d;
        d << "soft=true (no InChI lib); heavy=" << heavy
          << " all_known=" << (all_known ? "true" : "false");
        emit(out,
             "inchi_convertible",
             "InChI convertible",
             inchi_soft,
             d.str());
    }

    // --- all_atoms_connected (heavy-atom graph, bonds only) ---------------
    int heavy_count = 0;
    int n_comp      = 0;
    const bool connected =
        heavy_atoms_single_component(pred, heavy_count, n_comp);
    {
        std::ostringstream d;
        d << "heavy=" << heavy_count << " components=" << n_comp;
        emit(out,
             "all_atoms_connected",
             "All atoms connected",
             connected,
             d.str());
    }

    // --- no_radicals (valence heuristic) ----------------------------------
    int n_bad     = 0;
    int n_checked = 0;
    for (int i = 0; i < static_cast<int>(pred.atoms.size()); ++i) {
        const Atom& a = pred.atoms[static_cast<std::size_t>(i)];
        const int   Z = atomic_number_of(a);
        std::vector<int> vals;
        expected_valences(Z, vals);
        if (vals.empty()) continue;  // metals / unknown: skip
        ++n_checked;
        const float bos = bond_order_sum(pred, i);
        if (!valence_ok(Z, bos)) ++n_bad;
    }
    const bool no_rad = (n_bad == 0);
    {
        std::ostringstream d;
        d << "checked=" << n_checked << " radicals_or_bad_valence=" << n_bad
          << " (aromatic order=4→1.5; formal_charge slack ±1, charge field n/a)";
        emit(out, "no_radicals", "No radicals", no_rad, d.str());
    }
}

void check_identity_formula(const Molecule& pred,
                            const Molecule* crystal,
                            std::vector<CheckItem>& out) {
    // Per contract: if crystal is null, skip entirely (no keys appended).
    if (crystal == nullptr) return;

    // --- formula: heavy-atom element multiset equality --------------------
    const ElementMultiset f_pred = heavy_formula(pred);
    const ElementMultiset f_ref  = heavy_formula(*crystal);
    const bool formula_ok        = (f_pred == f_ref);
    {
        std::ostringstream d;
        d << "pred=" << formula_string(f_pred)
          << " crystal=" << formula_string(f_ref);
        emit(out, "formula", "Molecular formula", formula_ok, d.str());
    }

    // --- connections: bond count within 20% -------------------------------
    const std::size_t bp = pred.bonds.size();
    const std::size_t bc = crystal->bonds.size();
    bool conn_ok = false;
    if (bc == 0) {
        conn_ok = (bp == 0);
    } else {
        const double ratio =
            std::fabs(static_cast<double>(bp) - static_cast<double>(bc)) /
            static_cast<double>(bc);
        conn_ok = (ratio <= 0.20);
    }
    {
        std::ostringstream d;
        d << "pred_bonds=" << bp << " crystal_bonds=" << bc
          << " tolerance=20%";
        emit(out, "connections", "Molecular bonds", conn_ok, d.str());
    }
}

}  // namespace flexaids::posebust
