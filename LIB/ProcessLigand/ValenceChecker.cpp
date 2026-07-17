// ValenceChecker.cpp — Valence validation and implicit H computation
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "ValenceChecker.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace bonmol {
namespace valence {

// ---------------------------------------------------------------------------
// Expected valences per element and formal charge
// ---------------------------------------------------------------------------

std::vector<int> expected_valences(Element elem, int formal_charge) {
    // Normal valence = standard_valence - formal_charge (for cations) or
    // standard_valence + |formal_charge| (for anions) — depends on element.
    // We store the expected BOND ORDER SUM.
    //
    // Simplified model: for each element, list acceptable bond-order sums
    // after accounting for the most common formal charges.

    switch (elem) {
        case Element::H:
            return {1};
        case Element::B:
            if (formal_charge == -1) return {4};
            return {3};
        case Element::C:
            if (formal_charge ==  1) return {3}; // carbocation
            if (formal_charge == -1) return {3}; // carbanion
            return {4};
        case Element::N:
            if (formal_charge ==  1) return {4}; // ammonium
            if (formal_charge == -1) return {2}; // amide anion
            return {3, 5};
        case Element::O:
            if (formal_charge ==  1) return {3}; // oxonium
            if (formal_charge == -1) return {1}; // oxide/phenolate
            return {2, 3}; // 3 = aromatic O donation into 5-membered ring (furan, oxazole)
        case Element::F:
            return {1};
        case Element::Si:
            return {4};
        case Element::P:
            if (formal_charge ==  1) return {4};
            return {3, 5};
        case Element::S:
            if (formal_charge ==  1) return {3};
            if (formal_charge == -1) return {1};
            return {2, 3, 4, 6}; // 3 = aromatic S (thiophene, thiadiazole)
        case Element::Cl:
            if (formal_charge == -1) return {0};
            return {1, 3, 5, 7};
        case Element::Se:
            return {2, 4, 6};
        case Element::Br:
            if (formal_charge == -1) return {0};
            return {1, 3, 5};
        case Element::I:
            if (formal_charge == -1) return {0};
            return {1, 3, 5, 7};
        case Element::Fe:
        case Element::Zn:
        case Element::Cu:
        case Element::Ni:
        case Element::Ca:
        case Element::Mg:
        case Element::Na:
        case Element::K:
            // Metals: accept any valence 0-6
            return {0, 1, 2, 3, 4, 5, 6};
        default:
            return {0, 1, 2, 3, 4, 5, 6}; // unknown: permissive
    }
}

// ---------------------------------------------------------------------------
// Implicit H computation (for SDF/MOL2 atoms without explicit H)
// ---------------------------------------------------------------------------

int compute_implicit_h(const BonMol& mol, int atom_idx) {
    const Atom& a = mol.atoms[atom_idx];

    // Elements that don't get implicit H
    switch (a.element) {
        case Element::Fe: case Element::Zn: case Element::Cu:
        case Element::Ni: case Element::Ca: case Element::Mg:
        case Element::Na: case Element::K:  case Element::Unknown:
            return 0;
        default: break;
    }

    // Use the same effective bond-order model as validation. Aromatic bonds
    // count as 1.5, so benzene-like carbons with two aromatic bonds receive
    // one implicit H instead of two.
    float bos = mol.bond_order_sum(atom_idx);

    // Find smallest valid valence >= current bond-order sum.
    auto valences = expected_valences(a.element, a.formal_charge);
    std::sort(valences.begin(), valences.end());

    int target = -1;
    for (int v : valences) {
        if (static_cast<float>(v) + 0.05f >= bos) { target = v; break; }
    }
    if (target < 0) return 0; // over-valenced already

    int h = static_cast<int>(std::round(static_cast<float>(target) - bos));
    return std::max(0, h);
}

// ---------------------------------------------------------------------------
// Main valence check
// ---------------------------------------------------------------------------

ValenceCheckResult check_valence(BonMol& mol) {
    ValenceCheckResult result;
    result.valid = true;

    for (int i = 0; i < mol.num_atoms(); ++i) {
        Atom& a = mol.atoms[i];

        // H-stripped crystallographic SDF/MOL2 inputs still carry normal
        // heavy-atom bond orders. Infer implicit H before checking valence so
        // neutral O/N/C centers are validated against their completed valence.
        if (a.implicit_h_count == 0 &&
            a.element != Element::H) {
            a.implicit_h_count = compute_implicit_h(mol, i);
        }

        // Compute bond order sum (aromatic = 1.5)
        float bos = mol.bond_order_sum(i);

        // Add implicit H contribution
        float total_bos = bos + static_cast<float>(a.implicit_h_count);

        int bos_int = static_cast<int>(std::round(total_bos));

        auto valid_vals = expected_valences(a.element, a.formal_charge);

        bool ok = std::any_of(valid_vals.begin(), valid_vals.end(),
                              [&](int v){ return v == bos_int; });

        // ±0.6 float tolerance to absorb small rounding noise in the bond-order
        // sum. NOTE: this tolerance is NOT what handles aromatic heteroatoms.
        // An aromatic S/O in a 5-membered ring has BOS = 2 × 1.5 = 3.0, which is
        // a full 1.0 away from {2} — well outside ±0.6. That case is handled by
        // listing 3 as a valid valence in expected_valences() (see Element::O/S),
        // not here. This tolerance only catches near-integer BOS values.
        if (!ok) {
            // Try exact float check
            for (int v : valid_vals) {
                if (std::fabs(total_bos - static_cast<float>(v)) < 0.6f) {
                    ok = true;
                    break;
                }
            }
        }
        // Aromatic C: kekule-form or 1.5-order ring edges can report BOS≈4–5
        // before/after perception; accept when atom is marked aromatic and the
        // heavy graph is 6π-plausible (BOS in [3.5, 5.5]).
        if (!ok && a.is_aromatic && a.element == Element::C &&
            total_bos >= 3.4f && total_bos <= 5.6f) {
            ok = true;
        }
        // Aromatic N in 6π rings (pyridine, etc.)
        if (!ok && a.is_aromatic && a.element == Element::N &&
            total_bos >= 2.4f && total_bos <= 4.6f) {
            ok = true;
        }

        if (!ok) {
            // Find the closest expected valence for error message
            int expected = valid_vals.empty() ? 0 : valid_vals[0];
            int min_dist = std::numeric_limits<int>::max();
            for (int v : valid_vals) {
                int dist = std::abs(v - bos_int);
                if (dist < min_dist) { min_dist = dist; expected = v; }
            }

            std::ostringstream msg;
            msg << "atom " << i << " (" 
                << static_cast<int>(a.element)
                << ") has bond-order sum " << total_bos
                << " (formal charge " << a.formal_charge
                << "); expected one of {";
            for (size_t k = 0; k < valid_vals.size(); ++k) {
                if (k) msg << ",";
                msg << valid_vals[k];
            }
            msg << "}";

            ValenceError err;
            err.atom_idx       = i;
            err.element        = a.element;
            err.bond_order_sum = total_bos;
            err.expected_valence = expected;
            err.message        = msg.str();

            // S with extended valence (sulfone, sulfoxide) is a warning, not error
            if (a.element == Element::S && total_bos <= 6.5f) {
                result.warnings.push_back(err);
            } else if (a.element == Element::P && total_bos <= 5.5f) {
                result.warnings.push_back(err);
            } else if (a.element == Element::N && total_bos <= 5.5f &&
                       a.formal_charge >= 0) {
                // N with 5 bonds (nitro group)
                result.warnings.push_back(err);
            } else {
                result.errors.push_back(err);
                result.valid = false;
            }
        }

    }

    return result;
}

} // namespace valence
} // namespace bonmol
