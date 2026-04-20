// Aromaticity.cpp — Hückel aromaticity detection and Kekulé assignment
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0

#include "Aromaticity.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <unordered_set>

namespace bonmol {
namespace aromaticity {

// ---------------------------------------------------------------------------
// Hybridisation assignment
// ---------------------------------------------------------------------------

void assign_hybridisation(BonMol& mol) {
    for (int i = 0; i < mol.num_atoms(); ++i) {
        Atom& a = mol.atoms[i];
        if (a.hybrid != Hybridization::UNSET) continue; // already set (e.g. from MOL2)

        // For H: SP
        if (a.element == Element::H) { a.hybrid = Hybridization::SP; continue; }

        // Count total connections (explicit bonds + implicit H)
        [[maybe_unused]] int total_connections = mol.degree(i) + a.implicit_h_count;

        // Count double and triple bonds
        bool has_triple = false;
        int  double_bonds = 0;
        for (int bidx : mol.bond_adj[i]) {
            const Bond& b = mol.bonds[bidx];
            if (b.order == BondOrder::TRIPLE)  has_triple = true;
            if (b.order == BondOrder::DOUBLE)  ++double_bonds;
        }

        if (has_triple || (double_bonds >= 2))
            a.hybrid = Hybridization::SP;
        else if (double_bonds == 1 || a.is_aromatic)
            a.hybrid = Hybridization::SP2;
        else
            a.hybrid = Hybridization::SP3;
    }
}

// ---------------------------------------------------------------------------
// π electron count per atom in ring context
// ---------------------------------------------------------------------------

int pi_electron_count(const BonMol& mol, int atom_idx, const Ring& ring) {
    const Atom& a = mol.atoms[atom_idx];

    // Build set of ring atoms for quick membership test
    std::unordered_set<int> ring_set(ring.atom_indices.begin(), ring.atom_indices.end());

    switch (a.element) {
        case Element::C: {
            // sp2 C with a double bond to a ring neighbour → 1 π electron
            // sp2 C with exocyclic double bond → 1 π electron (e.g. fulvene)
            if (a.hybrid == Hybridization::SP2 || a.is_aromatic) return 1;
            return -1; // sp3 C breaks aromaticity
        }
        case Element::N: {
            // Pyridine-like N: sp2, lone pair in plane, 1 π electron in ring
            // Pyrrole-like N: sp3 (but in ring), lone pair into ring, 2 π electrons
            if (a.hybrid == Hybridization::SP2) return 1;
            if (a.hybrid == Hybridization::SP3) {
                // Check: has H or substituent; lone pair goes into ring
                // (pyrrole, indole NH pattern)
                return 2;
            }
            return -1;
        }
        case Element::O: {
            // Furan-like O: lone pair into ring, 2 π electrons
            if (a.hybrid == Hybridization::SP2 || a.hybrid == Hybridization::SP3)
                return 2;
            return -1;
        }
        case Element::S: {
            // Thiophene-like S: lone pair into ring, 2 π electrons
            return 2;
        }
        case Element::B: {
            // Borole-like B: empty p orbital, 0 π electrons (anti-aromatic contributor)
            if (a.hybrid == Hybridization::SP2) return 0;
            return -1;
        }
        case Element::P: {
            // Phosphole-like P: 2 π electrons from lone pair
            return 2;
        }
        default:
            return -1; // element cannot be part of aromatic ring
    }
}

// ---------------------------------------------------------------------------
// Hückel 4n+2 check for a single ring
// ---------------------------------------------------------------------------

static bool is_huckel_aromatic(BonMol& mol, const Ring& ring) {
    // First pass: assign hybridisation for ring atoms if not yet done
    for (int ai : ring.atom_indices) {
        if (mol.atoms[ai].hybrid == Hybridization::UNSET)
            assign_hybridisation(mol);
    }

    int pi_total = 0;
    for (int ai : ring.atom_indices) {
        int pi = pi_electron_count(mol, ai, ring);
        if (pi < 0) return false; // atom cannot be in aromatic ring
        pi_total += pi;
    }

    // 4n+2 rule: 2, 6, 10, 14, ...
    if (pi_total < 2) return false;
    int n = (pi_total - 2) % 4;
    return n == 0;
}

// ---------------------------------------------------------------------------
// Kekulé assignment for a single ring via backtracking matching
// ---------------------------------------------------------------------------

bool kekulize_ring(BonMol& mol, const Ring& ring) {
    // We need to assign alternating single/double bonds around the ring.
    // Atoms that donate 2 π electrons (O, S, pyrrole-N) must have SINGLE bonds
    // to both neighbours in the ring (their lone pair fills the π system).
    // Atoms that donate 1 π electron (C, pyridine-N) participate in one double bond.

    int sz = ring.size;
    std::vector<int>& ai = const_cast<Ring&>(ring).atom_indices; // local alias

    // Build list of ring bond indices
    std::vector<int> rbonds(sz);
    for (int i = 0; i < sz; ++i) {
        int a = ai[i];
        int b = ai[(i + 1) % sz];
        int bidx = mol.find_bond(a, b);
        if (bidx < 0) return false; // missing bond — cannot Kekulize
        rbonds[i] = bidx;
    }

    // Try two starting patterns (offset 0 and offset 1) for even-sized rings
    for (int start = 0; start < 2; ++start) {
        bool valid = true;
        for (int i = 0; i < sz; ++i) {
            // Alternate: bond i gets DOUBLE if (i + start) is even
            BondOrder new_order = ((i + start) % 2 == 0) ? BondOrder::DOUBLE
                                                           : BondOrder::SINGLE;
            mol.bonds[rbonds[i]].order       = new_order;
            mol.bonds[rbonds[i]].is_aromatic = true; // retain aromatic flag
        }

        // Validate: check each ring atom's total bond order doesn't exceed valence
        for (int aidx : ai) {
            float bos = mol.bond_order_sum(aidx);
            const Atom& a = mol.atoms[aidx];
            int expected = 0;
            switch (a.element) {
                case Element::C:  expected = 4; break;
                case Element::N:  expected = (a.formal_charge == 1) ? 4 : 3; break;
                case Element::O:  expected = 2; break;
                case Element::S:  expected = 2; break;
                case Element::B:  expected = 3; break;
                case Element::P:  expected = 3; break;
                default:          expected = 4; break;
            }
            if (bos > expected + 0.5f) { valid = false; break; }
        }
        if (valid) return true;
    }

    // Kekulization failed — revert to aromatic bond orders
    for (int bidx : rbonds) {
        mol.bonds[bidx].order       = BondOrder::AROMATIC;
        mol.bonds[bidx].is_aromatic = true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

AromaticityResult assign_aromaticity(BonMol& mol) {
    AromaticityResult result{};

    // Ensure hybridisation is assigned for all atoms
    assign_hybridisation(mol);

    // Mark atoms listed as aromatic by the SMILES parser
    // (they already have is_aromatic = true from parsing)

    // For each ring, test Hückel aromaticity
    for (Ring& ring : mol.rings) {
        if (ring.size < 4 || ring.size > 8) {
            ring.is_aromatic = false;
            continue;
        }

        if (is_huckel_aromatic(mol, ring)) {
            ring.is_aromatic = true;
            ++result.num_aromatic_rings;

            // Mark ring atoms and bonds as aromatic
            int sz = ring.size;
            for (int i = 0; i < sz; ++i) {
                int ai = ring.atom_indices[i];
                int aj = ring.atom_indices[(i + 1) % sz];

                if (!mol.atoms[ai].is_aromatic) {
                    mol.atoms[ai].is_aromatic = true;
                    mol.atoms[ai].hybrid      = Hybridization::SP2;
                    ++result.num_aromatic_atoms;
                }

                int bidx = mol.find_bond(ai, aj);
                if (bidx >= 0 && !mol.bonds[bidx].is_aromatic) {
                    mol.bonds[bidx].is_aromatic = true;
                    mol.bonds[bidx].order       = BondOrder::AROMATIC;
                    ++result.num_aromatic_bonds;
                }
            }
        } else {
            ring.is_aromatic = false;
        }
    }

    // Kekulize aromatic rings
    result.kekulized = true;
    for (const Ring& ring : mol.rings) {
        if (!ring.is_aromatic) continue;
        if (!kekulize_ring(mol, ring)) {
            result.kekulized = false;
            // Leave as AROMATIC bond order — downstream SYBYL typer handles it
        }
    }

    // Propagate aromaticity: atoms that are aromatic but not yet sp2
    for (int i = 0; i < mol.num_atoms(); ++i) {
        if (mol.atoms[i].is_aromatic &&
            mol.atoms[i].hybrid != Hybridization::SP2) {
            mol.atoms[i].hybrid = Hybridization::SP2;
        }
    }

    // Update aromatic atom count (may have been set via SMILES parsing)
    result.num_aromatic_atoms = static_cast<int>(
        std::count_if(mol.atoms.begin(), mol.atoms.end(),
                      [](const Atom& a){ return a.is_aromatic; }));
    result.num_aromatic_bonds = static_cast<int>(
        std::count_if(mol.bonds.begin(), mol.bonds.end(),
                      [](const Bond& b){ return b.is_aromatic; }));

    return result;
}

} // namespace aromaticity
} // namespace bonmol
