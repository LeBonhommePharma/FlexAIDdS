// RotatableBonds.h — Identify rotatable bonds for FlexAID dihedral gene encoding
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// A bond is rotatable if ALL of:
//   1. It is a single bond (BondOrder::SINGLE)
//   2. It is NOT in a ring (Bond::in_ring == false)
//   3. Both endpoint atoms have degree >= 2 (non-terminal)
//   4. It is NOT an amide bond (C(=O)-N)
//   5. It is NOT an S-S disulfide bond
//   6. It is NOT adjacent to a triple bond (propargylic bonds are often locked)
//   7. It is NOT a C-N bond where C is aromatic (aniline-like bonds treated as partial double)
//   8. It is NOT a conjugated C-N bond (urea, guanidine, vinylogous amide, N.am typed)
//
// Sets Bond::is_rotatable for all bonds in mol.
// Returns a list of bond indices that are rotatable.

#pragma once

#include "BonMol.h"
#include <vector>

namespace bonmol {
namespace rotatable_bonds {

struct RotatableBondsResult {
    std::vector<int> bond_indices; // indices into mol.bonds[]
    int count;
};

/// Identify rotatable bonds and set Bond::is_rotatable.
/// Requires ring perception to have been run (Bond::in_ring must be set).
RotatableBondsResult identify_rotatable_bonds(BonMol& mol);

/// Check whether bond at index bidx is an amide bond (C(=O)-N pattern).
bool is_amide_bond(const BonMol& mol, int bidx);

/// Check whether bond at index bidx is a disulfide bond (S-S).
bool is_disulfide_bond(const BonMol& mol, int bidx);

/// Conjugated / partial-double C–N (amide, urea, guanidine, vinylogous amide,
/// aromatic C–N, or MOL2/SYBYL N.am typed nitrogen). These must never be GA rotors.
bool is_conjugated_cn_bond(const BonMol& mol, int bidx);

/// True if either endpoint of the bond is adjacent to a triple bond.
bool is_triple_adjacent_bond(const BonMol& mol, int bidx);

} // namespace rotatable_bonds
} // namespace bonmol
