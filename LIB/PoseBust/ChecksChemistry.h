// ChecksChemistry.h — Clean-room chemistry plausibility checks (PoseBusters-compatible names)
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Independent reimplementation of chemistry-level plausibility tests used by
// PoseBusters-style docking validation. Check *keys* match the public PoseBusters
// vocabulary for interoperability; algorithms and code are original (no RDKit /
// posebusters source was copied).
//
// Depends on Types.h (Molecule, CheckItem, Vec3, Atom, Bond).

#pragma once

#include "Types.h"

#include <vector>

namespace flexaids::posebust {

/// Loading checks.
/// Appends CheckItems:
///   - key "mol_pred_loaded"  — predicted ligand pointer non-null and has ≥1 atom
///   - key "mol_cond_loaded"  — protein / condition pointer non-null and has ≥1 atom
void check_loading(const Molecule* pred,
                   const Molecule* protein,
                   std::vector<CheckItem>& out);

/// Native chemistry sanity (PoseBusters key names; native algorithms).
/// Appends CheckItems:
///   - "sanitization"
///       Finite coords, known elements, valid bond indices, no NaN/Inf.
///       (Named for interoperability; no RDKit dependency.)
///   - "inchi_convertible"
///       Real conversion via system `inchi-1` (IUPAC InChI CLI) when available
///       (FLEXAIDDS_INCHI_BIN / PATH / Homebrew). Requires connected graph,
///       known elements, sanitization OK; fails closed on conversion error.
///   - "all_atoms_connected"
///       Single connected component on the heavy-atom graph (bonds only).
///   - "no_radicals"
///       Over-valence only (bos > max expected + slack). Under-valence is
///       allowed for heavy-only poses (implicit H). Aromatic MDL order 4 → 1.5.
void check_chemistry_sanity(const Molecule& pred, std::vector<CheckItem>& out);

/// Identity vs. crystal reference (skipped entirely when crystal is null).
/// Appends CheckItems when crystal is non-null:
///   - "mol_true_loaded"
///   - "molecular_formula"  — heavy-atom element multiset equality
///   - "molecular_bonds"    — total bond count within 20% of the crystal
void check_identity_formula(const Molecule& pred,
                            const Molecule* crystal,
                            std::vector<CheckItem>& out);

/// Stereochemistry + chirality (PoseBusters key names; clean-room geometry).
/// When crystal is non-null and atom counts match, compares local stereo signs.
/// Otherwise emits vacuous pass with detail (no stereo inventory without ref).
/// Keys:
///   - "double_bond_stereochemistry"
///   - "tetrahedral_chirality"
void check_stereochemistry(const Molecule& pred,
                           const Molecule* crystal,
                           std::vector<CheckItem>& out);

/// Soft internal energy ratio vs covalent ideal geometry.
/// Key: "internal_energy" — pass when mean relative bond strain ≤ 0.25
/// (or ≤ 2× crystal strain when crystal provided).
void check_internal_energy(const Molecule& pred,
                           const Molecule* crystal,
                           std::vector<CheckItem>& out);

}  // namespace flexaids::posebust
