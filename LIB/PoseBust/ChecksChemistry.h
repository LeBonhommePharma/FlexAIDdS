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
///   - "passes_rdkit_sanity_checks"
///       Finite coords, known elements, valid bond indices, no NaN/Inf.
///       (Named for interoperability; no RDKit dependency.)
///   - "inchi_convertible"
///       Soft placeholder without an InChI library: true iff heavy-atom count > 0
///       and every atom has a known element. Documented as soft / non-strict.
///   - "all_atoms_connected"
///       Single connected component on the heavy-atom graph (bonds only).
///   - "no_radicals"
///       Valence heuristic via bond-order sum (aromatic MDL order 4 → 1.5).
///       Bases: C 4, N 3, O 2, H 1, S 2/6, P 3/5, F/Cl/Br/I 1;
///       formal-charge slack of ±1 (Types.h has no formal_charge field → treated as 0;
///       residual ±1 tolerance retained for charged/aromatic edge cases).
void check_chemistry_sanity(const Molecule& pred, std::vector<CheckItem>& out);

/// Identity vs. crystal reference (skipped entirely when crystal is null).
/// Appends CheckItems when crystal is non-null:
///   - "formula"      — heavy-atom element multiset equality
///   - "connections"  — total bond count within 20% of the crystal
void check_identity_formula(const Molecule& pred,
                            const Molecule* crystal,
                            std::vector<CheckItem>& out);

}  // namespace flexaids::posebust
