// SybylTyper.h — SYBYL atom-type assignment and 256-type encoding
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// Assigns SYBYL atom types (as used by FlexAID's 40-type system) from the
// molecular graph. The mapping is compatible with Mol2Reader.cpp's existing
// SYBYL→FlexAID mapping:
//
//   C.3→1   C.2→2   C.ar→3  C.1→0(special)
//   N.3→4   N.2→5   N.ar→6  N.am→7  N.pl3→8  N.4→9
//   O.3→10  O.2→11  O.co2→12
//   F→13    Cl→14   Br→15
//   S.3→16  S.2→17  S.O→18  S.O2→19
//   P.3→20
//   I→21
//   H→22
//   Fe→30
//
// After SYBYL type assignment, the 256-type encoding from atom_typing_256.h
// is applied via encode_from_sybyl().
//
// Also assigns H-bond donor/acceptor flags used in the 256-type encoding.

#pragma once

#include "BonMol.h"

namespace bonmol {
namespace sybyl {

/// SYBYL type strings (for display/debug)
const char* sybyl_type_name(int sybyl_type);

/// Assign SYBYL type (Atom::sybyl_type) for all atoms in mol.
/// Also sets Atom::type_256, Atom::is_hbond_donor, Atom::is_hbond_acceptor.
/// Requires hybridisation and aromaticity to be set (call assign_aromaticity first).
void assign_sybyl_types(BonMol& mol);

/// Assign SYBYL type for a single atom. Returns the FlexAID numeric SYBYL type (1-30).
int assign_sybyl_type_single(const BonMol& mol, int atom_idx);

/// Determine if an atom is an H-bond donor (has polar H or lone-pair NH).
bool is_hbond_donor(const BonMol& mol, int atom_idx);

/// Determine if an atom is an H-bond acceptor (N/O with lone pair).
bool is_hbond_acceptor(const BonMol& mol, int atom_idx);

/// Encode FlexAID SYBYL type + H-bond roles into 8-bit type.
/// Replicates atom_typing_256.h encode_from_sybyl() for inline use.
uint8_t encode_256(int sybyl_type, float partial_charge, bool is_donor,
                   bool is_acceptor);

} // namespace sybyl
} // namespace bonmol
