// SybylTyper.cpp — SYBYL atom-type assignment and 256-type encoding
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
//
// SYBYL type → FlexAID numeric type mapping mirrors Mol2Reader.cpp.
// 256-type encoding mirrors atom_typing_256.h encode_from_sybyl().

#include "SybylTyper.h"
#include "../atom_typing_256.h"

#include <cmath>
#include <algorithm>

namespace bonmol {
namespace sybyl {

// ---------------------------------------------------------------------------
// SYBYL type names (for display)
// ---------------------------------------------------------------------------

const char* sybyl_type_name(int sybyl_type) {
    switch (sybyl_type) {
        case  0: return "C.1";
        case  1: return "C.3";
        case  2: return "C.2";
        case  3: return "C.ar";
        case  4: return "N.3";
        case  5: return "N.2";
        case  6: return "N.ar";
        case  7: return "N.am";
        case  8: return "N.pl3";
        case  9: return "N.4";
        case 10: return "O.3";
        case 11: return "O.2";
        case 12: return "O.co2";
        case 13: return "F";
        case 14: return "Cl";
        case 15: return "Br";
        case 16: return "S.3";
        case 17: return "S.2";
        case 18: return "S.O";
        case 19: return "S.O2";
        case 20: return "P.3";
        case 21: return "I";
        case 22: return "H";
        case 30: return "Fe";
        default: return "X";
    }
}

// ---------------------------------------------------------------------------
// Helper: check if atom is in an aromatic ring
// ---------------------------------------------------------------------------

static bool in_aromatic_ring(const BonMol& mol, int atom_idx) {
    for (const Ring& r : mol.rings) {
        if (!r.is_aromatic) continue;
        for (int ai : r.atom_indices)
            if (ai == atom_idx) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Helper: count double bonds to O (carbonyl/sulfonyl checks)
// ---------------------------------------------------------------------------

static int count_double_bonds_to(const BonMol& mol, int atom_idx, Element target) {
    int cnt = 0;
    for (int bidx : mol.bond_adj[atom_idx]) {
        const Bond& b = mol.bonds[bidx];
        if (b.order == BondOrder::DOUBLE) {
            int nb = (b.atom_i == atom_idx) ? b.atom_j : b.atom_i;
            if (mol.atoms[nb].element == target) ++cnt;
        }
    }
    return cnt;
}

// ---------------------------------------------------------------------------
// Helper: check if atom is part of a carboxylate/carboxamide group
// ---------------------------------------------------------------------------

static bool is_carboxylate_oxygen(const BonMol& mol, int atom_idx) {
    // O.co2: oxygen in –COO– (carboxylate, ester, or carbamate)
    if (mol.atoms[atom_idx].element != Element::O) return false;
    // The O should be connected to a C that also has a =O or another O
    for (int nb : mol.adjacency[atom_idx]) {
        if (mol.atoms[nb].element != Element::C) continue;
        // Count oxygens attached to this C
        int o_count = 0;
        bool has_double_o = false;
        for (int nb2 : mol.adjacency[nb]) {
            if (mol.atoms[nb2].element == Element::O) {
                ++o_count;
                int bidx = mol.find_bond(nb, nb2);
                if (bidx >= 0 && mol.bonds[bidx].order == BondOrder::DOUBLE)
                    has_double_o = true;
            }
        }
        if (o_count >= 2 && has_double_o) return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// H-bond donor detection
// ---------------------------------------------------------------------------

bool is_hbond_donor(const BonMol& mol, int atom_idx) {
    const Atom& a = mol.atoms[atom_idx];
    // N or O with implicit or explicit H
    if (a.element == Element::N || a.element == Element::O) {
        if (a.implicit_h_count > 0) return true;
        // Check for explicit H neighbours
        for (int nb : mol.adjacency[atom_idx]) {
            if (mol.atoms[nb].element == Element::H) return true;
        }
    }
    return false;
}

// ---------------------------------------------------------------------------
// H-bond acceptor detection
// ---------------------------------------------------------------------------

bool is_hbond_acceptor(const BonMol& mol, int atom_idx) {
    const Atom& a = mol.atoms[atom_idx];
    // N or O with available lone pair. Donor-only amide/pyrrole-like N and
    // quaternary N are excluded from acceptor roles.
    if (a.element == Element::O) return true;
    if (a.element == Element::N) {
        bool has_h = a.implicit_h_count > 0;
        for (int nb : mol.adjacency[atom_idx]) {
            if (mol.atoms[nb].element == Element::H) {
                has_h = true;
                break;
            }
        }
        if (a.sybyl_type == 9) return false; // N.4
        if (a.sybyl_type == 7) return false; // N.am
        if (a.sybyl_type == 6 && has_h) return false; // pyrrole-like N.ar
        if (a.partial_charge >= 0.3f) return false;
        return true;
    }
    if (a.element == Element::F) return true;
    return false;
}

// ---------------------------------------------------------------------------
// 256-type encoding (mirrors atom_typing_256.h encode_from_sybyl)
// Bits 0-5: base type (6 bits, 64 classes — no Solvent fallback)
// Bit    6: H-bond acceptor role
// Bit    7: H-bond donor role
// ---------------------------------------------------------------------------

// Internal: SYBYL type → base type (6-bit, 0-63)
static uint8_t sybyl_to_base(int sybyl_type) {
    // ProcessLigand uses the Mol2Reader-era numeric SYBYL IDs shown in
    // sybyl_type_name(); map them onto atom_typing_256.h base classes.
    switch (sybyl_type) {
        case  0: return atom256::C_sp;       // C.1
        case  1: return atom256::C_sp3;      // C.3
        case  2: return atom256::C_sp2;      // C.2
        case  3: return atom256::C_ar;       // C.ar
        case  4: return atom256::N_sp3;      // N.3
        case  5: return atom256::N_sp2;      // N.2
        case  6: return atom256::N_ar;       // N.ar
        case  7: return atom256::N_am;       // N.am
        case  8: return atom256::N_pl3;      // N.pl3
        case  9: return atom256::N_quat;     // N.4
        case 10: return atom256::O_sp3;      // O.3
        case 11: return atom256::O_sp2;      // O.2
        case 12: return atom256::O_co2;      // O.co2
        case 13: return atom256::HAL_F;      // F
        case 14: return atom256::HAL_Cl;     // Cl
        case 15: return atom256::HAL_Br;     // Br
        case 16: return atom256::S_sp3;      // S.3
        case 17: return atom256::S_sp2;      // S.2
        case 18: return atom256::S_oxide;    // S.O
        case 19: return atom256::S_dioxide;  // S.O2
        case 20: return atom256::P_sp3;      // P.3
        case 21: return atom256::HAL_I;      // I
        case 30: return atom256::Metal_Fe;   // Fe
        default: return atom256::Dummy;
    }
}

uint8_t encode_256(int sybyl_type, float partial_charge, bool is_donor,
                   bool is_acceptor) {
    (void)partial_charge;
    uint8_t base = sybyl_to_base(sybyl_type) & 0x3F; // bits 0-5
    return atom256::encode_roles(base, is_donor, is_acceptor);
}

// ---------------------------------------------------------------------------
// Single-atom SYBYL type assignment
// ---------------------------------------------------------------------------

int assign_sybyl_type_single(const BonMol& mol, int atom_idx) {
    const Atom& a = mol.atoms[atom_idx];

    switch (a.element) {

        // ---- Hydrogen ----
        case Element::H:
            return 22;

        // ---- Iron ----
        case Element::Fe:
            return 30;

        // ---- Halogens ----
        case Element::F:  return 13;
        case Element::Cl: return 14;
        case Element::Br: return 15;
        case Element::I:  return 21;

        // ---- Carbon ----
        case Element::C: {
            if (in_aromatic_ring(mol, atom_idx) || a.is_aromatic) return 3; // C.ar
            switch (a.hybrid) {
                case Hybridization::SP:  return 0;  // C.1
                case Hybridization::SP2: return 2;  // C.2
                case Hybridization::SP3: return 1;  // C.3
                default:                 return 1;  // default C.3
            }
        }

        // ---- Nitrogen ----
        case Element::N: {
            if (in_aromatic_ring(mol, atom_idx) || a.is_aromatic) return 6; // N.ar

            // N.4: quaternary N (total bond order == 4)
            float bos = mol.bond_order_sum(atom_idx) + a.implicit_h_count;
            if (bos >= 3.9f && a.formal_charge >= 1) return 9; // N.4

            // N.am: amide N — bonded to C=O
            for (int nb : mol.adjacency[atom_idx]) {
                if (mol.atoms[nb].element == Element::C) {
                    if (count_double_bonds_to(mol, nb, Element::O) >= 1) return 7; // N.am
                }
            }

            // N.pl3: planar N not in ring, not amide (e.g. guanidinium, urea)
            // Detected by sp2 hybridisation with no double bonds in molecule
            if (a.hybrid == Hybridization::SP2) {
                // Check for double bond to N itself
                for (int bidx : mol.bond_adj[atom_idx]) {
                    if (mol.bonds[bidx].order == BondOrder::DOUBLE) return 5; // N.2
                }
                return 8; // N.pl3
            }

            // N.2: sp2 N with double bond
            if (a.hybrid == Hybridization::SP) return 5; // treat sp N as N.2

            return 4; // N.3 (default sp3 N)
        }

        // ---- Oxygen ----
        case Element::O: {
            if (is_carboxylate_oxygen(mol, atom_idx)) return 12; // O.co2
            if (a.hybrid == Hybridization::SP2 ||
                count_double_bonds_to(mol, atom_idx, Element::C) >= 1 ||
                count_double_bonds_to(mol, atom_idx, Element::S) >= 1 ||
                count_double_bonds_to(mol, atom_idx, Element::P) >= 1 ||
                count_double_bonds_to(mol, atom_idx, Element::N) >= 1) {
                return 11; // O.2
            }
            return 10; // O.3
        }

        // ---- Sulphur ----
        case Element::S: {
            if (in_aromatic_ring(mol, atom_idx) || a.is_aromatic) return 16; // S.3 (thiophene S)
            int dbo = count_double_bonds_to(mol, atom_idx, Element::O);
            if (dbo >= 2) return 19; // S.O2
            if (dbo == 1) return 18; // S.O
            if (a.hybrid == Hybridization::SP2) return 17; // S.2
            return 16; // S.3
        }

        // ---- Phosphorus ----
        case Element::P:
            return 20; // P.3

        default:
            return 1; // fallback C.3
    }
}

// ---------------------------------------------------------------------------
// Main entry point: assign types for all atoms
// ---------------------------------------------------------------------------

void assign_sybyl_types(BonMol& mol) {
    for (int i = 0; i < mol.num_atoms(); ++i) {
        Atom& a = mol.atoms[i];

        // Assign SYBYL type
        a.sybyl_type = assign_sybyl_type_single(mol, i);

        // Assign H-bond flags
        a.is_hbond_donor    = is_hbond_donor(mol, i);
        a.is_hbond_acceptor = is_hbond_acceptor(mol, i);

        // Encode 256-type
        a.type_256 = encode_256(a.sybyl_type, a.partial_charge,
                                a.is_hbond_donor, a.is_hbond_acceptor);
    }
}

} // namespace sybyl
} // namespace bonmol
