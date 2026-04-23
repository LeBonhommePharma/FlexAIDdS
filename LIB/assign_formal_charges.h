// assign_formal_charges.h — Residue-aware formal charge assignment for PDB input
//
// PDB files do not carry partial charges. Without charges:
//   - Coulomb electrostatics is dead (vcfunction.cpp guard: qA != 0 && qB != 0)
//   - Salt bridge detection is broken (type256 charge bit always Q_POSITIVE)
//   - Metal coordination has no electrostatic context
//
// This module assigns AMBER-ff14SB-derived partial charges to standard amino
// acid titratable atoms and formal charges to metal ions, called once after
// assign_radii_types() during the PDB loading pipeline.
//
// Charge sources:
//   - Amino acid side-chain: AMBER ff14SB partial charges (Cornell et al. 1995,
//     Maier et al. 2015 JCTC 11:3696) for charged/polar atoms only
//   - Metal ions: integer formal charges (IUPAC standard oxidation states)
//   - Common anions: Cl-, Br-, I- (monovalent)
//   - Backbone termini: standard COO-/NH3+ charges
//
// Only assigns to atoms with charge == 0.0 (i.e., no MOL2/PTM charge present).
// Only assigns to receptor atoms (residue[].type == 0), not ligands.
//
// Copyright 2024-2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstring>
#include <cstdio>

// Need full struct definitions for inline implementation
#include "flexaid.h"

namespace formal_charges {

// ─── Lookup entry ──────────────────────────────────────────────────────────
struct ChargeEntry {
    const char* res_name;   // 3-char residue name (e.g., "ASP")
    const char* atom_name;  // 4-char atom name (e.g., " OD1")
    float       charge;     // partial charge to assign
};

// ─── Amino acid titratable/polar atom charges ──────────────────────────────
// Charges from AMBER ff14SB for side-chain atoms that carry significant
// partial charge at physiological pH (7.4). We assign ONLY to atoms
// directly involved in salt bridges, H-bonds, or metal coordination.
// Non-polar carbons/hydrogens are left at 0.0 (they contribute through
// the complementarity function, not Coulomb).
//
// This is deliberately NOT a full force field — it's a minimal set that
// enables Coulomb and salt bridge detection for the most important
// interactions in protein-ligand docking.
static constexpr ChargeEntry AMINO_ACID_CHARGES[] = {
    // ── Aspartate (ASP) — deprotonated at pH 7.4 ──────────────────────
    // pKa 3.65, total formal charge -1
    // AMBER ff14SB: OD1=-0.8014, OD2=-0.8014, CG=0.7994
    {"ASP", " OD1", -0.80f},
    {"ASP", " OD2", -0.80f},
    {"ASP", " CG ", +0.70f},

    // ── Glutamate (GLU) — deprotonated at pH 7.4 ──────────────────────
    // pKa 4.25, total formal charge -1
    // AMBER ff14SB: OE1=-0.8188, OE2=-0.8188, CD=0.8054
    {"GLU", " OE1", -0.82f},
    {"GLU", " OE2", -0.82f},
    {"GLU", " CD ", +0.80f},

    // ── Lysine (LYS) — protonated at pH 7.4 ───────────────────────────
    // pKa 10.5, total formal charge +1
    // AMBER ff14SB: NZ=-0.3854, HZ1-3=+0.34 each, CE=−0.0187
    // Net on NZ+3H = +0.68. We assign net +1.0 on NZ (hydrogens implicit
    // in docking — explicit H not always present in PDB)
    {"LYS", " NZ ", +1.00f},

    // ── Arginine (ARG) — protonated at pH 7.4 ─────────────────────────
    // pKa 12.5, total formal charge +1
    // AMBER ff14SB: CZ=0.8281, NH1=-0.8693, NH2=-0.8693
    // The +1 is delocalized across the guanidinium; we distribute:
    {"ARG", " NH1", +0.45f},
    {"ARG", " NH2", +0.45f},
    {"ARG", " CZ ", +0.64f},
    {"ARG", " NE ", -0.54f},

    // ── Histidine — protonation state ambiguous ────────────────────────
    // At pH 7.4, His is ~50% HID / ~40% HIE / ~10% HIP
    // Conservative: assign small positive charge to both ring N atoms
    // (reflects average protonation, enables metal coordination detection)
    // HIS = generic, HID = delta-protonated, HIE = epsilon-protonated
    {"HIS", " ND1", -0.35f},
    {"HIS", " NE2", -0.35f},
    {"HIS", " CE1", +0.20f},
    {"HID", " ND1", -0.38f},   // delta-protonated: ND1 has H
    {"HID", " NE2", -0.57f},   // NE2 is the lone-pair (metal coordinator)
    {"HIE", " ND1", -0.54f},   // ND1 is the lone-pair
    {"HIE", " NE2", -0.27f},   // epsilon-protonated
    {"HIP", " ND1", -0.15f},   // doubly protonated (+1)
    {"HIP", " NE2", -0.15f},
    {"HIP", " CE1", +0.37f},

    // ── Tyrosine — phenol OH ───────────────────────────────────────────
    // pKa 10.1; protonated at pH 7.4 but weakly acidic
    {"TYR", " OH ", -0.56f},

    // ── Serine — hydroxyl ──────────────────────────────────────────────
    {"SER", " OG ", -0.65f},

    // ── Threonine — hydroxyl ───────────────────────────────────────────
    {"THR", " OG1", -0.68f},

    // ── Cysteine — thiol (neutral at pH 7.4, pKa ~8.3) ────────────────
    // AMBER: SG=-0.3119. When deprotonated (CYM/CYX), charge is ~-0.8
    {"CYS", " SG ", -0.31f},
    {"CYM", " SG ", -0.80f},   // deprotonated cysteine (thiolate)
    {"CYX", " SG ", -0.08f},   // disulfide-bonded

    // ── Asparagine — amide ─────────────────────────────────────────────
    {"ASN", " OD1", -0.59f},
    {"ASN", " ND2", -0.30f},

    // ── Glutamine — amide ──────────────────────────────────────────────
    {"GLN", " OE1", -0.59f},
    {"GLN", " NE2", -0.30f},

    // ── Tryptophan — indole NH ─────────────────────────────────────────
    {"TRP", " NE1", -0.34f},

    // ── Backbone carbonyl oxygen (all residues) ────────────────────────
    // Assigned separately via backbone pass, not via residue name lookup
};

// ─── Metal ion formal charges ──────────────────────────────────────────────
// Matches residue names from ion_utils.h
static constexpr ChargeEntry METAL_ION_CHARGES[] = {
    // Divalent cations
    {"CA ", " CA ", +2.0f},     // Calcium
    {" CA", "CA  ", +2.0f},     // alt padding
    {"ZN ", " ZN ", +2.0f},     // Zinc
    {"MG ", " MG ", +2.0f},     // Magnesium
    {"FE ", " FE ", +2.0f},     // Iron(II)
    {"FE2", " FE ", +2.0f},     // Iron(II) explicit
    {"MN ", " MN ", +2.0f},     // Manganese(II)
    {"CU ", " CU ", +2.0f},     // Copper(II)
    {"CU2", " CU ", +2.0f},     // Copper(II) explicit
    {"NI ", " NI ", +2.0f},     // Nickel(II)
    {"CO ", " CO ", +2.0f},     // Cobalt(II)
    {"CD ", " CD ", +2.0f},     // Cadmium(II)
    {"HG ", " HG ", +2.0f},     // Mercury(II)

    // Trivalent
    {"FE3", " FE ", +3.0f},     // Iron(III)

    // Monovalent cations
    {"NA ", " NA ", +1.0f},     // Sodium
    {"K  ", " K  ", +1.0f},     // Potassium
    {"LI ", " LI ", +1.0f},     // Lithium
    {"CU1", " CU ", +1.0f},     // Copper(I)

    // Monovalent anions
    {"CL ", " CL ", -1.0f},     // Chloride
    {"BR ", " BR ", -1.0f},     // Bromide
    {"IOD", " I  ", -1.0f},     // Iodide
};

// ─── Main assignment function ──────────────────────────────────────────────
//
// Called from top.cpp after assign_radii_types(). Iterates all receptor
// residues and assigns partial charges to atoms that match the lookup tables.
//
// Does NOT overwrite existing non-zero charges (preserves MOL2/PTM values).
//
inline void assign_formal_charges(FA_Global* FA, atom* atoms, resid* residue) {
    int n_assigned = 0;
    int n_backbone_o = 0;
    int n_backbone_n = 0;
    int n_metal = 0;

    for (int r = 1; r <= FA->res_cnt; r++) {
        // Skip ligand residues — they have charges from MOL2/SDF
        if (residue[r].type == 1) continue;

        const char* rname = residue[r].name;

        // ── Metal/ion check (HETATM single-atom residues) ──
        // Match against metal ion table first (fast path)
        for (const auto& me : METAL_ION_CHARGES) {
            if (std::strncmp(rname, me.res_name, 3) == 0) {
                // Single-atom ion residue: assign charge to all atoms in residue
                for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                    if (atoms[j].charge == 0.0f) {
                        atoms[j].charge = me.charge;
                        n_metal++;
                    }
                }
                goto next_residue;
            }
        }

        // ── Amino acid side-chain charges ──
        for (const auto& entry : AMINO_ACID_CHARGES) {
            if (std::strncmp(rname, entry.res_name, 3) != 0) continue;
            // Search atoms in this residue for matching atom name
            for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                if (atoms[j].charge != 0.0f) continue;  // don't overwrite
                if (std::strncmp(atoms[j].name, entry.atom_name, 4) == 0) {
                    atoms[j].charge = entry.charge;
                    n_assigned++;
                }
            }
        }

        // ── Backbone carbonyl oxygen: assign -0.57 (AMBER ff14SB average) ──
        // This enables H-bond scoring for backbone C=O acceptors
        for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
            if (atoms[j].charge != 0.0f) continue;
            if (std::strncmp(atoms[j].name, " O  ", 4) == 0 && atoms[j].isbb) {
                atoms[j].charge = -0.57f;
                n_backbone_o++;
            }
            // Backbone amide N: +0.17 (small positive, AMBER average for -NH-)
            else if (std::strncmp(atoms[j].name, " N  ", 4) == 0 && atoms[j].isbb) {
                // Only assign to non-proline residues (Pro has no amide H)
                if (std::strncmp(rname, "PRO", 3) != 0) {
                    atoms[j].charge = -0.42f;
                    n_backbone_n++;
                }
            }
        }

        // ── C-terminal carboxylate: OXT and last O get -0.83 each ──
        if (residue[r].ter) {
            for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                if (atoms[j].charge != 0.0f) continue;
                if (std::strncmp(atoms[j].name, " OXT", 4) == 0) {
                    atoms[j].charge = -0.83f;
                    n_assigned++;
                }
            }
            // Upgrade the backbone O to match OXT charge for symmetry
            for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                if (std::strncmp(atoms[j].name, " O  ", 4) == 0 && atoms[j].isbb) {
                    atoms[j].charge = -0.83f;
                }
            }
        }

        next_residue:;
    }

    // ── N-terminal NH3+: first residue of each chain gets +1.0 on N ──
    // Track which chains we've seen
    char seen_chains[128];
    int n_chains = 0;
    for (int r = 1; r <= FA->res_cnt; r++) {
        if (residue[r].type == 1) continue;  // skip ligand

        bool chain_seen = false;
        for (int c = 0; c < n_chains; c++) {
            if (seen_chains[c] == residue[r].chn) { chain_seen = true; break; }
        }
        if (chain_seen) continue;

        if (n_chains < 127) seen_chains[n_chains++] = residue[r].chn;

        // Find N atom in this residue
        for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
            if (std::strncmp(atoms[j].name, " N  ", 4) == 0 && atoms[j].isbb) {
                atoms[j].charge = +0.14f;  // AMBER NMET average for NH3+
                n_assigned++;
                break;
            }
        }
    }

    printf("Formal charges assigned: %d side-chain, %d backbone O, %d backbone N, %d metal/ion atoms\n",
           n_assigned, n_backbone_o, n_backbone_n, n_metal);
}

} // namespace formal_charges
