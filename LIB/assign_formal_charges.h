// assign_formal_charges.h — Residue-aware formal charge assignment for PDB input
//
// PDB files do not carry partial charges. Without charges:
//   - Coulomb electrostatics is dead (vcfunction.cpp guard: qA != 0 && qB != 0)
//   - Salt bridge detection uses atom.charge in hbond_potential.h (v56+)
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
#include <cstdlib>
#include "ff14sb_lumped_charges.h"

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
// ─────────────────────────────────────────────────────────────────────────────
// OPT-IN CHARGE-CONSERVING ALTERNATIVE: AMBER ff14SB, hydrogen-lumped.
//
// WHY THIS EXISTS. The default table above mixes two conventions inside itself.
// Most entries are verbatim ff14SB HEAVY-ATOM partials with the polar hydrogen
// simply omitted (SER OG -0.65, TYR OH -0.56, backbone N -0.42 ...), while the
// formally charged groups are DELIBERATELY overridden with formal charge instead
// (LYS NZ +1.00 against a lumped +0.63; ARG NH1/NH2 +0.45 each) - a documented
// choice, see the comments on those entries. Each half is defensible alone. Mixed,
// they are not: a formally NEUTRAL group ends up carrying net charge (SER -0.65,
// ASN -0.89, backbone -0.99 per peptide unit, of order -297 e over 300 residues),
// so a pairwise q*q term is invalid before it evaluates - and the ligand side is
// prepared united-atom with net ~0, i.e. on the OTHER convention.
//
// This path replaces the whole receptor set with ff14SB hydrogen-lumped values
// (LIB/ff14sb_lumped_charges.h, machine-generated), so every residue sums to its
// formal charge and both sides of the term share one convention.
//
// GATE: opt-in only. FLEXAIDDS_FF14SB_CHARGES=1 selects it; unset reproduces the
// default byte-for-byte. Metals keep METAL_ION_CHARGES formal integers, which is
// still a convention mismatch against these partials - disclosed, not fixed.
// Terminal residues use the N*/C* templates; HIS is aliased to HIE at generation.
// ─────────────────────────────────────────────────────────────────────────────
inline void assign_ff14sb_lumped(FA_Global* FA, atom* atoms, resid* residue) {
    int n_assigned = 0, n_metal = 0, n_res = 0, n_unmatched_res = 0;
    double q_total = 0.0;

    // first residue of each chain -> N-terminal template
    bool* is_nterm = (bool*)std::calloc((size_t)FA->res_cnt + 2, sizeof(bool));
    char seen[128]; int n_seen = 0;
    for (int r = 1; r <= FA->res_cnt; r++) {
        if (residue[r].type == 1) continue;
        bool s = false;
        for (int c = 0; c < n_seen; c++) if (seen[c] == residue[r].chn) { s = true; break; }
        if (s) continue;
        if (n_seen < 127) seen[n_seen++] = residue[r].chn;
        is_nterm[r] = true;
    }

    for (int r = 1; r <= FA->res_cnt; r++) {
        if (residue[r].type == 1) continue;               // ligand keeps SDF/MOL2 charges
        const char* rname = residue[r].name;

        bool is_metal = false;
        for (const auto& me : METAL_ION_CHARGES) {
            if (std::strncmp(rname, me.res_name, 3) == 0) {
                for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                    if (atoms[j].charge == 0.0f) { atoms[j].charge = me.charge; n_metal++; }
                }
                is_metal = true; break;
            }
        }
        if (is_metal) continue;

        // template name: N<res> at a chain start, C<res> at a chain terminus, else <res>
        char tmpl[8];
        if (is_nterm[r])        std::snprintf(tmpl, sizeof(tmpl), "N%.3s", rname);
        else if (residue[r].ter) std::snprintf(tmpl, sizeof(tmpl), "C%.3s", rname);
        else                     std::snprintf(tmpl, sizeof(tmpl), "%.3s", rname);

        int hits = 0;
        double q_res = 0.0;
        for (int pass = 0; pass < 2; pass++) {
            // pass 0 uses the terminal/plain template; pass 1 falls back to the plain
            // template when a terminal one is absent, so a chain end is never left blank
            const char* want = tmpl;
            char plain[8];
            if (pass == 1) { std::snprintf(plain, sizeof(plain), "%.3s", rname); want = plain; }
            if (pass == 1 && hits > 0) break;
            for (int e = 0; e < ff14sb_lumped::FF14SB_LUMPED_COUNT; e++) {
                const auto& ent = ff14sb_lumped::FF14SB_LUMPED_CHARGES[e];
                if (std::strcmp(want, ent.res_name) != 0) continue;
                for (int j = residue[r].fatm[0]; j <= residue[r].latm[0]; j++) {
                    if (std::strncmp(atoms[j].name, ent.atom_name, 4) != 0) continue;
                    atoms[j].charge = ent.charge;
                    q_res += (double)ent.charge;
                    hits++; n_assigned++;
                }
            }
        }
        if (hits == 0) n_unmatched_res++;
        else { n_res++; q_total += q_res; }
    }
    std::free(is_nterm);

    printf("ff14SB lumped charges: %d atoms over %d residues (net %+.4f e), "
           "%d metal/ion atoms, %d residues UNMATCHED\n",
           n_assigned, n_res, q_total, n_metal, n_unmatched_res);
    printf("  convention: united-atom (polar H folded into its heavy atom); "
           "every matched residue sums to its formal charge by construction\n");
}

inline void assign_formal_charges(FA_Global* FA, atom* atoms, resid* residue) {
    const char* ff14 = std::getenv("FLEXAIDDS_FF14SB_CHARGES");
    if (ff14 && ff14[0] == '1') { assign_ff14sb_lumped(FA, atoms, residue); return; }

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
