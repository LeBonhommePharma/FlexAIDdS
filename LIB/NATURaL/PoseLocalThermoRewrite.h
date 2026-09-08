// PoseLocalThermoRewrite.h — experimental docking-pose rewrite of local SS ΔH/ΔS
//
// Physics (same algebra for every polymer class):
//   docking poses rewrite *local* secondary-structure ΔH and ΔS on annotated
//   decision elements; ΔG = ΔH − TΔS; mixture over top-k poses.
//
// Class-specific increment tables are never shared:
//   RNA hairpin  ≠  DNA duplex/hairpin  ≠  protein α-helix  ≠  protein β-sheet.
// RNA uses Turner-shaped defaults; DNA uses SantaLucia-shaped defaults.
//
// This is thermodynamic rewrite glue. It is NOT k_fold / Arrhenius coupling
// in RibosomeElongation (that remains a later, optional PR).
//
// Status: EXPERIMENTAL, default OFF, not claim-ready. Draft increments are
// tunable placeholders. Do not feed Astex / FlexADS claim contracts. Do not
// silently mutate validated StatMech thermo. Optional DualAssembly hook stores
// a diagnostic mixture only; G_natural stays unset unless a labelled caller
// later opts in (not this PR).
//
// Lineage:
//   NATURAL = Native Assembly of Transcriptionally / Translationally Unified
//             Receptor and Ligand.
//   Transcriptional → RNA/DNA + ligand (Zhao–Zhang–Chen J. Chem. Phys. 135,
//     245101 (2011), doi:10.1063/1.3671644).
//   Translational → protein + ligand (Zhao et al. J. Phys. Chem. B 2011, 115,
//     3987 — ribosome/protein master equation only).
//   2014-07-25 LP Morency R2 Perl NATURAL: clustering, population inheritance,
//     pause at nt repeats, docking poses rewrite helix ΔH/ΔS.
//   2026 DualAssembly (LIB/NATURaL/) is the broader C++ home.
//
// Units (Turner / SantaLucia nearest-neighbour convention):
//   ΔH  kcal mol⁻¹
//   ΔS  cal mol⁻¹ K⁻¹  (entropy units)
//   ΔG  kcal mol⁻¹     = ΔH − T·(ΔS / 1000)
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace natural {

// ─── polymer / secondary-structure class ─────────────────────────────────────
enum class SSClass {
    RNA_STEM_LOOP,   // Turner-shaped defaults
    DNA_STEM_LOOP,   // SantaLucia-shaped defaults — NEVER alias the RNA table
    PROTEIN_HELIX,
    PROTEIN_SHEET,
    PROTEIN_OTHER    // default: no rewrite unless annotated with a filled table
};

// Contact kinds that a pose may contribute. Each kind is legal for exactly one
// SSClass; RNA kinds never read the DNA or protein tables.
enum class ContactKind {
    RNA_STEM_STACK = 0,
    RNA_LOOP_HBOND,
    DNA_STEM_STACK,
    DNA_LOOP_HBOND,
    DNA_INTERCALATION,
    PROTEIN_HELIX_BB_HBOND,
    PROTEIN_HELIX_PACKING,
    PROTEIN_HELIX_CLASH,
    PROTEIN_SHEET_BRIDGE_HBOND,
    PROTEIN_SHEET_HYDROPHOBIC_PACK,
    PROTEIN_SHEET_SHEAR,
    COUNT
};

inline constexpr int kContactKindCount = static_cast<int>(ContactKind::COUNT);

// Matches statmech::kB_kcal. Duplicated so this module has no GA / StatMech
// / hardware-dispatch dependency.
inline constexpr double kPoseRewrite_kB_kcal_mol_K = 0.001987206;
inline constexpr double kCalPerKcal = 1000.0;

struct ThermoIncrement {
    double dH_kcal = 0.0;            // kcal mol⁻¹ per unit geometry weight
    double dS_cal_per_mol_K = 0.0;   // cal mol⁻¹ K⁻¹ per unit geometry weight
};

// Dense per-kind table. Unused slots stay {0,0} so a DNA kind looked up in an
// RNA table cannot leak Turner increments.
struct ThermoIncrementTable {
    std::array<ThermoIncrement, kContactKindCount> increments{};
};

struct DecisionElement {
    SSClass     ss = SSClass::PROTEIN_OTHER;
    std::string label;
    // RNA/DNA stem-loop: inclusive 0-based nucleotide index range.
    int na_start = -1;
    int na_end   = -1;
    // Protein: inclusive 0-based residue range of the annotated helix/sheet.
    int res_start = -1;
    int res_end   = -1;
    // Optional β-sheet partner strand (inclusive). -1/-1 → partner not required.
    int sheet_partner_start = -1;
    int sheet_partner_end   = -1;
};

struct PoseContact {
    ContactKind kind = ContactKind::RNA_STEM_STACK;
    double      geometry_weight = 1.0;  // w in the draft increment tables
    int         nt_index = -1;          // RNA/DNA
    int         res_index = -1;         // protein
    int         partner_res_index = -1; // sheet partner (optional)
    std::string element_label;          // if non-empty, prefer this labelled element
};

struct PoseView {
    double score_kcal = 0.0;  // CF/proxy; used only to pick top-k (lower is better)
    std::vector<PoseContact> contacts;
};

struct PoseThermoRewriteConfig {
    ThermoIncrementTable rna;
    ThermoIncrementTable dna;
    ThermoIncrementTable protein_helix;
    ThermoIncrementTable protein_sheet;
    ThermoIncrementTable protein_other;  // empty in the draft defaults
    double temperature_K = 298.15;
    int    top_k = 10;
    bool   require_decision_element_only = true;
};

struct LocalThermoPatch {
    std::string label;
    SSClass     ss = SSClass::PROTEIN_OTHER;
    double      dH_kcal = 0.0;
    double      dS_cal_per_mol_K = 0.0;
    double      dG_kcal = 0.0;
};

struct LocalThermoMixture {
    double dH_kcal = 0.0;
    double dS_cal_per_mol_K = 0.0;
    double dG_kcal = 0.0;
    int    n_poses_used = 0;
    std::vector<double>           pose_weights;  // Boltzmann p_i for used poses
    std::vector<LocalThermoPatch> per_element;
    bool empty() const noexcept { return n_poses_used == 0 && per_element.empty(); }
};

const char* ss_class_name(SSClass ss) noexcept;
const char* contact_kind_name(ContactKind kind) noexcept;
SSClass     ss_class_for_kind(ContactKind kind) noexcept;

inline double delta_G_kcal(double dH_kcal, double dS_cal_per_mol_K, double T_K) noexcept
{
    return dH_kcal - T_K * (dS_cal_per_mol_K / kCalPerKcal);
}

PoseThermoRewriteConfig default_pose_thermo_rewrite_config();

const ThermoIncrementTable& table_for_ss(const PoseThermoRewriteConfig& cfg,
                                         SSClass ss) noexcept;

// Pure function: rewrite local SS ΔH/ΔS from labelled pose contacts.
// Returns an empty mixture (finite zeros) when nothing matches.
LocalThermoMixture rewrite_local_thermo_from_poses(
    const std::vector<PoseView>&         poses,
    const std::vector<DecisionElement>&  elements,
    const PoseThermoRewriteConfig&       cfg);

} // namespace natural
