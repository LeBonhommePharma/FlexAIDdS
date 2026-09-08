// PoseLocalThermoRewrite.h — experimental docking-pose rewrite of local SS ΔH/ΔS
//
// Physics (same algebra for every polymer class):
//   docking poses rewrite *local* secondary-structure ΔH and ΔS on annotated
//   decision elements; ΔG = ΔH − TΔS; mixture over top-k poses.
//
// Class-specific increment tables are never shared:
//   RNA hairpin  ≠  DNA duplex/hairpin  ≠  protein α-helix  ≠  protein β-sheet.
// Filled defaults are peer-reviewed NN / helix–coil means (Xia 1998 RNA,
// SantaLucia 2004 DNA, Scholtz 1991 helix ΔH). Motifs without a calorimetric
// consensus stay {0,0} experimental gaps — never invent numbers, never alias
// RNA Turner parameters onto DNA or protein.
//
// This is thermodynamic rewrite glue. It is NOT k_fold / Arrhenius coupling
// in RibosomeElongation (that remains a later, optional PR).
//
// Status: EXPERIMENTAL, default OFF, not claim-ready. Do not feed Astex /
// FlexADS claim contracts. Do not silently mutate validated StatMech thermo.
// Optional DualAssembly hook stores a diagnostic mixture only; G_natural stays
// unset unless a labelled caller later opts in (not this PR).
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
    RNA_STEM_LOOP,   // Xia 1998 INN-HB WC-stack mean; loop H-bond is an experimental gap
    DNA_STEM_LOOP,   // SantaLucia 2004 NN-stack mean — NEVER alias the RNA table
    PROTEIN_HELIX,   // Scholtz 1991 calorimetric helix-formation ΔH; packing/clash gaps
    PROTEIN_SHEET,   // all contact kinds experimental gaps until a per-bridge consensus
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

// Unset slot: motifs without a calorimetric / NN consensus. Burgundy harvest
// fills these; do not invent placeholder ΔH/ΔS.
inline constexpr ThermoIncrement kExperimentalGapIncrement{0.0, 0.0};

template <std::size_t N>
inline constexpr double pose_local_mean(const std::array<double, N>& values) noexcept
{
    double sum = 0.0;
    for (double v : values)
        sum += v;
    return sum / static_cast<double>(N);
}

// Xia, SantaLucia, Turner et al., Biochemistry 37:14719 (1998), doi:10.1021/bi9809425
// INN-HB unique Watson–Crick *propagation* stacks (exclude initiation, AU-end,
// symmetry). Values restated as “1998 Model” in Zuber et al., Nucleic Acids Res.
// 50:5251 (2022) Table 1A, doi:10.1093/nar/gkac261 (PMC9122537).
// Order: GC/CG, CC/GG, GA/CU, CG/GC, AC/UG, CA/GU, AG/UC, UA/AU, AU/UA, AA/UU.
// Sequence-specific NN lookup is future work; the default is this unweighted mean.
inline constexpr double kXia1998T37_K = 310.15;
inline constexpr std::array<double, 10> kXia1998RnaWcStack_dH_kcal{
    -14.88, -13.39, -12.44, -10.64, -11.40,
    -10.44, -10.48,  -7.69,  -9.38,  -6.82};
inline constexpr std::array<double, 10> kXia1998RnaWcStack_dG37_kcal{
    -3.42, -3.26, -2.35, -2.36, -2.08,
    -2.11, -2.08, -1.33, -1.10, -0.93};

inline constexpr double kXia1998RnaWcStackMean_dH_kcal =
    pose_local_mean(kXia1998RnaWcStack_dH_kcal);

inline constexpr double xia1998_rna_wc_stack_mean_dS_cal() noexcept
{
    double sum = 0.0;
    for (std::size_t i = 0; i < kXia1998RnaWcStack_dH_kcal.size(); ++i) {
        sum += (kXia1998RnaWcStack_dH_kcal[i] - kXia1998RnaWcStack_dG37_kcal[i]) /
               kXia1998T37_K * kCalPerKcal;
    }
    return sum / static_cast<double>(kXia1998RnaWcStack_dH_kcal.size());
}

inline constexpr double kXia1998RnaWcStackMean_dS_cal =
    xia1998_rna_wc_stack_mean_dS_cal();

// SantaLucia & Hicks, Annu. Rev. Biophys. Biomol. Struct. 33:415 (2004)
// doi:10.1146/annurev.biophys.32.110601.141800 — Table 1, 1 M NaCl,
// 10 WC propagation stacks (exclude Initiation, Terminal AT, Symmetry).
// Order: AA/TT, AT/TA, TA/AT, CA/GT, GT/CA, CT/GA, GA/CT, CG/GC, GC/CG, GG/CC.
inline constexpr std::array<double, 10> kSantaLucia2004DnaNn_dH_kcal{
    -7.6, -7.2, -7.2, -8.5, -8.4, -7.8, -8.2, -10.6, -9.8, -8.0};
inline constexpr std::array<double, 10> kSantaLucia2004DnaNn_dS_cal{
    -21.3, -20.4, -21.3, -22.7, -22.4, -21.0, -22.2, -27.2, -24.4, -19.9};

inline constexpr double kSantaLucia2004DnaNnMean_dH_kcal =
    pose_local_mean(kSantaLucia2004DnaNn_dH_kcal);
inline constexpr double kSantaLucia2004DnaNnMean_dS_cal =
    pose_local_mean(kSantaLucia2004DnaNn_dS_cal);

// Scholtz, Marqusee, Baldwin et al., PNAS 88:2854 (1991), doi:10.1073/pnas.88.7.2854
// Calorimetric helix *unfolding* ΔH ≈ +1.3 kcal mol⁻¹ residue⁻¹ (best estimate;
// lower limit 0.9 if ΔCp = 0). Formation increment used here is the negative:
// ΔH = −1.3 kcal mol⁻¹ residue⁻¹. Companion helix–coil CD fit: Scholtz, Qian,
// Baldwin, Biopolymers 31:1463 (1991), doi:10.1002/bip.360311304 (ΔH° ≈ −0.955
// kcal mol⁻¹ residue⁻¹). Prefer the calorimetric ΔH; per-residue ΔS is left as
// an experimental gap (do not invent ΔS from s without a published consensus).
inline constexpr double kScholtz1991HelixFormation_dH_kcal = -1.3;

static_assert(kXia1998RnaWcStackMean_dH_kcal != kSantaLucia2004DnaNnMean_dH_kcal,
              "RNA Xia 1998 WC-stack mean must not equal DNA SantaLucia 2004 mean");
static_assert(kScholtz1991HelixFormation_dH_kcal < 0.0,
              "helix formation ΔH is exothermic (Scholtz 1991 calorimetry)");

// Dense per-kind table. Unused slots stay {0,0} so a DNA kind looked up in an
// RNA table cannot leak Xia/Turner increments.
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
    double      geometry_weight = 1.0;  // w in ΔH += coeff_H · w, ΔS += coeff_S · w
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
    ThermoIncrementTable protein_other;  // remains empty (experimental / unused)
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
