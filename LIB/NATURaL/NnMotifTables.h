// NnMotifTables.h — peer-reviewed NN / helix–coil motif ΔH/ΔS for PoseLocalThermoRewrite
//
// Ligand pose weight w scales the matched published stack or loop term.
// RNA ≠ DNA ≠ helix ≠ sheet. Never alias Xia/Turner onto SantaLucia or protein.
//
// DNA Table 2: SantaLucia Jr (1998) PNAS 95:1460, doi:10.1073/pnas.95.4.1460
//   (PMC19045). Unified 1 M NaCl. Review: SantaLucia & Hicks 2004.
// RNA WC stacks: Xia et al. (1998) Biochemistry 37:14719, doi:10.1021/bi9809425
//   NNDB: Turner & Mathews (2010) NAR, doi:10.1093/nar/gkp892
//   https://rna.urmc.rochester.edu/NNDB/
// RNA hairpin ΔH°: Lu, Turner, Mathews (2006) NAR 34:4912, doi:10.1093/nar/gkl472
//
// TSV copies (do not invent numbers): LIB/NATURaL/data/
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include <array>
#include <cctype>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>

namespace natural {

struct MotifDeltaHS {
    double dH_kcal = 0.0;
    double dS_cal_per_mol_K = 0.0;
};

inline constexpr double kNnT37_K = 310.15;

inline constexpr double nn_delta_S_cal(double dH_kcal, double dG_kcal, double T_K) noexcept
{
    return (dH_kcal - dG_kcal) / T_K * 1000.0;
}

template <std::size_t N>
inline constexpr double motif_mean_dH(const std::array<MotifDeltaHS, N>& rows,
                                      std::size_t count) noexcept
{
    double sum = 0.0;
    for (std::size_t i = 0; i < count; ++i)
        sum += rows[i].dH_kcal;
    return sum / static_cast<double>(count);
}

template <std::size_t N>
inline constexpr double motif_mean_dS(const std::array<MotifDeltaHS, N>& rows,
                                      std::size_t count) noexcept
{
    double sum = 0.0;
    for (std::size_t i = 0; i < count; ++i)
        sum += rows[i].dS_cal_per_mol_K;
    return sum / static_cast<double>(count);
}

// ── DNA: SantaLucia 1998 PNAS Table 2 (embed literally) ─────────────────────
// Order matches the published table: 10 WC stacks, then Init G·C, Init A·T, Sym.
inline constexpr std::array<const char*, 13> kSantaLucia1998Motif{
    "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA", "CT/GA", "GA/CT",
    "CG/GC", "GC/CG", "GG/CC", "INIT_GC", "INIT_AT", "SYM"};

inline constexpr std::array<MotifDeltaHS, 13> kSantaLucia1998PnasTable2{{
    {-7.9, -22.2},
    {-7.2, -20.4},
    {-7.2, -21.3},
    {-8.5, -22.7},
    {-8.4, -22.4},
    {-7.8, -21.0},
    {-8.2, -22.2},
    {-10.6, -27.2},
    {-9.8, -24.4},
    {-8.0, -19.9},
    {0.1, -2.8},   // Init. w/term. G·C
    {2.3, 4.1},    // Init. w/term. A·T
    {0.0, -1.4},   // Symmetry correction
}};

inline constexpr std::size_t kSantaLucia1998WcStackCount = 10;

inline constexpr double kSantaLucia1998DnaWcStackMean_dH_kcal =
    motif_mean_dH(kSantaLucia1998PnasTable2, kSantaLucia1998WcStackCount);
inline constexpr double kSantaLucia1998DnaWcStackMean_dS_cal =
    motif_mean_dS(kSantaLucia1998PnasTable2, kSantaLucia1998WcStackCount);

// ── RNA: Xia 1998 INN-HB unique WC propagation stacks ───────────────────────
// ΔS from published ΔH and ΔG°37 at 310.15 K. Order: Zuber 2022 Table 1A 1998 Model.
inline constexpr std::array<const char*, 10> kXia1998RnaMotif{
    "GC/CG", "CC/GG", "GA/CU", "CG/GC", "AC/UG",
    "CA/GU", "AG/UC", "UA/AU", "AU/UA", "AA/UU"};

inline constexpr std::array<double, 10> kXia1998RnaWcStack_dH_kcal{
    -14.88, -13.39, -12.44, -10.64, -11.40,
    -10.44, -10.48,  -7.69,  -9.38,  -6.82};
inline constexpr std::array<double, 10> kXia1998RnaWcStack_dG37_kcal{
    -3.42, -3.26, -2.35, -2.36, -2.08,
    -2.11, -2.08, -1.33, -1.10, -0.93};

inline constexpr MotifDeltaHS xia1998_stack_at(std::size_t i) noexcept
{
    const double dH = kXia1998RnaWcStack_dH_kcal[i];
    const double dS = nn_delta_S_cal(dH, kXia1998RnaWcStack_dG37_kcal[i], kNnT37_K);
    return {dH, dS};
}

inline constexpr std::array<MotifDeltaHS, 10> kXia1998RnaWcStacks{{
    xia1998_stack_at(0), xia1998_stack_at(1), xia1998_stack_at(2), xia1998_stack_at(3),
    xia1998_stack_at(4), xia1998_stack_at(5), xia1998_stack_at(6), xia1998_stack_at(7),
    xia1998_stack_at(8), xia1998_stack_at(9)}};

inline constexpr double kXia1998RnaWcStackMean_dH_kcal =
    motif_mean_dH(kXia1998RnaWcStacks, 10);
inline constexpr double kXia1998RnaWcStackMean_dS_cal =
    motif_mean_dS(kXia1998RnaWcStacks, 10);

// ── RNA hairpin initiation ΔH° (Lu 2006 Table 1) ────────────────────────────
// Default unmatched loop size is n = 4 (canonical tetraloop published term).
inline constexpr int kLu2006DefaultHairpinN = 4;

inline constexpr double lu2006_hairpin_initiation_dH_kcal(int unpaired_n) noexcept
{
    if (unpaired_n < 3)
        return 0.0;
    switch (unpaired_n) {
    case 3: return 1.3;
    case 4: return 4.8;
    case 5: return 3.6;
    case 6: return -2.9;
    case 7: return 1.3;
    case 8: return -2.9;
    default: return 5.0;  // n >= 9
    }
}

inline constexpr MotifDeltaHS kLu2006DefaultHairpin{
    lu2006_hairpin_initiation_dH_kcal(kLu2006DefaultHairpinN),
    0.0};  // ΔS not tabulated; extra length cost is entropic (Lu 2006)

// ── protein helix / sheet (class-specific; not NN nucleic-acid tables) ───────
// Scholtz 1991 PNAS doi:10.1073/pnas.88.7.2854 — calorimetric helix unfolding
// ΔH ≈ +1.3 kcal mol⁻¹ residue⁻¹ → formation ΔH = −1.3.
inline constexpr double kScholtz1991HelixFormation_dH_kcal = -1.3;
// Zavrtanik, Lah, Hadži, Biophys. J. 125:305 (2026), doi:10.1016/j.bpj.2025.11.2689
// PubMed 41318999. Helix→coil ΔS_BB = 5.2 ± 0.3 cal mol⁻¹ K⁻¹ per peptide unit
// → formation ΔS = −5.2 (backbone baseline).
inline constexpr double kZavrtanik2026HelixFormation_dS_cal = -5.2;

// Meier & Seelig, J. Am. Chem. Soc. (2008) doi:10.1021/ja077231r — membrane
// coil⇄β, length-dependent ΔH_fold ≈ −0.2 to −0.6 kcal mol⁻¹ residue⁻¹,
// TΔS_fold ≈ −0.1 to −0.5 kcal mol⁻¹ residue⁻¹. Default mid, labelled experimental:
// ΔH = −0.4 ; ΔS ≈ −1.0 e.u. at 298 K from TΔS ≈ −0.3.
// Caveat: Deechongkit et al. Nature 430:101 (2004) doi:10.1038/nature02611 —
// β H-bond energetics are context-dependent. No universal NN sheet table.
inline constexpr double kMeierSeelig2008SheetFold_dH_kcal = -0.4;
inline constexpr double kMeierSeelig2008SheetFold_dS_cal = -1.0;
inline constexpr bool kMeierSeelig2008SheetFoldIsExperimental = true;

static_assert(kXia1998RnaWcStackMean_dH_kcal != kSantaLucia1998DnaWcStackMean_dH_kcal,
              "RNA Xia 1998 mean must not equal DNA SantaLucia 1998 Table 2 mean");
static_assert(kScholtz1991HelixFormation_dH_kcal < 0.0);
static_assert(kZavrtanik2026HelixFormation_dS_cal < 0.0);
static_assert(kMeierSeelig2008SheetFold_dH_kcal != kScholtz1991HelixFormation_dH_kcal);

inline std::string normalize_nn_motif(std::string_view raw)
{
    std::string out;
    out.reserve(raw.size());
    for (unsigned char ch : raw) {
        if (ch == ' ' || ch == '\t')
            continue;
        out.push_back(static_cast<char>(std::toupper(ch)));
    }
    return out;
}

inline std::optional<MotifDeltaHS> lookup_santalucia1998_dna(std::string_view motif)
{
    if (motif.empty())
        return MotifDeltaHS{kSantaLucia1998DnaWcStackMean_dH_kcal,
                            kSantaLucia1998DnaWcStackMean_dS_cal};
    const std::string key = normalize_nn_motif(motif);
    for (std::size_t i = 0; i < kSantaLucia1998Motif.size(); ++i) {
        if (key == kSantaLucia1998Motif[i])
            return kSantaLucia1998PnasTable2[i];
    }
    return std::nullopt;
}

inline std::optional<MotifDeltaHS> lookup_xia1998_rna(std::string_view motif)
{
    if (motif.empty())
        return MotifDeltaHS{kXia1998RnaWcStackMean_dH_kcal, kXia1998RnaWcStackMean_dS_cal};
    const std::string key = normalize_nn_motif(motif);
    for (std::size_t i = 0; i < kXia1998RnaMotif.size(); ++i) {
        if (key == kXia1998RnaMotif[i])
            return kXia1998RnaWcStacks[i];
    }
    return std::nullopt;
}

// loop_unpaired < 0 → Lu 2006 n=4 default. 0..2 → not a hairpin (nullopt).
inline std::optional<MotifDeltaHS> lookup_lu2006_rna_hairpin(int loop_unpaired)
{
    const int n = (loop_unpaired < 0) ? kLu2006DefaultHairpinN : loop_unpaired;
    if (n < 3)
        return std::nullopt;
    return MotifDeltaHS{lu2006_hairpin_initiation_dH_kcal(n), 0.0};
}

} // namespace natural
