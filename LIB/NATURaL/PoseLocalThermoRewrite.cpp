// PoseLocalThermoRewrite.cpp — see header for physics, citations, and status.
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#include "PoseLocalThermoRewrite.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>

namespace natural {

namespace {

bool in_closed_range(int value, int lo, int hi) noexcept
{
    if (lo < 0 || hi < 0 || hi < lo) return false;
    return value >= lo && value <= hi;
}

bool na_range_set(const DecisionElement& el) noexcept
{
    return el.na_start >= 0 && el.na_end >= el.na_start;
}

bool res_range_set(const DecisionElement& el) noexcept
{
    return el.res_start >= 0 && el.res_end >= el.res_start;
}

bool partner_range_set(const DecisionElement& el) noexcept
{
    return el.sheet_partner_start >= 0 &&
           el.sheet_partner_end >= el.sheet_partner_start;
}

bool contact_hits_element(const PoseContact& c, const DecisionElement& el) noexcept
{
    if (!c.element_label.empty())
        return c.element_label == el.label;

    switch (el.ss) {
    case SSClass::RNA_STEM_LOOP:
    case SSClass::DNA_STEM_LOOP:
        if (!na_range_set(el))
            return false;
        return in_closed_range(c.nt_index, el.na_start, el.na_end);
    case SSClass::PROTEIN_HELIX:
    case SSClass::PROTEIN_OTHER:
        if (!res_range_set(el))
            return false;
        return in_closed_range(c.res_index, el.res_start, el.res_end);
    case SSClass::PROTEIN_SHEET: {
        if (!res_range_set(el))
            return false;
        const bool on_primary = in_closed_range(c.res_index, el.res_start, el.res_end);
        const bool on_partner_as_primary =
            partner_range_set(el) &&
            in_closed_range(c.res_index, el.sheet_partner_start, el.sheet_partner_end);
        if (!on_primary && !on_partner_as_primary)
            return false;
        if (!partner_range_set(el))
            return true;
        const bool partner_ok =
            in_closed_range(c.partner_res_index, el.sheet_partner_start, el.sheet_partner_end) ||
            in_closed_range(c.partner_res_index, el.res_start, el.res_end) ||
            c.partner_res_index < 0;
        return partner_ok;
    }
    }
    return false;
}

const DecisionElement* first_matching_element(
    const PoseContact& c,
    const std::vector<DecisionElement>& elements) noexcept
{
    const SSClass need = ss_class_for_kind(c.kind);
    for (const DecisionElement& el : elements) {
        if (el.ss != need)
            continue;
        if (contact_hits_element(c, el))
            return &el;
    }
    return nullptr;
}

ThermoIncrement lookup(const ThermoIncrementTable& table, ContactKind kind) noexcept
{
    const int i = static_cast<int>(kind);
    if (i < 0 || i >= kContactKindCount)
        return {};
    return table.increments[static_cast<std::size_t>(i)];
}

ThermoIncrement from_motif(const MotifDeltaHS& m) noexcept
{
    return ThermoIncrement{m.dH_kcal, m.dS_cal_per_mol_K};
}

// Sequence-specific NN / loop terms override the class-table mean when the
// contact carries a motif or unpaired-loop annotation. Unknown motifs skip.
std::optional<ThermoIncrement> increment_for_contact(
    const PoseContact& c,
    const DecisionElement& el,
    const PoseThermoRewriteConfig& cfg)
{
    ThermoIncrement inc = lookup(table_for_ss(cfg, el.ss), c.kind);
    if (c.kind == ContactKind::RNA_STEM_STACK) {
        const auto m = lookup_xia1998_rna(c.motif);
        if (!m.has_value())
            return std::nullopt;
        inc = from_motif(*m);
    } else if (c.kind == ContactKind::DNA_STEM_STACK) {
        const auto m = lookup_santalucia1998_dna(c.motif);
        if (!m.has_value())
            return std::nullopt;
        inc = from_motif(*m);
    } else if (c.kind == ContactKind::RNA_LOOP_HBOND) {
        const auto m = lookup_lu2006_rna_hairpin(c.loop_unpaired);
        if (!m.has_value())
            return std::nullopt;
        inc = from_motif(*m);
    }
    return inc;
}

void set_inc(ThermoIncrementTable& table, ContactKind kind,
             double dH_kcal, double dS_cal_per_mol_K) noexcept
{
    const int i = static_cast<int>(kind);
    if (i < 0 || i >= kContactKindCount)
        return;
    table.increments[static_cast<std::size_t>(i)] = ThermoIncrement{dH_kcal, dS_cal_per_mol_K};
}

double log_sum_exp(const std::vector<double>& values) noexcept
{
    if (values.empty())
        return -std::numeric_limits<double>::infinity();
    double xmax = values[0];
    for (std::size_t i = 1; i < values.size(); ++i)
        if (values[i] > xmax) xmax = values[i];
    if (!std::isfinite(xmax))
        return xmax;
    double sum = 0.0;
    for (double v : values)
        sum += std::exp(v - xmax);
    if (!(sum > 0.0) || !std::isfinite(sum))
        return xmax;
    return xmax + std::log(sum);
}

struct PoseAccum {
    std::size_t original_index = 0;
    double score_kcal = 0.0;
    double dH_kcal = 0.0;
    double dS_cal_per_mol_K = 0.0;
    double dG_kcal = 0.0;
    std::vector<LocalThermoPatch> per_element;
};

LocalThermoPatch& patch_for(std::vector<LocalThermoPatch>& patches,
                            const DecisionElement& el)
{
    for (LocalThermoPatch& p : patches) {
        if (p.label == el.label && p.ss == el.ss)
            return p;
    }
    LocalThermoPatch created;
    created.label = el.label;
    created.ss = el.ss;
    patches.push_back(created);
    return patches.back();
}

bool accumulate_pose(const PoseView& pose,
                     const std::vector<DecisionElement>& elements,
                     const PoseThermoRewriteConfig& cfg,
                     PoseAccum& out)
{
    out.score_kcal = pose.score_kcal;
    out.dH_kcal = 0.0;
    out.dS_cal_per_mol_K = 0.0;
    out.per_element.clear();

    bool any = false;
    for (const PoseContact& c : pose.contacts) {
        if (!std::isfinite(c.geometry_weight))
            continue;
        const DecisionElement* el = first_matching_element(c, elements);
        if (el == nullptr)
            continue;  // require_decision_element_only: unmatched contacts drop
        const std::optional<ThermoIncrement> maybe = increment_for_contact(c, *el, cfg);
        if (!maybe.has_value())
            continue;
        const ThermoIncrement inc = *maybe;
        if (el->ss == SSClass::PROTEIN_OTHER &&
            inc.dH_kcal == 0.0 && inc.dS_cal_per_mol_K == 0.0)
            continue;
        const double w = c.geometry_weight;
        const double dH = inc.dH_kcal * w;
        const double dS = inc.dS_cal_per_mol_K * w;
        if (!std::isfinite(dH) || !std::isfinite(dS))
            continue;
        LocalThermoPatch& patch = patch_for(out.per_element, *el);
        patch.dH_kcal += dH;
        patch.dS_cal_per_mol_K += dS;
        out.dH_kcal += dH;
        out.dS_cal_per_mol_K += dS;
        any = true;
    }
    out.dG_kcal = delta_G_kcal(out.dH_kcal, out.dS_cal_per_mol_K, cfg.temperature_K);
    for (LocalThermoPatch& p : out.per_element)
        p.dG_kcal = delta_G_kcal(p.dH_kcal, p.dS_cal_per_mol_K, cfg.temperature_K);
    return any && std::isfinite(out.dG_kcal) && std::isfinite(out.dH_kcal) &&
           std::isfinite(out.dS_cal_per_mol_K);
}

} // namespace

const char* ss_class_name(SSClass ss) noexcept
{
    switch (ss) {
    case SSClass::RNA_STEM_LOOP:  return "RNA_STEM_LOOP";
    case SSClass::DNA_STEM_LOOP:  return "DNA_STEM_LOOP";
    case SSClass::PROTEIN_HELIX:  return "PROTEIN_HELIX";
    case SSClass::PROTEIN_SHEET:  return "PROTEIN_SHEET";
    case SSClass::PROTEIN_OTHER:  return "PROTEIN_OTHER";
    }
    return "UNKNOWN";
}

const char* contact_kind_name(ContactKind kind) noexcept
{
    switch (kind) {
    case ContactKind::RNA_STEM_STACK:                 return "RNA_STEM_STACK";
    case ContactKind::RNA_LOOP_HBOND:                 return "RNA_LOOP_HBOND";
    case ContactKind::DNA_STEM_STACK:                 return "DNA_STEM_STACK";
    case ContactKind::DNA_LOOP_HBOND:                 return "DNA_LOOP_HBOND";
    case ContactKind::DNA_INTERCALATION:              return "DNA_INTERCALATION";
    case ContactKind::PROTEIN_HELIX_BB_HBOND:         return "PROTEIN_HELIX_BB_HBOND";
    case ContactKind::PROTEIN_HELIX_PACKING:          return "PROTEIN_HELIX_PACKING";
    case ContactKind::PROTEIN_HELIX_CLASH:            return "PROTEIN_HELIX_CLASH";
    case ContactKind::PROTEIN_SHEET_BRIDGE_HBOND:     return "PROTEIN_SHEET_BRIDGE_HBOND";
    case ContactKind::PROTEIN_SHEET_HYDROPHOBIC_PACK: return "PROTEIN_SHEET_HYDROPHOBIC_PACK";
    case ContactKind::PROTEIN_SHEET_SHEAR:            return "PROTEIN_SHEET_SHEAR";
    case ContactKind::COUNT:                          return "COUNT";
    }
    return "UNKNOWN";
}

SSClass ss_class_for_kind(ContactKind kind) noexcept
{
    switch (kind) {
    case ContactKind::RNA_STEM_STACK:
    case ContactKind::RNA_LOOP_HBOND:
        return SSClass::RNA_STEM_LOOP;
    case ContactKind::DNA_STEM_STACK:
    case ContactKind::DNA_LOOP_HBOND:
    case ContactKind::DNA_INTERCALATION:
        return SSClass::DNA_STEM_LOOP;
    case ContactKind::PROTEIN_HELIX_BB_HBOND:
    case ContactKind::PROTEIN_HELIX_PACKING:
    case ContactKind::PROTEIN_HELIX_CLASH:
        return SSClass::PROTEIN_HELIX;
    case ContactKind::PROTEIN_SHEET_BRIDGE_HBOND:
    case ContactKind::PROTEIN_SHEET_HYDROPHOBIC_PACK:
    case ContactKind::PROTEIN_SHEET_SHEAR:
        return SSClass::PROTEIN_SHEET;
    case ContactKind::COUNT:
        break;
    }
    return SSClass::PROTEIN_OTHER;
}

PoseThermoRewriteConfig default_pose_thermo_rewrite_config()
{
    PoseThermoRewriteConfig cfg;

    // RNA stem stack: Xia 1998 INN-HB mean of 10 unique WC propagation stacks.
    // doi:10.1021/bi9809425 ; Zuber 2022 NAR Table 1A "1998 Model"
    // doi:10.1093/nar/gkac261. PoseContact::motif selects a published dimer.
    set_inc(cfg.rna, ContactKind::RNA_STEM_STACK,
            kXia1998RnaWcStackMean_dH_kcal, kXia1998RnaWcStackMean_dS_cal);

    // RNA loop: Lu, Turner, Mathews, NAR 34:4912 (2006), doi:10.1093/nar/gkl472
    // Table 1 hairpin-loop *initiation* ΔH°(n). Default n=4. ΔS is not tabulated
    // (extra length cost is treated as entropic). This is initiation enthalpy,
    // not a ligand–loop H-bond increment. PoseContact::loop_unpaired selects n.
    set_inc(cfg.rna, ContactKind::RNA_LOOP_HBOND,
            kLu2006DefaultHairpin.dH_kcal, kLu2006DefaultHairpin.dS_cal_per_mol_K);

    // DNA stem stack: SantaLucia 1998 PNAS Table 2 mean of 10 WC propagation
    // stacks (1 M NaCl). doi:10.1073/pnas.95.4.1460 (PMC19045).
    // NEVER a copy of the RNA table. PoseContact::motif selects a Table 2 dimer.
    set_inc(cfg.dna, ContactKind::DNA_STEM_STACK,
            kSantaLucia1998DnaWcStackMean_dH_kcal, kSantaLucia1998DnaWcStackMean_dS_cal);

    // TODO(burgundy): DNA ligand–loop H-bond ΔH/ΔS. SantaLucia hairpin
    // initiation tables are loop-initiation ΔG, not this contact increment.
    set_inc(cfg.dna, ContactKind::DNA_LOOP_HBOND,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);

    // TODO(burgundy): DNA intercalation ΔH/ΔS. Ligand-specific (e.g. Chaires
    // reviews); no generic calorimetric consensus for a geometry-weighted
    // increment. Leave unset — do not invent a mean intercalator.
    set_inc(cfg.dna, ContactKind::DNA_INTERCALATION,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);

    // Protein α-helix backbone: Scholtz 1991 PNAS calorimetric helix-formation
    // ΔH (doi:10.1073/pnas.88.7.2854) + Zavrtanik, Lah, Hadži, Biophys. J.
    // 125:305 (2026) helix→coil ΔS_BB = 5.2 ± 0.3 e.u. → formation −5.2
    // (doi:10.1016/j.bpj.2025.11.2689, PubMed 41318999).
    set_inc(cfg.protein_helix, ContactKind::PROTEIN_HELIX_BB_HBOND,
            kScholtz1991HelixFormation_dH_kcal,
            kZavrtanik2026HelixFormation_dS_cal);

    // TODO(burgundy): helix side-chain packing ΔH/ΔS — no calorimetric
    // per-contact consensus distinct from the backbone H-bond term above.
    set_inc(cfg.protein_helix, ContactKind::PROTEIN_HELIX_PACKING,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);

    // TODO(burgundy): helix steric clash ΔH/ΔS — no calorimetric consensus.
    set_inc(cfg.protein_helix, ContactKind::PROTEIN_HELIX_CLASH,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);

    // Protein β-sheet bridge: labelled experimental midpoint from Meier &
    // Seelig, J. Am. Chem. Soc. 130:1017 (2008), doi:10.1021/ja077231r
    // (membrane coil⇄β, ΔH_fold ≈ −0.2 to −0.6 kcal mol⁻¹ residue⁻¹).
    // Caveat: Deechongkit et al. Nature 430:101 (2004) doi:10.1038/nature02611
    // — β H-bond energetics are context-dependent. Not a universal NN table.
    // Maynard 1998 whole-hairpin folding is the wrong scale — do not reuse.
    set_inc(cfg.protein_sheet, ContactKind::PROTEIN_SHEET_BRIDGE_HBOND,
            kMeierSeelig2008SheetFold_dH_kcal,
            kMeierSeelig2008SheetFold_dS_cal);
    // TODO(burgundy): sheet hydrophobic packing ΔH/ΔS.
    set_inc(cfg.protein_sheet, ContactKind::PROTEIN_SHEET_HYDROPHOBIC_PACK,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);
    // TODO(burgundy): sheet shear / mismatch ΔH/ΔS.
    set_inc(cfg.protein_sheet, ContactKind::PROTEIN_SHEET_SHEAR,
            kExperimentalGapIncrement.dH_kcal,
            kExperimentalGapIncrement.dS_cal_per_mol_K);

    return cfg;
}

const ThermoIncrementTable& table_for_ss(const PoseThermoRewriteConfig& cfg,
                                         SSClass ss) noexcept
{
    switch (ss) {
    case SSClass::RNA_STEM_LOOP: return cfg.rna;
    case SSClass::DNA_STEM_LOOP: return cfg.dna;
    case SSClass::PROTEIN_HELIX: return cfg.protein_helix;
    case SSClass::PROTEIN_SHEET: return cfg.protein_sheet;
    case SSClass::PROTEIN_OTHER: return cfg.protein_other;
    }
    return cfg.protein_other;
}

LocalThermoMixture rewrite_local_thermo_from_poses(
    const std::vector<PoseView>&        poses,
    const std::vector<DecisionElement>& elements,
    const PoseThermoRewriteConfig&      cfg)
{
    LocalThermoMixture mix;
    if (poses.empty() || !std::isfinite(cfg.temperature_K) || cfg.temperature_K <= 0.0)
        return mix;
    if (cfg.require_decision_element_only && elements.empty())
        return mix;

    std::vector<PoseAccum> scored;
    scored.reserve(poses.size());
    for (std::size_t i = 0; i < poses.size(); ++i) {
        if (!std::isfinite(poses[i].score_kcal))
            continue;
        PoseAccum acc;
        acc.original_index = i;
        if (!accumulate_pose(poses[i], elements, cfg, acc))
            continue;
        scored.push_back(std::move(acc));
    }
    if (scored.empty())
        return mix;

    std::stable_sort(scored.begin(), scored.end(),
                     [](const PoseAccum& a, const PoseAccum& b) {
                         if (a.score_kcal != b.score_kcal)
                             return a.score_kcal < b.score_kcal;
                         return a.original_index < b.original_index;
                     });
    const int k = (cfg.top_k <= 0) ? static_cast<int>(scored.size())
                                   : std::min(cfg.top_k, static_cast<int>(scored.size()));
    scored.resize(static_cast<std::size_t>(k));

    const double kT = kPoseRewrite_kB_kcal_mol_K * cfg.temperature_K;
    std::vector<double> log_w;
    log_w.reserve(scored.size());
    if (!(kT > 0.0) || !std::isfinite(kT)) {
        for (std::size_t i = 0; i < scored.size(); ++i)
            log_w.push_back(0.0);
    } else {
        for (const PoseAccum& acc : scored)
            log_w.push_back(-acc.dG_kcal / kT);
    }
    const double lnZ = log_sum_exp(log_w);
    mix.pose_weights.resize(scored.size(), 0.0);
    if (!std::isfinite(lnZ))
        return mix;

    for (std::size_t i = 0; i < scored.size(); ++i) {
        const double p = std::exp(log_w[i] - lnZ);
        mix.pose_weights[i] = std::isfinite(p) ? p : 0.0;
    }

    mix.n_poses_used = static_cast<int>(scored.size());
    std::vector<LocalThermoPatch> mixed_elements;
    for (std::size_t i = 0; i < scored.size(); ++i) {
        const double p = mix.pose_weights[i];
        mix.dH_kcal += p * scored[i].dH_kcal;
        mix.dS_cal_per_mol_K += p * scored[i].dS_cal_per_mol_K;
        for (const LocalThermoPatch& patch : scored[i].per_element) {
            LocalThermoPatch* dest = nullptr;
            for (LocalThermoPatch& m : mixed_elements) {
                if (m.label == patch.label && m.ss == patch.ss) {
                    dest = &m;
                    break;
                }
            }
            if (dest == nullptr) {
                LocalThermoPatch created;
                created.label = patch.label;
                created.ss = patch.ss;
                mixed_elements.push_back(created);
                dest = &mixed_elements.back();
            }
            dest->dH_kcal += p * patch.dH_kcal;
            dest->dS_cal_per_mol_K += p * patch.dS_cal_per_mol_K;
        }
    }
    mix.dG_kcal = delta_G_kcal(mix.dH_kcal, mix.dS_cal_per_mol_K, cfg.temperature_K);
    for (LocalThermoPatch& p : mixed_elements)
        p.dG_kcal = delta_G_kcal(p.dH_kcal, p.dS_cal_per_mol_K, cfg.temperature_K);
    mix.per_element = std::move(mixed_elements);

    if (!std::isfinite(mix.dH_kcal)) mix.dH_kcal = 0.0;
    if (!std::isfinite(mix.dS_cal_per_mol_K)) mix.dS_cal_per_mol_K = 0.0;
    if (!std::isfinite(mix.dG_kcal)) mix.dG_kcal = 0.0;
    return mix;
}

} // namespace natural
