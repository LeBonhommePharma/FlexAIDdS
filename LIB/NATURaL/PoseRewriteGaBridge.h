// PoseRewriteGaBridge.h — opt-in DualAssembly / NATURAL GA-coord plumbing
//
// Wires elected/GA Cartesian snapshots into PoseHelixThermoRewrite and
// PoseLocalThermoRewrite. Both rewrites stay DEFAULT OFF. Missing coords
// fail closed (no fabricated contacts, no mutation of validated DualAssembly
// ΔG). This header has no FlexAID GA / FA_Global dependency.
//
// Flags: DualAssemblyConfig::enable_pose_helix_rewrite /
//        DualAssemblyConfig::enable_pose_local_thermo_rewrite
//        NATURaLConfig::enable_pose_helix_rewrite /
//        NATURaLConfig::pose_local_thermo_rewrite
// Env:   FLEXAIDDS_POSE_HELIX_THERMO_REWRITE=1
//        FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE=1
//
// Copyright 2026 Le Bonhomme Pharma. SPDX-License-Identifier: Apache-2.0
#pragma once

#include "PoseHelixThermoRewrite.h"
#include "PoseLocalThermoRewrite.h"

#include <cmath>
#include <vector>

namespace natural {

// True when a GA backend emitted at least one receptor site and one ligand
// pose with Cartesian atoms. Empty snapshots are treated as missing coords.
inline bool ga_atom_coords_present(
    const std::vector<ReceptorNtCoord>& nts,
    const std::vector<LigandPose>& poses) noexcept
{
    if (nts.empty() || poses.empty()) return false;
    for (const LigandPose& p : poses) {
        if (!p.atoms.empty()) return true;
    }
    return false;
}

namespace pose_rewrite_ga_detail {

inline bool finite3(double x, double y, double z) noexcept
{
    return std::isfinite(x) && std::isfinite(y) && std::isfinite(z);
}

inline double dist_A(const ReceptorNtCoord& nt, const LigandAtomCoord& lig) noexcept
{
    const double dx = nt.x - lig.x;
    const double dy = nt.y - lig.y;
    const double dz = nt.z - lig.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

inline double vec_norm(double x, double y, double z) noexcept
{
    return std::sqrt(x * x + y * y + z * z);
}

// Folded into [0, 90] so parallel and antiparallel stacks both count.
// Degenerate normals → 0° (distance-only), matching PoseHelixThermoRewrite.
inline double stack_angle_deg(const ReceptorNtCoord& nt,
                              const LigandAtomCoord& lig) noexcept
{
    constexpr double kPi = 3.14159265358979323846;
    const double n1 = vec_norm(nt.nx, nt.ny, nt.nz);
    const double n2 = vec_norm(lig.nx, lig.ny, lig.nz);
    if (n1 < 1e-12 || n2 < 1e-12) return 0.0;
    double cosang = (nt.nx * lig.nx + nt.ny * lig.ny + nt.nz * lig.nz) / (n1 * n2);
    if (cosang > 1.0) cosang = 1.0;
    if (cosang < -1.0) cosang = -1.0;
    return std::acos(std::fabs(cosang)) * 180.0 / kPi;
}

inline bool index_hits_element(const ReceptorNtCoord& nt,
                               const DecisionElement& el) noexcept
{
    if (nt.nt_index < 0) return false;
    switch (el.ss) {
    case SSClass::RNA_STEM_LOOP:
    case SSClass::DNA_STEM_LOOP:
        return el.na_start >= 0 && el.na_end >= el.na_start &&
               nt.nt_index >= el.na_start && nt.nt_index <= el.na_end;
    case SSClass::PROTEIN_HELIX:
    case SSClass::PROTEIN_SHEET:
    case SSClass::PROTEIN_OTHER:
        return el.res_start >= 0 && el.res_end >= el.res_start &&
               nt.nt_index >= el.res_start && nt.nt_index <= el.res_end;
    }
    return false;
}

inline ContactKind stack_kind_for(SSClass ss) noexcept
{
    switch (ss) {
    case SSClass::RNA_STEM_LOOP: return ContactKind::RNA_STEM_STACK;
    case SSClass::DNA_STEM_LOOP: return ContactKind::DNA_STEM_STACK;
    case SSClass::PROTEIN_HELIX: return ContactKind::PROTEIN_HELIX_PACKING;
    case SSClass::PROTEIN_SHEET: return ContactKind::PROTEIN_SHEET_HYDROPHOBIC_PACK;
    case SSClass::PROTEIN_OTHER: break;
    }
    return ContactKind::COUNT;
}

inline ContactKind hbond_kind_for(SSClass ss) noexcept
{
    switch (ss) {
    case SSClass::RNA_STEM_LOOP: return ContactKind::RNA_LOOP_HBOND;
    case SSClass::DNA_STEM_LOOP: return ContactKind::DNA_LOOP_HBOND;
    case SSClass::PROTEIN_HELIX: return ContactKind::PROTEIN_HELIX_BB_HBOND;
    case SSClass::PROTEIN_SHEET: return ContactKind::PROTEIN_SHEET_BRIDGE_HBOND;
    case SSClass::PROTEIN_OTHER: break;
    }
    return ContactKind::COUNT;
}

inline PoseContact make_contact(ContactKind kind,
                                const DecisionElement& el,
                                const ReceptorNtCoord& nt)
{
    PoseContact c;
    c.kind = kind;
    c.geometry_weight = 1.0;
    c.element_label = el.label;
    if (el.ss == SSClass::RNA_STEM_LOOP || el.ss == SSClass::DNA_STEM_LOOP)
        c.nt_index = nt.nt_index;
    else
        c.res_index = nt.nt_index;
    return c;
}

} // namespace pose_rewrite_ga_detail

// Map GA atom coords onto PoseLocal PoseView contacts using PoseHelix
// geometry cutoffs (H-bond / aromatic stack). Empty input → empty output.
inline std::vector<PoseView> pose_views_from_ga_atom_coords(
    const std::vector<ReceptorNtCoord>& nts,
    const std::vector<LigandPose>& poses,
    const std::vector<DecisionElement>& elements,
    const PoseHelixRewriteConfig& geom = {})
{
    using namespace pose_rewrite_ga_detail;
    std::vector<PoseView> views;
    if (nts.empty() || poses.empty() || elements.empty()) return views;

    for (const LigandPose& pose : poses) {
        if (pose.atoms.empty()) continue;
        PoseView view;
        view.score_kcal = pose.score_kcal;
        for (const DecisionElement& el : elements) {
            if (el.ss == SSClass::PROTEIN_OTHER) continue;
            const ContactKind stack_k = stack_kind_for(el.ss);
            const ContactKind hbond_k = hbond_kind_for(el.ss);
            for (const ReceptorNtCoord& nt : nts) {
                if (!finite3(nt.x, nt.y, nt.z)) continue;
                if (!index_hits_element(nt, el)) continue;
                bool stack = false;
                bool hbond = false;
                for (const LigandAtomCoord& lig : pose.atoms) {
                    if (!finite3(lig.x, lig.y, lig.z)) continue;
                    const double d = dist_A(nt, lig);
                    if (nt.is_hbond_partner && lig.is_hbond_partner &&
                        d <= geom.hbond_max_A) {
                        hbond = true;
                    }
                    if (nt.is_aromatic && lig.is_aromatic &&
                        d <= geom.stack_max_A &&
                        stack_angle_deg(nt, lig) <= geom.stack_angle_max_deg) {
                        stack = true;
                    }
                }
                if (stack && stack_k != ContactKind::COUNT)
                    view.contacts.push_back(make_contact(stack_k, el, nt));
                if (hbond && hbond_k != ContactKind::COUNT)
                    view.contacts.push_back(make_contact(hbond_k, el, nt));
            }
        }
        if (!view.contacts.empty()) views.push_back(std::move(view));
    }
    return views;
}

// PoseHelix sidecar: empty + applied=false unless flag ON, annotated helices,
// and GA atom coords are present. Never invents coords.
inline PoseHelixRewriteResult pose_helix_from_ga_or_closed(
    bool enabled,
    const std::vector<HelixSegment>& helices,
    const std::vector<ReceptorNtCoord>& nts,
    const std::vector<LigandPose>& poses,
    const PoseHelixRewriteConfig& cfg = {})
{
    PoseHelixRewriteResult empty;
    empty.applied = false;
    if (!enabled) return empty;
    if (helices.empty() || !ga_atom_coords_present(nts, poses)) return empty;
    PoseHelixRewriteResult r = rewrite_helices_from_poses(helices, nts, poses, cfg);
    r.applied = true;
    return r;
}

inline void accumulate_helix_ensemble(const PoseHelixRewriteResult& r,
                                      double T_K,
                                      double& dH_kcal,
                                      double& dS_cal_per_mol_K,
                                      double& dG_kcal) noexcept
{
    dH_kcal = 0.0;
    dS_cal_per_mol_K = 0.0;
    dG_kcal = 0.0;
    if (!r.applied) return;
    for (const HelixThermoPatch& e : r.ensemble) {
        dH_kcal += e.delta_H_kcal;
        dS_cal_per_mol_K += e.delta_S_cal_per_K;
    }
    dG_kcal = helix_thermo_delta_g_kcal(dH_kcal, dS_cal_per_mol_K, T_K);
}

// PoseLocal sidecar from GA coords (preferred) or labelled PoseView fallback.
// Flag OFF or missing both sources → empty mixture (fail closed).
inline LocalThermoMixture pose_local_from_ga_or_closed(
    bool enabled,
    const std::vector<DecisionElement>& elements,
    const std::vector<PoseView>& labelled_poses,
    const std::vector<ReceptorNtCoord>& nts,
    const std::vector<LigandPose>& poses,
    const PoseThermoRewriteConfig& cfg,
    const PoseHelixRewriteConfig& geom = {})
{
    LocalThermoMixture empty;
    if (!enabled) return empty;
    if (elements.empty()) return empty;

    std::vector<PoseView> views;
    if (ga_atom_coords_present(nts, poses)) {
        views = pose_views_from_ga_atom_coords(nts, poses, elements, geom);
    } else {
        views = labelled_poses;
    }
    if (views.empty()) return empty;
    return rewrite_local_thermo_from_poses(views, elements, cfg);
}

} // namespace natural
