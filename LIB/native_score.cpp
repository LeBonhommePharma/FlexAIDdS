// =============================================================================
// native_score.cpp — IdealPopulation native-pose CF scoring diagnostic
//
// Implementation of score_native_pose().  See native_score.h for full docs.
//
// Design notes:
//   • atoms[fa..la].coor[] are in FlexAID's PROCESSED internal frame (same
//     coordinate system as the cleft grid, MC_st0r matrix, and vcfunction).
//     Scoring directly on these coordinates produces the physically meaningful
//     CF value — the expected ≈−154 for 2GBP's glucose pose.
//
//   • FLEXAIDDS_RMSDST intentionally ignored for coordinate input (see header).
//     Raw SDF/PDB crystallographic coordinates are in the un-centred frame and
//     caused +500M steric-clash explosions in the previous implementation.
//
//   • IC round-trip (buildic → nearest_grid_point → ic2cf) is intentionally
//     avoided.  The IC reconstruction has an axis-frame bug: ic2cf places GPA0
//     at a cleft-grid point (≠ crystal position), which changes the direction of
//     the GPA0→virtual-atom bond axis.  All subsequent atoms' dihedrals were
//     computed relative to the crystal GPA0 axis via buildic(), so reconstructing
//     them against the (slightly different) grid-point axis places every atom
//     wrong → intra-molecular clashes → vcfunction() fires the penalty path →
//     wal = 520M instead of the expected ≈−154.
//
//   • Direct vcfunction() call is the correct approach.  ic2cf() itself does
//     exactly this after buildcc() updates coor[]; since coor[] already holds
//     the crystal pose in the right frame, buildcc() is not needed.
//
//   • FA->ori is saved and restored because ic2cf/buildcc can drift it.
//
// Copyright 2026 Le Bonhomme Pharma.  Licensed under Apache-2.0.
// =============================================================================

#include "native_score.h"

#include <cstdio>
#include <vector>
#include <utility>

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

void score_native_pose(FA_Global* FA, VC_Global* VC, atom* atoms,
                       resid* residue, gridpoint* /*cleftgrid*/)
{
    const int lig_res = FA->res_cnt;
    const int fa      = residue[lig_res].fatm[0];
    const int la      = residue[lig_res].latm[0];
    const int n_lig   = la - fa + 1;

    if (n_lig <= 0 || FA->npar <= 0) {
        fprintf(stderr, "[NATIVE_CF] SKIP: no ligand atoms (n_lig=%d npar=%d)\n",
                n_lig, FA->npar);
        return;
    }

    // ── 1. Save ligand coor[] and FA->ori ─────────────────────────────────────
    //    atoms[fa..la].coor[] are already in FlexAID's processed internal frame.
    //    Save them so we can restore after vcfunction() in case it has side-effects.
    std::vector<float> orig_coor(static_cast<size_t>(n_lig) * 3);
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        orig_coor[static_cast<size_t>(k+0)] = atoms[i].coor[0];
        orig_coor[static_cast<size_t>(k+1)] = atoms[i].coor[1];
        orig_coor[static_cast<size_t>(k+2)] = atoms[i].coor[2];
    }
    const float ori_save[3] = {FA->ori[0], FA->ori[1], FA->ori[2]};

    // ── 2. Evaluate CF directly at the pre-loaded crystal coordinates ─────────
    //    vcfunction() reads atoms[].coor[] and populates FA->optres[i].cf.
    //    ic2cf() does exactly this after buildcc(); since coor[] already holds
    //    the crystal pose in the correct frame, no buildcc() is needed.
    std::vector<std::pair<int,int>> intraclashes;
    bool error = false;
    const double penalty = vcfunction(FA, VC, atoms, residue, intraclashes, &error);

    cfstr cf{};
    if (error) {
        // Intra-molecular clash in the crystal pose — almost certainly an
        // atom-type assignment or topology issue (real crystal poses don't clash).
        // Report as wal so DatasetRunner sees a clearly-wrong large positive value.
        cf.wal = penalty;
        fprintf(stderr,
            "[NATIVE_CF] WARN: vcfunction intra-clash (penalty=%.4f) "
            "— check atom-type assignment or topology for this ligand\n",
            penalty);
    } else {
        for (int i = 0; i < FA->num_optres; ++i) {
            cf.com += FA->optres[i].cf.com;
            cf.wal += FA->optres[i].cf.wal;
            cf.sas += FA->optres[i].cf.sas;
            cf.con += FA->optres[i].cf.con;
        }
    }

    // ── 3. Restore FA->ori and ligand coor[] ─────────────────────────────────
    FA->ori[0] = ori_save[0];
    FA->ori[1] = ori_save[1];
    FA->ori[2] = ori_save[2];
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        atoms[i].coor[0] = orig_coor[static_cast<size_t>(k+0)];
        atoms[i].coor[1] = orig_coor[static_cast<size_t>(k+1)];
        atoms[i].coor[2] = orig_coor[static_cast<size_t>(k+2)];
    }

    // ── 4. Emit [NATIVE_CF] line for DatasetRunner parsing ───────────────────
    //    Machine-parseable format:
    //      [NATIVE_CF] cf=<total> breakdown=com:<v>,wal:<v>,sas:<v>,con:<v>
    const double cf_total = get_cf_evalue(&cf);
    fprintf(stderr,
        "[NATIVE_CF] cf=%.6f breakdown=com:%.4f,wal:%.4f,sas:%.4f,con:%.4f\n",
        cf_total,
        static_cast<double>(cf.com),
        static_cast<double>(cf.wal),
        static_cast<double>(cf.sas),
        static_cast<double>(cf.con));
}
