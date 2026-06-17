// =============================================================================
// native_score.cpp — IdealPopulation native-pose CF scoring diagnostic
//
// Implementation of score_native_pose().  See native_score.h for full docs.
//
// Design notes:
//   • atoms[fa..la].coor[] at call time hold the BLINDED ligand pose (same
//     centroid as crystal, randomly rotated) — NOT the crystal pose.
//     DatasetRunner blinds the input SDF to prevent the GA seed from injecting
//     the crystal answer.  Scoring the blinded pose against the receptor causes
//     random ligand-receptor overlaps (30M+ wall energy for 1LPZ, etc.).
//
//   • Correct approach: override coor[] with the CRYSTAL pose coordinates
//     loaded from FLEXAIDDS_RMSDST (set by DatasetRunner to the unblinded
//     crystal SDF).  The crystal SDF is in the same PDB coordinate frame as
//     the receptor atoms → direct vcfunction() gives the physically meaningful
//     native CF.
//
//   • Why FLEXAIDDS_RMSDST was previously ignored: an earlier attempt loaded
//     raw SDF coordinates via an IC round-trip (buildic → nearest_grid_point
//     → ic2cf).  That round-trip has an axis-frame bug: ic2cf places GPA0 at
//     a cleft-grid point (≠ crystal position), which shifts the GPA0→virtual-
//     atom bond axis.  All subsequent atoms' dihedrals were computed relative
//     to the crystal GPA0 axis via buildic(), so reconstructing them against
//     the (slightly different) grid-point axis places every atom wrong →
//     intra-molecular clashes → 500M+ penalty.  Direct coordinate injection
//     (no IC round-trip) avoids this entirely.
//
//   • FA->ori is saved and restored as a precaution.
//
// Copyright 2026 Le Bonhomme Pharma.  Licensed under Apache-2.0.
// =============================================================================

#include "native_score.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <utility>

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

// Parse atom-block XYZ from an SDF file into coor_out[n_atoms * 3].
// Returns true on success; n_atoms must match the SDF atom count.
static bool load_crystal_coor_from_sdf(const char* sdf_path, int n_atoms,
                                        float* coor_out)
{
    FILE* f = fopen(sdf_path, "r");
    if (!f) return false;

    char buf[256];
    // Skip 3 header lines
    for (int i = 0; i < 3; ++i)
        if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }

    // Counts line: natoms nbonds ...
    if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }
    int sdf_natoms = 0, sdf_nbonds = 0;
    sscanf(buf, "%d %d", &sdf_natoms, &sdf_nbonds);
    if (sdf_natoms != n_atoms) {
        fprintf(stderr,
            "[NATIVE_CF] WARN: crystal SDF has %d atoms but ligand has %d — "
            "skipping crystal override\n", sdf_natoms, n_atoms);
        fclose(f);
        return false;
    }

    // Atom block: x y z element ...
    for (int i = 0; i < n_atoms; ++i) {
        if (!fgets(buf, sizeof(buf), f)) { fclose(f); return false; }
        float x = 0.f, y = 0.f, z = 0.f;
        if (sscanf(buf, "%f %f %f", &x, &y, &z) != 3) { fclose(f); return false; }
        coor_out[i * 3 + 0] = x;
        coor_out[i * 3 + 1] = y;
        coor_out[i * 3 + 2] = z;
    }
    fclose(f);
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

void score_native_pose(FA_Global* FA, VC_Global* VC, atom* atoms,
                       resid* residue, gridpoint* cleftgrid)
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

    // ── 0. [NATIVE-SEED-RMSD] diagnostic ──────────────────────────────────────
    //   Reconstructs the gen-0 native-seed pose via the SAME IC→Cartesian path
    //   the GA uses: ic2cf(FA->opt_par).  FA->opt_par holds the reference IC
    //   (orientation/dihedral genes = input-ligand atoms[].ang/.dih, grid gene
    //   = 0 ⇒ GPA0 snapped to cleftgrid[0]) — identical to what native_direct_seed
    //   writes into chrom[0].genes[].to_ic.  RMSD vs the crystal SDF therefore
    //   measures the pure IC-reconstruction floor of the native seed, with NO
    //   5° int32 quantization (the seed evaluates to_ic directly, not via
    //   genetoic).  Answers: is a residual native RMSD a reconstruction-frame
    //   floor (quantization-class) or GA displacement to a false CF minimum?
    {
        const char* rmsdst = getenv("FLEXAIDDS_RMSDST");
        if (rmsdst && rmsdst[0] != '\0') {
            std::vector<float> save(static_cast<size_t>(n_lig) * 3);
            for (int i = fa; i <= la; ++i) {
                const int k = (i - fa) * 3;
                save[static_cast<size_t>(k+0)] = atoms[i].coor[0];
                save[static_cast<size_t>(k+1)] = atoms[i].coor[1];
                save[static_cast<size_t>(k+2)] = atoms[i].coor[2];
            }
            const float ori0[3] = {FA->ori[0], FA->ori[1], FA->ori[2]};

            // Reconstruct native seed via the GA's IC→Cartesian path.
            ic2cf(FA, VC, atoms, residue, cleftgrid, FA->npar, FA->opt_par);

            std::vector<float> cry(static_cast<size_t>(n_lig) * 3);
            if (load_crystal_coor_from_sdf(rmsdst, n_lig, cry.data())) {
                double ss = 0.0;
                for (int i = fa; i <= la; ++i) {
                    const int k = (i - fa) * 3;
                    const double dx = atoms[i].coor[0] - cry[static_cast<size_t>(k+0)];
                    const double dy = atoms[i].coor[1] - cry[static_cast<size_t>(k+1)];
                    const double dz = atoms[i].coor[2] - cry[static_cast<size_t>(k+2)];
                    ss += dx*dx + dy*dy + dz*dz;
                }
                const double rmsd = std::sqrt(ss / static_cast<double>(n_lig));
                fprintf(stderr,
                    "[NATIVE-SEED-RMSD] %s round-trip RMSD = %.2f A "
                    "(%d genes, %d atoms, gene0=%.3f)\n",
                    rmsdst, rmsd, FA->npar, n_lig,
                    FA->npar > 0 ? FA->opt_par[0] : 0.0);

                // ── Fix 5: preflight gate (loud, non-aborting) ─────────────
                // A native-seed IC round-trip worse than 1.0 A means the oracle
                // seed cannot faithfully encode the crystal pose, so the whole
                // run's oracle mode is suspect. Make it impossible to miss in the
                // logs; LP decides whether to promote this to a hard abort.
                if (rmsd > 1.0) {
                    fprintf(stderr,
                        "[WARN] NATIVE-SEED-RMSD = %.2f A exceeds 1.0 A threshold.\n"
                        "       IC round-trip fidelity is poor — oracle seed may not "
                        "faithfully encode crystal pose.\n"
                        "       Check grid resolution, cleft centroid, and IC "
                        "parameter bounds before trusting this run.\n",
                        rmsd);
                }
            } else {
                fprintf(stderr,
                    "[NATIVE-SEED-RMSD] SKIP: could not load crystal SDF %s "
                    "(atom-count mismatch?)\n", rmsdst);
            }

            // Restore blinded coor[] and FA->ori for the CF scoring below.
            FA->ori[0] = ori0[0]; FA->ori[1] = ori0[1]; FA->ori[2] = ori0[2];
            for (int i = fa; i <= la; ++i) {
                const int k = (i - fa) * 3;
                atoms[i].coor[0] = save[static_cast<size_t>(k+0)];
                atoms[i].coor[1] = save[static_cast<size_t>(k+1)];
                atoms[i].coor[2] = save[static_cast<size_t>(k+2)];
            }
        }
    }

    // ── 1. Save ligand coor[] and FA->ori ─────────────────────────────────────
    std::vector<float> orig_coor(static_cast<size_t>(n_lig) * 3);
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        orig_coor[static_cast<size_t>(k+0)] = atoms[i].coor[0];
        orig_coor[static_cast<size_t>(k+1)] = atoms[i].coor[1];
        orig_coor[static_cast<size_t>(k+2)] = atoms[i].coor[2];
    }
    const float ori_save[3] = {FA->ori[0], FA->ori[1], FA->ori[2]};

    // ── 2. Load crystal pose from FLEXAIDDS_RMSDST ────────────────────────────
    //    DatasetRunner sets FLEXAIDDS_RMSDST to the unblinded crystal SDF.
    //    Override coor[] so vcfunction() scores the actual crystal pose.
    bool crystal_loaded = false;
    {
        const char* rmsdst = getenv("FLEXAIDDS_RMSDST");
        if (rmsdst && rmsdst[0] != '\0') {
            std::vector<float> crystal_coor(static_cast<size_t>(n_lig) * 3);
            if (load_crystal_coor_from_sdf(rmsdst, n_lig, crystal_coor.data())) {
                for (int i = fa; i <= la; ++i) {
                    const int k = (i - fa) * 3;
                    atoms[i].coor[0] = crystal_coor[static_cast<size_t>(k+0)];
                    atoms[i].coor[1] = crystal_coor[static_cast<size_t>(k+1)];
                    atoms[i].coor[2] = crystal_coor[static_cast<size_t>(k+2)];
                }
                crystal_loaded = true;
            }
        }
    }
    if (!crystal_loaded) {
        fprintf(stderr,
            "[NATIVE_CF] WARN: no crystal SDF loaded (FLEXAIDDS_RMSDST not set "
            "or atom-count mismatch) — scoring blinded pose\n");
    }

    // ── 3. Evaluate CF at crystal coordinates ─────────────────────────────────
    std::vector<std::pair<int,int>> intraclashes;
    bool error = false;
    const double penalty = vcfunction(FA, VC, atoms, residue, intraclashes, &error);

    cfstr cf{};
    if (error) {
        // Crystal pose caused an intra-clash — almost certainly a data-prep
        // issue (atom-count mismatch, wrong SDF source, unusual geometry).
        // Report as wal so DatasetRunner flags this run as problematic.
        cf.wal = penalty;
        fprintf(stderr,
            "[NATIVE_CF] WARN: vcfunction intra-clash (penalty=%.4f) "
            "— crystal SDF: %s\n",
            penalty, crystal_loaded ? getenv("FLEXAIDDS_RMSDST") : "(blinded pose)");

        // ── Diagnostic: identify which pairs cause the clash ──────────────────
        double permea = FA->permeability;
        int** rb = residue[lig_res].bonded;
        constexpr double KWALL_D = 1.0e6;
        double diag_total = 0.0;
        int diag_n = 0;
        fprintf(stderr, "[NATIVE_CF] CLASH_DIAG: n_lig=%d fatm=%d permeability=%.3f\n",
                n_lig, fa, permea);
        for (int ii = fa; ii <= la; ++ii) {
            for (int jj = ii+1; jj <= la; ++jj) {
                float dx = atoms[ii].coor[0]-atoms[jj].coor[0];
                float dy = atoms[ii].coor[1]-atoms[jj].coor[1];
                float dz = atoms[ii].coor[2]-atoms[jj].coor[2];
                double dist = std::sqrt(static_cast<double>(dx*dx+dy*dy+dz*dz));
                double classdist = permea*(atoms[ii].radius+atoms[jj].radius);
                if (dist < classdist) {
                    int ib = ii - fa, jb = jj - fa;
                    int bonded_val = (rb != NULL) ? rb[ib][jb] : -999;
                    double ew = KWALL_D*(std::pow(dist,-12.0)-std::pow(classdist,-12.0));
                    diag_total += ew;
                    if (diag_n < 20) {
                        fprintf(stderr,
                            "[NATIVE_CF] PAIR ii=%d(%s,t=%d,r=%.2f) jj=%d(%s,t=%d,r=%.2f)"
                            " dist=%.4f cd=%.4f bonded=%d ew=%.2f\n",
                            atoms[ii].number, atoms[ii].name, atoms[ii].type, atoms[ii].radius,
                            atoms[jj].number, atoms[jj].name, atoms[jj].type, atoms[jj].radius,
                            dist, classdist, bonded_val, ew);
                    }
                    ++diag_n;
                }
            }
        }
        fprintf(stderr, "[NATIVE_CF] CLASH_DIAG: %d intra-lig pairs within classdist, intra-lig wall=%.2f\n",
                diag_n, diag_total);

        // ── Ligand-receptor contact diagnostic ───────────────────────────────
        struct LRContact { int lig, rec; double dist, cd, ew; };
        std::vector<LRContact> worst;
        int lr_n = 0;
        double lr_total = 0.0;
        for (int ii = fa; ii <= la; ++ii) {
            for (int jj = 0; jj < fa; ++jj) {
                float dx = atoms[ii].coor[0] - atoms[jj].coor[0];
                float dy = atoms[ii].coor[1] - atoms[jj].coor[1];
                float dz = atoms[ii].coor[2] - atoms[jj].coor[2];
                double dist = std::sqrt(static_cast<double>(dx*dx + dy*dy + dz*dz));
                double cd = permea * (atoms[ii].radius + atoms[jj].radius);
                if (dist < cd) {
                    double ew = KWALL_D * (std::pow(dist, -12.0) - std::pow(cd, -12.0));
                    lr_total += ew;
                    ++lr_n;
                    worst.push_back({ii, jj, dist, cd, ew});
                }
            }
        }
        std::sort(worst.begin(), worst.end(),
            [](const LRContact& a, const LRContact& b){ return a.ew > b.ew; });
        int print_lr = std::min(static_cast<int>(worst.size()), 20);
        fprintf(stderr, "[NATIVE_CF] LIG-REC_DIAG: %d lig-rec pairs within classdist, lig-rec wall=%.2e\n",
                lr_n, lr_total);
        for (int k = 0; k < print_lr; ++k) {
            const auto& c = worst[static_cast<size_t>(k)];
            fprintf(stderr,
                "[NATIVE_CF] LIG-REC ii=%d(%s,t=%d,r=%.2f) jj=%d(%s,t=%d,r=%.2f)"
                " dist=%.4f cd=%.4f ew=%.2e\n",
                atoms[c.lig].number, atoms[c.lig].name, atoms[c.lig].type, atoms[c.lig].radius,
                atoms[c.rec].number, atoms[c.rec].name, atoms[c.rec].type, atoms[c.rec].radius,
                c.dist, c.cd, c.ew);
        }
    } else {
        for (int i = 0; i < FA->num_optres; ++i) {
            cf.com += FA->optres[i].cf.com;
            cf.wal += FA->optres[i].cf.wal;
            cf.sas += FA->optres[i].cf.sas;
            cf.con += FA->optres[i].cf.con;
            cf.hbond += FA->optres[i].cf.hbond;
        }
    }

    // ── 4. Restore FA->ori and ligand coor[] ─────────────────────────────────
    FA->ori[0] = ori_save[0];
    FA->ori[1] = ori_save[1];
    FA->ori[2] = ori_save[2];
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        atoms[i].coor[0] = orig_coor[static_cast<size_t>(k+0)];
        atoms[i].coor[1] = orig_coor[static_cast<size_t>(k+1)];
        atoms[i].coor[2] = orig_coor[static_cast<size_t>(k+2)];
    }

    // ── 5. Emit [NATIVE_CF] line for DatasetRunner parsing ───────────────────
    const double cf_total = get_cf_evalue(&cf);
    fprintf(stderr,
        "[NATIVE_CF] cf=%.6f breakdown=com:%.4f,wal:%.4f,sas:%.4f,con:%.4f,hbond:%.4f\n",
        cf_total,
        static_cast<double>(cf.com),
        static_cast<double>(cf.wal),
        static_cast<double>(cf.sas),
        static_cast<double>(cf.con),
        static_cast<double>(cf.hbond));
}
