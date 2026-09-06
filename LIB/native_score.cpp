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
#include "ensemble_pipeline.h"
#include "fileio.h"  // Terminate() for strict frame-chart gate
// Atom-type-pair contact-surface vector. Gate FLEXAIDDS_CONTACT_PROFILE,
// DEFAULT OFF: unset, nothing is snapshotted and no sidecar is written, so
// the [NATIVE_CF] diagnostic is byte-identical to HEAD.
#include "contact_profile.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>
#include <cstdlib>
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

                // Layer 1 frame-chart gate (ensemble_pipeline): identity-class
                // Cartesian⇄gene chart. Soft warn at 1.0 A; hard fail at 0.1 A
                // only when FLEXAIDDS_FRAME_CHART_STRICT=1 (CI / product gate).
                {
                    const bool strict = ensemble::frame_chart_strict_enabled();
                    const auto st = ensemble::classify_frame_chart_rmsd(rmsd, strict);
                    const char* st_s =
                        (st == ensemble::FrameChartStatus::Ok)   ? "ok" :
                        (st == ensemble::FrameChartStatus::Warn) ? "warn" : "fail";
                    fprintf(stderr,
                        "[FRAME_CHART] status=%s rmsd=%.3f strict=%d "
                        "warn_A=%.1f strict_A=%.1f\n",
                        st_s, rmsd, strict ? 1 : 0,
                        ensemble::kFrameChartWarnRmsdA,
                        ensemble::kFrameChartStrictRmsdA);
                    if (st == ensemble::FrameChartStatus::Warn) {
                        fprintf(stderr,
                            "[WARN] NATIVE-SEED-RMSD = %.2f A exceeds %.1f A threshold.\n"
                            "       IC round-trip fidelity is poor — oracle seed may not "
                            "faithfully encode crystal pose.\n"
                            "       Check grid resolution, cleft centroid, and IC "
                            "parameter bounds before trusting this run.\n",
                            rmsd, ensemble::kFrameChartWarnRmsdA);
                    }
                    if (st == ensemble::FrameChartStatus::Fail) {
                        fprintf(stderr,
                            "[FRAME_CHART] FATAL: native-seed RMSD %.3f A > strict "
                            "%.1f A (FLEXAIDDS_FRAME_CHART_STRICT). "
                            "Gene chart cannot represent the crystal pose.\n",
                            rmsd, ensemble::kFrameChartStrictRmsdA);
                        Terminate(2);
                    }
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

    // ── FLEXAIDDS_CONTACT_PROFILE (default OFF) ──────────────────────────────
    // Freeze the NATIVE complex's atom-type-pair contact-surface vector right
    // after the only vcfunction() call in this function, so nothing downstream
    // can substitute a different profile. The sidecar is written at the end of
    // the function, once cf_total is known. THIS PROFILE IS THE ORACLE: any
    // comparison of a docked pose against it consumes the answer and is valid
    // for benchmark analysis and diagnosis ONLY — never as a scoring or
    // selection term. See scripts/contact_profile_tanimoto.py.
    if (flexaids::contact_profile::enabled()) {
        flexaids::contact_profile::snapshot(FA->contributions, FA->ntypes);
    }

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
            cf.elec += FA->optres[i].cf.elec;
            cf.gist_desolv += FA->optres[i].cf.gist_desolv;
            cf.metal_coord += FA->optres[i].cf.metal_coord;
            cf.hbond += FA->optres[i].cf.hbond;
            cf.entropy += FA->optres[i].cf.entropy;
            cf.pb_clash += FA->optres[i].cf.pb_clash;
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
    const double cf_total = get_cf_evalue(&cf, FA);
    fprintf(stderr,
        "[NATIVE_CF] cf=%.6f breakdown=com:%.4f,wal:%.4f,sas:%.4f,con:%.4f,"
        "elec:%.4f,hbond:%.4f,gist_desolv:%.4f,metal_coord:%.4f,entropy:%.4f,"
        "pb_clash:%.4f\n",
        cf_total,
        static_cast<double>(cf.com),
        static_cast<double>(cf.wal),
        static_cast<double>(cf.sas),
        static_cast<double>(cf.con),
        static_cast<double>(cf.elec),
        static_cast<double>(cf.hbond),
        static_cast<double>(cf.gist_desolv),
        static_cast<double>(cf.metal_coord),
        static_cast<double>(cf.entropy),
        static_cast<double>(cf.pb_clash));

    // ── Write the ORACLE contact-profile sidecar (default OFF) ──────────────
    // Path is derived from FA->rrgfile, which top.cpp sets to end_strfile — the
    // same output prefix cluster.cpp builds <prefix>_<j>.pdb from — so the
    // native profile lands beside the pose profiles of the very same run and
    // joins to them by directory. No new environment variable: the single gate
    // FLEXAIDDS_CONTACT_PROFILE governs this too. Reached only when
    // FLEXAIDDS_SCORE_NATIVE is also on, since that is what calls this function.
    if (flexaids::contact_profile::enabled()) {
        const std::string prefix = (FA->rrgfile[0] != '\0')
                                 ? std::string(FA->rrgfile)
                                 : std::string("flexaid");
        const std::string cp_path = prefix + "_native.cprof.csv";
        const bool cp_ok = flexaids::contact_profile::write_csv(
            cp_path, "native_crystal_pose", cp_path, -1, "total", cf_total);
        if (cp_ok) {
            std::fprintf(stderr,
                "[CPROF] wrote ORACLE native profile %s (ntypes=%d, %lld contacts)"
                " -- oracle metric, benchmark analysis only, never a scoring term\n",
                cp_path.c_str(),
                flexaids::contact_profile::last().ntypes,
                flexaids::contact_profile::last().ncontacts);
        } else {
            std::fprintf(stderr,
                "[CPROF] WARNING: no native profile written to %s (ntypes=%d)\n",
                cp_path.c_str(),
                flexaids::contact_profile::last().ntypes);
        }
        std::fflush(stderr);
    }
}


// =============================================================================
// DUMP_POP refstructure loader — audit RMSD without pose seeding
// =============================================================================
// cluster.cpp DUMP_POP (.pop.tsv) and write_rrd require FA->refstructure==1 and
// atoms[].coor_ref. Classic RMSDST (read_rmsdst PDB match) is not set on the
// modern JSON/direct path. This loader fills coor_ref from the crystal SDF
// (same atom order as the docked ligand) when FLEXAIDDS_DUMP_POP is on.
// Does NOT modify atoms[].coor (blinded GA coords stay blinded).
bool load_dump_pop_refstructure(FA_Global* FA, atom* atoms, resid* residue)
{
    if (FA == nullptr || atoms == nullptr || residue == nullptr)
        return false;
    if (FA->refstructure == 1)
        return true;

    const char* dump_env = std::getenv("FLEXAIDDS_DUMP_POP");
    if (!dump_env || dump_env[0] == '\0' || std::atoi(dump_env) == 0)
        return false;

    const char* rmsdst = std::getenv("FLEXAIDDS_RMSDST");
    if (!rmsdst || rmsdst[0] == '\0') {
        std::fprintf(stderr,
            "[DUMP_POP] WARN: FLEXAIDDS_DUMP_POP set but FLEXAIDDS_RMSDST empty "
            "— refstructure stays 0; .pop.tsv will not be written\n");
        return false;
    }

    const int lig_res = FA->res_cnt;
    if (lig_res < 1)
        return false;
    const int fa = residue[lig_res].fatm[0];
    const int la = residue[lig_res].latm[0];
    const int n_lig = la - fa + 1;
    if (n_lig <= 0) {
        std::fprintf(stderr, "[DUMP_POP] WARN: no ligand atoms for refstructure\n");
        return false;
    }

    std::vector<float> crystal(static_cast<size_t>(n_lig) * 3u);
    if (!load_crystal_coor_from_sdf(rmsdst, n_lig, crystal.data())) {
        std::fprintf(stderr,
            "[DUMP_POP] WARN: could not load crystal SDF %s "
            "(open fail or atom-count mismatch vs n_lig=%d)\n",
            rmsdst, n_lig);
        return false;
    }

    int loaded = 0;
    for (int i = fa; i <= la; ++i) {
        if (atoms[i].coor_ref == nullptr) {
            atoms[i].coor_ref = static_cast<float*>(std::malloc(3 * sizeof(float)));
            if (atoms[i].coor_ref == nullptr) {
                std::fprintf(stderr,
                    "[DUMP_POP] ERROR: malloc coor_ref failed at atom %d\n", i);
                return false;
            }
        }
        const int k = (i - fa) * 3;
        atoms[i].coor_ref[0] = crystal[static_cast<size_t>(k + 0)];
        atoms[i].coor_ref[1] = crystal[static_cast<size_t>(k + 1)];
        atoms[i].coor_ref[2] = crystal[static_cast<size_t>(k + 2)];
        ++loaded;
    }

    FA->refstructure = 1;
    std::fprintf(stdout,
        "[DUMP_POP] refstructure=1 from FLEXAIDDS_RMSDST=%s (%d ligand atoms; "
        "audit only — no GA seed)\n",
        rmsdst, loaded);
    return true;
}
