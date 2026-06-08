// =============================================================================
// native_score.cpp — IdealPopulation native-pose CF scoring diagnostic
//
// Implementation of score_native_pose().  See native_score.h for full docs.
//
// Design notes:
//   • Uses a LOCAL ic-vector (ref_icv), never overwrites FA->opt_par, so the
//     GA's initial population is untouched.
//   • Saves and restores atoms[fa..la].coor[] around the buildic() call so
//     the subsequent ic2cf(FA->opt_par) in top.cpp sees the original state.
//   • Translational gene (map_par[i].typ == -1): resolved by nearest-grid-
//     point search in Cartesian space (more robust than IC-distance matching
//     because Cartesian distance has no angle-wrap issues).
//   • SDF loading is by atom index (header gives natoms; lines 5…natoms+4
//     hold x y z).  PDB loading delegates to the existing read_rmsdst().
//
// Copyright 2026 Le Bonhomme Pharma.  Licensed under Apache-2.0.
// =============================================================================

#include "native_score.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <limits>
#include <vector>

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Read x/y/z from an SDF atom block into atoms[fa..fa+min(natoms,n_lig)-1].coor[].
/// Returns the number of atoms whose coordinates were successfully loaded.
static int load_sdf_coords(const char* sdf_path,
                            atom* atoms, int fa, int la)
{
    FILE* fp = fopen(sdf_path, "r");
    if (!fp) return 0;

    char buf[256];

    // SDF: line 1 = molecule name, line 2 = program/timestamp, line 3 = comment
    for (int h = 0; h < 3; ++h)
        if (!fgets(buf, (int)sizeof(buf), fp)) { fclose(fp); return 0; }

    // Counts line (V2000): "aaabbblllfffcccsssxxxrrrpppiiimmmvvvvvv"
    //   first token = number of atoms
    if (!fgets(buf, (int)sizeof(buf), fp)) { fclose(fp); return 0; }
    int natoms = 0;
    if (sscanf(buf, "%d", &natoms) != 1 || natoms <= 0) { fclose(fp); return 0; }

    const int n_lig   = la - fa + 1;
    const int n_read  = (natoms < n_lig) ? natoms : n_lig;
    int matched = 0;

    for (int ai = 0; ai < n_read; ++ai) {
        if (!fgets(buf, (int)sizeof(buf), fp)) break;
        float x = 0.0f, y = 0.0f, z = 0.0f;
        // Atom line: "   x.xxxx   y.yyyy   z.zzzz ..."
        if (sscanf(buf, "%f %f %f", &x, &y, &z) == 3) {
            atoms[fa + ai].coor[0] = x;
            atoms[fa + ai].coor[1] = y;
            atoms[fa + ai].coor[2] = z;
            ++matched;
        }
    }
    fclose(fp);
    return matched;
}

/// Return the index (1-based in cleftgrid) of the grid point nearest to (rx,ry,rz).
/// Starts at index 1 — index 0 is the pseudo-origin (FA->ori).
static int nearest_grid_point(const gridpoint* cleftgrid, int num_grd,
                               float rx, float ry, float rz)
{
    int   best   = 1;
    float best_d2 = std::numeric_limits<float>::max();

    for (int g = 1; g < num_grd; ++g) {
        float dx = cleftgrid[g].coor[0] - rx;
        float dy = cleftgrid[g].coor[1] - ry;
        float dz = cleftgrid[g].coor[2] - rz;
        float d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < best_d2) { best_d2 = d2; best = g; }
    }
    return best;
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

    // ── 1. Save original ligand coor[] ───────────────────────────────────────
    std::vector<float> orig_coor(static_cast<size_t>(n_lig) * 3);
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        orig_coor[static_cast<size_t>(k+0)] = atoms[i].coor[0];
        orig_coor[static_cast<size_t>(k+1)] = atoms[i].coor[1];
        orig_coor[static_cast<size_t>(k+2)] = atoms[i].coor[2];
    }

    // ── 2. Optionally load reference coordinates from FLEXAIDDS_RMSDST ───────
    // If not set: use pre-loaded coor[] (crystal coords when DatasetRunner
    // passes the unblinded crystal SDF as the ligand for the native-score pass).
    const char* rmsdst_path = std::getenv("FLEXAIDDS_RMSDST");
    bool ref_overwritten = false;

    if (rmsdst_path && rmsdst_path[0] != '\0') {
        // Detect format by extension (.sdf / .mol → SDF; everything else → PDB)
        const char* ext = std::strrchr(rmsdst_path, '.');
        const bool is_sdf = ext &&
            (strcmp(ext, ".sdf") == 0 || strcmp(ext, ".SDF") == 0 ||
             strcmp(ext, ".mol") == 0 || strcmp(ext, ".MOL") == 0);

        if (is_sdf) {
            const int n = load_sdf_coords(rmsdst_path, atoms, fa, la);
            if (n > 0) {
                ref_overwritten = true;
            } else {
                fprintf(stderr,
                    "[NATIVE_CF] WARNING: SDF loader read 0 atoms from '%s'"
                    " — using pre-loaded coor[]\n", rmsdst_path);
            }
        } else {
            // PDB / mol2: delegate to read_rmsdst(), then copy coor_ref→coor
            const int n = read_rmsdst(FA, atoms, residue,
                                      const_cast<char*>(rmsdst_path));
            if (n > 0) {
                for (int i = fa; i <= la; ++i) {
                    if (atoms[i].coor_ref) {
                        atoms[i].coor[0] = atoms[i].coor_ref[0];
                        atoms[i].coor[1] = atoms[i].coor_ref[1];
                        atoms[i].coor[2] = atoms[i].coor_ref[2];
                    }
                }
                ref_overwritten = true;
            } else {
                fprintf(stderr,
                    "[NATIVE_CF] WARNING: read_rmsdst('%s') matched 0 atoms"
                    " — using pre-loaded coor[]\n", rmsdst_path);
            }
        }
    }
    (void)ref_overwritten; // only used for the warning path

    // ── 3. Compute ICs from reference coor[] ─────────────────────────────────
    // buildic() writes atoms[i].dis / .ang / .dih from atoms[i].coor[].
    // It operates only on atoms whose .recs == 'm' (moveable ligand atoms).
    buildic(FA, atoms, residue, lig_res);

    // ── 4. Build IC vector for the native pose ───────────────────────────────
    std::vector<double> ref_icv(static_cast<size_t>(FA->npar), 0.0);
    for (int i = 0; i < FA->npar; ++i) {
        const int typ = FA->map_par[i].typ;
        const int atm = FA->map_par[i].atm;

        if (typ == -1) {
            // Translational gene: grid index — find nearest cleft-grid point
            // to the reference GPA0 position in Cartesian space.
            ref_icv[static_cast<size_t>(i)] = static_cast<double>(
                nearest_grid_point(cleftgrid, FA->num_grd,
                                   atoms[atm].coor[0],
                                   atoms[atm].coor[1],
                                   atoms[atm].coor[2]));
        } else if (typ == 0) {
            ref_icv[static_cast<size_t>(i)] = static_cast<double>(atoms[atm].dis);
        } else if (typ == 1) {
            ref_icv[static_cast<size_t>(i)] = static_cast<double>(atoms[atm].ang);
        } else if (typ == 2) {
            ref_icv[static_cast<size_t>(i)] = static_cast<double>(atoms[atm].dih);
        }
        // typ 3 (normal mode) and typ 4 (rotamer side-chain): leave as 0.0.
        // Rotamers at 0 → first rotamer; normal-mode amplitude 0 → no perturbation.
    }

    // ── 5. Restore original coor[] ───────────────────────────────────────────
    // ic2cf() will rebuild coor[] from ref_icv internally, so we restore here
    // so that the subsequent ic2cf(FA->opt_par) call in top.cpp (initial-pose
    // evaluation for the INI.pdb) starts from the correct blinded/input state.
    for (int i = fa; i <= la; ++i) {
        const int k = (i - fa) * 3;
        atoms[i].coor[0] = orig_coor[static_cast<size_t>(k+0)];
        atoms[i].coor[1] = orig_coor[static_cast<size_t>(k+1)];
        atoms[i].coor[2] = orig_coor[static_cast<size_t>(k+2)];
    }

    // ── 6. Evaluate CF at native IC vector ───────────────────────────────────
    // Passes ref_icv (NOT FA->opt_par) so the GA's initial state is untouched.
    cfstr cf = ic2cf(FA, VC, atoms, residue, cleftgrid,
                     FA->npar, ref_icv.data());
    const double cf_total = get_cf_evalue(&cf);

    // ── 7. Emit [NATIVE_CF] line for DatasetRunner parsing ───────────────────
    // Format is intentionally machine-parseable:
    //   [NATIVE_CF] cf=<total> breakdown=com:<v>,wal:<v>,sas:<v>,con:<v>
    fprintf(stderr,
        "[NATIVE_CF] cf=%.6f breakdown=com:%.4f,wal:%.4f,sas:%.4f,con:%.4f\n",
        cf_total,
        static_cast<double>(cf.com),
        static_cast<double>(cf.wal),
        static_cast<double>(cf.sas),
        static_cast<double>(cf.con));
}
