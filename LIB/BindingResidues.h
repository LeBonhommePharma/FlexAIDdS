// BindingResidues.h — Identify key binding-site residues from MIF scores
//
// Given MIF energies on a cleftgrid, finds the protein residues that contribute
// most to favorable binding interactions. Uses SpatialGrid for fast atom lookup.
//
// Header-only, C++20, Apache-2.0 © 2026 Le Bonhomme Pharma

#pragma once

#include "flexaid.h"
#include "metal_coordination.h"   // metal_coord::is_metal_type / is_coord_donor_type
#include "MIFGrid.h"
#include "CavityDetect/SpatialGrid.h"

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <numeric>

namespace binding_residues {

// ── Result: one residue's contribution to binding ──────────────────────────

struct ResidueContribution {
    int    res_index;        // internal residue index (into resid[])
    char   name[4];          // 3-letter residue name (e.g. "ASP")
    int    number;           // PDB residue number
    char   chain;            // chain ID
    float  mif_score;        // summed MIF energy from nearby grid points (more negative = more favorable)
    int    contact_count;    // number of favorable grid points near this residue
    float  min_distance;     // closest distance from any residue atom to a favorable grid point
};

// ── Identify key binding-site residues ─────────────────────────────────────
//
// Algorithm:
//   1. Select top-K% most favorable grid points by MIF energy
//   2. For each favorable grid point, find nearby protein atoms (within cutoff)
//   3. Map atoms → residues, accumulate MIF score per residue
//   4. Sort residues by total MIF score (most favorable first)
//
// Returns a sorted vector of ResidueContribution (most favorable first).

inline std::vector<ResidueContribution> identify_key_residues(
    const gridpoint* cleftgrid, int num_grd,
    const float* mif_energies,
    const atom* atoms, int atm_cnt,
    const resid* residues,
    const cavity_detect::SpatialGrid& spatial_grid,
    float top_k_percent = 30.0f,
    float contact_cutoff = 4.5f)
{
    if (num_grd <= 1 || !mif_energies) return {};

    // Step 1: Find favorable grid points (top-K% by MIF energy)
    const int n_points = num_grd - 1;  // index 0 is unused
    std::vector<int> sorted_indices(static_cast<size_t>(n_points));
    std::iota(sorted_indices.begin(), sorted_indices.end(), 1);

    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&](int a, int b) { return mif_energies[a] < mif_energies[b]; });

    int k = std::max(1, static_cast<int>(
        static_cast<float>(n_points) * top_k_percent / 100.0f));
    if (k > n_points) k = n_points;

    // Step 2 & 3: For each favorable grid point, find nearby atoms → residues
    // Use a map indexed by residue index for accumulation
    struct ResAccum {
        float  score = 0.0f;
        int    count = 0;
        float  min_dist = 1e10f;
    };
    std::vector<ResAccum> res_accum(static_cast<size_t>(atm_cnt > 0 ? atm_cnt : 1));

    // Track which residue indices we've seen
    std::vector<bool> res_seen(static_cast<size_t>(atm_cnt > 0 ? atm_cnt : 1), false);
    int max_res_idx = 0;

    const float cutoff_sq = contact_cutoff * contact_cutoff;

    for (int ki = 0; ki < k; ++ki) {
        int gp = sorted_indices[static_cast<size_t>(ki)];
        float e = mif_energies[gp];
        if (e >= 0.0f) continue;  // only consider favorable (negative) energies

        float gx = cleftgrid[gp].coor[0];
        float gy = cleftgrid[gp].coor[1];
        float gz = cleftgrid[gp].coor[2];

        // Query nearby atoms using SpatialGrid (returns indices within cell_size range)
        float qcoord[3] = {gx, gy, gz};
        auto neighbors = spatial_grid.query_neighbors(qcoord);

        for (std::size_t nb_idx : neighbors) {
            int atom_idx = static_cast<int>(nb_idx);
            if (atom_idx < 0 || atom_idx >= atm_cnt) continue;

            float dx = atoms[atom_idx].coor[0] - gx;
            float dy = atoms[atom_idx].coor[1] - gy;
            float dz = atoms[atom_idx].coor[2] - gz;
            float dist_sq = dx * dx + dy * dy + dz * dz;

            if (dist_sq > cutoff_sq) continue;

            int ri = atoms[atom_idx].ofres;
            if (ri < 0 || ri >= atm_cnt) continue;

            auto uri = static_cast<size_t>(ri);
            if (uri >= res_accum.size()) {
                res_accum.resize(uri + 1);
                res_seen.resize(uri + 1, false);
            }

            float dist = std::sqrt(dist_sq);
            // Weight contribution by inverse distance (closer atoms contribute more)
            float weight = 1.0f / (1.0f + dist);
            res_accum[uri].score += e * weight;
            res_accum[uri].count += 1;
            if (dist < res_accum[uri].min_dist) {
                res_accum[uri].min_dist = dist;
            }
            res_seen[uri] = true;
            if (ri > max_res_idx) max_res_idx = ri;
        }
    }

    // Step 4: Collect and sort results
    std::vector<ResidueContribution> results;
    for (int ri = 0; ri <= max_res_idx; ++ri) {
        auto uri = static_cast<size_t>(ri);
        if (uri >= res_seen.size() || !res_seen[uri]) continue;

        ResidueContribution rc{};
        rc.res_index = ri;
        rc.mif_score = res_accum[uri].score;
        rc.contact_count = res_accum[uri].count;
        rc.min_distance = res_accum[uri].min_dist;

        // Copy residue info
        std::strncpy(rc.name, residues[ri].name, 3);
        rc.name[3] = '\0';
        rc.number = residues[ri].number;
        rc.chain = residues[ri].chn;

        results.push_back(rc);
    }

    // Sort by MIF score (most negative = most favorable = first)
    std::sort(results.begin(), results.end(),
              [](const ResidueContribution& a, const ResidueContribution& b) {
                  return a.mif_score < b.mif_score;
              });

    return results;
}

// ── Convenience: identify from FA_Global (uses pre-computed MIF) ───────────

inline std::vector<ResidueContribution> identify_key_residues_from_fa(
    const FA_Global* FA,
    const gridpoint* cleftgrid,
    const atom* atoms,
    const resid* residues,
    float top_k_percent = 30.0f,
    float contact_cutoff = 4.5f)
{
    if (!FA->mif_energies || FA->mif_count == 0) return {};

    std::vector<atom> protein_atoms(atoms, atoms + FA->atm_cnt_real);
    cavity_detect::SpatialGrid sg;
    sg.build(protein_atoms);

    return identify_key_residues(
        cleftgrid, FA->num_grd, FA->mif_energies,
        atoms, FA->atm_cnt_real, residues, sg,
        top_k_percent, contact_cutoff);
}

// ── Print summary to stdout ────────────────────────────────────────────────

inline void print_key_residues(const std::vector<ResidueContribution>& residues,
                                int max_display = 20) {
    printf("─── Key Binding-Site Residues (by MIF score) ───\n");
    printf("%-5s %-4s %5s %6s  %8s  %5s  %7s\n",
           "Rank", "Name", "Num", "Chain", "MIF_Score", "Contacts", "MinDist");

    int n = std::min(static_cast<int>(residues.size()), max_display);
    for (int i = 0; i < n; ++i) {
        const auto& r = residues[static_cast<size_t>(i)];
        printf("%-5d %-4s %5d %4c    %8.2f  %5d   %6.2f\n",
               i + 1, r.name, r.number, r.chain,
               r.mif_score, r.contact_count, r.min_distance);
    }
    if (static_cast<int>(residues.size()) > max_display) {
        printf("... and %d more residues\n",
               static_cast<int>(residues.size()) - max_display);
    }
}

// ── Auto-flex: add key binding residues to FA->flex_res[] ──────────────────
//
// Adds the top-N most favorable binding residues as flexible side-chains.
// Skips GLY, ALA, PRO (no rotameric freedom), ligand residues, and
// residues already in flex_res[]. Returns count of residues added.
//
// Call AFTER compute_mif_and_reflig() and BEFORE build_rotamers()/add2_optimiz_vec("SC").

inline int add_key_residues_as_flexible(
    FA_Global* FA,
    const gridpoint* cleftgrid,
    const atom* atoms,
    const resid* residues,
    int max_auto_flex = 5,
    float top_k_percent = 30.0f,
    float contact_cutoff = 4.5f)
{
    if (!FA->mif_energies || FA->mif_count == 0 || max_auto_flex <= 0) return 0;

    // Identify key residues
    auto key_res = identify_key_residues_from_fa(
        FA, cleftgrid, atoms, residues, top_k_percent, contact_cutoff);

    if (key_res.empty()) return 0;

    // Residues to skip (no side-chain rotamers)
    auto is_skip_residue = [](const char* name) {
        return std::strcmp(name, "GLY") == 0 ||
               std::strcmp(name, "ALA") == 0 ||
               std::strcmp(name, "PRO") == 0;
    };

    // ─── Exclude residues whose SIDE CHAIN coordinates a metal ────────────────
    // autoflex ranks by MIF score, and the interaction field peaks at the
    // catalytic centre, so it preferentially selects metal-coordinating residues
    // -- precisely the ones that must NOT move in a holo structure. A holo
    // coordination sphere is not a degree of freedom; leaving it selectable makes
    // the arm measure whether the GA can break a metal site rather than induced fit.
    //
    // MEASURED, 1JD0 (carbonic anhydrase XII + acetazolamide): HIS 94, one of the
    // three Zn ligands, was ranked FIRST (MIF -901.57, 566 contacts), and flexing
    // it moved the Zn-NE2 bond from the crystal/rigid 2.064 A to 2.140 A. The
    // rotamer library is metal-blind -- rotobs.lst carries only
    // (RESNAME, chi1, chi2, ...) with no coordination geometry -- so nothing
    // downstream can recover the constraint. 32 of the 85 Astex Diverse targets
    // carry at least one metal (ZN 13, CA 9, MG 7, K 3, NA 2, MN 2, HG 1, LI 1),
    // so this is not an edge case.
    //
    // Both predicates are the ENGINE'S OWN, from metal_coordination.h, so this
    // exclusion and the metal_coord scoring term cannot disagree about what counts
    // as a metal or a donor:
    //   is_metal_type()       -- the same test vcfunction.cpp:852 applies
    //   is_coord_donor_type() -- N (6-12 except N.4), O (13-16), S (17-18), P.3
    // Backbone N/CA/C/O/OXT are exempt: chi rotation cannot move them, so a
    // backbone-to-metal contact carries no risk and excluding on it would be
    // over-aggressive.
    std::vector<const float*> metal_xyz;
    for (int i = 1; i <= FA->atm_cnt; ++i) {
        if (metal_coord::is_metal_type(atoms[i].type)) {
            metal_xyz.push_back(atoms[i].coor);
        }
    }

    auto atom_is_backbone = [](const char* nm) {
        char t[8]; int k = 0;
        for (const char* p = nm; *p && k < 7; ++p) { if (*p != ' ') t[k++] = *p; }
        t[k] = '\0';
        return std::strcmp(t, "N") == 0 || std::strcmp(t, "CA") == 0 ||
               std::strcmp(t, "C") == 0 || std::strcmp(t, "O") == 0 ||
               std::strcmp(t, "OXT") == 0;
    };

    // Returns the offending metal distance (Angstrom) if the side chain coordinates
    // a metal, or -1.0 otherwise, so the caller can report WHY it excluded.
    auto sidechain_metal_distance = [&](int res_index) -> float {
        if (metal_xyz.empty()) return -1.0f;
        const resid& R = residues[res_index];
        if (!R.fatm || !R.latm) return -1.0f;
        const float cut2 = 2.8f * 2.8f;
        float best = -1.0f;
        for (int i = R.fatm[0]; i <= R.latm[0]; ++i) {
            if (atom_is_backbone(atoms[i].name)) continue;
            if (!metal_coord::is_coord_donor_type(atoms[i].type)) continue;
            for (const float* m : metal_xyz) {
                const float dx = atoms[i].coor[0] - m[0];
                const float dy = atoms[i].coor[1] - m[1];
                const float dz = atoms[i].coor[2] - m[2];
                const float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 < cut2) {
                    const float d = std::sqrt(d2);
                    if (best < 0.0f || d < best) best = d;
                }
            }
        }
        return best;
    };

    // Check if residue is already in flex_res[]
    auto already_flexible = [&](int res_index) {
        for (int i = 0; i < FA->nflxsc; ++i) {
            if (FA->flex_res[i].inum == res_index) return true;
        }
        return false;
    };

    // Ensure flex_res array is allocated
    if (!FA->flex_res) {
        FA->MIN_FLEX_RESIDUE = std::max(FA->MIN_FLEX_RESIDUE, max_auto_flex + 5);
        FA->flex_res = static_cast<flxsc*>(
            calloc(static_cast<size_t>(FA->MIN_FLEX_RESIDUE), sizeof(flxsc)));
        if (!FA->flex_res) return 0;
    }

    // Slot policy: a METAL-excluded candidate either yields its slot to the next
    // candidate (backfill, constant gene count) or consumes it and leaves it empty
    // (shrink, gene count varies with metal content). Only metal exclusions are
    // affected -- GLY/ALA/PRO and ligand skips still backfill in both modes -- so the
    // two settings touch exactly one MECHANISM. NOTE, measured 2026-08-30: the SWITCH is
    // single-mechanism but the OUTCOME is NOT single-variable -- it changes both WHICH
    // residues are flexed and HOW MANY genes the chromosome carries (shrink used 1-2 fewer
    // genes on 9 of 12 metal-bearing Astex targets). At a fixed evaluation budget shrink
    // therefore searches a smaller space with more evaluations per gene, so a shrink win
    // cannot be attributed to residue quality without a gene-matched control arm.
    const bool metal_shrink = (FA->autoflex_metal_shrink != 0);
    int added = 0;
    int slots_used = 0;
    for (const auto& rc : key_res) {
        if (slots_used >= max_auto_flex) break;
        if (is_skip_residue(rc.name)) continue;
        if (residues[rc.res_index].type != 0) continue;  // skip ligand residues
        if (already_flexible(rc.res_index)) continue;

        // Metal-coordinating side chains are not degrees of freedom. Report the
        // exclusion rather than skipping silently -- a silently rigid or silently
        // dropped residue is exactly the failure mode this project keeps paying for.
        const float mdist = sidechain_metal_distance(rc.res_index);
        if (mdist >= 0.0f) {
            printf("AUTOFLEX: %s %d:%c EXCLUDED, side chain coordinates a metal "
                   "(%.2f A < 2.8 A cutoff; MIF=%.2f would have ranked it) [%s]\n",
                   rc.name, rc.number, rc.chain, mdist, rc.mif_score,
                   metal_shrink ? "SHRINK: slot left empty"
                                : "BACKFILL: next candidate takes the slot");
            if (metal_shrink) ++slots_used;
            continue;
        }

        // Grow flex_res if needed
        if (FA->nflxsc >= FA->MIN_FLEX_RESIDUE) {
            FA->MIN_FLEX_RESIDUE += 5;
            FA->flex_res = static_cast<flxsc*>(
                realloc(FA->flex_res,
                        static_cast<size_t>(FA->MIN_FLEX_RESIDUE) * sizeof(flxsc)));
            if (!FA->flex_res) return added;
            std::memset(&FA->flex_res[FA->MIN_FLEX_RESIDUE - 5], 0, 5 * sizeof(flxsc));
        }

        // Add to flex_res
        flxsc& fr = FA->flex_res[FA->nflxsc];
        std::strncpy(fr.name, rc.name, 3);
        fr.name[3] = '\0';
        fr.chn = rc.chain;
        fr.num = rc.number;
        fr.inum = rc.res_index;
        set_intprob(&fr);

        FA->nflxsc++;
        added++;
        slots_used++;

        printf("AUTOFLEX: %s %d:%c added as flexible (MIF=%.2f, %d contacts)\n",
               rc.name, rc.number, rc.chain, rc.mif_score, rc.contact_count);
    }

    return added;
}

} // namespace binding_residues
