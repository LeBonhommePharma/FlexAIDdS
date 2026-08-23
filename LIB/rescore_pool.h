// =============================================================================
// rescore_pool.h — offline pool rescoring with the exact production CF
//
// PURPOSE (near-miss-autopsy follow-up, 2026-08-23):
//   The campaign's deficit is scoring-proxy blindness, not sampling. Testing
//   scoring-term candidates requires scoring EXISTING emitted poses without
//   re-docking. This header exposes the pose-file → coordinate-slot loader;
//   rescore_pool_mode() (rescore_pool.cpp) drives full-pool evaluation.
//
// USAGE (top.cpp hook, before the GA):
//   FLEXAIDDS_RESCORE_POOL=<dir>   pool root; <dir>/<TARGET>/*.pdb if the
//                                  subdirectory exists, else <dir>/*.pdb
//   FLEXAIDDS_RESCORE_OUT=<csv>    optional per-pose CSV output path
//   The mode exits the process after scoring (the GA never runs).
//
// FLEXIBILITY: each pose file carries the FULL complex as emitted (receptor
// atoms included, in whatever side-chain/rotamer state that evaluation used).
// Coordinates are restored for every atom serial found in the file and mapped
// onto atoms[] slots by PDB serial — so both target-side optimisable DoF
// (FA->optres) and ligand torsion states are honoured exactly as docked.
// Files whose atom coverage does not match FA->atm_cnt are refused (fail-closed)
// rather than silently scored against stale coordinates.
//
// Apache-2.0 © 2026 Le Bonhomme Pharma
// =============================================================================

#pragma once

// NOTE: intentionally includes flexaid.h + Vcontacts.h (house pattern, cf.
// native_score.h) so rescore_pool_mode's signature uses canonical typedefs.
#include "flexaid.h"
#include "Vcontacts.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

namespace flexaids {

/// Parse a pose/complex PDB and copy coordinates into coor_out[slot] for every
/// ATOM/HETATM record whose serial exists in serial_to_slot.
///
/// @param pdb_path       pose file to read
/// @param serial_to_slot PDB atom serial -> atoms[] slot index
/// @param coor_out       caller buffer, 3 * n_slots floats; only matched slots
///                       are written (caller must initialise unmatched slots)
/// @param out_matched    optional; receives number of records matched
/// @param out_skipped    optional; receives number of atom records skipped
///                       (unmapped serials or malformed lines)
/// @returns true if the file was opened and parsed to EOF
inline bool load_complex_coor_from_pdb(const char* pdb_path,
                                       const std::unordered_map<int,int>& serial_to_slot,
                                       float* coor_out,
                                       int* out_matched = nullptr,
                                       int* out_skipped = nullptr,
                                       unsigned char* match_mask = nullptr)
{
    if (out_matched) *out_matched = 0;
    if (out_skipped) *out_skipped = 0;
    FILE* f = std::fopen(pdb_path, "r");
    if (!f) return false;

    char buf[256];
    int matched = 0, skipped = 0;
    while (std::fgets(buf, sizeof(buf), f)) {
        if (std::strncmp(buf, "ATOM  ", 6) != 0 &&
            std::strncmp(buf, "HETATM", 6) != 0)
            continue;
        if (std::strlen(buf) < 54) { ++skipped; continue; }
        char serial_buf[8];
        std::memcpy(serial_buf, buf + 6, 5);
        serial_buf[5] = '\0';
        char* endp = nullptr;
        const long serial = std::strtol(serial_buf, &endp, 10);
        if (endp == serial_buf) { ++skipped; continue; }
        auto it = serial_to_slot.find(static_cast<int>(serial));
        if (it == serial_to_slot.end()) {
            ++skipped;
            // Prep may legitimately drop atoms the writer still emits (e.g.
            // united-atom non-polar H); report under debug only.
            if (std::getenv("FLEXAIDDS_RESCORE_DEBUG"))
                std::fprintf(stderr, "[RESCORE] unmapped serial %ld: %.30s\n",
                             serial, buf);
            continue;
        }
        float x, y, z;
        if (std::sscanf(buf + 30, "%f%f%f", &x, &y, &z) != 3) { ++skipped; continue; }
        const int slot = it->second;
        coor_out[slot * 3 + 0] = x;
        coor_out[slot * 3 + 1] = y;
        coor_out[slot * 3 + 2] = z;
        if (match_mask) match_mask[slot] = 1;
        ++matched;
    }
    std::fclose(f);
    if (out_matched) *out_matched = matched;
    if (out_skipped) *out_skipped = skipped;
    return true;
}

} // namespace flexaids

/// Score every pose in the FLEXAIDDS_RESCORE_POOL pool with the exact
/// production CF (full-complex coordinate restore; both-side flexibility
/// honoured), emit [RESCORE] stderr lines + optional CSV, then return.
/// Caller (top.cpp) exits before the GA when the env gate is set.
void rescore_pool_mode(FA_Global* FA, VC_Global* VC, atom* atoms,
                       resid* residue, gridpoint* cleftgrid);
