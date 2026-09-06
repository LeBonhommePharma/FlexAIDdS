// =============================================================================
// flexed_receptor.h — write the receptor AS SCORED, and let the validator
//                     choose which receptor frame it judges the pose in
//
//   gates: FLEXAIDDS_WRITE_FLEXED_RECEPTOR  (engine, DEFAULT OFF)
//          FLEXAIDDS_PB_RECEPTOR=crystal|flexed (harness, DEFAULT crystal)
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma
//
// THE PROBLEM
// -----------
// PoseBusters is handed entry.receptor_path — the CRYSTAL receptor — for every
// arm (DatasetRunner.cpp, the validate_elected_pose call). In the RIGID arm
// that is exactly right: the receptor the engine scored against and the file
// on disk are the same coordinates.
//
// In the FLEXIBLE arm they are not. The GA moved side chains, so a reported
// `volume_overlap_with_protein` failure is UNATTRIBUTABLE: it may be a genuine
// clash, or it may be an artifact of measuring an induced-fit pose against an
// unflexed reference in which the side chain is still sitting where the ligand
// now is. The campaign could not separate the two, because
// FLEXAIDDS_SCORED_ONLY=1 suppressed the flexed receptor coordinates from the
// emitted pose file, so the frame the engine actually scored in was never
// written to disk at all.
//
// WHAT THIS ADDS
// --------------
// 1. FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1 makes cluster.cpp write, next to every
//    emitted pose, the receptor AS SCORED: every residue at its live rotamer
//    index residue[].rot — the same field write_pdb.cpp:140 reads — with the
//    docked ligand removed. Cofactors, metals and waters are kept, because the
//    engine scored against them too.
// 2. FLEXAIDDS_PB_RECEPTOR selects which of those two receptors the validator
//    is handed. It defaults to "crystal", so today's validity outcomes are
//    reproduced EXACTLY and no published number moves.
//
// NEITHER GATE TOUCHES SCORING. No CF channel, no REMARK on the pose, no
// existing artifact. With both unset the engine and the harness are
// bit-identical to HEAD: the writer is never called and the validator receives
// the same entry.receptor_path it receives today.
//
// WHY A SUBDIRECTORY, NOT A SIBLING FILE
// --------------------------------------
// The companion is written to <dir-of-pose>/flexed_receptor/<stem>_receptor.pdb
// and NOT beside the pose. Offline analysis in scripts/ enumerates poses with
// globs as loose as "*.pdb", "*_*.pdb", "{pdb_id}_*.pdb" and "*_[0-9]*.pdb"
// (scripts/audit_native_cf.py, scripts/backfill_inline_rmsd.py,
// scripts/ranking_bias_audit.py, scripts/e10_election_vs_scoring.py). ANY name
// ending in .pdb in that directory would be silently counted as a pose by at
// least one of them. A subdirectory is invisible to every one of those globs,
// and to the harness's own scans, which either match exact
// <prefix>_<rank>.pdb names, require "_mode_"/"_cluster_" in the filename, or
// test is_regular_file() first.
//
// HOW TO READ A FLEXED-FRAME VERDICT (important)
// ----------------------------------------------
// Scoring the pose against the flexed receptor answers "did this pose clash in
// the state it was scored in". It CANNOT answer "is the receptor state itself
// physical": if the search evicted a side chain into a strained or clashing
// rotamer, the ligand-protein overlap disappears because the protein moved out
// of the way, and the flexed-frame verdict gets BETTER for the wrong reason.
// The flexed-frame number is therefore only interpretable next to a
// receptor-internal charge — that is what FLEXAIDDS_RECEPTOR_STRAIN exists for.
// Report the pair, never the flexed-frame number alone.
// =============================================================================

#pragma once

#include <cctype>
#include <cstdlib>
#include <string>

#include "EnvFlags.h"  // flexaids::env_bool — one parser for FLEXAIDDS_* switches

namespace flexaids {
namespace flexed_receptor {

/// Subdirectory (relative to the directory holding the pose) that receives the
/// as-scored receptor companions. Deliberately NOT a sibling file — see the
/// header comment.
inline const char* subdir_name() noexcept { return "flexed_receptor"; }

/// Engine gate: FLEXAIDDS_WRITE_FLEXED_RECEPTOR. DEFAULT OFF.
/// Unset / empty / 0 / false / no / off -> false (EnvFlags parsing).
inline bool write_enabled() noexcept
{
    return flexaids::env_bool("FLEXAIDDS_WRITE_FLEXED_RECEPTOR", false);
}

/// Raw value of FLEXAIDDS_PB_RECEPTOR, lowercased and whitespace-trimmed.
/// Empty when the variable is unset or empty.
inline std::string pb_receptor_raw()
{
    const char* s = std::getenv("FLEXAIDDS_PB_RECEPTOR");
    if (s == nullptr) return std::string();
    std::string v(s);
    std::size_t b = 0;
    while (b < v.size() && std::isspace(static_cast<unsigned char>(v[b]))) ++b;
    std::size_t e = v.size();
    while (e > b && std::isspace(static_cast<unsigned char>(v[e - 1]))) --e;
    v = v.substr(b, e - b);
    for (std::size_t i = 0; i < v.size(); ++i)
        v[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(v[i])));
    return v;
}

/// True only when FLEXAIDDS_PB_RECEPTOR names a mode this build understands.
/// An unset variable is "recognised" (it selects the documented default).
inline bool pb_receptor_recognised()
{
    const std::string v = pb_receptor_raw();
    return v.empty() || v == "crystal" || v == "flexed";
}

/// Harness gate: which receptor the validator is handed.
/// "crystal" (DEFAULT, and the value returned for any unrecognised input, so a
/// typo can never silently move a validity number) or "flexed".
inline std::string pb_receptor_mode()
{
    const std::string v = pb_receptor_raw();
    if (v == "flexed") return "flexed";
    return "crystal";
}

/// Directory that holds the as-scored receptor companion for `pose_pdb_path`.
/// Never empty; a bare filename yields "./flexed_receptor".
inline std::string companion_dir(const std::string& pose_pdb_path)
{
    const std::size_t slash = pose_pdb_path.find_last_of("/\\");
    const std::string dir = (slash == std::string::npos)
                                ? std::string(".")
                                : pose_pdb_path.substr(0, slash);
    return dir + "/" + subdir_name();
}

/// Full path of the as-scored receptor companion for `pose_pdb_path`:
///     .../<dir>/flexed_receptor/<pose-stem>_receptor.pdb
/// The ".pdb" suffix of the pose is stripped before "_receptor.pdb" is added,
/// so <prefix>_3.pdb -> <dir>/flexed_receptor/<prefix>_3_receptor.pdb.
/// Pure string arithmetic: identical in the engine (writer) and the harness
/// (reader), which is what makes the two halves agree without a manifest.
inline std::string companion_path(const std::string& pose_pdb_path)
{
    const std::size_t slash = pose_pdb_path.find_last_of("/\\");
    std::string base = (slash == std::string::npos)
                           ? pose_pdb_path
                           : pose_pdb_path.substr(slash + 1);
    if (base.size() > 4 && base.compare(base.size() - 4, 4, ".pdb") == 0)
        base = base.substr(0, base.size() - 4);
    return companion_dir(pose_pdb_path) + "/" + base + "_receptor.pdb";
}

}  // namespace flexed_receptor
}  // namespace flexaids
