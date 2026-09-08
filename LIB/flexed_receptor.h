// =============================================================================
// flexed_receptor.h — write the receptor AS SCORED, and let the validator
//                     choose which receptor frame it judges the pose in
//
//   gates: FLEXAIDDS_WRITE_FLEXED_RECEPTOR  (engine, DEFAULT OFF)
//          FLEXAIDDS_PB_RECEPTOR=crystal|flexed|scored (harness, DEFAULT crystal)
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
// THE WATER SET — why "flexed" was not enough, and what "scored" is for
// ---------------------------------------------------------------------
// The header above frames the crystal/as-scored difference as a SIDE-CHAIN
// difference, and concludes that on a RIGID arm the two frames are the same
// conformation and selecting the companion is a structural no-op. The
// conformation claim is true. The no-op conclusion is FALSE, and it hid a
// measured defect for an entire campaign.
//
// The receptor the engine scores is not the receptor file on disk. Every Astex
// arm runs with protein.remove_water=true, keep_structural_waters=true,
// structural_water_bfactor_max=20, so modify_pdb() (modify_pdb.cpp:157-166)
// drops every water with B > 20 while writing the temp receptor that read_pdb()
// then parses into atoms[]/residue[]. The validator, meanwhile, was handed
// entry.receptor_path — the cache receptor, with EVERY water still in it.
//
// So on a rigid arm the two frames differ by exactly the discarded waters, and
// PoseBusters was scoring ligand-water clashes against waters the engine had
// already removed. Measured over 84 Astex targets per arm, water checks were
// 45 of 70 failures (guard) and 50 of 67 (seed2): minimum_distance_to_waters
// 38/41 and volume_overlap_with_waters 7/9. Re-validating offline against a
// B<=20-filtered receptor flipped 6 targets per arm from fail to pass. The cut
// is structure-dependent and large: 1GPK keeps 0 of 529 waters, 2BSM 0 of 247,
// 1XOZ 36 of 352, 1K3U 123 of 712, 1JD0 156 of 504.
//
// The companion already carries the right water set, because it serialises the
// in-memory atoms[] the filtered temp file produced — verified on the one
// rigid-arm companion that exists on this box: 1JD0_1_receptor.pdb has 156 HOH
// atom lines and 4325 atom lines total, against 504 HOH and 4673 total in
// cache_v2/astex_diverse/1JD0/1JD0_apo.pdb, and 504-156 = 4673-4325 = 348. The
// companion IS the apo file minus precisely the waters the engine discarded.
//
// WATERS ARE NOT THE ONLY COMPOSITION DIFFERENCE. modify_pdb() also strips
// hydrogens and keeps only alternate conformation 'A' (see its own banner:
// "Hydrogens are removed", "'A' alternate conformation ONLY is chosen"), and
// that is measurable: 1GPK's cache receptor carries 91 atoms with an altloc
// other than 'A', and the companion carries none of them. On 1GPK the companion
// is 4162 atoms against the apo's 4782 = 529 waters + 91 alt-conformers + 4162.
// Both files then agree atom-for-atom on the remainder (C 2683, N 702, O 755,
// S 22). 1JD0 has no altlocs and no hydrogens, which is why its arithmetic is
// water-only.
//
// So the honest statement of what "scored" gives you is the receptor
// COMPOSITION as modify_pdb produced it — water-filtered by B-factor, altloc-A
// only, hydrogen-free — not "the same file minus some waters". On a structure
// with alternate conformations, attributing a validity flip to waters alone
// would be wrong; the attributable quantity is the composition, and the water
// term is merely the largest part of it on most Astex targets.
//
// Mode "scored" is therefore not a new mechanism — it selects the same
// companion file "flexed" selects. It exists because:
//   * the two modes answer different questions. "flexed" asks "did this pose
//     clash in the side-chain state it was scored in", which is a FLEXIBLE-arm
//     question and needs FLEXAIDDS_AUTOFLEX_MAX>0 to mean anything. "scored"
//     asks "was this pose judged against the receptor composition the engine
//     actually saw", which is a question on EVERY arm, rigid included.
//   * the announce text for "flexed" tells a rigid-arm user the mode is a
//     no-op and to go set AUTOFLEX_MAX. That advice is correct for attributing
//     a side-chain overlap and wrong for the water confound, and a user
//     chasing water failures has no reason to reach for a flag called
//     "flexed".
//   * pb_receptor_used must distinguish them in the record. A rigid arm
//     validated against the engine's water set is not a flexed-side-chain arm,
//     and labelling it "flexed" would make the two uncombinable in analysis.
// What "scored" does NOT do: it does not filter anything itself. modify_pdb.cpp
// remains the single source of truth for receptor composition; a second
// B-factor comparison in the harness is the defect that let this project's
// Cartesian and torsional entropy bases drift apart, and it is not repeated
// here.
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

/// The mode set this build understands, as one string, for every human-facing
/// message. Warning text that lists the modes by hand drifts the moment a mode
/// is added — that is exactly how "crystal|flexed" survived the addition of a
/// third mode in the first draft of this change.
inline const char* pb_receptor_modes_csv() noexcept { return "crystal|flexed|scored"; }

/// True only when FLEXAIDDS_PB_RECEPTOR names a mode this build understands.
/// An unset variable is "recognised" (it selects the documented default).
inline bool pb_receptor_recognised()
{
    const std::string v = pb_receptor_raw();
    return v.empty() || v == "crystal" || v == "flexed" || v == "scored";
}

/// Harness gate: which receptor the validator is handed.
/// "crystal" (DEFAULT, and the value returned for any unrecognised input, so a
/// typo can never silently move a validity number), "flexed", or "scored".
inline std::string pb_receptor_mode()
{
    const std::string v = pb_receptor_raw();
    if (v == "flexed") return "flexed";
    if (v == "scored") return "scored";
    return "crystal";
}

/// True for the modes whose receptor is the ENGINE'S as-scored companion rather
/// than the crystal input file. Both "flexed" and "scored" select the same
/// file — they differ in intent and in what gets recorded, not in the bytes
/// they read (see THE WATER SET, above).
///
/// This predicate exists so the selection site, the announce block and the
/// provenance writer cannot disagree about which modes need a companion. A
/// fourth mode added to pb_receptor_mode() without being added here is a
/// half-wired mode that silently validates against the crystal receptor, which
/// is the failure shape this whole file exists to prevent.
inline bool pb_receptor_uses_companion(const std::string& mode)
{
    return mode == "flexed" || mode == "scored";
}

/// The REMARK key carrying the number of residues sitting off their input
/// rotamer, written by the companion writer (cluster.cpp) and read back by the
/// harness. ONE literal, shared by writer and reader, including the trailing
/// space: a companion whose key the reader cannot find is indistinguishable
/// from a companion with zero moved side chains, and those two mean opposite
/// things when you are deciding whether a validity delta came from the water
/// set or from side-chain motion.
inline const char* remark_key_n_res_off_rotamer() noexcept
{
    return "REMARK n_residues_off_input_rotamer ";
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
