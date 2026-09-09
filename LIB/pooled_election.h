// =============================================================================
// pooled_election.h — choose the OBJECTIVE the pooled cluster-head election
//                     optimises, without touching the SCOPE it elects over
//
//   gate: FLEXAIDDS_POOLED_ELECTION=default|mincf   (harness, DEFAULT default)
//
// Apache-2.0 (c) 2026 Le Bonhomme Pharma
//
// READ THIS FIRST: CROSS-RESTART POOLING IS ALREADY THE DEFAULT
// -------------------------------------------------------------
// This file was written to add "pooled-restart election". It does not add it,
// because the engine has been pooling restarts before electing all along, and
// that was MEASURED on the stored arms rather than assumed:
//
//   * DatasetRunner.cpp: n_restarts = max(1, protocol_cfg_.restarts), and
//     ProtocolConfig.h declares `int restarts{5}` — multi-restart is the
//     DEFAULT, not an opt-in.
//   * the restart loop pushes EVERY restart prefix into `all_prefixes`
//     unconditionally (one push_back per ri, no success predicate guarding it),
//     and the cached/skip path rebuilds the same list from the r<N>/
//     subdirectories already on disk.
//   * all three call sites of select_pose_freq_gated_pooled() are handed
//     `all_prefixes`, and its candidate loop iterates over every prefix in it.
//   * the `[POOL] <ID> pooling N/M restart(s) for Fix B selection` receipt
//     read `pooling 3/3` on 254 of 254 target receipts across the three
//     84-target Astex arms guard_arm_20260907_062124 / seed2_arm_20260908_044750
//     / seed3_arm_20260908_123008.
//   * the `[ELECTED-POSE] ... restart=` receipt shows the elected pose came
//     from a restart OTHER than r0 in 165 of 252 cells (65.5%): guard 48/84,
//     seed2 55/84, seed3 62/84.
//
// So there is no pooling mechanism to add and no pooling default to flip. The
// single-prefix elector select_pose_freq_gated() — the one that really does see
// restart 0 only — has NO CALLERS; it survives as dead code, named in a comment
// at the election site that describes behaviour the call does not have. That
// comment is the most likely origin of the standing belief that "the engine
// elects from restart 0 only", and it is corrected in the commit before this
// one.
//
// WHAT THE MOTIVATING MEASUREMENT ACTUALLY DECOMPOSES INTO
// -------------------------------------------------------
// Lane B (SELECTION_RECOVERY_REPORT.md) re-elected over 31,870 stored poses,
// 84 targets x 3 arms, and measured, in targets with a sub-2 A top-1 pose out
// of 84 (symmetry-corrected RMSD):
//
//     the arm's own elected pose                        40.0
//     min-CF over the pooled restarts                   45.3   (+5.3)
//     oracle (best pose in the written pool)            67.3
//     random pick                                        4.3
//
// and decomposed the +5.3 as +1.3 for "elect min-CF WITHIN restart 0" and +4.0
// for "then pool the three restarts". That +4.0 is a real measured step, but it
// is measured between two COUNTERFACTUAL electors — an r0-restricted min-CF
// rule and a pooled min-CF rule — and the engine is neither of them. Against
// the engine's actual baseline, which already pools, the entire +5.3 is a
// change of OBJECTIVE (frequency-gated composite -> plain min-CF) and exactly
// 0.0 of it is a change of SCOPE.
//
// Pooling is still load-bearing for the SIZE of that gain: confining the same
// objective change to restart 0 recovers only +1.3 of the +5.3. It is
// load-bearing and already switched on. Both halves of that sentence matter,
// and shipping a flag that claims to "enable pooling" would have made the
// second half unsayable.
//
// WHAT THIS GATE THEREFORE DOES
// -----------------------------
// FLEXAIDDS_POOLED_ELECTION=mincf replaces the election OBJECTIVE with plain
// min-CF over the already-pooled candidate set — Lane B's elector (b), the
// 45.3 row. It does not widen, narrow or reorder the candidate set. Set to
// anything else, or left unset, the elector is bit-for-bit today's.
//
// The min-CF rule is deliberately the whole of it:
//   * every finite-CF cluster head from every pooled restart prefix competes;
//   * the frequency gate (prefer Frequency>1 heads), seed-anchored elitism,
//     macro-cluster consensus demotion and the Softbeta free-energy branch are
//     all bypassed — they are the objective being replaced;
//   * `<prefix>_INI.pdb` seeds CANNOT win, because they are not in the
//     candidate set at all: the enumerator matches digit-only rank suffixes,
//     so "INI" never enumerates, and the seed list is built separately and
//     downstream of the short-circuit. A crystal-seeded pose winning a
//     min-CF election would be an oracle leak, so this is a property worth
//     stating rather than assuming.
//   * ties resolve to the first candidate in enumeration order (prefix order,
//     then rank), which is deterministic.
// Sentinel poses need no special case: the engine flags them with CF >= 1e3,
// and a min-CF argmin cannot select a large positive CF while any scored head
// is negative. Lane B excluded them explicitly; here they are excluded by the
// arithmetic, which is why the two agree.
//
// WHAT IT DOES NOT DO
// -------------------
// No CF channel moves, no pose coordinates change, no REMARK is written, and
// the engine never reads this variable — it is read by the HARNESS only. The
// gate re-elects from poses already on disk, so it is meaningful on a cached
// re-score and needs no re-dock. It is NOT a benchmark endpoint: Lane B's
// per-arm McNemar on this contrast is p = 0.09-0.39, so on a single arm the
// gain is not separable from seed noise, and the case for it rests on 3/3 sign
// consistency plus the mechanism being understood. Ship it as an arm, not as a
// default.
//
// WHY A MODE STRING AND NOT A BOOLEAN
// -----------------------------------
// Because "pooled election" has already been misread once, and a boolean has
// no room to be read wrongly in a recorded field. `pooled_election_used`
// distinguishes "mincf" from "mincf_single_restart_pool" — a one-restart pool
// where the objective changed but pooling contributed nothing. That is a NULL
// CONTROL, not an error, and it must not be aggregated with genuine
// multi-restart pooled elections. A boolean would have made the two
// indistinguishable in the record.
// =============================================================================

#pragma once

#include <cctype>
#include <cstdlib>
#include <string>

namespace flexaids {
namespace pooled_election {

/// Raw value of FLEXAIDDS_POOLED_ELECTION, lowercased and whitespace-trimmed.
/// Empty when the variable is unset or empty.
inline std::string raw()
{
    const char* s = std::getenv("FLEXAIDDS_POOLED_ELECTION");
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

/// raw() reduced to characters that are safe to paste into a JSON string
/// value, and capped at 32 characters. The unrecognised-value note quotes what
/// the user actually typed back at them, and that note is written into
/// validator_provenance.json; a value containing a quote or a backslash would
/// otherwise turn a diagnostic into an unparseable provenance file, which is
/// the one file that must survive a garbage input intact.
inline std::string raw_sanitized()
{
    const std::string v = raw();
    std::string out;
    for (std::size_t i = 0; i < v.size() && out.size() < 32; ++i) {
        const unsigned char c = static_cast<unsigned char>(v[i]);
        const bool ok = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
                        c == '_' || c == '-' || c == '.' || c == '+';
        out.push_back(ok ? static_cast<char>(c) : '?');
    }
    if (out.size() < v.size()) out += "...";
    return out;
}

/// The mode set this build understands, as one string, for every human-facing
/// message. Warning text that lists the modes by hand drifts the moment a mode
/// is added.
inline const char* modes_csv() noexcept { return "default|mincf"; }

/// True only when FLEXAIDDS_POOLED_ELECTION names a mode this build
/// understands. An unset variable is "recognised" (it selects the documented
/// default).
inline bool recognised()
{
    const std::string v = raw();
    return v.empty() || v == "default" || v == "mincf";
}

/// Which election objective was asked for. "default" (DEFAULT, and the value
/// returned for any unrecognised input, so a typo can never silently move an
/// elected pose) or "mincf".
inline std::string mode()
{
    const std::string v = raw();
    if (v == "mincf") return "mincf";
    return "default";
}

/// mode(), snapshotted on first use.
///
/// This exists for a correctness reason, not for speed. The elector is called
/// from three sites per target — the one that fills the CSV `best_score`, the
/// one that elects the pose RMSD and PoseBusters are computed on, and the
/// no-crystal-reference fallback — and `best_score` and `rmsd_to_crystal` are
/// only interpretable if all three describe ONE pose. Re-reading getenv at each
/// site would make that invariant depend on the environment not changing
/// mid-run. Snapshotting makes the three sites provably agree.
inline const std::string& mode_snapshot()
{
    static const std::string m = mode();
    return m;
}

/// True for the modes whose objective is plain min-CF over the pooled
/// candidate set rather than today's frequency-gated composite.
///
/// This predicate exists so the elector, the receipt line and the provenance
/// writer cannot disagree about which modes change the objective. A third mode
/// added to mode() without being added here is a half-wired mode that silently
/// elects with the default objective while the record claims otherwise — the
/// failure shape that made this file necessary in the first place.
inline bool elects_pooled_mincf(const std::string& m)
{
    return m == "mincf";
}

}  // namespace pooled_election
}  // namespace flexaids
