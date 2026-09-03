#ifndef FLEXAIDS_VERSION_INFO_H
#define FLEXAIDS_VERSION_INFO_H

// Build identity for the FlexAID∆S engine.
//
// WHY THIS EXISTS
// ---------------
// On 2026-08-28 a three-pass docking comparison was invalidated and could not
// be reconstructed.  `build/FlexAIDdS` was relinked seven seconds before pass
// 3; the earlier passes' binaries were overwritten and their hashes are gone.
// argv was `./build/FlexAIDdS` throughout -- the *path* was stable, the *file*
// behind it was not, and nothing in any of the three output trees recorded
// which engine produced them.
//
// The engine previously embedded its commit only in a `REMARK FLEXAID.commit=`
// line inside emitted pose files, so a binary that had not yet run anything was
// unidentifiable, and a binary whose poses had been post-processed was
// unidentifiable in retrospect.  `--version` closes that gap: the build
// identity is recoverable from the binary alone, months later, with no inputs
// and no prior run.
//
// The output is deliberately the same shape as the frozen `KIND` and `DONE`
// sidecars (docs/run-uniformity/): one `key=value` per line, LF-terminated,
// UTF-8, no quoting, no comments, no blank lines, order not significant,
// readers ignore unknown keys.  A third format for the same job is exactly the
// defect that effort exists to remove.

namespace flexaids {
namespace version {

// Writes the build identity to stdout, one `key=value` per line.
//
// Contract, in order of importance:
//   * Never claims a clean tree it did not observe.  `git_dirty` is 0 ONLY
//     when `git status --porcelain` actually ran and returned empty.
//   * No side effects: no file access, no RNG, no allocation of engine state.
//     Safe to call as the very first statement of main().
//   * stdout only, and nothing but key=value lines -- no banner, no colour.
void print_build_identity();

}  // namespace version
}  // namespace flexaids

#endif  // FLEXAIDS_VERSION_INFO_H
