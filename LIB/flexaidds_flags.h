// flexaidds_flags.h — unified runtime/compile gate registry + convenience overlay
//
// Streamlines *usage* of FLEXAIDDS_* / FLEXAID_* / FLEXAIDS_* knobs without
// removing any option, env var, or -D flag. Individual getenv / compile
// defines keep working. Mutually exclusive pairs keep both names in the
// registry; resolve_once() auto-disables the loser and emits a one-time
// stderr warning.
//
// Convenience overlay (UNION with individual env vars, not a replacement):
//   FLEXAIDDS_FLAGS=hoist,epoch,fastpath,rng-stream-fix
// Tokens are case-insensitive; '-' and '_' are equivalent.
//   FLEXAIDDS_FLAGS_DUMP=1  → dump the full registry once at resolve.
//
// Query RIGID_FASTPATH from C++:
//   flexaidds::flags::active("FLEXAIDDS_RIGID_FASTPATH")
//   flexaidds::flags::rigid_fastpath()
//
// Copyright 2026 Le Bonhomme Pharma
// SPDX-License-Identifier: Apache-2.0
#ifndef FLEXAIDDS_FLAGS_H
#define FLEXAIDDS_FLAGS_H

#include <cstdio>

namespace flexaidds {
namespace flags {

/// Idempotent, thread-safe resolve of env + compile defs + overlay + exclusions.
void resolve_once();

/// Drop the resolved snapshot so the next resolve_once() re-reads the environment.
/// Tests that setenv/unsetenv must call this (or rely on requested/active after it).
void reset_for_tests();

/// Raw env / compile / FLEXAIDDS_FLAGS request, before mutual-exclusion.
bool requested(const char* name);

/// After implication + exclusion (what should actually be treated as on).
bool active(const char* name);

/// Env/overlay value string for a known gate, or "" if unset / unknown.
const char* value(const char* name);

/// Disable reason when requested && !active; otherwise "".
const char* reason(const char* name);

/// Print every known gate: name, requested, active, value, reason-if-disabled.
void dump(FILE* out);

/// Push the resolved registry into the process environment so existing
/// `getenv("FLEXAIDDS_*")` call sites honour the overlay and exclusions
/// without being rewritten. Does not delete any flag from the API:
///   * active but unset  → setenv(name, "1") (overlay / implications)
///   * requested but inactive → unsetenv(name) (mutual-exclusion loser)
/// Compile-time -D gates are never written. Individual env vars that are
/// already set and still active are left untouched.
void apply_to_environ();

/// Convenience: same as active("FLEXAIDDS_RIGID_FASTPATH"). Default OFF.
inline bool rigid_fastpath() { return active("FLEXAIDDS_RIGID_FASTPATH"); }

inline bool hoist_receptor_index() {
    return active("FLEXAIDDS_HOIST_RECEPTOR_INDEX");
}

inline bool contacts_epoch() { return active("FLEXAIDDS_CONTACTS_EPOCH"); }

}  // namespace flags
}  // namespace flexaidds

#endif  // FLEXAIDDS_FLAGS_H
