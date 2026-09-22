#ifndef FLEXAIDS_TEMP_DIR_H
#define FLEXAIDS_TEMP_DIR_H

#include <string>

// Single, never-throwing resolution point for the process temporary directory.
//
// WHY THIS EXISTS
// ---------------
// std::filesystem::temp_directory_path() has two overloads. The throwing one
// aborts the process with
//
//     Fatal error: filesystem error: in temp_directory_path: <path>
//
// when the directory named by $TMPDIR no longer exists. On 2026-09-20 that cost
// nine cells of an 85-target Astex campaign: $TMPDIR pointed into an agent
// session workspace that its harness deleted after a few hours idle, and every
// invocation started after the sweep died at startup with no poses written and
// a non-zero exit the driver reported only as "Incomplete docking". The same
// nine targets had docked normally in an earlier campaign; nothing was wrong
// with the engine, the ligands or the receptors.
//
// A path that *resolves* is not a path that *works*: the failure mode here is a
// directory that existed when the process list was built and was gone by the
// time the cell ran. Resolution must therefore validate, not merely parse.
//
// CONTRACT
// --------
//  * never throws, never aborts, always returns a usable directory
//  * validates that the candidate exists, is a directory, and is WRITABLE, by
//    creating and removing a probe file -- access() lies under some sandboxes
//  * resolves once per process and caches; later deletion of the chosen
//    directory cannot change the answer mid-run, which keeps a long docking run
//    self-consistent
//  * emits exactly one diagnostic line per rejected candidate, naming the
//    source (environment variable or fallback) and the reason
//
// FALLBACK CHAIN, first writable candidate wins:
//    1. $FLEXAIDDS_TMPDIR   explicit override, for campaign drivers that must
//                           pin a durable location outside an ephemeral tree
//    2. std::filesystem::temp_directory_path(error_code)   honours TMPDIR/TMP/TEMP
//    3. /tmp                POSIX fallback
//    4. the current working directory
//
// FOR LOG MONITORS -- read this before writing a grep against engine stderr.
// Every line this file emits is prefixed "[TEMPDIR] " and is INFORMATIONAL:
// it means a candidate was rejected and the chain moved on, i.e. the guard
// working. The only line that indicates the process could not proceed is
// "[TEMPDIR] FATAL-SOFT", and even that returns rather than aborting.
// A genuine abort is the C++ runtime's own "Fatal error: filesystem error".
// MEASURED 2026-09-20: a monitor grepping for "temp_directory_path" matched
// BOTH the abort and the healthy recovery, because the recovery message named
// the function it had just survived. Key on "Fatal error:" or "[TEMPDIR] FATAL-SOFT",
// never on the function name.
//
// If every candidate fails the writability probe the function returns "." and
// flexaids_temp_dir_ok() reports false, so a caller that needs to refuse can.

namespace flexaids {

// Resolved, validated, writable temporary directory. No trailing slash.
const std::string& temp_dir();

// False when every candidate failed and temp_dir() fell back to ".".
bool temp_dir_ok();

// Human-readable provenance of the resolved directory, e.g.
// "FLEXAIDDS_TMPDIR" or "std::temp_directory_path" or "/tmp". Intended for
// run receipts so a campaign records where its scratch actually lived.
const std::string& temp_dir_source();

// Testing seam: drop the cached resolution so the next temp_dir() re-resolves.
// Not for production use; the cache is deliberate.
void temp_dir_reset_for_testing();

}  // namespace flexaids

#endif  // FLEXAIDS_TEMP_DIR_H
