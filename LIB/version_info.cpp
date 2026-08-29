#include "version_info.h"

#include <cstdio>

// ── Fallbacks ────────────────────────────────────────────────────────────
//
// Every fallback below is the value that is honest when the build system did
// NOT supply the real one.  None of them is optimistic.  A build that reaches
// these defaults reports itself as unidentifiable, which is recoverable; a
// build that guesses "clean" is not, because a false clean claim is
// indistinguishable from a true one downstream.
//
// This is the same failure the packaging work found on
// chore/installability-matrix (682f8090): `git status` failing silently left
// FLEXAIDS_GIT_DIRTY=0, so every pose from a tarball build asserted a clean
// tree with nothing behind it.

#ifndef FLEXAIDS_VERSION_STR
#define FLEXAIDS_VERSION_STR "unknown"
#endif

#ifndef FLEXAIDS_GIT_COMMIT
#define FLEXAIDS_GIT_COMMIT "unknown"
#endif

// 0 = clean, 1 = dirty, 2 = unknown.  The default is 2 and MUST stay 2.
#ifndef FLEXAIDS_GIT_DIRTY
#define FLEXAIDS_GIT_DIRTY 2
#endif

#ifndef FLEXAIDS_SRC_PROVENANCE
#define FLEXAIDS_SRC_PROVENANCE "unknown"
#endif

#ifndef FLEXAIDS_BUILD_TYPE
#define FLEXAIDS_BUILD_TYPE "unknown"
#endif

#ifndef FLEXAIDS_COMPILER
#define FLEXAIDS_COMPILER "unknown"
#endif

#ifndef FLEXAIDS_BUILT_UTC
#define FLEXAIDS_BUILT_UTC "unknown"
#endif

namespace flexaids {
namespace version {

namespace {

// Collapse anything that is not a defined tri-state value onto 2 (unknown).
//
// The macro is an integer substituted by the build system.  If it ever arrives
// malformed -- a stale cache, a hand-passed -D, a packaging script that writes
// an empty string -- the only safe reading is "I do not know".  Clamping up to
// 2 rather than down to 0 is the whole point: the failure mode being defended
// against is a binary that claims clean without evidence.
constexpr int normalized_git_dirty()
{
    constexpr int raw = static_cast<int>(FLEXAIDS_GIT_DIRTY);
    return (raw == 0 || raw == 1) ? raw : 2;
}

}  // namespace

void print_build_identity()
{
    // One key=value per line, LF-terminated.  Field order is for humans only;
    // per docs/run-uniformity/KIND_SCHEMA.md order is not significant and
    // readers must ignore unknown keys.
    //
    // Values are printed with %s from string literals fixed at compile time,
    // so there is no formatting that can fail and no input that can influence
    // this output.

    std::printf("name=FlexAIDdS\n");
    std::printf("version=%s\n", FLEXAIDS_VERSION_STR);
    std::printf("git_commit=%s\n", FLEXAIDS_GIT_COMMIT);

    // Stated explicitly rather than left to be inferred from length, because
    // an 8-hex-character full SHA is not distinguishable from an abbreviated
    // one by inspection.  `short` here means `git rev-parse --short HEAD`, the
    // SAME form stamped into `REMARK FLEXAID.commit=` by
    // BindingMode::output_BindingMode -- so a pose file and a `--version`
    // output can be compared directly, which is the point.
    std::printf("git_commit_form=short\n");

    std::printf("git_dirty=%d\n", normalized_git_dirty());
    std::printf("git_dirty_meaning=0=clean 1=dirty 2=unknown\n");

    // How the commit above was obtained: git | archive | override | unknown.
    // Distinguishes "read from a live checkout" from "asserted by a packager",
    // which are different strengths of evidence for the same field.
    std::printf("src_provenance=%s\n", FLEXAIDS_SRC_PROVENANCE);

    std::printf("build_type=%s\n", FLEXAIDS_BUILD_TYPE);
    std::printf("compiler=%s\n", FLEXAIDS_COMPILER);
    std::printf("built_utc=%s\n", FLEXAIDS_BUILT_UTC);
}

}  // namespace version
}  // namespace flexaids
