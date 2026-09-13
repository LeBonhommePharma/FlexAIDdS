# ═══════════════════════════════════════════════════════════════════════════
#  Eigen — ONE declared version for every machine, CI leg and container.
# ═══════════════════════════════════════════════════════════════════════════
#
# ── Why this file exists ──────────────────────────────────────────────────
#
# Eigen resolution used to be a four-way first-match-wins chain whose winner
# depended on what happened to be installed on the machine:
#
#   1. vendored LIB/vendor/eigen  — never fired.  .gitmodules declares the
#      submodule, but the index carries NO mode-160000 gitlink, so
#      `git submodule update --init` populates nothing and the directory
#      never appears.  Measured: `git ls-files -s | awk '$1=="160000"'`
#      returned no entries for the whole repository.
#   2. find_package(Eigen3 <floor>) — a FLOOR, not a pin.  ANY Eigen at or
#      above it satisfied the constraint, silently.
#   3. pkg_check_modules(eigen3)  — no version constraint at all.
#   4. FetchContent               — the ONLY pinned branch, and it was
#      pinned to the previous major series.  The one branch that could not
#      fire was also the only branch that was reproducible.
#
# Measured consequence: macOS resolved branch 2 to brew Eigen (a 5.x
# release) while every Linux CI leg installed `libeigen3-dev` (a 3.x
# package).  Two different MAJOR versions of the self-adjoint eigensolver,
# both satisfying the same floor, with nothing recording which one produced
# any given result.
#
# That is a CORRECTNESS problem, not untidiness.  tENCoM classifies a mode
# as a zero mode by the SIGN of a mathematically-zero eigenvalue.  The
# torsional cutoff is exactly 0.0, so the dropped-mode COUNT tracks that
# sign, and the sign is noise-signed even within a single solver: on one
# machine, one solver, one code path, with only the helix length varying,
# the minimum eigenvalue of a torsional basis was measured at +2.2e-14,
# +1.0e-14, -1.0e-13, +1.4e-13 and +3.3e-14 for n = 12, 16, 20, 24, 28.
# The n=20 case is negative and drops a mode that should not be dropped.
# tests/test_svib_invariants.cpp asserts that count with EXPECT_EQ and no
# tolerance, so a different eigensolver can turn that assertion red.
#
# S_vib itself did NOT move in those five cases (delta exactly 0 —
# compute_vibrational_entropy filters those modes regardless), so the
# exposure is the integer mode count and the silent assumption violation,
# not the entropy value.  Do not overstate it.
#
# Two independent fixes are needed and NEITHER substitutes for the other.
# Pinning makes the solver IDENTICAL everywhere (this file).  Making the
# classification solver-INDEPENDENT (the exact-0.0 cutoff, owned elsewhere)
# is the one that survives a future version bump, a different BLAS, or a
# container build.  This file does not make that fix unnecessary.
#
# ── The single source of truth ────────────────────────────────────────────
#
# The declared version and its tarball hash move TOGETHER and appear
# NOWHERE else in the tree.  The system-version comparison, the download
# URL, the integrity hash, every status message, the divergence gate and
# the provenance stamp are all derived from these two lines.  CI reads the
# version out of THIS file too (scripts/ci/eigen_pin.sh), so no workflow
# carries a version literal either.
#
# Version rationale: this is the version that produced every benchmark arm
# currently stored under ~/flexaidds_results — those were all built by
# locally compiled engines, and this machine's system Eigen is that
# version.  Pinning it makes CI comparable with the machine that generated
# the historical baseline instead of invalidating it.  It is also, as of
# this commit, the latest Eigen release, so "pin the newest" and "keep the
# stored arms comparable" happen to name the same version; no trade-off had
# to be resolved.

set(FLEXAIDS_EIGEN_VERSION "5.0.1")
set(FLEXAIDS_EIGEN_SHA256  "e9c326dc8c05cd1e044c71f30f1b2e34a6161a3b6ecf445d56b53ff1669e3dec")

# Derived, not declared.
set(FLEXAIDS_EIGEN_URL
    "https://gitlab.com/libeigen/eigen/-/archive/${FLEXAIDS_EIGEN_VERSION}/eigen-${FLEXAIDS_EIGEN_VERSION}.tar.gz")

# A system Eigen whose version is not the pinned one is a hard configure
# error, not a silent skip (see flexaids_resolve_eigen below).  Set this to
# ignore the system copy entirely and build against the pinned tarball.
# This is the documented way past that error for a developer who cannot
# uninstall their distro Eigen.
option(FLEXAIDS_EIGEN_IGNORE_SYSTEM
       "Ignore any system Eigen and build against the pinned release tarball" OFF)

# Escape hatch for deliberately bisecting a solver difference.  It cannot
# be used to launder a mismatch into a stored result: provenance records
# the version ACTUALLY compiled against, read back from the headers the
# compiler will see, whether or not this is ON.
option(FLEXAIDS_ALLOW_EIGEN_MISMATCH
       "Downgrade the Eigen version-divergence gate from an error to a warning (bisection only)" OFF)

# Offline / air-gapped builds: point this at a local copy of the pinned
# tarball.  It is hash-checked against FLEXAIDS_EIGEN_SHA256 exactly like a
# fresh download, so it is a transport shortcut, not a trust shortcut.
set(FLEXAIDS_EIGEN_TARBALL "" CACHE FILEPATH
    "Local path to the pinned Eigen release tarball (skips the network download)")

# ── Eigen's macro triple is NOT its release version ──────────────────────
#
# Eigen's own Eigen/Version header states that the "WORLD" version will
# forever remain 3 for the Eigen3 library, and that Eigen moved to Semantic
# Versioning as of its 5.0.0 release.  Measured from both release trees,
# for a release numbered A.B.C:
#
#   A >= 5 -> WORLD=3  MAJOR=A  MINOR=B  PATCH=C
#             and EIGEN_VERSION_STRING exists, carrying "A.B.C" directly
#   A <  5 -> WORLD=A  MAJOR=B  MINOR=C
#             no PATCH macro exists, no EIGEN_VERSION_STRING exists, and
#             there is no Eigen/Version header at all -- the macros live in
#             Eigen/src/Core/util/Macros.h
#
# So for any Eigen >= 5 the triple WORLD.MAJOR.MINOR reads "3.A.B", which
# matches no Eigen release and sorts BELOW every 3.x release.  Anything
# recording or comparing an Eigen version must prefer EIGEN_VERSION_STRING
# where it exists and fall back to the triple only below 5.  This function
# encodes that convention once so no caller has to remember it.
#
# Sets in the caller's scope: <out>_WORLD <out>_MAJOR <out>_MINOR <out>_PATCH
function(flexaids_eigen_expected_macros version out)
    string(REPLACE "." ";" _parts "${version}")
    list(LENGTH _parts _n)
    if(_n LESS 3)
        message(FATAL_ERROR
            "FLEXAIDS_EIGEN_VERSION must be MAJOR.MINOR.PATCH; got '${version}'")
    endif()
    list(GET _parts 0 _a)
    list(GET _parts 1 _b)
    list(GET _parts 2 _c)
    if(_a GREATER_EQUAL 5)
        set(${out}_WORLD 3       PARENT_SCOPE)
        set(${out}_MAJOR "${_a}" PARENT_SCOPE)
        set(${out}_MINOR "${_b}" PARENT_SCOPE)
        set(${out}_PATCH "${_c}" PARENT_SCOPE)
    else()
        set(${out}_WORLD "${_a}" PARENT_SCOPE)
        set(${out}_MAJOR "${_b}" PARENT_SCOPE)
        set(${out}_MINOR "${_c}" PARENT_SCOPE)
        set(${out}_PATCH ""      PARENT_SCOPE)
    endif()
endfunction()

# ── Detect the version of an Eigen tree by reading its headers ───────────
#
# Header parsing rather than try_run: it needs no compiler, executes
# nothing, works when cross-compiling, and it REPORTS the version it found
# instead of merely succeeding or failing.  Both candidate header locations
# are checked, so this works on 3.x and 5.x trees alike.
#
# Sets <out> to the canonical release string (or "" if this is not a
# recognisable Eigen tree), and <out>_TRIPLE to the literal
# EIGEN_WORLD_VERSION.EIGEN_MAJOR_VERSION.EIGEN_MINOR_VERSION macro triple.
# Both are recorded in provenance: the release string is the one humans and
# comparisons should use, the triple is what the preprocessor actually says.
function(flexaids_eigen_detect_version inc_dir out)
    set(_vstring "")
    set(_w "")
    set(_m "")
    set(_n "")
    foreach(_hdr "${inc_dir}/Eigen/Version"
                 "${inc_dir}/Eigen/src/Core/util/Macros.h")
        if(NOT EXISTS "${_hdr}")
            continue()
        endif()
        file(STRINGS "${_hdr}" _lines
             REGEX "^[ \t]*#[ \t]*define[ \t]+EIGEN_(WORLD|MAJOR|MINOR)_VERSION|^[ \t]*#[ \t]*define[ \t]+EIGEN_VERSION_STRING")
        foreach(_l IN LISTS _lines)
            if(_l MATCHES "EIGEN_VERSION_STRING[ \t]+\"([^\"]+)\"")
                set(_vstring "${CMAKE_MATCH_1}")
            elseif(_l MATCHES "EIGEN_WORLD_VERSION[ \t]+([0-9]+)")
                set(_w "${CMAKE_MATCH_1}")
            elseif(_l MATCHES "EIGEN_MAJOR_VERSION[ \t]+([0-9]+)")
                set(_m "${CMAKE_MATCH_1}")
            elseif(_l MATCHES "EIGEN_MINOR_VERSION[ \t]+([0-9]+)")
                set(_n "${CMAKE_MATCH_1}")
            endif()
        endforeach()
    endforeach()

    if(NOT _w STREQUAL "" AND NOT _m STREQUAL "" AND NOT _n STREQUAL "")
        set(${out}_TRIPLE "${_w}.${_m}.${_n}" PARENT_SCOPE)
    else()
        set(${out}_TRIPLE "unknown" PARENT_SCOPE)
    endif()

    if(NOT _vstring STREQUAL "")
        set(${out} "${_vstring}" PARENT_SCOPE)
    elseif(NOT _w STREQUAL "" AND NOT _m STREQUAL "" AND NOT _n STREQUAL "")
        set(${out} "${_w}.${_m}.${_n}" PARENT_SCOPE)
    else()
        set(${out} "" PARENT_SCOPE)
    endif()
endfunction()

# ── Resolve Eigen to the declared version ────────────────────────────────
#
# Three branches, in this order.  Every one of them yields the declared
# version or stops the build:
#
#   1. vendored LIB/vendor/eigen — used when the submodule gitlink has been
#      initialised.  Version-checked like everything else, so a stale
#      submodule fails loudly instead of quietly changing the solver.
#   2. system find_package — DISCOVERY ONLY, with no version constraint,
#      followed by an equality test against the pin.  A non-matching
#      system Eigen is a FATAL configure error, not a silent skip.  The
#      reason it is not a silent skip: a developer carrying a 3.x
#      libeigen3-dev should be told that their machine disagrees with the
#      pin, because that is the exact condition that produced two major
#      versions across the matrix in the first place.  An EXACT constraint
#      on find_package would have quietly fallen through to branch 3 and
#      built correctly while leaving the developer's machine still
#      misconfigured for every other tool that reads the system Eigen.
#   3. pinned release tarball, verified against FLEXAIDS_EIGEN_SHA256.
#
# The pkg-config branch is deliberately GONE.  pkg_check_modules cannot be
# relied on to enforce an exact version across the pkg-config
# implementations in play, and it was the last remaining unversioned source
# of Eigen.  Branch 3 covers every platform it used to cover.
#
# Branch 3 has no precondition beyond a reachable network (or a local
# tarball via FLEXAIDS_EIGEN_TARBALL), so a PINNED branch always fires when
# no system Eigen is present.  That is what makes removing the distro,
# Homebrew and Chocolatey Eigen packages from CI safe rather than hopeful.
function(flexaids_resolve_eigen)
    set(_vendor "${CMAKE_SOURCE_DIR}/LIB/vendor/eigen")
    set(_inc "")
    set(_how "")

    if(EXISTS "${_vendor}/Eigen/Dense")
        set(_inc "${_vendor}")
        set(_how "vendored submodule (LIB/vendor/eigen)")
    else()
        if(NOT FLEXAIDS_EIGEN_IGNORE_SYSTEM)
            find_package(Eigen3 QUIET NO_MODULE)
            if(Eigen3_FOUND)
                if(Eigen3_VERSION VERSION_EQUAL "${FLEXAIDS_EIGEN_VERSION}")
                    get_target_property(_sys_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
                    list(GET _sys_inc 0 _inc)
                    set(_how "system cmake config (${Eigen3_DIR})")
                else()
                    message(FATAL_ERROR
                        "System Eigen does not match the pinned version.\n"
                        "  pinned (cmake/FlexAIDEigen.cmake): ${FLEXAIDS_EIGEN_VERSION}\n"
                        "  found on this machine            : ${Eigen3_VERSION}\n"
                        "  found via                        : ${Eigen3_DIR}\n"
                        "\n"
                        "This is a hard error rather than a silent fallback because a "
                        "mismatched system Eigen is exactly the condition that let macOS "
                        "and Linux CI compile two different major versions of the "
                        "self-adjoint eigensolver. tENCoM's zero-mode classification "
                        "depends on the sign of a numerically-zero eigenvalue, so the "
                        "solver is part of what produces a dS_vib mode count.\n"
                        "\n"
                        "Pick one:\n"
                        "  * install the pinned version   (macOS: brew install eigen)\n"
                        "  * build against the pinned tarball and leave the system copy "
                        "alone:  -DFLEXAIDS_EIGEN_IGNORE_SYSTEM=ON\n"
                        "  * initialise the vendored submodule:  git submodule update "
                        "--init --recursive LIB/vendor/eigen")
                endif()
            endif()
        endif()

        if(_inc STREQUAL "")
            set(_deps "${CMAKE_BINARY_DIR}/_deps")
            set(_root "${_deps}/eigen-${FLEXAIDS_EIGEN_VERSION}")
            set(_src  "${_root}/eigen-${FLEXAIDS_EIGEN_VERSION}")

            if(NOT EXISTS "${_src}/Eigen/Dense")
                file(MAKE_DIRECTORY "${_root}")
                if(FLEXAIDS_EIGEN_TARBALL AND EXISTS "${FLEXAIDS_EIGEN_TARBALL}")
                    set(_tar "${FLEXAIDS_EIGEN_TARBALL}")
                    file(SHA256 "${_tar}" _got)
                    if(NOT _got STREQUAL "${FLEXAIDS_EIGEN_SHA256}")
                        message(FATAL_ERROR
                            "FLEXAIDS_EIGEN_TARBALL hash mismatch.\n"
                            "  file     : ${_tar}\n"
                            "  expected : ${FLEXAIDS_EIGEN_SHA256}\n"
                            "  found    : ${_got}")
                    endif()
                else()
                    set(_tar "${_deps}/eigen-${FLEXAIDS_EIGEN_VERSION}.tar.gz")
                    message(STATUS
                        "Eigen3: downloading pinned release ${FLEXAIDS_EIGEN_VERSION}")
                    file(DOWNLOAD "${FLEXAIDS_EIGEN_URL}" "${_tar}"
                         EXPECTED_HASH "SHA256=${FLEXAIDS_EIGEN_SHA256}"
                         TLS_VERIFY ON
                         STATUS _dl_status)
                    list(GET _dl_status 0 _dl_rc)
                    if(NOT _dl_rc EQUAL 0)
                        list(GET _dl_status 1 _dl_msg)
                        message(FATAL_ERROR
                            "Eigen3 ${FLEXAIDS_EIGEN_VERSION} download failed: ${_dl_msg}\n"
                            "  url: ${FLEXAIDS_EIGEN_URL}\n"
                            "If this machine has no network access, fetch that tarball "
                            "elsewhere and pass -DFLEXAIDS_EIGEN_TARBALL=/path/to/it.")
                    endif()
                endif()
                file(ARCHIVE_EXTRACT INPUT "${_tar}" DESTINATION "${_root}")
            endif()

            set(_inc "${_src}")
            set(_how "pinned release tarball (hash-verified)")
        endif()
    endif()

    if(_inc STREQUAL "" OR NOT EXISTS "${_inc}/Eigen/Dense")
        message(FATAL_ERROR
            "Eigen3 could not be resolved to the declared version "
            "${FLEXAIDS_EIGEN_VERSION}.")
    endif()

    # ── The divergence gate ──────────────────────────────────────────────
    #
    # Runs on whatever branch fired, so there is ONE chokepoint rather than
    # one check per branch -- a stale vendored submodule and a shadowing
    # include path are caught by the same code that catches a bad download.
    # Two independent steps, because they can disagree and the
    # disagreement is precisely the bug worth catching:
    #
    #   (a) read the headers CMake resolved, and report the version found;
    #   (b) compile a probe with the real compiler and the real include
    #       path, and static_assert the macros.  This is what catches a
    #       second Eigen earlier on the include search path (CPATH,
    #       CXXFLAGS, a toolchain file, an SDK copy) shadowing (a), which
    #       header parsing alone cannot see.
    flexaids_eigen_detect_version("${_inc}" _detected)
    if(_detected STREQUAL "")
        message(FATAL_ERROR
            "Eigen3: '${_inc}' does not look like an Eigen tree — no version "
            "macros found in Eigen/Version or Eigen/src/Core/util/Macros.h.")
    endif()

    if(NOT _detected STREQUAL "${FLEXAIDS_EIGEN_VERSION}")
        set(_msg
            "Eigen version divergence.\n"
            "  declared (cmake/FlexAIDEigen.cmake): ${FLEXAIDS_EIGEN_VERSION}\n"
            "  actually resolved                  : ${_detected}\n"
            "  source                             : ${_how}\n"
            "  include dir                        : ${_inc}\n"
            "\n"
            "This machine would compile a DIFFERENT self-adjoint eigensolver "
            "than the declared one. tENCoM's zero-mode classification depends "
            "on the sign of a numerically-zero eigenvalue, so dS_vib mode "
            "counts from this build are not guaranteed comparable with any "
            "stored arm.\n"
            "\n"
            "Fixes, in order of preference:\n"
            "  * if the vendored submodule is stale:  git submodule update "
            "--init --recursive LIB/vendor/eigen\n"
            "  * ignore a mismatched system copy:  "
            "-DFLEXAIDS_EIGEN_IGNORE_SYSTEM=ON\n"
            "  * install the declared version  (macOS: brew install eigen)\n"
            "\n"
            "To bisect a solver difference on purpose, reconfigure with "
            "-DFLEXAIDS_ALLOW_EIGEN_MISMATCH=ON. Provenance will still record "
            "${_detected}, not the declared version.")
        if(FLEXAIDS_ALLOW_EIGEN_MISMATCH)
            message(WARNING ${_msg})
        else()
            message(FATAL_ERROR ${_msg})
        endif()
    endif()

    flexaids_eigen_expected_macros("${FLEXAIDS_EIGEN_VERSION}" _exp)
    set(_probe "
#include <Eigen/Core>
static_assert(EIGEN_WORLD_VERSION == ${_exp_WORLD}, \"Eigen WORLD version is not the pinned one\");
static_assert(EIGEN_MAJOR_VERSION == ${_exp_MAJOR}, \"Eigen MAJOR version is not the pinned one\");
static_assert(EIGEN_MINOR_VERSION == ${_exp_MINOR}, \"Eigen MINOR version is not the pinned one\");
")
    if(NOT _exp_PATCH STREQUAL "")
        string(APPEND _probe
            "static_assert(EIGEN_PATCH_VERSION == ${_exp_PATCH}, \"Eigen PATCH version is not the pinned one\");\n")
    endif()
    string(APPEND _probe "int main() { return 0; }\n")

    include(CheckCXXSourceCompiles)
    set(CMAKE_REQUIRED_INCLUDES "${_inc}")
    set(CMAKE_REQUIRED_QUIET TRUE)
    check_cxx_source_compiles("${_probe}" FLEXAIDS_EIGEN_MACROS_AGREE)
    unset(CMAKE_REQUIRED_INCLUDES)
    unset(CMAKE_REQUIRED_QUIET)

    if(NOT FLEXAIDS_EIGEN_MACROS_AGREE)
        set(_msg2
            "Eigen version gate: the compiler does not see the Eigen that CMake "
            "resolved.\n"
            "  resolved include dir : ${_inc}\n"
            "  resolved version     : ${_detected}\n"
            "  expected macros      : WORLD=${_exp_WORLD} MAJOR=${_exp_MAJOR} "
            "MINOR=${_exp_MINOR} PATCH=${_exp_PATCH}\n"
            "\n"
            "The usual cause is a second Eigen earlier on the include search "
            "path (CPATH, CXXFLAGS, a toolchain file, or an SDK copy) "
            "shadowing the resolved one.\n"
            "\n"
            "Note: Eigen's macro triple is not its release version for Eigen "
            ">= 5, because WORLD is frozen at 3. The expected values above "
            "already account for that.")
        if(FLEXAIDS_ALLOW_EIGEN_MISMATCH)
            message(WARNING ${_msg2})
        else()
            message(FATAL_ERROR ${_msg2})
        endif()
    endif()

    if(NOT TARGET Eigen3::Eigen)
        add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
        set_target_properties(Eigen3::Eigen PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${_inc}")
    endif()

    message(STATUS
        "Eigen3 ${_detected} (macros ${_detected_TRIPLE}) via ${_how} — "
        "pin ${FLEXAIDS_EIGEN_VERSION}, gate OK")

    # Consumed by the provenance stamps in the top-level CMakeLists.txt.
    # These are the DETECTED values, not the declared ones, so what gets
    # recorded beside a result is a measurement of what was compiled rather
    # than a restatement of intent.
    set(FLEXAIDS_EIGEN_VERSION_DETECTED "${_detected}"        PARENT_SCOPE)
    set(FLEXAIDS_EIGEN_MACRO_TRIPLE     "${_detected_TRIPLE}" PARENT_SCOPE)
    set(FLEXAIDS_EIGEN_SOURCE           "${_how}"             PARENT_SCOPE)
    set(Eigen3_FOUND                    TRUE                  PARENT_SCOPE)
endfunction()
