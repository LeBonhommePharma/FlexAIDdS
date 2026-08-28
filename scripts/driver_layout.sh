#!/usr/bin/env bash
# driver_layout.sh — one place where a benchmark driver decides its run directory.
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
#
# ── WHAT THIS IS ────────────────────────────────────────────────────────────────
# Sourceable library. Exports exactly one public function:
#
#     new_run_dir <kind> <tier> <name>      # kind: test|arm   tier: T|A
#
# It creates the run directory, writes a KIND file into it, and echoes the
# absolute path on stdout. Everything else it prints goes to stderr, so
#
#     O=$(new_run_dir arm A 15_dofscale) || exit 1
#
# is always safe.
#
# ── WHY IT EXISTS ───────────────────────────────────────────────────────────────
# Across ~/flexaidds_results there are 60 result-bearing trees in 53 distinct
# structural shapes. Coverage of an artifact tracks who writes it, not how much
# it matters: engine-written files land ~85% of the time, shell-driver-written
# files land 11-40% of the time. 11% is what "the driver author remembers to do"
# achieves, and no amount of care raises it, because each driver re-decides the
# same five things from scratch.
#
# Those five decisions, taken as the union of what
#   run_arm15_dofscale.sh:53,83,107,123,134
#   run_t13_twotarget.sh:145,195,221,238,254
#   run_cachefix_2targets.sh:48,64,83,95,106
# each open-code, are:
#
#   1. Timestamp + batch naming     all three: date +%Y%m%d_%H%M%S, $R/<slug>_$STAMP
#   2. Run-dir naming + subdirs     all three mkdir "$O/bin" "$O/run" — but name it
#                                   three different ways: arm15_dofscale, run_t13_s374,
#                                   13_intragenes. Three conventions, three parsers.
#   3. Collision policy             arm15:82 and t13:194 mv aside, NEVER delete;
#                                   cachefix has no policy at all.
#   4. Disk floor + SKIPPED         arm15/cachefix hardcode 6 GB; t13 uses
#                                   ${FLEXAIDDS_DISK_FLOOR_GB:-2}. Same intent, two
#                                   floors, one env knob only one driver honours.
#   5. Engine identity              all three shasum the binary into provenance.txt,
#                                   but only t13:205 verifies the staged copy still
#                                   matches. --redock writes no receipt at all, which
#                                   is why a mid-sequence relink on 2026-08-28 left the
#                                   engine identity unrecoverable: inodes gone, and the
#                                   project builds with Unix Makefiles, so there is no
#                                   build ledger to fall back on.
#
# engine_sha is therefore mandatory here. new_run_dir REFUSES to create a
# directory it cannot stamp with the identity of the binary that will run in it.
# A run directory without engine identity is not a cheap run directory, it is a
# run directory whose results cannot be attributed later.
#
# ── COMPOSITION ─────────────────────────────────────────────────────────────────
# Independent of driver_preamble.sh. If the preamble has been sourced its
# FLEXAIDDS_RESULTS is picked up automatically; if not, set FLEXAIDDS_RUN_ROOT.
# This file is a library: it never calls exit, never sets shell options, and does
# nothing at source time. Sourcing it cannot change the caller's behaviour.
#
# ── ADOPTION ────────────────────────────────────────────────────────────────────
# NOT wired into any existing driver. Adoption is Phase 3 and forward-only:
# new drivers call it, existing drivers are left exactly as they are.
#
# ── ENVIRONMENT ─────────────────────────────────────────────────────────────────
#   FLEXAIDDS_ENGINE_BIN      REQUIRED. Absolute path to the FlexAIDdS binary that
#                             will be staged and invoked. Hashed into KIND.
#   FLEXAIDDS_SEAT            REQUIRED. One of: grok science claude-code dispatch lp
#   FLEXAIDDS_RUN_ROOT        Batch parent. Defaults to $FLEXAIDDS_RESULTS.
#   FLEXAIDDS_BATCH_DIR       Existing batch to place this run inside. When unset,
#                             a new batch <root>/<name>_<STAMP> is created.
#   FLEXAIDDS_DISK_FLOOR_GB   Default 6 (the stricter of the two floors in use).
#
# ── RETURN CODES ────────────────────────────────────────────────────────────────
#   0  directory created, KIND written, path on stdout
#   2  refused: bad arguments, missing/unhashable engine, unusable root
#   3  created + KIND + SKIPPED sentinel, below disk floor. Path still on stdout.
#      "never started" is a different fact from "died", and only the driver that
#      declined to start can record it.

# ── internals (prefixed _dl_, not part of the contract) ─────────────────────────

_dl_err() { printf '[driver_layout] %s\n' "$*" >&2; }

_dl_sha256() {                      # $1 file -> 64 hex on stdout
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" 2>/dev/null | cut -d' ' -f1
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" 2>/dev/null | cut -d' ' -f1
    else
        return 1
    fi
}

_dl_utc() { date -u +%Y-%m-%dT%H:%M:%SZ; }

_dl_free_gb() {                     # $1 path — free GB on the volume HOLDING $1
    # Deliberately not `df /`. On macOS 10.15+ the root filesystem is the sealed,
    # read-only system volume; results live on /System/Volumes/Data, which is why
    # all three existing drivers hardcode that path. Asking about the run root
    # itself is correct on every platform and hardcodes nothing.
    #
    # BSD/macOS first (df -g), then GNU (df -BG). Fails OPEN with a sentinel
    # rather than blocking a launch on an unparseable df — refusing to run
    # because we could not read the disk would be a worse failure than running.
    local g
    g=$(df -g "$1" 2>/dev/null | awk 'NR==2{print $4}')
    case "$g" in
        ''|*[!0-9]*)
            g=$(df -BG "$1" 2>/dev/null | awk 'NR==2{gsub(/G$/,"",$4); print $4}')
            ;;
    esac
    case "$g" in
        ''|*[!0-9]*) echo 999999 ;;
        *)           echo "$g" ;;
    esac
}

_dl_is_slug() {                     # lowercase alnum, _ and - ; must not start with - or _
    case "$1" in
        ''|-*|_*)        return 1 ;;
        *[!a-z0-9_-]*)   return 1 ;;
        *)               return 0 ;;
    esac
}

# ── public ──────────────────────────────────────────────────────────────────────

# new_run_dir <kind> <tier> <name>
#
#   kind  test | arm
#   tier  T | A     ('-' derives the canonical tier from kind)
#   name  slug, no prefix. The directory basename is <prefix><name>, where
#         prefix is 'arm' for kind=arm and 't' for kind=test — which reproduces
#         the two names the existing drivers already got right:
#           new_run_dir arm A 15_dofscale   -> .../arm15_dofscale
#           new_run_dir test T 13_twotarget -> .../t13_twotarget
#
# Creates <dir>, <dir>/bin, <dir>/run and <dir>/KIND. Echoes <dir>.
new_run_dir() {
    if [ $# -ne 3 ]; then
        _dl_err "usage: new_run_dir <test|arm> <T|A|-> <name>"
        return 2
    fi

    local kind="$1" tier="$2" name="$3"
    local prefix canon

    case "$kind" in
        test) prefix=t   ; canon=T ;;
        arm)  prefix=arm ; canon=A ;;
        *)    _dl_err "kind must be 'test' or 'arm', got '$kind'"; return 2 ;;
    esac

    [ "$tier" = "-" ] && tier="$canon"
    case "$tier" in
        T|A) : ;;
        *)   _dl_err "tier must be 'T', 'A' or '-', got '$tier'"; return 2 ;;
    esac

    # kind and tier are 1:1 in the convention, so a mismatch is almost always a
    # typo rather than an intention. Refuse by default; let a deliberate case
    # through loudly. If this override is never used, the field is redundant and
    # one of the two should be dropped from the schema.
    if [ "$tier" != "$canon" ] && [ "${FLEXAIDDS_KIND_ALLOW_TIER_MISMATCH:-0}" != "1" ]; then
        _dl_err "kind=$kind implies tier=$canon, got tier=$tier"
        _dl_err "set FLEXAIDDS_KIND_ALLOW_TIER_MISMATCH=1 if that is deliberate"
        return 2
    fi

    if ! _dl_is_slug "$name"; then
        _dl_err "name must match [a-z0-9][a-z0-9_-]*, got '$name'"
        return 2
    fi

    # ── seat ────────────────────────────────────────────────────────────────
    local seat="${FLEXAIDDS_SEAT:-}"
    case "$seat" in
        grok|science|claude-code|dispatch|lp) : ;;
        '') _dl_err "FLEXAIDDS_SEAT unset (grok|science|claude-code|dispatch|lp)"; return 2 ;;
        *)  _dl_err "FLEXAIDDS_SEAT='$seat' is not a known seat"; return 2 ;;
    esac

    # ── engine identity: mandatory, hashed before anything is created ───────
    local engine="${FLEXAIDDS_ENGINE_BIN:-}"
    if [ -z "$engine" ]; then
        _dl_err "FLEXAIDDS_ENGINE_BIN unset — refusing to create an unattributable run dir"
        return 2
    fi
    case "$engine" in
        /*) : ;;
        *)  _dl_err "FLEXAIDDS_ENGINE_BIN must be absolute, got '$engine'"; return 2 ;;
    esac
    if [ ! -f "$engine" ]; then
        _dl_err "FLEXAIDDS_ENGINE_BIN is not a regular file: $engine"
        return 2
    fi
    local engine_sha
    engine_sha=$(_dl_sha256 "$engine")
    if [ ${#engine_sha} -ne 64 ]; then
        _dl_err "could not sha256 $engine (need shasum or sha256sum)"
        return 2
    fi

    # ── root and batch ──────────────────────────────────────────────────────
    local root="${FLEXAIDDS_RUN_ROOT:-${FLEXAIDDS_RESULTS:-}}"
    if [ -z "$root" ]; then
        _dl_err "set FLEXAIDDS_RUN_ROOT (or source driver_preamble.sh for FLEXAIDDS_RESULTS)"
        return 2
    fi
    case "$root" in
        /*) : ;;
        *)  _dl_err "run root must be absolute, got '$root'"; return 2 ;;
    esac
    if [ ! -d "$root" ]; then
        _dl_err "run root is not a directory: $root"
        return 2
    fi

    local batch="${FLEXAIDDS_BATCH_DIR:-}"
    if [ -n "$batch" ]; then
        case "$batch" in
            /*) : ;;
            *)  _dl_err "FLEXAIDDS_BATCH_DIR must be absolute, got '$batch'"; return 2 ;;
        esac
    else
        batch="$root/${name}_$(date +%Y%m%d_%H%M%S)"
    fi
    if ! mkdir -p "$batch" 2>/dev/null; then
        _dl_err "cannot create batch dir: $batch"
        return 2
    fi

    local dir="$batch/${prefix}${name}"

    # ── collision: move aside, never delete ─────────────────────────────────
    # A run directory is evidence. arm15:82 and t13:194 both already say NEVER
    # delete; cachefix says nothing, which is how a re-run silently overwrites
    # the run it was supposed to be compared against.
    if [ -e "$dir" ]; then
        local aside="${dir}_aside_$(date +%s)"
        if ! mv "$dir" "$aside" 2>/dev/null; then
            _dl_err "refusing to touch existing $dir (could not move it aside)"
            return 2
        fi
        _dl_err "existing run dir moved aside -> $aside"
    fi

    if ! mkdir -p "$dir/bin" "$dir/run" 2>/dev/null; then
        _dl_err "cannot create run dir: $dir"
        return 2
    fi

    # ── KIND ────────────────────────────────────────────────────────────────
    # Written atomically: a half-written KIND is worse than none, because a
    # backfill would then decline to replace it.
    #
    # status starts at 'unknown' and is sealed by the driver at end-of-run. It is
    # deliberately NOT 'ok' here: a directory that has produced nothing yet has
    # not earned that claim, and a KIND that reads 'ok' from the moment of mkdir
    # would make `find -name KIND` enumerate intentions rather than results.
    local tmp="$dir/.KIND.$$"
    {
        echo "kind=$kind"
        echo "tier=$tier"
        echo "status=unknown"
        echo "name=$name"
        echo "created=$(_dl_utc)"
        echo "by=$seat"
        echo "engine_sha=$engine_sha"
    } > "$tmp" || { _dl_err "cannot write KIND in $dir"; return 2; }
    mv "$tmp" "$dir/KIND" || { _dl_err "cannot finalise KIND in $dir"; return 2; }

    # ── disk floor ──────────────────────────────────────────────────────────
    local floor="${FLEXAIDDS_DISK_FLOOR_GB:-6}"
    local free; free=$(_dl_free_gb "$dir")
    if [ "$free" -lt "$floor" ]; then
        echo "only ${free} GB free, floor ${floor} — never started, SKIPPED not killed" \
            > "$dir/SKIPPED"
        _dl_err "$dir: ${free} GB < floor ${floor} GB — SKIPPED sentinel written"
        echo "$dir"
        return 3
    fi

    echo "$dir"
    return 0
}
