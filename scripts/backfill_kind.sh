#!/usr/bin/env bash
# backfill_kind.sh — stamp a KIND sidecar onto existing result trees.
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
#
# ── CONTRACT ────────────────────────────────────────────────────────────────────
# DRY RUN BY DEFAULT. Writing requires an explicit --apply.
#
# This script is a SIDECAR WRITER. Per run tree it creates exactly one file,
# named KIND, and it can only ever create — never truncate, never append, never
# replace. It does not write RUN_RECEIPT.json, poses, .mcf, result.csv, DONE, or
# anything else, anywhere, ever.
#
# That is enforced structurally, not by comment:
#
#   * every filesystem write in this script goes through _bk_create_kind()
#   * _bk_create_kind() builds its target as "$tree/KIND" and asserts the
#     basename is literally KIND before touching anything
#   * the write is `set -C` (noclobber) — the open is O_EXCL, so an existing
#     KIND makes the create fail at the syscall, not at a check we might skip
#   * there is no other `>` or `>>` to a file path in the file; every other
#     redirection goes to stdout or stderr
#
# ── LIVE-TREE SAFETY ────────────────────────────────────────────────────────────
# A benchmark batch that is being written right now must not be walked. Three
# independent exclusions, any one of which skips a whole BATCH (not just a tree):
#
#   1. process    a runner whose command line mentions the batch path is alive
#   2. mtime      anything under the batch changed inside the quarantine window
#   3. explicit   --exclude <substring>, repeatable
#
# Batch-level rather than tree-level, because sibling directories of a live run
# are exactly the ones a runner is about to create.
#
# ── WHY status= EXISTS ──────────────────────────────────────────────────────────
# Stamping KIND onto a void tree would launder it into looking like a real run.
# `find -name KIND` is meant to enumerate runs; if enumeration silently implied
# validity, the first void directory to be backfilled would poison every census
# built on it. So enumerable and valid are kept as separate claims: KIND makes a
# tree enumerable, status= says what it is worth.
#
# Backfill defaults to status=unknown and only promotes on positive evidence:
#
#   ok        DONE exists and records rc=0
#   partial   DONE exists and records a non-zero rc, or results exist with no DONE
#   void      an explicit VOID/void marker is present
#   unknown   everything else — including SKIPPED (never started is not void)
#
# engine_sha and by are written as the literal 'unknown' when they cannot be
# recovered. That is the point: an unrecoverable engine identity is the finding,
# and a backfilled record is identifiable precisely because it carries them.
#
# ── USAGE ───────────────────────────────────────────────────────────────────────
#   bash scripts/backfill_kind.sh                        # dry run, every tree
#   bash scripts/backfill_kind.sh --scope campaign       # only <campaign>_<TS>/ trees
#   bash scripts/backfill_kind.sh --exclude gan2vsq5 --exclude ga1jd0
#   bash scripts/backfill_kind.sh --apply                # actually create sidecars
#
#   --root <dir>                default $FLEXAIDDS_RESULTS, else ~/flexaidds_results
#   --scope all|campaign        default all. LP has not ruled on scope; this is the flag.
#   --quarantine-minutes <n>    default 120
#   --exclude <substring>       repeatable
#   --apply                     leave dry-run mode
#
# Exit 0 always unless arguments are bad (2). A tree that cannot be classified is
# reported and skipped, never guessed at.

set -uo pipefail

APPLY=0
SCOPE=all
QUARANTINE_MIN=120
EXCLUDES=""
ROOT="${FLEXAIDDS_RESULTS:-$HOME/flexaidds_results}"

while [ $# -gt 0 ]; do
    case "$1" in
        --apply)              APPLY=1 ;;
        --scope)              SCOPE="${2:-}"; shift ;;
        --root)               ROOT="${2:-}"; shift ;;
        --quarantine-minutes) QUARANTINE_MIN="${2:-}"; shift ;;
        --exclude)            EXCLUDES="$EXCLUDES ${2:-}"; shift ;;
        -h|--help)            sed -n '1,70p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
    shift
done

case "$SCOPE" in all|campaign) : ;; *) echo "--scope must be all|campaign" >&2; exit 2 ;; esac
case "$ROOT" in /*) : ;; *) echo "--root must be absolute" >&2; exit 2 ;; esac
[ -d "$ROOT" ] || { echo "--root is not a directory: $ROOT" >&2; exit 2; }

n_seen=0; n_would=0; n_wrote=0; n_have=0; n_skip_live=0; n_skip_excl=0; n_unclass=0

note() { printf '%s\n' "$*" >&2; }

# ── the one and only writer ─────────────────────────────────────────────────────
# create-only, single named file, basename asserted, O_EXCL via noclobber.
_bk_create_kind() {                 # $1 tree  $2 body
    local tree="$1" body="$2" target
    target="$tree/KIND"
    [ "${target##*/}" = "KIND" ] || { note "REFUSE: target is not a KIND file: $target"; return 1; }
    [ -e "$target" ] && { note "REFUSE: KIND already exists: $target"; return 1; }
    ( set -C; printf '%s\n' "$body" > "$target" ) 2>/dev/null || {
        note "REFUSE: could not create (exists or unwritable): $target"; return 1; }
    return 0
}

# ── exclusions ──────────────────────────────────────────────────────────────────
_bk_excluded() {                    # $1 path
    local p="$1" e
    [ -z "$EXCLUDES" ] && return 1
    for e in $EXCLUDES; do
        case "$p" in *"$e"*) return 0 ;; esac
    done
    return 1
}

_bk_batch_live() {                  # $1 batch dir — process OR recent mtime
    local b="$1"
    if pgrep -f "$b" >/dev/null 2>&1; then
        return 0
    fi
    # any file touched inside the quarantine window means a writer is mid-flight
    if [ -n "$(find "$b" -mmin "-$QUARANTINE_MIN" -print -quit 2>/dev/null)" ]; then
        return 0
    fi
    return 1
}

# ── read-only inference ─────────────────────────────────────────────────────────
_bk_first_match() {                 # $1 file  $2 key= -> value on stdout, empty if none
    [ -r "$1" ] || return 0
    grep -m1 -o "$2[0-9a-zA-Z._/-]*" "$1" 2>/dev/null | head -1 | sed "s|^$2||"
}

_bk_infer_kind() {                  # $1 basename -> "kind tier" or empty
    case "$1" in
        arm[0-9]*|arm_*)   echo "arm A" ;;
        run_t[0-9]*|t[0-9]*) echo "test T" ;;
        S[0-9]_*|s[0-9]_*) echo "arm A" ;;
        *) echo "" ;;
    esac
}

_bk_infer_status() {                # $1 tree
    local t="$1"
    if [ -e "$t/VOID" ] || [ -e "$t/void" ]; then echo void; return; fi
    if [ -r "$t/DONE" ]; then
        if grep -q 'rc=0' "$t/DONE" 2>/dev/null; then echo ok; else echo partial; fi
        return
    fi
    if [ -e "$t/SKIPPED" ]; then echo unknown; return; fi
    if [ -n "$(find "$t" -name result.csv -print -quit 2>/dev/null)" ]; then echo partial; return; fi
    echo unknown
}

_bk_infer_engine_sha() {            # $1 tree
    local t="$1" v=""
    v=$(_bk_first_match "$t/provenance.txt" "engine_sha256=")
    if [ -z "$v" ] && [ -r "$t/RUN_RECEIPT.json" ]; then
        v=$(grep -m1 -o '"binary_sha256"[^,]*' "$t/RUN_RECEIPT.json" 2>/dev/null \
            | grep -o '[0-9a-f]\{64\}' | head -1)
    fi
    [ ${#v} -eq 64 ] && echo "$v" || echo unknown
}

_bk_infer_by() {                    # $1 tree
    local t="$1"
    if grep -qi 'claude-science' "$t/claim.log" 2>/dev/null; then echo science; return; fi
    if grep -qi 'claude-science' "$t/provenance.txt" 2>/dev/null; then echo science; return; fi
    if grep -qi 'grok' "$t/provenance.txt" 2>/dev/null; then echo grok; return; fi
    echo unknown
}

_bk_created() {                     # $1 tree — prefer the _YYYYmmdd_HHMMSS in the batch name
    local t="$1" stamp
    stamp=$(printf '%s' "$t" | grep -o '_[0-9]\{8\}_[0-9]\{6\}' | tail -1)
    if [ -n "$stamp" ]; then
        printf '%s-%s-%sT%s:%s:%sZ' \
            "$(echo "$stamp" | cut -c2-5)"   "$(echo "$stamp" | cut -c6-7)" \
            "$(echo "$stamp" | cut -c8-9)"   "$(echo "$stamp" | cut -c11-12)" \
            "$(echo "$stamp" | cut -c13-14)" "$(echo "$stamp" | cut -c15-16)"
        return
    fi
    date -u -r "$t" +%Y-%m-%dT%H:%M:%SZ 2>/dev/null || echo unknown
}

# ── walk ────────────────────────────────────────────────────────────────────────
# A "tree" is a directory containing a run/ subdir — the shape every driver agrees
# on even when it agrees on nothing else.
note "root=$ROOT scope=$SCOPE quarantine=${QUARANTINE_MIN}m apply=$APPLY"
note ""

for batch in "$ROOT"/*/; do
    batch="${batch%/}"
    [ -d "$batch" ] || continue

    if _bk_excluded "$batch"; then
        n_skip_excl=$((n_skip_excl+1)); note "SKIP excluded : $(basename "$batch")"; continue
    fi
    if _bk_batch_live "$batch"; then
        n_skip_live=$((n_skip_live+1)); note "SKIP LIVE     : $(basename "$batch") (process or mtime inside quarantine)"; continue
    fi

    for tree in "$batch" "$batch"/*/; do
        tree="${tree%/}"
        [ -d "$tree/run" ] || continue
        if [ "$SCOPE" = campaign ]; then
            case "$(basename "$batch")" in *_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9][0-9][0-9]) : ;; *) continue ;; esac
        fi

        n_seen=$((n_seen+1))
        rel="${tree#"$ROOT"/}"

        if [ -e "$tree/KIND" ]; then
            n_have=$((n_have+1)); note "HAVE KIND     : $rel"; continue
        fi

        kt=$(_bk_infer_kind "$(basename "$tree")")
        if [ -z "$kt" ]; then
            n_unclass=$((n_unclass+1)); note "UNCLASSIFIED  : $rel (name matches no kind rule — not guessing)"; continue
        fi
        kind=${kt%% *}; tier=${kt##* }

        body="kind=$kind
tier=$tier
status=$(_bk_infer_status "$tree")
name=$(basename "$tree")
created=$(_bk_created "$tree")
by=$(_bk_infer_by "$tree")
engine_sha=$(_bk_infer_engine_sha "$tree")"

        if [ "$APPLY" -eq 1 ]; then
            if _bk_create_kind "$tree" "$body"; then
                n_wrote=$((n_wrote+1)); note "WROTE         : $rel"
            fi
        else
            n_would=$((n_would+1))
            note "WOULD WRITE   : $rel"
            printf '%s\n' "$body" | sed 's/^/                  /' >&2
        fi
    done
done

note ""
note "seen=$n_seen  already_have=$n_have  unclassified=$n_unclass"
note "skipped_live=$n_skip_live  skipped_excluded=$n_skip_excl"
if [ "$APPLY" -eq 1 ]; then
    note "wrote=$n_wrote"
else
    note "would_write=$n_would   (dry run — pass --apply to create sidecars)"
fi
exit 0
