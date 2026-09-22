#!/usr/bin/env bash
#
# campaign_preflight.sh -- refuse to launch when a dependency will not
# outlive the run.
#
# WHY THIS EXISTS
# ---------------
# On 2026-09-20 nine cells of an 85-target Astex campaign died at startup
# because $TMPDIR pointed into an agent session workspace and that harness
# swept the directory after a few hours idle. Nothing about the engine, the
# ligands or the receptors was wrong; a dependency simply stopped existing
# halfway through. Every check below is one that would have caught a real
# failure in this project, and none is hypothetical.
#
# The engine no longer aborts on a vanished temp dir (LIB/temp_dir.{h,cpp}),
# but degrading gracefully into /tmp is not the same as launching a 20-hour
# campaign whose scratch, output root or cache sits on borrowed ground. This
# script is the gate before the launcher, not a substitute for the fix.
#
# THE RULE THAT COST A DAY: RESOLVE BY EVAL, NEVER BY sed
# -------------------------------------------------------
# A preflight written earlier today sed-substituted a launcher's text to find
# the engine path. It expanded $R but not $BB, tested the literal string
# "$BB/FlexAIDdS", found no such file, and reported the engine ABSENT while
# the engine was sitting right there. A check that reports on its own string
# handling is worse than no check, because it spends the credibility of a
# real check on noise. Path values are therefore expanded with eval, and
# launcher variables are recovered by eval'ing assignment lines in order --
# never by pattern-substituting text.
#
# NO ABSOLUTE PATHS ARE HARDCODED. Everything comes from arguments or the
# environment, so this runs on any machine.
#
# Exit: 0 all checks passed, 1 one or more failed, 2 usage error.
#
# Written for bash 3.2.57 (the macOS system bash): no ${var,,}, no mapfile,
# no declare -A, no wait -n, and PIPESTATUS rather than $? after a pipeline.

set -u

PROG="$(basename "$0")"

# ---------------------------------------------------------------------------
# Tunables (override from the environment)
# ---------------------------------------------------------------------------
: "${PREFLIGHT_MIN_FREE_GB:=20}"
: "${PREFLIGHT_EXPECT_TARGETS:=}"
: "${PREFLIGHT_TARGET_NAME_LEN:=4}"

# Engine data files that must sit beside the binary.
ENGINE_DATA_FILES="MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def"

# An "agent session workspace" is ephemeral: a harness sweeps it after a few
# hours idle. Both markers must be present, so that a durable directory that
# merely contains the word "workspaces" is not condemned.
EPHEMERAL_MARKER_A="/.claude-science/"
EPHEMERAL_MARKER_B="/workspaces/"

N_PASS=0
N_FAIL=0
N_SKIP=0
REQUIRE_ALL=0
FAILED_NAMES=""

# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
_pass() { N_PASS=$((N_PASS + 1)); printf 'PASS  %-22s %s\n' "$1" "$2"; }
_fail() {
    N_FAIL=$((N_FAIL + 1))
    FAILED_NAMES="$FAILED_NAMES $1"
    printf 'FAIL  %-22s %s\n' "$1" "$2"
}
_skip() {
    N_SKIP=$((N_SKIP + 1))
    printf 'SKIP  %-22s %s\n' "$1" "$2"
    if [ "$REQUIRE_ALL" -eq 1 ]; then
        N_FAIL=$((N_FAIL + 1))
        FAILED_NAMES="$FAILED_NAMES $1(skipped)"
    fi
}

usage() {
    cat <<EOF
$PROG -- refuse to launch a campaign whose dependencies will not outlive it.

Usage:
  $PROG [options]

Options (all paths may contain \$VARIABLES; they are expanded with eval):
  --engine PATH        engine binary (data files must sit beside it)
  --output-root PATH   campaign output root
  --cache PATH         target cache directory
  --expect-targets N   exact number of ${PREFLIGHT_TARGET_NAME_LEN}-char target dirs the cache must hold
  --min-free-gb N      free-space floor on the output root (default ${PREFLIGHT_MIN_FREE_GB})
  --env-file PATH      file of KEY=VALUE lines, sourced before resolution
  --launcher PATH      launcher script; its simple assignments are eval'd in
                       order to recover \$R, \$BB, \$ENGINE and friends
  --require-all        treat a SKIP as a failure
  -h, --help           this text

Environment equivalents: FLEXAIDDS_ENGINE, CAMPAIGN_OUTPUT_ROOT, CAMPAIGN_CACHE,
PREFLIGHT_EXPECT_TARGETS, PREFLIGHT_MIN_FREE_GB.

Exit: 0 pass, 1 failure, 2 usage.
EOF
}

# ---------------------------------------------------------------------------
# Resolution by eval -- the whole point of the script's contract
# ---------------------------------------------------------------------------

# resolve_path VALUE
# Expand $VAR / ${VAR} references using the CURRENT environment. Command
# substitution is refused rather than executed: a preflight must not run
# arbitrary code found in a path string.
resolve_path() {
    _rp_raw="$1"
    case "$_rp_raw" in
        *'$('*|*'`'*)
            printf '%s' "$_rp_raw"
            return 1
            ;;
    esac
    eval "_rp_out=\"$_rp_raw\"" 2>/dev/null || {
        printf '%s' "$_rp_raw"
        return 1
    }
    printf '%s' "$_rp_out"
    return 0
}

# load_launcher_vars PATH
# Recover a launcher's variables WITHOUT running the launcher. Only simple
# assignment lines are eval'd, in file order, so that BB=$R/build resolves
# against an R that was eval'd a moment earlier -- which plain text
# substitution cannot do. Lines carrying command substitution, pipes,
# redirection or command chaining are skipped and counted, never executed.
load_launcher_vars() {
    _ll_file="$1"
    _ll_taken=0
    _ll_skipped=0
    _ll_dropped=0
    if [ ! -f "$_ll_file" ]; then
        printf 'launcher not found: %s\n' "$_ll_file" >&2
        return 1
    fi
    while IFS= read -r _ll_line || [ -n "$_ll_line" ]; do
        # ORDER IS LOAD-BEARING. `case` takes the first matching pattern, and
        # [A-Za-z_]*=* matches "export FOO=bar" too -- the leading * happily
        # swallows "xport FOO". With the plain-assignment pattern first, the
        # export branch is unreachable, the "export " prefix survives into
        # the name, the identifier check then sees a space and drops the line
        # WITHOUT counting it. Every exported launcher variable vanished in
        # silence. The export pattern must be tested first.
        case "$_ll_line" in
            'export '[A-Za-z_]*=*) _ll_line="${_ll_line#export }" ;;
            [A-Za-z_]*=*) : ;;
            *) continue ;;
        esac
        case "$_ll_line" in
            *'$('*|*'`'*|*';'*|*'|'*|*'&'*|*'>'*|*'<'*) 
                _ll_skipped=$((_ll_skipped + 1))
                continue
                ;;
        esac
        _ll_name="${_ll_line%%=*}"
        case "$_ll_name" in
            # Not an identifier after all. Counted, never silent: an
            # uncounted drop is exactly how the export bug stayed hidden.
            *[!A-Za-z0-9_]*) _ll_dropped=$((_ll_dropped + 1)); continue ;;
        esac
        if eval "$_ll_line" 2>/dev/null; then
            eval "export $_ll_name" 2>/dev/null
            _ll_taken=$((_ll_taken + 1))
        else
            _ll_dropped=$((_ll_dropped + 1))
        fi
    done < "$_ll_file"
    printf 'launcher %s: %d assignments eval'"'"'d, %d skipped (command substitution or chaining), %d unparsable\n' \
        "$(basename "$_ll_file")" "$_ll_taken" "$_ll_skipped" "$_ll_dropped"
    return 0
}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

# check_ephemeral LABEL PATH
# A path under an agent session workspace is swept after hours of idling. The
# nine dead cells of 2026-09-20 are exactly this.
check_ephemeral() {
    _ce_label="$1"
    _ce_path="$2"
    if [ -z "$_ce_path" ]; then
        _skip "ephemeral:$_ce_label" "not set"
        return 0
    fi
    case "$_ce_path" in
        *"$EPHEMERAL_MARKER_A"*)
            case "$_ce_path" in
                *"$EPHEMERAL_MARKER_B"*)
                    _fail "ephemeral:$_ce_label" \
                        "lives in an agent session workspace (swept when idle): $_ce_path"
                    return 1
                    ;;
            esac
            ;;
    esac
    _pass "ephemeral:$_ce_label" "durable: $_ce_path"
    return 0
}

# check_exists_dir LABEL PATH
check_exists_dir() {
    _cd_label="$1"
    _cd_path="$2"
    if [ -z "$_cd_path" ]; then
        _skip "exists:$_cd_label" "not set"
        return 0
    fi
    if [ -d "$_cd_path" ]; then
        _pass "exists:$_cd_label" "$_cd_path"
        return 0
    fi
    _fail "exists:$_cd_label" "directory missing: $_cd_path"
    return 1
}

# check_engine PATH -- binary present and executable, data files beside it.
check_engine() {
    _ceg_bin="$1"
    if [ -z "$_ceg_bin" ]; then
        _skip "engine" "no --engine given"
        return 0
    fi
    if [ ! -e "$_ceg_bin" ]; then
        _fail "engine" "binary not found: $_ceg_bin"
        return 1
    fi
    if [ ! -x "$_ceg_bin" ]; then
        _fail "engine" "binary not executable: $_ceg_bin"
        return 1
    fi
    _pass "engine" "executable: $_ceg_bin"

    _ceg_dir="$(dirname "$_ceg_bin")"
    _ceg_missing=""
    for _ceg_f in $ENGINE_DATA_FILES; do
        if [ ! -f "$_ceg_dir/$_ceg_f" ]; then
            _ceg_missing="$_ceg_missing $_ceg_f"
        fi
    done
    if [ -n "$_ceg_missing" ]; then
        _fail "engine-data" "missing beside binary:$_ceg_missing (in $_ceg_dir)"
        return 1
    fi
    _pass "engine-data" "all present in $_ceg_dir"
    return 0
}

# check_cache PATH EXPECTED -- exactly EXPECTED N-char target dirs, zero strays.
check_cache() {
    _cc_dir="$1"
    _cc_expect="$2"
    if [ -z "$_cc_dir" ]; then
        _skip "cache" "no --cache given"
        return 0
    fi
    if [ ! -d "$_cc_dir" ]; then
        _fail "cache" "cache directory missing: $_cc_dir"
        return 1
    fi

    _cc_len="$PREFLIGHT_TARGET_NAME_LEN"
    _cc_re="^[0-9A-Za-z]{$_cc_len}\$"
    _cc_targets=0
    _cc_strays=""
    for _cc_e in "$_cc_dir"/*; do
        [ -e "$_cc_e" ] || continue
        _cc_b="$(basename "$_cc_e")"
        if [ -d "$_cc_e" ] && printf '%s' "$_cc_b" | grep -Eq "$_cc_re"; then
            _cc_targets=$((_cc_targets + 1))
        else
            _cc_strays="$_cc_strays $_cc_b"
        fi
    done

    if [ -n "$_cc_strays" ]; then
        _fail "cache-strays" "non-target entries in cache:$_cc_strays"
    else
        _pass "cache-strays" "no strays in $_cc_dir"
    fi

    if [ -z "$_cc_expect" ]; then
        _skip "cache-count" "found $_cc_targets target dirs; no --expect-targets to compare"
        return 0
    fi
    if [ "$_cc_targets" -eq "$_cc_expect" ]; then
        _pass "cache-count" "$_cc_targets target dirs (expected $_cc_expect)"
        return 0
    fi
    _fail "cache-count" "cache holds $_cc_targets target dirs, expected $_cc_expect"
    return 1
}

# check_disk PATH MIN_GB -- headroom on the filesystem holding PATH.
# df -P keeps the output on one line even when the device name is long;
# plain df wraps and shifts the column, which would silently misread.
check_disk() {
    _cdk_path="$1"
    _cdk_min="$2"
    if [ -z "$_cdk_path" ]; then
        _skip "disk" "no path to measure"
        return 0
    fi
    if [ ! -e "$_cdk_path" ]; then
        _fail "disk" "cannot measure free space, path missing: $_cdk_path"
        return 1
    fi
    _cdk_avail_k="$(df -Pk "$_cdk_path" 2>/dev/null | awk 'NR==2 {print $4}')"
    if [ -z "$_cdk_avail_k" ]; then
        _fail "disk" "could not parse df output for $_cdk_path"
        return 1
    fi
    _cdk_avail_gb=$((_cdk_avail_k / 1048576))
    if [ "$_cdk_avail_gb" -lt "$_cdk_min" ]; then
        _fail "disk" "${_cdk_avail_gb} GB free on $_cdk_path, floor is ${_cdk_min} GB"
        return 1
    fi
    _pass "disk" "${_cdk_avail_gb} GB free on $_cdk_path (floor ${_cdk_min} GB)"
    return 0
}

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

ARG_ENGINE="${FLEXAIDDS_ENGINE:-}"
ARG_OUTPUT_ROOT="${CAMPAIGN_OUTPUT_ROOT:-}"
ARG_CACHE="${CAMPAIGN_CACHE:-}"
ARG_ENV_FILE=""
ARG_LAUNCHER=""

while [ $# -gt 0 ]; do
    case "$1" in
        --engine)         ARG_ENGINE="${2:-}"; shift 2 ;;
        --output-root)    ARG_OUTPUT_ROOT="${2:-}"; shift 2 ;;
        --cache)          ARG_CACHE="${2:-}"; shift 2 ;;
        --expect-targets) PREFLIGHT_EXPECT_TARGETS="${2:-}"; shift 2 ;;
        --min-free-gb)    PREFLIGHT_MIN_FREE_GB="${2:-}"; shift 2 ;;
        --env-file)       ARG_ENV_FILE="${2:-}"; shift 2 ;;
        --launcher)       ARG_LAUNCHER="${2:-}"; shift 2 ;;
        --require-all)    REQUIRE_ALL=1; shift ;;
        -h|--help)        usage; exit 0 ;;
        *) printf '%s: unknown argument: %s\n' "$PROG" "$1" >&2; usage >&2; exit 2 ;;
    esac
done

echo "======================================================================"
echo "CAMPAIGN PREFLIGHT"
echo "======================================================================"

# Sourcing happens BEFORE resolution so that a path written as "$BB/FlexAIDdS"
# resolves against the launcher's own BB.
if [ -n "$ARG_ENV_FILE" ]; then
    if [ -f "$ARG_ENV_FILE" ]; then
        # shellcheck disable=SC1090
        . "$ARG_ENV_FILE"
        echo "env-file sourced: $ARG_ENV_FILE"
    else
        printf '%s: --env-file not found: %s\n' "$PROG" "$ARG_ENV_FILE" >&2
        exit 2
    fi
fi

if [ -n "$ARG_LAUNCHER" ]; then
    load_launcher_vars "$ARG_LAUNCHER" || exit 2
fi

# Late defaults: a launcher or env-file may define these.
[ -z "$ARG_ENGINE" ] && ARG_ENGINE="${ENGINE:-}"
[ -z "$ARG_OUTPUT_ROOT" ] && ARG_OUTPUT_ROOT="${OUTROOT:-}"
[ -z "$ARG_CACHE" ] && ARG_CACHE="${CACHE:-}"

RES_ENGINE="$(resolve_path "$ARG_ENGINE")" || \
    echo "NOTE  engine path contains command substitution; left unexpanded"
RES_OUTPUT_ROOT="$(resolve_path "$ARG_OUTPUT_ROOT")" || \
    echo "NOTE  output root contains command substitution; left unexpanded"
RES_CACHE="$(resolve_path "$ARG_CACHE")" || \
    echo "NOTE  cache path contains command substitution; left unexpanded"

RES_ENGINE_DIR=""
[ -n "$RES_ENGINE" ] && RES_ENGINE_DIR="$(dirname "$RES_ENGINE")"

echo "engine      : ${RES_ENGINE:-<unset>}"
echo "output root : ${RES_OUTPUT_ROOT:-<unset>}"
echo "cache       : ${RES_CACHE:-<unset>}"
echo "TMPDIR      : ${TMPDIR:-<unset>}"
echo "----------------------------------------------------------------------"

# --- ephemerality: the 2026-09-20 failure mode --------------------------------
check_ephemeral "TMPDIR"    "${TMPDIR:-}"
check_ephemeral "TMP"       "${TMP:-}"
check_ephemeral "TEMP"      "${TEMP:-}"
check_ephemeral "FLEXAIDDS_TMPDIR" "${FLEXAIDDS_TMPDIR:-}"
check_ephemeral "output"    "$RES_OUTPUT_ROOT"
check_ephemeral "engine"    "$RES_ENGINE_DIR"
check_ephemeral "cache"     "$RES_CACHE"

# --- existence and shape ------------------------------------------------------
check_exists_dir "output" "$RES_OUTPUT_ROOT"
check_engine "$RES_ENGINE"
check_cache "$RES_CACHE" "$PREFLIGHT_EXPECT_TARGETS"

# --- headroom -----------------------------------------------------------------
if [ -n "$RES_OUTPUT_ROOT" ]; then
    check_disk "$RES_OUTPUT_ROOT" "$PREFLIGHT_MIN_FREE_GB"
else
    _skip "disk" "no output root to measure"
fi

echo "----------------------------------------------------------------------"
printf 'preflight: %d passed, %d failed, %d skipped\n' "$N_PASS" "$N_FAIL" "$N_SKIP"
if [ "$N_FAIL" -gt 0 ]; then
    printf 'REFUSING TO LAUNCH. Failed:%s\n' "$FAILED_NAMES"
    exit 1
fi
echo "OK to launch."
exit 0
