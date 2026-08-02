#!/usr/bin/env bash
# =============================================================================
# bench_flexaid_fast.sh — P0 build-level speedup + parity harness for the
#                         classic FlexAID docking target.
#
# WHAT IT DOES
#   1. Configures + builds the `FlexAID` target TWICE from the current git HEAD:
#        - baseline : -DBUILD_FLEXAID_FAST=OFF -DFLEXAID_PGO=off  (plain -O3 -ffast-math)
#        - fast     : -DBUILD_FLEXAID_FAST=ON                     (LTO/IPO + -march=native
#                                                                  + -DNDEBUG + strip)
#      and records each build's wall-clock time.
#   2. Runs ONE reference dock with each binary (fixed FLEXAID_SEED) and records
#      the docking wall-clock time + speedup.
#   3. Parity check: extracts best-CF (lowest `REMARK CF=` across output poses)
#      and the top-pose file from each run and compares them. Under the repo's
#      drift policy the two builds must stay rank/pose-equivalent, NOT bit-identical
#      — so the check passes when |ΔCF| <= FLEXAID_BENCH_CF_TOL (default 0.50).
#
#   This script NEVER fabricates a number. If a phase cannot run (no reference
#   input, binary missing, no CF in output) it says so and reports N/A.
#
# REFERENCE INPUT (a real receptor + ligand is required to run the dock phase)
#   The classic `FlexAID` binary does NOT support --redock (that is FlexAIDdS-only),
#   so it needs an on-disk receptor + cognate ligand. Resolution order:
#     1. $FLEXAID_BENCH_RECEPTOR + $FLEXAID_BENCH_LIGAND         (explicit override)
#     2. Astex 1G9V staged by scripts/reproduce_astex85.sh:
#          $ASTEX_DIR/1G9V/1G9V_apo.pdb  +  $ASTEX_DIR/1G9V/1G9V_ligand.{sdf,mol2}
#          (ASTEX_DIR defaults to ~/FlexAIDdS_reviewer_benchmark/astex_diverse)
#     3. none found -> build timings are still reported; the dock+parity phase is
#        SKIPPED with instructions to stage inputs first.
#
#   To stage the 1G9V reference set without hand-collecting files, run the Astex
#   reproducer once (it downloads + preps apo receptors + cognate ligands):
#       bash scripts/reproduce_astex85.sh          # full 85-target run, or
#       # stop it after the "Astex download/prep" step to just stage inputs.
#   then point this script at the staged files via $ASTEX_DIR or the explicit
#   $FLEXAID_BENCH_RECEPTOR / $FLEXAID_BENCH_LIGAND overrides.
#
# USAGE
#   bash scripts/bench_flexaid_fast.sh [--build-only] [--jobs N] [--keep]
#
# ENV OVERRIDES
#   FLEXAID_BENCH_RECEPTOR   explicit receptor file (.pdb/.cif)
#   FLEXAID_BENCH_LIGAND     explicit ligand file (.sdf/.mol2/.mol/.pdb)
#   ASTEX_DIR                Astex staging dir (default ~/FlexAIDdS_reviewer_benchmark/astex_diverse)
#   FLEXAID_SEED             fixed RNG seed for the dock (default 42)
#   FLEXAID_BENCH_CF_TOL     |ΔCF| tolerance for parity PASS (default 0.50)
#   FLEXAID_BENCH_EXTRA_ARGS extra args passed verbatim to the FlexAID dock
#   FLEXAID_BENCH_TIMEOUT    per-dock timeout in seconds (default 3600)
#   FLEXAID_BENCH_OUT        results/log dir (default <repo>/WRK/bench_flexaid_fast)
# =============================================================================
set -euo pipefail

# ── repo-relative paths (no machine-specific absolutes; AGENTS.md hygiene) ────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
if command -v git >/dev/null 2>&1; then
    _top="$(git -C "${REPO_ROOT}" rev-parse --show-toplevel 2>/dev/null || true)"
    [[ -n "${_top}" ]] && REPO_ROOT="${_top}"
fi

# ── colour helpers ───────────────────────────────────────────────────────────
if [[ -t 1 ]]; then
    RED='\033[0;31m'; GRN='\033[0;32m'; YLW='\033[0;33m'; CYN='\033[0;36m'; BLD='\033[1m'; RST='\033[0m'
else
    RED=''; GRN=''; YLW=''; CYN=''; BLD=''; RST=''
fi
info(){ printf "${CYN}[INFO]${RST} %s\n" "$*"; }
ok(){   printf "${GRN}[OK]${RST}   %s\n" "$*"; }
warn(){ printf "${YLW}[WARN]${RST} %s\n" "$*"; }
die(){  printf "${RED}[FATAL]${RST} %s\n" "$*" >&2; exit 1; }
banner(){ printf "\n${BLD}${CYN}== %s ==${RST}\n" "$*"; }

# ── args ─────────────────────────────────────────────────────────────────────
BUILD_ONLY=0
KEEP=0
JOBS="$( (command -v nproc >/dev/null && nproc) || sysctl -n hw.ncpu 2>/dev/null || echo 4)"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-only) BUILD_ONLY=1; shift;;
        --keep)       KEEP=1; shift;;
        --jobs)       JOBS="${2:?}"; shift 2;;
        -h|--help)    sed -n '2,60p' "${BASH_SOURCE[0]}"; exit 0;;
        *) die "unknown argument: $1";;
    esac
done

# ── config ───────────────────────────────────────────────────────────────────
OUT_DIR="${FLEXAID_BENCH_OUT:-${REPO_ROOT}/WRK/bench_flexaid_fast}"
BASELINE_BUILD="${OUT_DIR}/build-baseline"
FAST_BUILD="${OUT_DIR}/build-fast"
SEED="${FLEXAID_SEED:-42}"
CF_TOL="${FLEXAID_BENCH_CF_TOL:-0.50}"
DOCK_TIMEOUT="${FLEXAID_BENCH_TIMEOUT:-3600}"
ASTEX_DIR="${ASTEX_DIR:-${HOME}/FlexAIDdS_reviewer_benchmark/astex_diverse}"
REPORT="${OUT_DIR}/bench_report.txt"
mkdir -p "${OUT_DIR}"
: > "${REPORT}"

log(){ printf '%s\n' "$*" | tee -a "${REPORT}"; }

# ── portable wall-clock timer ────────────────────────────────────────────────
now(){ date +%s.%N; }
elapsed(){ awk -v a="$1" -v b="$2" 'BEGIN{printf "%.2f", b-a}'; }

# ── build one variant ────────────────────────────────────────────────────────
# $1 = build dir, $2 = human label, $3.. = extra cmake -D flags
build_variant(){
    local dir="$1" label="$2"; shift 2
    banner "Building FlexAID (${label})"
    local t0 t1
    t0="$(now)"
    cmake -S "${REPO_ROOT}" -B "${dir}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DFLEXAIDS_BUILD_CORE=ON \
        "$@" > "${dir}.cmake.log" 2>&1 \
        || { tail -40 "${dir}.cmake.log"; die "cmake configure failed (${label}) — see ${dir}.cmake.log"; }
    cmake --build "${dir}" --target FlexAID -j "${JOBS}" > "${dir}.build.log" 2>&1 \
        || { tail -60 "${dir}.build.log"; die "build failed (${label}) — see ${dir}.build.log"; }
    t1="$(now)"
    local secs; secs="$(elapsed "${t0}" "${t1}")"
    local bin="${dir}/FlexAID"
    [[ -x "${bin}" ]] || die "FlexAID binary not produced (${label}): ${bin}"
    printf '%s' "${secs}"
}

banner "FlexAID-Fast P0 build + parity harness"
log "repo_root        : ${REPO_ROOT}"
log "git_head         : $(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
log "jobs             : ${JOBS}"
log "seed             : ${SEED}"
log "cf_tolerance     : ${CF_TOL}"
log ""

BASE_SECS="$(build_variant "${BASELINE_BUILD}" "baseline (-O3 -ffast-math, BUILD_FLEXAID_FAST=OFF)" \
                -DBUILD_FLEXAID_FAST=OFF -DFLEXAID_PGO=off)"
ok "baseline build: ${BASE_SECS}s"
FAST_SECS="$(build_variant "${FAST_BUILD}" "fast (LTO/IPO + native + NDEBUG, BUILD_FLEXAID_FAST=ON)" \
                -DBUILD_FLEXAID_FAST=ON)"
ok "fast build:     ${FAST_SECS}s"

BASE_BIN="${BASELINE_BUILD}/FlexAID"
FAST_BIN="${FAST_BUILD}/FlexAID"
BASE_SZ="$(stat -c%s "${BASE_BIN}" 2>/dev/null || stat -f%z "${BASE_BIN}" 2>/dev/null || echo '?')"
FAST_SZ="$(stat -c%s "${FAST_BIN}" 2>/dev/null || stat -f%z "${FAST_BIN}" 2>/dev/null || echo '?')"

log ""
log "─── BUILD RESULTS ────────────────────────────────────────────"
log "baseline build wall-time : ${BASE_SECS}s   binary=${BASE_SZ} bytes"
log "fast     build wall-time : ${FAST_SECS}s   binary=${FAST_SZ} bytes"
log "(fast build is expected to be SLOWER to compile — LTO trades build time for run time)"

if [[ "${BUILD_ONLY}" -eq 1 ]]; then
    ok "--build-only: skipping dock + parity phase"
    log ""
    log "Both FlexAID variants built successfully. Report: ${REPORT}"
    exit 0
fi

# ── resolve reference input ──────────────────────────────────────────────────
RECEPTOR="${FLEXAID_BENCH_RECEPTOR:-}"
LIGAND="${FLEXAID_BENCH_LIGAND:-}"
if [[ -z "${RECEPTOR}" || -z "${LIGAND}" ]]; then
    _apo="${ASTEX_DIR}/1G9V/1G9V_apo.pdb"
    if [[ -f "${_apo}" ]]; then
        RECEPTOR="${_apo}"
        for ext in sdf mol2 mol pdb; do
            for cand in "${ASTEX_DIR}/1G9V/1G9V_ligand.${ext}" "${ASTEX_DIR}/1G9V/1G9V_lig.${ext}"; do
                [[ -f "${cand}" ]] && { LIGAND="${cand}"; break 2; }
            done
        done
    fi
fi

if [[ -z "${RECEPTOR}" || ! -f "${RECEPTOR}" || -z "${LIGAND}" || ! -f "${LIGAND}" ]]; then
    banner "Reference dock: SKIPPED (no local input set found)"
    warn "No runnable receptor + ligand pair was found."
    log ""
    log "─── DOCK + PARITY: SKIPPED ───────────────────────────────────"
    log "The classic FlexAID binary needs an on-disk receptor + cognate ligand"
    log "(it does not support --redock, which is FlexAIDdS-only, and RCSB download"
    log " is unavailable in this environment)."
    log ""
    log "To stage the 1G9V reference set, run the Astex reproducer once:"
    log "    bash scripts/reproduce_astex85.sh"
    log "then re-run this script (it auto-detects \$ASTEX_DIR/1G9V/1G9V_apo.pdb +"
    log "cognate ligand), or pass an explicit pair:"
    log "    FLEXAID_BENCH_RECEPTOR=path/to/receptor.pdb \\"
    log "    FLEXAID_BENCH_LIGAND=path/to/ligand.sdf \\"
    log "    bash scripts/bench_flexaid_fast.sh"
    ok "Build phase complete; dock/parity deferred. Report: ${REPORT}"
    exit 0
fi

log ""
log "reference receptor : ${RECEPTOR}"
log "reference ligand   : ${LIGAND}"

# ── run one dock ─────────────────────────────────────────────────────────────
# echoes wall-seconds to stdout; writes poses under $out_prefix*
run_dock(){
    local bin="$1" out_prefix="$2" logf="$3"
    local t0 t1
    t0="$(now)"
    # shellcheck disable=SC2086
    FLEXAID_SEED="${SEED}" timeout "${DOCK_TIMEOUT}" \
        "${bin}" "${RECEPTOR}" "${LIGAND}" -o "${out_prefix}" \
        ${FLEXAID_BENCH_EXTRA_ARGS:-} > "${logf}" 2>&1 || return $?
    t1="$(now)"
    elapsed "${t0}" "${t1}"
}

# ── extract best CF (lowest REMARK CF=) + owning pose file from an output dir ──
# prints "<best_cf> <top_pose_file>" or "NA NA" if none found. No fabrication.
extract_best_cf(){
    local dir="$1"
    local best_cf="" best_file=""
    while IFS= read -r f; do
        # skip the pre-optimisation initial pose
        [[ "${f}" == *_INI.pdb ]] && continue
        local cf
        cf="$(grep -aoE 'REMARK CF=[ ]*-?[0-9]+\.[0-9]+' "${f}" 2>/dev/null | head -1 | grep -oE '\-?[0-9]+\.[0-9]+' || true)"
        [[ -z "${cf}" ]] && continue
        if [[ -z "${best_cf}" ]] || awk -v a="${cf}" -v b="${best_cf}" 'BEGIN{exit !(a<b)}'; then
            best_cf="${cf}"; best_file="${f}"
        fi
    done < <(find "${dir}" -maxdepth 1 -name '*.pdb' 2>/dev/null)
    if [[ -z "${best_cf}" ]]; then printf 'NA NA'; else printf '%s %s' "${best_cf}" "$(basename "${best_file}")"; fi
}

banner "Reference dock — baseline"
BASE_OUT="${OUT_DIR}/dock-baseline"; mkdir -p "${BASE_OUT}"
if BASE_DOCK_SECS="$(run_dock "${BASE_BIN}" "${BASE_OUT}/1G9V" "${BASE_OUT}/dock.log")"; then
    ok "baseline dock: ${BASE_DOCK_SECS}s"
else
    rc=$?; BASE_DOCK_SECS="FAILED(rc=${rc})"; warn "baseline dock did not complete (rc=${rc}); see ${BASE_OUT}/dock.log"
fi

banner "Reference dock — fast"
FAST_OUT="${OUT_DIR}/dock-fast"; mkdir -p "${FAST_OUT}"
if FAST_DOCK_SECS="$(run_dock "${FAST_BIN}" "${FAST_OUT}/1G9V" "${FAST_OUT}/dock.log")"; then
    ok "fast dock: ${FAST_DOCK_SECS}s"
else
    rc=$?; FAST_DOCK_SECS="FAILED(rc=${rc})"; warn "fast dock did not complete (rc=${rc}); see ${FAST_OUT}/dock.log"
fi

read -r BASE_CF BASE_POSE < <(extract_best_cf "${BASE_OUT}")
read -r FAST_CF FAST_POSE < <(extract_best_cf "${FAST_OUT}")

log ""
log "─── DOCK RESULTS ─────────────────────────────────────────────"
log "baseline dock wall-time : ${BASE_DOCK_SECS}s"
log "fast     dock wall-time : ${FAST_DOCK_SECS}s"
if [[ "${BASE_DOCK_SECS}" =~ ^[0-9.]+$ && "${FAST_DOCK_SECS}" =~ ^[0-9.]+$ ]]; then
    SPEEDUP="$(awk -v a="${BASE_DOCK_SECS}" -v b="${FAST_DOCK_SECS}" 'BEGIN{ if(b>0) printf "%.3f", a/b; else printf "NA"}')"
    log "dock speedup (baseline/fast) : ${SPEEDUP}x"
fi

log ""
log "─── PARITY (rank/pose-equivalence, drift allowed) ────────────"
log "baseline best-CF : ${BASE_CF}   top-pose : ${BASE_POSE}"
log "fast     best-CF : ${FAST_CF}   top-pose : ${FAST_POSE}"
if [[ "${BASE_CF}" == "NA" || "${FAST_CF}" == "NA" ]]; then
    warn "Could not extract best-CF from one or both runs — parity UNDETERMINED (no fabricated value)."
    log "parity : UNDETERMINED (missing REMARK CF= in output)"
else
    DCF="$(awk -v a="${BASE_CF}" -v b="${FAST_CF}" 'BEGIN{d=a-b; if(d<0)d=-d; printf "%.5f", d}')"
    if awk -v d="${DCF}" -v t="${CF_TOL}" 'BEGIN{exit !(d<=t)}'; then
        ok "parity PASS  |ΔCF|=${DCF} <= tol ${CF_TOL}"
        log "parity : PASS  |ΔCF|=${DCF}  (tol=${CF_TOL})"
    else
        warn "parity FAIL |ΔCF|=${DCF} > tol ${CF_TOL}"
        log "parity : FAIL  |ΔCF|=${DCF}  (tol=${CF_TOL})"
    fi
fi

log ""
log "Full report: ${REPORT}"
if [[ "${KEEP}" -eq 0 ]]; then
    info "(build dirs kept under ${OUT_DIR}; pass nothing to reuse, or delete manually to reclaim disk)"
fi
ok "bench_flexaid_fast.sh complete"
