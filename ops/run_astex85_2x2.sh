#!/usr/bin/env bash
# =============================================================================
# OWNERSHIP — this campaign is run by Claude Science, not by hand.
#
# Do not launch this script manually. Claude Science owns the run: it starts the
# campaign, watches it to completion, and is responsible for the results.
#
# Reporting is also Science's job, and goes through the flexaidds Python metrics
# pipeline against the DatasetRunner output under each cell directory — not
# through hand-read logs. The success rate this script prints in its own summary
# table is an at-a-glance sanity check during the run; the reported numbers are
# the ones the Python pipeline produces.
# =============================================================================
# run_astex85_2x2.sh — Astex-85 2x2 factorial campaign
#
# Replaces the three-arm design in run_astex85_twoarm.sh, which confounded its
# factors: arms A and B both retained waters while only B carried COM_FLOOR, and
# B additionally carried VCT_NORM and the (inert) thermo env, so no arm pair
# isolated a single variable.
#
# The two independent variables under test:
#
#   Factor 1  water retention : FLEXAIDDS_SMART_WATER=1  (retain bridging waters)
#                               vs FLEXAIDDS_STRIP_ALL_WATERS=1 (strip all)
#   Factor 2  com lower clamp : FLEXAIDDS_COM_FLOOR=500 vs unset (no floor)
#
#   cell          waters      COM_FLOOR
#   ------------  ----------  ---------
#   baseline      stripped    unset      current v50b behaviour
#   floor_only    stripped    500
#   water_only    retained    unset
#   full          retained    500
#
# Main effects and the interaction are then read directly off the four cells:
#   water effect = (water_only - baseline) and (full - floor_only)
#   floor effect = (floor_only - baseline) and (full - water_only)
#   interaction  = (full - floor_only) - (water_only - baseline)
#
# VCT_NORM is deliberately NOT a factor and stays unset in every cell: it is a
# third variable that would double the design, and the canary run showed
# VCT_NORM=1 gives 0/5 on 5 Astex targets (CF.com -77K to -369K).
#
# The thermo env (THERMO_SCORE / T_EFF / SOFTBETA_ELECTION) is also omitted: the
# thermodynamic gate is printf-only (gaboom.cpp:1341) and runs after the QuickSort
# at gaboom.cpp:705 that finalizes the elected pose, so it cannot move results.
#
# Cells run SEQUENTIALLY on purpose: benchmark_datasets instances share one cache
# directory and concurrent runs corrupt it.
#
# Reproducibility: OMP_NUM_THREADS=1 and serial restarts
# (FLEXAIDDS_PARALLEL_RESTARTS=0) make each worker bit-deterministic; the workers
# are independent processes on different targets, so cross-target parallelism does
# not affect any single target's result.
# =============================================================================
set -euo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
CACHE="${HOME}/.flexaidds/benchmarks"
ORACLE_DIR="${REPO}/benchmarks/astex_diverse/astex_diverse"
# WORKERS: 18 GB box. 25 concurrent FlexAIDdS (~1.5 GB each) = ~37 GB >> 18 GB → thrash/OOM.
# Capped to 6 (~9 GB peak) to fit physical RAM with headroom on a shared machine. — OPS 2026-07-23
WORKERS="${FLEXAIDDS_CAMPAIGN_WORKERS:-6}"
JOB_TIMEOUT=3600
RMSD_CUTOFF=2.0

# FIX (OPS 2026-07-23): caffeinate aborts in restricted/sandboxed contexts
# ("Failed to create PreventUserIdleSystemSleep assertion" → rc=134), killing the
# cell. Use it only if present AND able to create an assertion; otherwise run bare.
CAFF=()
if command -v caffeinate >/dev/null 2>&1 && caffeinate -i true >/dev/null 2>&1; then
    CAFF=(caffeinate -i)
fi
# FIX (OPS 2026-07-23): allow resume. Default keeps --force (clean run, as intended);
# FLEXAIDDS_CAMPAIGN_RESUME=1 drops it so an interrupted campaign skips finished targets.
FORCE_FLAG="--force"
[ "${FLEXAIDDS_CAMPAIGN_RESUME:-0}" = "1" ] && FORCE_FLAG=""

STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/campaign_2x2_${STAMP}"
LOG="${ROOT}/campaign.log"
CELLS=(baseline floor_only water_only full)

for f in "${ENGINE}" "${RUNNER}"; do
    [ -x "${f}" ] || { echo "[ABORT] missing or non-executable: ${f}" >&2; exit 1; }
done
[ -d "${ORACLE_DIR}" ] || { echo "[ABORT] oracle site dir not found: ${ORACLE_DIR}" >&2; exit 1; }

# ── Single-instance guard (PID file; flock is not available on macOS) ─────────
# Stale locks self-heal: if the recorded PID is gone, we take the lock over.
LOCKFILE="/tmp/flexaidds_campaign_2x2.pid"
if [ -f "${LOCKFILE}" ]; then
    other="$(cat "${LOCKFILE}" 2>/dev/null || true)"
    if [ -n "${other}" ] && kill -0 "${other}" 2>/dev/null; then
        echo "[ABORT] campaign already running as PID ${other} (${LOCKFILE})" >&2
        exit 1
    fi
    echo "[INFO] clearing stale lock from PID ${other:-unknown}"
fi
echo "$$" > "${LOCKFILE}"
trap 'rm -f "${LOCKFILE}"' EXIT

# Survive session teardown (v7 multi-worker SIGTERM kill): only HUP is ignored,
# so the campaign stays deliberately killable.
trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== Astex-85 2x2 factorial campaign ${STAMP} ==="
echo "engine : ${ENGINE}"
echo "sha256 : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit : $(cd "${REPO}" && git rev-parse HEAD)"
echo "root   : ${ROOT}"
echo "workers: ${WORKERS}"

# ── Reproducibility / protocol env common to every cell ───────────────────────
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0     # serial restarts → deterministic
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_ORACLE_SITE_DIR="${ORACLE_DIR}"

# ── Metrics ───────────────────────────────────────────────────────────────────
# The runner's own `success` column means "docking ran", NOT "RMSD < cutoff", so
# the success rate is recomputed here from rmsd_to_crystal. Prints "<n> <total>".
cell_metrics() {
    local out="$1"
    local csv
    csv="$(find "${out}" -maxdepth 2 -name 'astex_diverse_results.csv' -print -quit)"
    if [ -z "${csv}" ]; then
        echo "NA NA"
        return 0
    fi
    RMSD_CUTOFF="${RMSD_CUTOFF}" python3 - "${csv}" <<'PY'
import csv, os, sys
cut = float(os.environ["RMSD_CUTOFF"])
hits = total = 0
with open(sys.argv[1], newline="") as fh:
    for row in csv.DictReader(fh):
        total += 1
        try:
            if float(row.get("rmsd_to_crystal", "")) < cut:
                hits += 1
        except (TypeError, ValueError):
            pass
print(hits, total)
PY
}

# ── Mechanism logger ──────────────────────────────────────────────────────────
# Walk the cell's elected poses and pull the CF.com channel plus the single
# largest per-pair contact contribution out of the pose REMARKs, so the campaign
# reads as a direct mechanism test ("did com blow up, and on which contact
# type?") rather than an RMSD-only inference. Healthy com ~ -130; blow-up >~ -1000.
com_summary() {
    local out="$1"
    local csv="${out}/com_mechanism.csv"
    echo "pdb_id,cf_com,cf_sas,cf_wal,cf_total,dominant_contact_type,dominant_contact_value" > "${csv}"
    local pose
    for pose in "${out}"/*/elected_pose.pdb; do
        [ -f "${pose}" ] || continue
        local pdb com sas wal tot dom dtype dval
        pdb="$(basename "$(dirname "${pose}")")"
        com="$(awk '/^REMARK CF\.com=/{split($0,a,"="); print a[2]; exit}' "${pose}")"
        sas="$(awk '/^REMARK CF\.sas=/{split($0,a,"="); print a[2]; exit}' "${pose}")"
        wal="$(awk '/^REMARK CF\.wal=/{split($0,a,"="); print a[2]; exit}' "${pose}")"
        tot="$(awk -F= '/^REMARK CF=/{print $2; exit}' "${pose}")"
        # FIX (OPS 2026-07-23): the previous pattern /^REMARK contact / never matched.
        # Real REMARK format is:  REMARK   1:  <A>-<B> with energy of <VALUE>
        # under the "Most important POSITIVE" header. Take the first (rank-1 = largest
        # magnitude) contact: dtype=the "A-B" atom-type-index pair, dval=its energy.
        # FIX (OPS 2026-07-23): the previous pattern /^REMARK contact / never matched.
        # Real format under the "Most important POSITIVE" header:
        #   REMARK   1:  <A>-<B> with energy of <VALUE>
        # The A-B pair may be spaced ("1- 4") or not ("1-14"); normalise by stripping
        # the "REMARK  N:" prefix and all spaces from the "A-B" token. Take rank-1
        # (largest-magnitude contact). dtype = "A-B" indices, dval = energy.
        dom="$(awk '/Most important POSITIVE/{p=1; next}
                    p && /with energy of/{
                        e=$NF;
                        s=$0; sub(/^REMARK[[:space:]]+[0-9]+:[[:space:]]*/,"",s);
                        sub(/[[:space:]]+with energy of.*/,"",s); gsub(/[[:space:]]/,"",s);
                        print s, e; exit
                    }' "${pose}")"
        dtype="$(awk '{print $1}' <<< "${dom}")"
        dval="$(awk '{print $2}' <<< "${dom}")"
        echo "${pdb},${com:-NA},${sas:-NA},${wal:-NA},${tot:-NA},${dtype:-NA},${dval:-NA}" >> "${csv}"
    done
    echo "  [MECHANISM] wrote $(( $(wc -l < "${csv}") - 1 )) rows → ${csv}"
}

# run_cell <name> — env for the cell is set by the caller.
# Returns the runner's exit status; it is NOT swallowed.
run_cell() {
    local name="$1"
    local out="${ROOT}/${name}"
    mkdir -p "${out}"
    echo ""
    echo "=== CELL ${name} starting $(date -u +%FT%TZ) ==="
    env | grep -E '^FLEXAIDDS_(SMART_WATER|STRIP_ALL_WATERS|COM_FLOOR|VCT_NORM)=' | sort || true

    local rc=0
    "${CAFF[@]}" "${RUNNER}" \
        --benchmark astex \
        --output "${out}" \
        --cache  "${CACHE}" \
        --threads "${WORKERS}" \
        --omp-threads 1 \
        --job-timeout-seconds "${JOB_TIMEOUT}" \
        ${FORCE_FLAG} \
        >"${ROOT}/cell_${name}.log" 2>"${ROOT}/cell_${name}.err" || rc=$?

    echo "=== CELL ${name} finished rc=${rc} $(date -u +%FT%TZ) ==="
    if [ "${rc}" -ne 0 ]; then
        echo "  [ERROR] see ${ROOT}/cell_${name}.err"
        tail -20 "${ROOT}/cell_${name}.err" || true
        return "${rc}"
    fi

    com_summary "${out}"
    cell_metrics "${out}" > "${ROOT}/cell_${name}.metrics"
    echo "  [METRICS] $(cat "${ROOT}/cell_${name}.metrics") (hits total, RMSD < ${RMSD_CUTOFF} A)"
    return 0
}

clear_factor_env() {
    unset FLEXAIDDS_SMART_WATER FLEXAIDDS_STRIP_ALL_WATERS FLEXAIDDS_COM_FLOOR || true
    unset FLEXAIDDS_VCT_NORM FLEXAIDDS_THERMO_SCORE FLEXAIDDS_T_EFF FLEXAIDDS_SOFTBETA_ELECTION || true
}

# ── Cell 1/4 baseline: waters stripped, no COM_FLOOR (current v50b) ────────────
clear_factor_env
export FLEXAIDDS_STRIP_ALL_WATERS=1
run_cell baseline

# ── Cell 2/4 floor_only: waters stripped, COM_FLOOR=500 ───────────────────────
clear_factor_env
export FLEXAIDDS_STRIP_ALL_WATERS=1
export FLEXAIDDS_COM_FLOOR=500
run_cell floor_only

# ── Cell 3/4 water_only: waters retained, no COM_FLOOR ────────────────────────
clear_factor_env
export FLEXAIDDS_SMART_WATER=1
run_cell water_only

# ── Cell 4/4 full: waters retained, COM_FLOOR=500 ─────────────────────────────
clear_factor_env
export FLEXAIDDS_SMART_WATER=1
export FLEXAIDDS_COM_FLOOR=500
run_cell full

clear_factor_env

# ── Summary table ─────────────────────────────────────────────────────────────
echo ""
echo "=== SUMMARY (elected pose RMSD < ${RMSD_CUTOFF} A) ==="
printf '%-12s | %-8s | %-9s | %-7s\n' "cell" "waters" "com_floor" "N/85"
printf '%-12s-+-%-8s-+-%-9s-+-%-7s\n' "------------" "--------" "---------" "-------"
for cell in "${CELLS[@]}"; do
    case "${cell}" in
        baseline)   waters="strip"; floor="off" ;;
        floor_only) waters="strip"; floor="500" ;;
        water_only) waters="keep";  floor="off" ;;
        full)       waters="keep";  floor="500" ;;
    esac
    hits="NA"; total="NA"; pct="NA"
    if [ -s "${ROOT}/cell_${cell}.metrics" ]; then
        read -r hits total < "${ROOT}/cell_${cell}.metrics"
        if [ "${total}" != "NA" ] && [ "${total}" -gt 0 ] 2>/dev/null; then
            pct="$(awk -v h="${hits}" -v t="${total}" 'BEGIN{printf "%.1f", 100*h/t}')"
        fi
    fi
    printf '%-12s | %-8s | %-9s | %s/%s (%s%%)\n' "${cell}" "${waters}" "${floor}" "${hits}" "${total}" "${pct}"
done

echo ""
echo "=== CAMPAIGN COMPLETE $(date -u +%FT%TZ) ==="
touch "${ROOT}/.CAMPAIGN_DONE"
