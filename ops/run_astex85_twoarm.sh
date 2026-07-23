#!/usr/bin/env bash
# =============================================================================
# run_astex85_twoarm.sh — Astex-85 three-arm campaign
#
#   Arm A (baseline)   : min(CF) election, no thermodynamic ranking
#   Arm B (dG_eff)     : FLEXAIDDS_THERMO_SCORE=1, T_EFF=0.596,
#                        FLEXAIDDS_SOFTBETA_ELECTION=1
#   Arm C (smartwater) : Arm B thermodynamics + selective crystallographic
#                        water retention (FLEXAIDDS_SMART_WATER=1 →
#                        binding_site_water_radius 4.5 Å, H-bond required).
#                        Isolates the effect of stripping bulk/surface waters
#                        while keeping bridging waters as receptor atoms.
#
# Arms run SEQUENTIALLY on purpose: benchmark_datasets instances share one cache
# directory and concurrent runs corrupt it.
#
# Reproducibility: OMP_NUM_THREADS=1 and serial restarts
# (FLEXAIDDS_PARALLEL_RESTARTS=0) make each worker bit-deterministic; the 6
# workers are independent processes on different targets, so cross-target
# parallelism does not affect any single target's result.
# =============================================================================
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_benchmark_results/astex85_twoarm_${STAMP}"
CACHE="${HOME}/.flexaidds/benchmarks"
LOG="${ROOT}/campaign.log"

# Survive session teardown (v7 multi-worker SIGTERM kill). setsid puts us in our
# own process group; only HUP is ignored so the campaign stays deliberately killable.
trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== Astex-85 three-arm campaign ${STAMP} ==="
echo "engine : ${ENGINE}"
echo "sha256 : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit : $(cd "${REPO}" && git rev-parse HEAD)"
echo "root   : ${ROOT}"

# Single-instance guard.
if pgrep -f "benchmark_datasets --benchmark" >/dev/null 2>&1; then
    echo "[ABORT] a benchmark_datasets run is already active — refusing to share the cache"
    exit 1
fi

# ── Reproducibility / protocol env common to both arms ────────────────────────
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0     # serial restarts → deterministic
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_ORACLE_SITE_DIR="${REPO}/benchmarks/astex_diverse/astex_diverse"

run_arm() {
    local name="$1"; shift
    local out="${ROOT}/${name}"
    mkdir -p "${out}"
    echo ""
    echo "=== ARM ${name} starting $(date -u +%FT%TZ) ==="
    env | grep -E '^FLEXAIDDS_(THERMO_SCORE|T_EFF|SOFTBETA_ELECTION|SMART_WATER)=' | sort || true
    caffeinate -i "${RUNNER}" \
        --benchmark astex \
        --output "${out}" \
        --cache  "${CACHE}" \
        --threads 6 \
        --omp-threads 1 \
        --job-timeout-seconds 3600 \
        --force
    echo "=== ARM ${name} finished rc=$? $(date -u +%FT%TZ) ==="
}

# ── Arm A: baseline min(CF) ───────────────────────────────────────────────────
(
  unset FLEXAIDDS_THERMO_SCORE FLEXAIDDS_T_EFF FLEXAIDDS_SOFTBETA_ELECTION FLEXAIDDS_SMART_WATER
  run_arm "armA_mincf"
)

# ── Arm B: dG_eff, T=0.596, Softbeta election ─────────────────────────────────
(
  unset FLEXAIDDS_SMART_WATER
  export FLEXAIDDS_THERMO_SCORE=1
  export FLEXAIDDS_T_EFF=0.596
  export FLEXAIDDS_SOFTBETA_ELECTION=1
  run_arm "armB_dgeff"
)

# ── Arm C: smart water retention, same thermodynamics as Arm B ────────────────
(
  export FLEXAIDDS_THERMO_SCORE=1
  export FLEXAIDDS_T_EFF=0.596
  export FLEXAIDDS_SOFTBETA_ELECTION=1
  # binding_site_water_radius: 4.5 will be in dock_config default
  export FLEXAIDDS_SMART_WATER=1
  run_arm "armC_smartwater"
)

echo ""
echo "=== CAMPAIGN COMPLETE $(date -u +%FT%TZ) ==="
touch "${ROOT}/.CAMPAIGN_DONE"
