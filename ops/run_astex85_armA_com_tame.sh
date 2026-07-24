#!/usr/bin/env bash
# =============================================================================
# run_astex85_armA_com_tame.sh — Arm A canary WITH com-taming fixes ON
#
# Diff vs a bare Arm A min-CF run:
#   + FLEXAIDDS_VCT_NORM=1
#   + FLEXAIDDS_COM_FLOOR=500  (softplus floor; engine transition S=5)
#
# Use this for a single-arm hypothesis test before changing campaign defaults.
# Full multi-arm campaigns should only adopt these if this canary recovers
# seeded targets and stops CF.com ~ -3000 blow-ups.
# =============================================================================
set -uo pipefail

REPO="${FLEXAIDDS_ROOT:-${FLEXAIDDS_REPO:-$(cd "$(dirname "$0")/.." && pwd)}}"
ENGINE="${FLEXAIDDS_BINARY:-${REPO}/build/FlexAIDdS}"
RUNNER="${FLEXAIDDS_RUNNER:-${REPO}/build/benchmark_datasets}"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_benchmark_results/astex85_armA_com_tame_${STAMP}"
CACHE="${HOME}/.flexaidds/benchmarks"
LOG="${ROOT}/campaign.log"

trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== Astex-85 Arm A com-tame canary ${STAMP} ==="
echo "engine : ${ENGINE}"
echo "sha256 : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit : $(cd "${REPO}" && git rev-parse HEAD)"
echo "root   : ${ROOT}"
echo "env    : FLEXAIDDS_VCT_NORM=1 FLEXAIDDS_COM_FLOOR=500 (Arm A min-CF)"

if pgrep -f "benchmark_datasets --benchmark" >/dev/null 2>&1; then
    echo "[ABORT] a benchmark_datasets run is already active — refusing to share the cache"
    exit 1
fi

export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_ORACLE_SITE_DIR="${REPO}/benchmarks/astex_diverse/astex_diverse"

# Dormant fixes — the whole point of this canary
export FLEXAIDDS_VCT_NORM=1
export FLEXAIDDS_COM_FLOOR=500

# Pure Arm A ranking
unset FLEXAIDDS_THERMO_SCORE FLEXAIDDS_T_EFF FLEXAIDDS_SOFTBETA_ELECTION FLEXAIDDS_SMART_WATER

echo "active FLEXAIDDS_* :"
env | grep -E '^FLEXAIDDS_' | sort || true

mkdir -p "${ROOT}/armA_com_tame"
caffeinate -i "${RUNNER}" \
    --benchmark astex \
    --output "${ROOT}/armA_com_tame" \
    --cache  "${CACHE}" \
    --threads 6 \
    --omp-threads 1 \
    --job-timeout-seconds 3600 \
    --force
echo "=== finished rc=$? $(date -u +%FT%TZ) ==="
touch "${ROOT}/.CAMPAIGN_DONE"
