#!/usr/bin/env bash
# launch_comcap_fixed.sh — v_comcap_fixed full-85 Astex campaign
# FIXES vs v_comcap_softbeta:
#   1. ELECTION_ENTROPY=0 (soft-β OFF, reverts to CF rank-0)
#   2. Vcontacts.cpp clash cap restored on main (unbounded accumulation reverted)
# Still testing: COM_BURIAL_CAP=-130 + autonomous mode
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/v_comcap_fixed_${STAMP}"
LOG="${ROOT}/campaign.log"

trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== v_comcap_fixed ${STAMP} ==="
echo "engine   : ${ENGINE}"
echo "sha256   : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit   : $(cd "${REPO}" && git rev-parse HEAD)"
echo "runner   : ${RUNNER}"
echo "root     : ${ROOT}"
echo "env:"
echo "  FLEXAIDDS_COM_BURIAL_CAP=-130"
echo "  FLEXAIDDS_VCT_R0=4.0"
echo "  FLEXAIDDS_ELECTION_ENTROPY=0  [FIXED: was 1]"
echo "  FLEXAIDDS_SEED_ELITISM=0"
echo "  OMP_NUM_THREADS=1"
echo "  FLEXAIDDS_PARALLEL_RESTARTS=0"
echo "  mode: autonomous (seed_elitism=OFF)"
echo "  threads: 6 x omp-threads 1"

# ── Reproducibility
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"

# ── COM burial cap — hypothesis still under test
export FLEXAIDDS_COM_BURIAL_CAP=-130

# ── Voronoi contact R0
export FLEXAIDDS_VCT_R0=4.0

# ── Soft-β election OFF (fix #1: was causing election to pick largest cluster)
export FLEXAIDDS_ELECTION_ENTROPY=0

# ── No seed elitism
export FLEXAIDDS_SEED_ELITISM=0

echo ""
echo "active FLEXAIDDS_* :"
env | grep -E '^FLEXAIDDS_' | sort || true
echo ""

caffeinate -i "${RUNNER}" \
    --benchmark astex \
    --mode autonomous \
    --output "${ROOT}" \
    --threads 6 \
    --omp-threads 1 \
    --job-timeout-seconds 3600 \
    --force

RC=$?
echo "=== finished rc=${RC} $(date -u +%FT%TZ) ==="
[ "${RC}" -eq 0 ] && touch "${ROOT}/.CAMPAIGN_DONE"
exit "${RC}"
