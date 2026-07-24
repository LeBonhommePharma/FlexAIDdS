#!/usr/bin/env bash
# launch_comcap_softbeta.sh — v_comcap_softbeta full-85 Astex campaign
# COM_BURIAL_CAP=-130 + soft-β election + autonomous mode (seed_elitism=OFF)
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAID"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/v_comcap_softbeta_${STAMP}"
LOG="${ROOT}/campaign.log"

trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== v_comcap_softbeta ${STAMP} ==="
echo "engine   : ${ENGINE}"
echo "sha256   : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit   : $(cd "${REPO}" && git rev-parse HEAD)"
echo "runner   : ${RUNNER}"
echo "root     : ${ROOT}"
echo "env:"
echo "  FLEXAIDDS_COM_BURIAL_CAP=-130"
echo "  FLEXAIDDS_VCT_R0=4.0"
echo "  FLEXAIDDS_ELECTION_ENTROPY=1"
echo "  FLEXAIDDS_SEED_ELITISM=0"
echo "  OMP_NUM_THREADS=1"
echo "  FLEXAIDDS_PARALLEL_RESTARTS=0"
echo "  mode: autonomous (seed_elitism=OFF, blinding=ON)"
echo "  threads: 6 x omp-threads 1"

# ── Reproducibility
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"

# ── COM burial cap — main hypothesis for this campaign
export FLEXAIDDS_COM_BURIAL_CAP=-130

# ── Voronoi contact R0
export FLEXAIDDS_VCT_R0=4.0

# ── Soft-β entropy-aware election (default ON, explicit here)
export FLEXAIDDS_ELECTION_ENTROPY=1

# ── No seed elitism (blinding=ON, matches --mode autonomous)
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
