#!/usr/bin/env bash
# launch_seeded.sh — oracle-ceiling seeded comparison run
# Binary ba42e2bd (Vcontacts cap + COM_FLOOR S=5 fix)
# Matches v50b conditions for apples-to-apples comparison vs 76.5%
# Reduced to 4 workers out of consideration for shared resources
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/v_seeded_${STAMP}"
LOG="${ROOT}/campaign.log"

trap '' HUP
mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true
exec >>"${LOG}" 2>&1

echo "=== v_seeded ${STAMP} (oracle-ceiling, ba42e2bd) ==="
echo "engine   : ${ENGINE}"
echo "sha256   : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit   : $(cd "${REPO}" && git rev-parse HEAD)"
echo "runner   : ${RUNNER}"
echo "root     : ${ROOT}"

export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_VCT_R0=4.0
export FLEXAIDDS_ELECTION_ENTROPY=0
# NO COM_FLOOR — clean comparison against v50b

echo ""
echo "active FLEXAIDDS_* :"
env | grep -E '^FLEXAIDDS_' | sort || true
echo ""

caffeinate -i "${RUNNER}" \
    --benchmark astex \
    --mode oracle-ceiling \
    --output "${ROOT}" \
    --threads 4 \
    --omp-threads 1 \
    --job-timeout-seconds 3600 \
    --force

RC=$?
echo "=== finished rc=${RC} $(date -u +%FT%TZ) ==="
[ "${RC}" -eq 0 ] && touch "${ROOT}/.CAMPAIGN_DONE"
exit "${RC}"
