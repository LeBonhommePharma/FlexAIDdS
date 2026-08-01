#!/usr/bin/env bash
# launch_comfloor.sh — v_comfloor full-85 Astex campaign
# Uses FLEXAIDDS_COM_FLOOR=130 (softplus floor at -130, transition S=5 in binary).
# COM_BURIAL_CAP is a separate lever; this script only exercises COM_FLOOR.
set -uo pipefail

REPO="${FLEXAIDDS_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
ENGINE="${FLEXAIDDS_BINARY:-${REPO}/build/FlexAIDdS}"
RUNNER="${FLEXAIDDS_RUNNER:-${REPO}/build/benchmark_datasets}"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/v_comfloor_${STAMP}"
LOG="${ROOT}/campaign.log"

trap '' HUP
mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true
exec >>"${LOG}" 2>&1

echo "=== v_comfloor ${STAMP} ==="
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
export FLEXAIDDS_SEED_ELITISM=0

# Promoted reference config (COM floor + the P1 pocket/clash steer it enables,
# plus the determinism knobs).  Sourced rather than duplicated so this campaign
# cannot drift from the configuration that was actually measured.
#
# NOTE: this script previously set FLEXAIDDS_COM_FLOOR=130 alone.  That is only
# half the configuration -- the 2026-08-01 sweep shows the floor WITHOUT the P1
# weights, and P1 without the floor, both score docking_power_top1 = 0.000.
set -a
. "${REPO}/ops/reference_config.env"
set +a

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
