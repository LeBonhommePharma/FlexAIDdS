#!/usr/bin/env bash
# launch_autonomous_blind.sh — full-85 Astex campaign, autonomous blind.
#
# No oracle site, no native seeding, no seed elitism: the only inputs are the
# apo receptor and the ligand topology. Carries the 2026-07-24 fix set:
#   * dead VCT rows 16/20/21/27 aliased to live rows (sulfonamide S, ring O/S)
#   * retained water oxygen typed O.3 instead of C.1
#   * COM_FLOOR softplus transition scale decoupled from the floor depth
# Spread guard and the 1TW6 peptide fallback are already in the binary.
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_results/v_autonomous_${STAMP}"
LOG="${ROOT}/campaign.log"

trap '' HUP
mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true
exec >>"${LOG}" 2>&1

echo "=== v_autonomous ${STAMP} ==="
echo "engine   : ${ENGINE}"
echo "sha256   : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit   : $(cd "${REPO}" && git rev-parse HEAD)"
echo "branch   : $(cd "${REPO}" && git rev-parse --abbrev-ref HEAD)"
echo "runner   : ${RUNNER}"
echo "root     : ${ROOT}"

export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_VCT_R0=4.0
export FLEXAIDDS_ELECTION_ENTROPY=0   # soft-beta off: T=298 is ~500x wrong scale
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_COM_FLOOR=130        # softplus floor at -130, S=5 transition
# No FLEXAIDDS_ORACLE_SITE_DIR — site comes from blind cavity detection.

echo ""
echo "active FLEXAIDDS_* :"
env | grep -E '^FLEXAIDDS_|^OMP_NUM_THREADS|^FLEXAID_SEED' | sort || true
echo ""

caffeinate -i "${RUNNER}" \
    --benchmark astex \
    --mode autonomous \
    --output "${ROOT}" \
    --threads 4 \
    --omp-threads 1 \
    --job-timeout-seconds 3600 \
    --force

RC=$?
echo "=== finished rc=${RC} $(date -u +%FT%TZ) ==="
[ "${RC}" -eq 0 ] && touch "${ROOT}/.CAMPAIGN_DONE"
exit "${RC}"
