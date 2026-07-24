#!/usr/bin/env bash
# launch_coarse6.sh — 6-target coarse_init probe with per-target grid_step for 1K3U.
#
# Targets: 1G9V 1J3J 1K3U 1L7F 1M2Z 1N1M (Astex Diverse subset)
#
# Protocol (genuine, no-seeding):
#   mode UNSET  -> seed_elitism forced OFF by DatasetRunner (seed_echo=0 metric)
#   coarse_init is ALWAYS ON in this build (not gated on oracle-ceiling), so the
#   grid_step patch is exercised regardless of mode.
#   FLEXAIDDS_RESTARTS=8, seed_fraction=0.0, seed_elitism=false,
#   sas_weight=0.40, hbond_rank=on.
#   1K3U ONLY: FLEXAIDDS_COARSE_GRID_STEP=1.5 (finer scan for narrow/deep pocket).
#
# Because the grid-step env var is process-global, 1K3U runs as a SEPARATE
# invocation (subdir k3u/) from the other 5 (subdir main/) so the finer grid
# does not leak onto other targets and result.csv files do not clobber.
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RUNNER="${ROOT}/build/benchmark_datasets"
BINARY="${ROOT}/build/FlexAIDdS"

OUTDIR="${COARSE6_OUTDIR:?set COARSE6_OUTDIR}"
MAIN_OUT="${OUTDIR}/main"
K3U_OUT="${OUTDIR}/k3u"
mkdir -p "${MAIN_OUT}" "${K3U_OUT}"

[[ -x "${RUNNER}" ]] || { echo "FAIL: runner not executable: ${RUNNER}" >&2; exit 1; }
[[ -x "${BINARY}" ]] || { echo "FAIL: engine not executable: ${BINARY}" >&2; exit 1; }

# Shared genuine-protocol env (no seeding).
export FLEXAIDDS_BINARY="${BINARY}"
export FLEXAIDDS_RESTARTS=8
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_NATIVE_SEED_FRAC=0
export FLEXAIDDS_SAS_WEIGHT=0.40
export FLEXAIDDS_HBOND_RANK=1
export FLEXAIDDS_VCT_R0=4.0
export OMP_WAIT_POLICY=passive

POP=1000
GEN=2000
TEMP=298
TIMEOUT=10800

echo "=== COARSE6 LAUNCH $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "RUNNER=${RUNNER}"
echo "BINARY=${BINARY}"
echo "OUTDIR=${OUTDIR}"
echo "mode=UNSET (seed_elitism OFF) RESTARTS=8 pop=${POP} gen=${GEN} T=${TEMP}"
echo "sas_weight=0.40 hbond_rank=on native_seed_frac=0 vct_r0=4.0"
echo "runner_sha256=$(shasum -a 256 "${RUNNER}" | awk '{print $1}')"
echo "engine_sha256=$(shasum -a 256 "${BINARY}" | awk '{print $1}')"

# ---- Phase 1: the 4 standard-grid targets (1M2Z excluded — sentinel, blind placement fails) ---
echo "--- Phase 1: 1G9V,1J3J,1L7F,1N1M (grid_step default 3.0, vct_r0=4.0) ---"
"${RUNNER}" \
  --benchmark astex \
  --only-codes 1G9V,1J3J,1L7F,1N1M \
  --output "${MAIN_OUT}/" \
  --threads 4 \
  --omp-threads 2 \
  --ga-population "${POP}" \
  --ga-generations "${GEN}" \
  --temperature "${TEMP}" \
  --job-timeout-seconds "${TIMEOUT}"

# ---- Phase 2: 1K3U with finer coarse grid --------------------------------
echo "--- Phase 2: 1K3U (FLEXAIDDS_COARSE_GRID_STEP=1.5) ---"
FLEXAIDDS_COARSE_GRID_STEP=1.5 \
"${RUNNER}" \
  --benchmark astex \
  --only-codes 1K3U \
  --output "${K3U_OUT}/" \
  --threads 1 \
  --omp-threads 8 \
  --ga-population "${POP}" \
  --ga-generations "${GEN}" \
  --temperature "${TEMP}" \
  --job-timeout-seconds "${TIMEOUT}"

echo "=== COARSE6 DONE $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
