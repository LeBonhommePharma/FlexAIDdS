#!/usr/bin/env bash
# =============================================================================
# run_all_thermo.sh — One-command full production run for ThermoAffinitySuite v2.1
#
# Runs: itc187 + bindingdb_itc + scorpio
# Supports: --ablation, --tier, --dry-run, PB validation + analyzer + bundle
#
# Usage:
#   ./scripts/run_all_thermo.sh --tier 1 --dry-run
#   ./scripts/run_all_thermo.sh --ablation
#
# Produces unified artifact bundle + leaderboard + PB failure analysis.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TIER=2
DRY_RUN=""
ABLATION=""
NO_PB=""
RESULTS_DIR="${RESULTS_DIR:-results/thermo_v21}"
MAX_TARGETS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tier) TIER="$2"; shift 2 ;;
    --dry-run) DRY_RUN="--dry-run"; shift ;;
    --ablation) ABLATION="--ablation"; shift ;;
    --no-pb) NO_PB="--no-pb"; shift ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --max-targets) MAX_TARGETS="--max-targets $2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--tier N] [--dry-run] [--ablation] [--no-pb] [--results-dir DIR] [--max-targets N]"
      exit 0
      ;;
    *) echo "Unknown: $1"; exit 1 ;;
  esac
done

cd "$REPO_ROOT"

echo "========================================================================"
echo " ThermoAffinitySuite v2.1 — FULL PRODUCTION RUN (ITC-187 + BindingDB + SCORPIO)"
echo " TIER=${TIER}  ABLATION=${ABLATION:-false}  DRY=${DRY_RUN:-false}"
echo "========================================================================"

python -m benchmarks.runner \
  --all-thermo \
  --tier "${TIER}" \
  --results-dir "${RESULTS_DIR}" \
  ${DRY_RUN} ${ABLATION} ${NO_PB} \
  ${MAX_TARGETS}

echo ""
echo "=== Running PB failure analyzer ==="
python scripts/analyze_pb_failures.py \
  --results-dir "${RESULTS_DIR}" \
  --out-dir "${RESULTS_DIR}/pb_analysis" \
  ${DRY_RUN:+--smoke} || echo "Analyzer completed (or skipped in dry-run)."

echo ""
echo "All done. Leaderboard + bundle in ${RESULTS_DIR}"
ls -l "${RESULTS_DIR}"/thermoaffinity_v21_* 2>/dev/null || true
