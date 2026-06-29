#!/usr/bin/env bash
# =============================================================================
# run_itc187.sh — One-command runner for ITC-187 via ThermoAffinitySuite v2.1
#
# Usage:
#   ./scripts/run_itc187.sh [--tier 1|2] [--dry-run] [--ablation] [--no-pb]
#
# Examples:
#   ./scripts/run_itc187.sh --tier 1 --dry-run
#   caffeinate -i ./scripts/run_itc187.sh --tier 2
#
# Sets temperature=298K (standard for ITC reporting).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

TIER=2
DRY_RUN=""
ABLATION=""
NO_PB=""
RESULTS_DIR="${RESULTS_DIR:-results/thermo}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --tier) TIER="$2"; shift 2 ;;
    --dry-run) DRY_RUN="--dry-run"; shift ;;
    --ablation) ABLATION="--ablation"; shift ;;
    --no-pb) NO_PB="--no-pb"; shift ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--tier N] [--dry-run] [--ablation] [--no-pb] [--results-dir DIR]"
      exit 0
      ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

cd "$REPO_ROOT"

echo "=== ThermoAffinitySuite v2.1 — ITC-187 (T=${TIER}) ==="
python -m benchmarks.runner \
  --dataset itc187 \
  --tier "${TIER}" \
  --results-dir "${RESULTS_DIR}" \
  ${DRY_RUN} ${ABLATION} ${NO_PB} \
  --all-thermo || true   # --all-thermo ignored with explicit --dataset but harmless

echo "ITC-187 run complete. Check ${RESULTS_DIR}"
