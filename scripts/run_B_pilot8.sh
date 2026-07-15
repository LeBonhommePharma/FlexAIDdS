#!/usr/bin/env bash
# Arm B (FlexAID master, TEMPER 21 / LP-optimized entropy ranking) pilot8 → iCloud.
# Results always under $FLEXAIDDS_RESULTS (iCloud). Override TEMPER with FLEXAID_TEMPER.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
exec bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" B "$@"
