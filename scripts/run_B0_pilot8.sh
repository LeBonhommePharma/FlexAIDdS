#!/usr/bin/env bash
# Arm B0 (FlexAID master, TEMPER 0 / CF) pilot8 → iCloud.
# Deferred relative to A / B@T21 / C0; still must land on iCloud when run.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" B0 "$@"
