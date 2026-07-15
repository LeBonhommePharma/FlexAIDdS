#!/usr/bin/env bash
# Arm B (FlexAID master, TEMPER 298 / entropy ranking) pilot8 → iCloud.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" B "$@"
