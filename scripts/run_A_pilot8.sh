#!/usr/bin/env bash
# Arm A (FlexAID 2015-era) pilot8 → iCloud.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" A "$@"
