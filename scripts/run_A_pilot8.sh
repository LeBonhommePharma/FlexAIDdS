#!/usr/bin/env bash
# Arm A (FlexAID JCIM 2015) pilot8 → iCloud only.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
exec bash "$ROOT/scripts/run_flexaid_arm_pilot8.sh" A "$@"
