#!/usr/bin/env bash
# run_benchmark_ops_monitor.sh — one-shot automated ops + results monitor.
# Intended for scheduler / launchd / cron every 15–30 min.
#
#   bash scripts/run_benchmark_ops_monitor.sh
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
cd "$ROOT"

# shellcheck disable=SC1091
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"
# shellcheck disable=SC1091
source "$ROOT/scripts/use_icloud_benchmark_storage.sh"

export FLEXAIDDS_MONITOR_SCRATCH="${FLEXAIDDS_MONITOR_SCRATCH:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-19923fdc9045/implementer}"
mkdir -p "$FLEXAIDDS_MONITOR_SCRATCH"

exec python3 "$ROOT/scripts/benchmark_ops_monitor.py" \
  --scratch "$FLEXAIDDS_MONITOR_SCRATCH" \
  "$@"
