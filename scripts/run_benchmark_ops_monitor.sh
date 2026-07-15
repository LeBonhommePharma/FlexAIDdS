#!/usr/bin/env bash
# run_benchmark_ops_monitor.sh — one-shot automated ops + finished-run monitor.
# Intended for scheduler / launchd / cron every 15–30 min.
#
# Production anti-hang:
#   - Prefer local campaign trees (~/flexaidds_results)
#   - Never block the scheduler on CloudDocs walks
#   - Optional reap of hung iCloud find/md5 walkers (NOT dockers)
#
#   bash scripts/run_benchmark_ops_monitor.sh
#   bash scripts/run_benchmark_ops_monitor.sh --reap-walkers
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export FLEXAIDDS_ROOT="${FLEXAIDDS_ROOT:-$ROOT}"
cd "$ROOT"

# shellcheck disable=SC1091
[[ -f "${HOME}/.flexaidds_env" ]] && source "${HOME}/.flexaidds_env"
# Local-first paths (do not force iCloud RESULTS for live scan)
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"
# Keep ICLOUD for mirror metadata; live scan uses LOCAL campaigns first
export FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$FLEXAIDDS_ICLOUD/results}"
export FLEXAIDDS_QUEUE_ROOT="${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}"

REAP=0
for a in "$@"; do
  [[ "$a" == "--reap-walkers" ]] && REAP=1
done

if (( REAP )); then
  bash "$ROOT/scripts/reap_hung_icloud_walkers.sh" || true
fi

export FLEXAIDDS_MONITOR_SCRATCH="${FLEXAIDDS_MONITOR_SCRATCH:-$FLEXAIDDS_LOCAL_ROOT/logs/ops_monitor}"
mkdir -p "$FLEXAIDDS_MONITOR_SCRATCH"

# Hard wall-clock for the whole monitor (CloudDocs must never stick cron)
TIMEOUT_SEC="${FLEXAIDDS_MONITOR_TIMEOUT_SEC:-90}"
if command -v timeout >/dev/null 2>&1; then
  exec timeout "$TIMEOUT_SEC" python3 "$ROOT/scripts/benchmark_ops_monitor.py" \
    --scratch "$FLEXAIDDS_MONITOR_SCRATCH" \
    ${@//--reap-walkers/}
else
  # macOS without GNU timeout: background + kill
  python3 "$ROOT/scripts/benchmark_ops_monitor.py" \
    --scratch "$FLEXAIDDS_MONITOR_SCRATCH" \
    ${@//--reap-walkers/} &
  mon_pid=$!
  (
    sleep "$TIMEOUT_SEC"
    kill -TERM "$mon_pid" 2>/dev/null || true
  ) &
  wait "$mon_pid" 2>/dev/null || true
fi
