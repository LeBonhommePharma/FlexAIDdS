#!/usr/bin/env bash
# conditioned_full_launch.sh — ps-check then launch FULL astex_nonnative + posex benchmarks.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SCRATCH="${SCRATCH:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-1e3f1c67d740/implementer}"
STATUS="$SCRATCH/run_status.json"
LOG="$SCRATCH/conditioned_full_launch.log"
PYTHON_PKG="$REPO_ROOT/python"

mkdir -p "$SCRATCH"
exec > >(tee -a "$LOG") 2>&1

echo "=== conditioned_full_launch $(date -Iseconds) ==="

SIBLING_COUNT=$(pgrep -fl "astex_diverse|benchmark_datasets.*astex" 2>/dev/null | wc -l | tr -d ' ')
echo "Sibling astex_diverse/benchmark processes: $SIBLING_COUNT"

PYTHONPATH="$PYTHON_PKG${PYTHONPATH:+:$PYTHONPATH}" python3 -m flexaidds.dataset_runner.launch_queue \
  --scratch "$SCRATCH" \
  --repo-root "$REPO_ROOT" \
  --sibling-count "$SIBLING_COUNT" \
  --write-status "$STATUS"

if [[ "$SIBLING_COUNT" -gt 0 ]]; then
  echo "Siblings active — deferring to queue_large_benchmarks.sh (detached)"
  nohup env SCRATCH="$SCRATCH" MAX_WAIT_SECONDS="${MAX_WAIT_SECONDS:-7200}" \
    bash "$REPO_ROOT/benchmarks/m3pro/queue_large_benchmarks.sh" \
    >> "$SCRATCH/queue_large_benchmarks.log" 2>&1 &
  echo "queue_pid=$!"
  exit 0
fi

echo "No siblings — launching FULL benchmarks now"
PYTHONPATH="$PYTHON_PKG${PYTHONPATH:+:$PYTHONPATH}" python3 -m flexaidds.dataset_runner.launch_queue \
  --scratch "$SCRATCH" \
  --repo-root "$REPO_ROOT" \
  --sibling-count 0 \
  --write-status "$STATUS" \
  --execute

echo "=== conditioned_full_launch complete ==="