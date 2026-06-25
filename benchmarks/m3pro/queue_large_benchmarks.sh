#!/usr/bin/env bash
# queue_large_benchmarks.sh — Wait for Astex Diverse sibling runs, then launch
# astex_nonnative + PoseX cross-docking benchmarks (detached, with run_status.json).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SCRATCH="${SCRATCH:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-1e3f1c67d740/implementer}"
STATUS="$SCRATCH/run_status.json"
LOG="$SCRATCH/queue_large_benchmarks.log"
PYTHON_PKG="$REPO_ROOT/python"
MAX_WAIT_SECONDS="${MAX_WAIT_SECONDS:-7200}"

mkdir -p "$SCRATCH"
exec > >(tee -a "$LOG") 2>&1

echo "=== queue_large_benchmarks $(date -Iseconds) MAX_WAIT_SECONDS=$MAX_WAIT_SECONDS ==="

wait_for_astex_diverse_clear() {
  local max_wait_s="${1:-7200}"
  if [[ "$max_wait_s" -le 0 ]]; then
    echo "MAX_WAIT_SECONDS=0 — skipping sibling wait, proceeding to launch."
    return 0
  fi
  local elapsed=0
  local interval=60
  while (( elapsed < max_wait_s )); do
    if ! pgrep -fl "astex_diverse|benchmark_datasets.*astex" >/dev/null 2>&1; then
      echo "No active astex_diverse / astex benchmark processes detected."
      return 0
    fi
    echo "[$(date +%H:%M:%S)] Astex Diverse sibling still active — waiting ${interval}s..."
    pgrep -fl "astex_diverse|benchmark_datasets.*astex" | head -5 || true
    sleep "$interval"
    elapsed=$((elapsed + interval))
  done
  echo "WARNING: timed out waiting for sibling Astex Diverse runs; proceeding anyway."
  return 0
}

wait_for_astex_diverse_clear "$MAX_WAIT_SECONDS"

echo "=== Launching astex_nonnative (1113) + posex_cd (1312) ==="
PYTHONPATH="$PYTHON_PKG${PYTHONPATH:+:$PYTHONPATH}" python3 -m flexaidds.dataset_runner.launch_queue \
  --scratch "$SCRATCH" \
  --repo-root "$REPO_ROOT" \
  --sibling-count 0 \
  --write-status "$STATUS" \
  --execute

echo "=== queue_large_benchmarks complete $(date -Iseconds) ==="