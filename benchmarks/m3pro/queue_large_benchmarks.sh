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

wait_for_astex_nonnative_clear() {
  local max_wait_s="${1:-0}"
  local interval=120
  local elapsed=0
  while pgrep -fl "--benchmark astex_nonnative|astex_nonnative_298K" >/dev/null 2>&1; do
    if [[ "$max_wait_s" -gt 0 && "$elapsed" -ge "$max_wait_s" ]]; then
      echo "WARNING: astex_nonnative still running after ${max_wait_s}s; launching posex anyway."
      return 0
    fi
    echo "[$(date +%H:%M:%S)] astex_nonnative resume still active — waiting ${interval}s before posex_cd..."
    pgrep -fl "--benchmark astex_nonnative|astex_nonnative_298K" | head -3 || true
    sleep "$interval"
    elapsed=$((elapsed + interval))
  done
  echo "astex_nonnative clear — safe to launch posex_cd."
}

if pgrep -fl "--benchmark astex_nonnative|astex_nonnative_298K" >/dev/null 2>&1; then
  echo "Detected running astex_nonnative — will skip relaunch and wait before posex_cd."
  export SKIP_ASTEX_NONNATIVE=1
  wait_for_astex_nonnative_clear "${WAIT_FOR_ASTEX_NONNATIVE_SECONDS:-0}"
fi

echo "=== Launching posex_cd (1312); astex_nonnative skipped if already running ==="
LAUNCH_ARGS=(--scratch "$SCRATCH" --repo-root "$REPO_ROOT" --sibling-count 0 --write-status "$STATUS")
if [[ "${LAUNCH_DRY_RUN:-0}" == "1" ]]; then
  LAUNCH_ARGS+=(--execute --dry-run)
else
  LAUNCH_ARGS+=(--execute)
fi
PYTHONPATH="$PYTHON_PKG${PYTHONPATH:+:$PYTHONPATH}" python3 -m flexaidds.dataset_runner.launch_queue "${LAUNCH_ARGS[@]}"

echo "=== queue_large_benchmarks complete $(date -Iseconds) ==="