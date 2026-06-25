#!/usr/bin/env bash
# queue_large_benchmarks.sh — Wait for Astex Diverse sibling runs, then launch
# astex_nonnative + PoseX cross-docking benchmarks (detached, with run_status.json).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SCRATCH="${SCRATCH:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-1e3f1c67d740/implementer}"
LAUNCHER="$REPO_ROOT/.grok/skills/flexaid-docking/scripts/launch_full_benchmark.sh"
POSEX_JSON="$REPO_ROOT/benchmarks/datasets/posex_cd_1312.json"
LOG="$SCRATCH/queue_large_benchmarks.log"

mkdir -p "$SCRATCH"
exec > >(tee -a "$LOG") 2>&1

echo "=== queue_large_benchmarks $(date -Iseconds) ==="

wait_for_astex_diverse_clear() {
  local max_wait_s="${1:-7200}"
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

wait_for_astex_diverse_clear 7200

echo "=== Launching astex_nonnative (298 K dry-run prep via benchmark_datasets if launcher unavailable) ==="
if [[ -x "$LAUNCHER" || -f "$LAUNCHER" ]]; then
  bash "$LAUNCHER" astex_nonnative 298 "astex_nonnative_298K_vht" || {
    echo "launch_full_benchmark astex_nonnative failed (rc=$?); see $LOG"
  }
else
  echo "FATAL: launcher missing at $LAUNCHER"
  exit 1
fi

echo "=== Launching PoseX crossdock_json (1312 pairs) ==="
if [[ -f "$POSEX_JSON" ]] && command -v benchmark_datasets >/dev/null 2>&1; then
  POSEX_OUT="${FLEXAIDDS_RESULTS:-$HOME/flexaidds_results}/posex_cd_298K_$(date +%s)"
  mkdir -p "$POSEX_OUT"
  STATUS="$POSEX_OUT/run_status.json"
  cat > "$STATUS" <<EOF
{
  "status": "launched",
  "dataset": "crossdock_json:posex_cd_1312",
  "n_pairs": 1312,
  "start_time": "$(date -Iseconds)",
  "output_dir": "$POSEX_OUT"
}
EOF
  nohup benchmark_datasets \
    --benchmark "crossdock_json:$POSEX_JSON" \
    --output "$POSEX_OUT" \
    --threads 8 \
    >> "$POSEX_OUT/binary.log" 2>> "$POSEX_OUT/stderr.log" </dev/null &
  echo "{\"benchmark_runner_pid\": $!}" > "$POSEX_OUT/pid.json"
  echo "PoseX launched detached → $POSEX_OUT (status: $STATUS)"
else
  echo "SKIP PoseX: missing $POSEX_JSON or benchmark_datasets not on PATH"
fi

echo "=== queue_large_benchmarks complete $(date -Iseconds) ==="