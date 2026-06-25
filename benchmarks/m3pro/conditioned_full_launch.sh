#!/usr/bin/env bash
# conditioned_full_launch.sh — ps-check then launch FULL astex_nonnative + posex benchmarks.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SCRATCH="${SCRATCH:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-1e3f1c67d740/implementer}"
LAUNCHER="$REPO_ROOT/.grok/skills/flexaid-docking/scripts/launch_full_benchmark.sh"
POSEX_JSON="$REPO_ROOT/benchmarks/datasets/posex_cd_1312.json"
STATUS="$SCRATCH/run_status.json"
LOG="$SCRATCH/conditioned_full_launch.log"

mkdir -p "$SCRATCH"
exec > >(tee -a "$LOG") 2>&1

echo "=== conditioned_full_launch $(date -Iseconds) ==="

SIBLING_COUNT=$(pgrep -fl "astex_diverse|benchmark_datasets.*astex" 2>/dev/null | wc -l | tr -d ' ')
echo "Sibling astex_diverse/benchmark processes: $SIBLING_COUNT"

if [[ "$SIBLING_COUNT" -gt 0 ]]; then
  python3 - <<PY
import json, time
status = {
    "status": "waiting_for_astex_diverse_siblings",
    "sibling_process_count": int("$SIBLING_COUNT"),
    "datasets_pending": ["astex_nonnative_1113", "posex_cd_1312"],
    "checked_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    "next_action": "nohup benchmarks/m3pro/queue_large_benchmarks.sh",
}
open("$STATUS", "w").write(json.dumps(status, indent=2))
PY
  echo "Siblings active — deferring to queue_large_benchmarks.sh (detached)"
  nohup env SCRATCH="$SCRATCH" bash "$REPO_ROOT/benchmarks/m3pro/queue_large_benchmarks.sh" >> "$SCRATCH/queue_large_benchmarks.log" 2>&1 &
  echo "queue_pid=$!"
  exit 0
fi

echo "No siblings — launching FULL benchmarks now"
NN_OUT="${FLEXAIDDS_RESULTS:-$SCRATCH}/astex_nonnative_298K_full_$(date +%s)"
POSEX_OUT="${FLEXAIDDS_RESULTS:-$SCRATCH}/posex_cd_298K_full_$(date +%s)"
mkdir -p "$NN_OUT" "$POSEX_OUT"

bash "$LAUNCHER" astex_nonnative 298 "astex_nonnative_298K_vht_full" || true

export PATH="${FLEXAIDDS_BUILD:-/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS/build}:$PATH"
nohup benchmark_datasets \
  --benchmark "crossdock_json:$POSEX_JSON" \
  --output "$POSEX_OUT" \
  --threads 8 \
  >> "$POSEX_OUT/binary.log" 2>> "$POSEX_OUT/stderr.log" </dev/null &
POSEX_PID=$!

python3 - <<PY
import json, time
status = {
    "status": "launched_full",
    "astex_nonnative": {"n_entries": 1113, "launcher": "$LAUNCHER"},
    "posex_cd": {"n_pairs": 1312, "output_dir": "$POSEX_OUT", "pid": $POSEX_PID},
    "launched_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
}
open("$STATUS", "w").write(json.dumps(status, indent=2))
PY
echo "=== conditioned_full_launch complete ==="