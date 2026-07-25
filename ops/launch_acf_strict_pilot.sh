#!/usr/bin/env bash
# W1.1 serial ACF_STRICT on/off pilot — 9 targets, workers=2.
# DOES NOT use Projects/FlexAIDdS/build while v_autonomous_* is live.
# Prefer: worktree build_wave0/FlexAIDdS or baseline_engine copy AFTER baseline finishes.
#
# Usage (after baseline done):
#   export FLEXAIDDS_BINARY=/path/to/worktree/build_wave0/FlexAIDdS
#   bash ops/launch_acf_strict_pilot.sh off   # ACF_STRICT unset
#   bash ops/launch_acf_strict_pilot.sh on    # ACF_STRICT=1
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
ARM="${1:-}"
if [[ "$ARM" != "on" && "$ARM" != "off" ]]; then
  echo "Usage: $0 on|off" >&2
  exit 2
fi

# Refuse if live autonomous campaign still holds the shared box.
# macOS: match process command lines on the OUT path; exclude waiters/triggers
# that only *mention* the path in their own scripts.
BASELINE_OUT="${FLEXAIDDS_BASELINE_OUT:-$HOME/flexaidds_results/v_autonomous_20260724_160919}"
if ps -axo command= 2>/dev/null | grep -F -- "$BASELINE_OUT" \
    | grep -v grep | grep -v wait_baseline | grep -v wave_pilots \
    | grep -v manual_pilot | grep -v pilot_acf | grep -v pilot_w3 \
    | grep -E 'FlexAIDdS|benchmark_datasets' >/dev/null 2>&1; then
  echo "REFUSE: v_autonomous still running at $BASELINE_OUT — do not start pilot docks" >&2
  exit 92
fi

BIN="${FLEXAIDDS_BINARY:-}"
if [[ -z "$BIN" || ! -x "$BIN" ]]; then
  echo "error: set FLEXAIDDS_BINARY to isolated FlexAIDdS (not live mmap binary)" >&2
  exit 2
fi

STAMP=$(date +%Y%m%d_%H%M%S)
OUT="${HOME}/flexaidds_results/pilot_acf_strict_${ARM}_${STAMP}"
export OUT
mkdir -p "$OUT"
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_BINARY="$BIN"
if [[ "$ARM" == "on" ]]; then
  export FLEXAIDDS_ACF_STRICT=1
else
  unset FLEXAIDDS_ACF_STRICT || true
fi

# Full scoring env into receipt sidecar
python3 - <<PY
import json, os
from pathlib import Path
out = Path(os.environ["OUT"] if "OUT" in os.environ else "$OUT")
env = {k: v for k, v in os.environ.items() if k.startswith("FLEXAIDDS_") or k in ("OMP_NUM_THREADS",)}
(out / "SCORING_PROVENANCE.json").write_text(json.dumps({
    "arm": "$ARM",
    "binary": "$BIN",
    "scoring_env": env,
    "panel": "1G9V,1M2Z,1N1M,1J3J,1K3U,1L7F,1HNN,1HP0,1HQ2",
    "workers": 2,
}, indent=2) + "\n")
print("wrote", out / "SCORING_PROVENANCE.json")
PY

ONLY="1G9V,1M2Z,1N1M,1J3J,1K3U,1L7F,1HNN,1HP0,1HQ2"
RUNNER="${FLEXAIDDS_RUNNER:-$ROOT/build_wave0/benchmark_datasets}"
if [[ ! -x "$RUNNER" ]]; then
  echo "error: benchmark_datasets not found at $RUNNER — build in worktree build_wave0" >&2
  exit 2
fi

echo "=== ACF_STRICT pilot arm=$ARM out=$OUT binary=$BIN ==="
caffeinate -i "$RUNNER" \
  --benchmark astex \
  --mode autonomous \
  --output "$OUT" \
  --threads 2 \
  --omp-threads 1 \
  --only-codes "$ONLY" \
  --job-timeout-seconds 3600 \
  --force
echo "=== done arm=$ARM ==="
