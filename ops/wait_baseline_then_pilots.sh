#!/usr/bin/env bash
# When live v_autonomous finishes, run deferred Wave 1–3 pilots from worktree binaries.
# macOS: use ps argv matching (pgrep -a does not print full command lines here).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_ROOT="${1:-$HOME/flexaidds_results/workorders/wave_pilots_$(date +%Y%m%d_%H%M%S)}"
export OUT_ROOT
BASELINE="${FLEXAIDDS_BASELINE_OUT:-$HOME/flexaidds_results/v_autonomous_20260724_160919}"
POLL_SEC="${POLL_SEC:-180}"
LOG="$OUT_ROOT/waiter.log"
mkdir -p "$OUT_ROOT"
# Optional GOAL_SCRATCH for post-pilot copy; export so nested Python sees it.
if [[ -n "${GOAL_SCRATCH:-}" ]]; then
  export GOAL_SCRATCH
fi

baseline_live() {
  # True if any process command line references the baseline OUT path.
  # Uses ps (macOS pgrep -a is PID-only and breaks path greps).
  ps -axo command= 2>/dev/null | grep -F -- "$BASELINE" | grep -v grep | grep -v wait_baseline | grep -v wave_pilots >/dev/null 2>&1
}

{
  echo "Watching baseline: $BASELINE"
  echo "Pilot OUT: $OUT_ROOT"
  echo "POLL_SEC=$POLL_SEC"
} | tee "$LOG"

while baseline_live; do
  n=$(ls "$BASELINE"/*/result.csv 2>/dev/null | wc -l | tr -d ' ')
  echo "$(date -u +%Y-%m-%dT%H:%MZ) still_running result_csv=$n sleep=${POLL_SEC}s" | tee -a "$LOG"
  sleep "$POLL_SEC"
done

n=$(ls "$BASELINE"/*/result.csv 2>/dev/null | wc -l | tr -d ' ')
echo "$(date -u +%Y-%m-%dT%H:%MZ) baseline_exited result_csv=$n" | tee -a "$LOG"

export FLEXAIDDS_BINARY="${FLEXAIDDS_BINARY:-$ROOT/build_wave0/FlexAIDdS}"
export FLEXAIDDS_RUNNER="${FLEXAIDDS_RUNNER:-$ROOT/build_wave0/benchmark_datasets}"
PROBE="${PROBE_CF:-/Users/lp.more/Projects/FlexAIDdS/build/probe_cf}"

if [[ ! -x "$FLEXAIDDS_BINARY" ]]; then
  echo "error: missing binary $FLEXAIDDS_BINARY" | tee -a "$LOG"
  exit 2
fi

echo "=== W2 wall multi-panel ===" | tee -a "$LOG"
python3 "$ROOT/scripts/wall_coercive_oracle.py" \
  --repo "$ROOT" \
  --probe-cf "$PROBE" \
  --binary "$FLEXAIDDS_BINARY" \
  --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
  --panel 1G9V 1M2Z 1N1M 1J3J 1K3U 1L7F 1HNN 1HP0 \
  --out-dir "$OUT_ROOT/w2_wall_oracle" \
  2>&1 | tee "$OUT_ROOT/w2_wall_oracle.log" || true

echo "=== W1.2 elec oracle ===" | tee -a "$LOG"
python3 "$ROOT/scripts/elec_native_cf_oracle.py" \
  --repo "$ROOT" \
  --probe-cf "$PROBE" \
  --binary "$FLEXAIDDS_BINARY" \
  --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
  --out-dir "$OUT_ROOT/w1_elec" \
  2>&1 | tee "$OUT_ROOT/w1_elec.log" || true

echo "=== W1.1 ACF_STRICT dock pilots ===" | tee -a "$LOG"
if [[ -x "$FLEXAIDDS_RUNNER" ]]; then
  # launch script refuses if v_autonomous still in ps
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" off 2>&1 | tee "$OUT_ROOT/acf_off.log" || true
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" on 2>&1 | tee "$OUT_ROOT/acf_on.log" || true
else
  echo "WARN: no runner $FLEXAIDDS_RUNNER" | tee -a "$LOG"
fi

echo "=== W3 E10 snapshot on completed baseline ===" | tee -a "$LOG"
python3 "$ROOT/scripts/e10_election_vs_scoring.py" \
  --campaign-dir "$BASELINE" \
  --out-json "$OUT_ROOT/w3_baseline_e10.json" \
  --out-md "$OUT_ROOT/w3_baseline_e10.md" \
  2>&1 | tee -a "$LOG" || true

# Sync into goal scratch if provided
if [[ -n "${GOAL_SCRATCH:-}" && -d "${GOAL_SCRATCH}" ]]; then
  mkdir -p "$GOAL_SCRATCH/w2_wall_oracle" "$GOAL_SCRATCH/w1_elec" "$GOAL_SCRATCH/w1_acf_strict_pilot" "$GOAL_SCRATCH/w3_sampling"
  cp -f "$OUT_ROOT/w2_wall_oracle/"* "$GOAL_SCRATCH/w2_wall_oracle/" 2>/dev/null || true
  cp -f "$OUT_ROOT/w1_elec/"* "$GOAL_SCRATCH/w1_elec/" 2>/dev/null || true
  cp -f "$OUT_ROOT/w3_baseline_e10."* "$GOAL_SCRATCH/w3_sampling/" 2>/dev/null || true
  cp -f "$OUT_ROOT/"*.log "$GOAL_SCRATCH/" 2>/dev/null || true
  # summarize pilots
  python3 - <<'PY' || true
import json, os
from pathlib import Path
root = Path(os.environ["OUT_ROOT"])
sc = Path(os.environ["GOAL_SCRATCH"])
lines = ["# Auto pilot completion summary\n"]
wj = root / "w2_wall_oracle" / "wall_oracle.json"
if wj.is_file():
    d = json.loads(wj.read_text())
    lines.append(f"- wall_pilot_pass: {d.get('wall_pilot_pass')} n_on={d.get('native_wins_on')}/{d.get('n_scored')}\n")
ej = root / "w1_elec" / "elec_oracle.json"
if ej.is_file():
    d = json.loads(ej.read_text())
    lines.append(f"- elec mass_invert: {d.get('mass_invert')} n_inv={d.get('n_inverted_by_elec')}\n")
(sc / "PILOTS_DONE.md").write_text("".join(lines))
print("".join(lines))
PY
fi

echo "DONE $(date -u +%Y-%m-%dT%H:%MZ) under $OUT_ROOT" | tee -a "$LOG"
