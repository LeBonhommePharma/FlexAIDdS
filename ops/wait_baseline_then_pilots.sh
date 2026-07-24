#!/usr/bin/env bash
# When live v_autonomous finishes, run deferred Wave 1–3 pilots from worktree binaries.
# Safe: refuses to start while autonomous OUT is still being written by benchmark_datasets.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCRATCH_DEFAULT="${HOME}/flexaidds_results/workorders/wave_pilots_$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${1:-$SCRATCH_DEFAULT}"
BASELINE="${FLEXAIDDS_BASELINE_OUT:-$HOME/flexaidds_results/v_autonomous_20260724_160919}"
POLL_SEC="${POLL_SEC:-120}"

echo "Watching baseline: $BASELINE"
echo "Pilot OUT root: $OUT_ROOT"
mkdir -p "$OUT_ROOT"

while pgrep -f 'v_autonomous_20260724_160919' >/dev/null 2>&1; do
  n=$(ls "$BASELINE"/*/result.csv 2>/dev/null | wc -l | tr -d ' ')
  echo "$(date -u +%H:%M:%SZ) still running; result.csv count=$n; sleep ${POLL_SEC}s"
  sleep "$POLL_SEC"
done
echo "Baseline process gone. result.csv=$(ls "$BASELINE"/*/result.csv 2>/dev/null | wc -l | tr -d ' ')"

export FLEXAIDDS_BINARY="${FLEXAIDDS_BINARY:-$ROOT/build_wave0/FlexAIDdS}"
export FLEXAIDDS_RUNNER="${FLEXAIDDS_RUNNER:-$ROOT/build_wave0/benchmark_datasets}"
PROBE="${PROBE_CF:-$ROOT/../build/probe_cf}"
[[ -x "$PROBE" ]] || PROBE="/Users/lp.more/Projects/FlexAIDdS/build/probe_cf"

if [[ ! -x "$FLEXAIDDS_BINARY" ]]; then
  echo "error: missing $FLEXAIDDS_BINARY — cmake --build build_wave0 --target FlexAIDdS" >&2
  exit 2
fi

# W2 wall multi-panel (score-only)
python3 "$ROOT/scripts/wall_coercive_oracle.py" \
  --repo "$ROOT" \
  --probe-cf "$PROBE" \
  --binary "$FLEXAIDDS_BINARY" \
  --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
  --panel 1G9V 1M2Z 1N1M 1J3J 1K3U 1L7F 1HNN 1HP0 \
  --out-dir "$OUT_ROOT/w2_wall_oracle" \
  | tee "$OUT_ROOT/w2_wall_oracle.log"

# W1.2 elec
python3 "$ROOT/scripts/elec_native_cf_oracle.py" \
  --repo "$ROOT" \
  --probe-cf "$PROBE" \
  --binary "$FLEXAIDDS_BINARY" \
  --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
  --out-dir "$OUT_ROOT/w1_elec" \
  | tee "$OUT_ROOT/w1_elec.log"

# W1.1 ACF_STRICT dock pilots (workers=2 each, sequential arms)
if [[ -x "$FLEXAIDDS_RUNNER" ]]; then
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" off
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" on
else
  echo "WARN: no benchmark_datasets at $FLEXAIDDS_RUNNER — skip dock pilots" | tee -a "$OUT_ROOT/WARN.txt"
fi

# W3 one-var: document control metrics from baseline only (no extra dock unless requested)
python3 "$ROOT/scripts/e10_election_vs_scoring.py" \
  --campaign-dir "$BASELINE" \
  --out-json "$OUT_ROOT/w3_baseline_e10.json" \
  --out-md "$OUT_ROOT/w3_baseline_e10.md"

echo "DONE pilots under $OUT_ROOT"
echo "If wall_oracle.json wall_pilot_pass=true, you may set FLEXAIDDS_WALL_PILOT_PASS=1 for later memetic."
