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
  # True only if a real dock engine/runner still holds the baseline OUT path.
  # Monitors/waiters that merely mention the path must NOT keep us blocked.
  # macOS: ps argv matching (pgrep -a is PID-only here).
  ps -axo command= 2>/dev/null | grep -F -- "$BASELINE" \
    | grep -v grep | grep -v wait_baseline | grep -v wave_pilots \
    | grep -v manual_pilot | grep -v force_pilots | grep -v baseline_wait \
    | grep -E 'FlexAIDdS|benchmark_datasets' >/dev/null 2>&1
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

# Score-only oracles: skip if GOAL_SCRATCH already has multi-panel results (session may have
# pre-run them while baseline was still live). Panel = 6 armA refs (1HNN/1HP0 lack falsemin).
WALL_PANEL=(1G9V 1M2Z 1N1M 1J3J 1K3U 1L7F)
if [[ -n "${GOAL_SCRATCH:-}" && -f "${GOAL_SCRATCH}/w2_wall_oracle/wall_oracle.json" ]]; then
  echo "=== W2 wall multi-panel SKIP (reuse GOAL_SCRATCH) ===" | tee -a "$LOG"
  mkdir -p "$OUT_ROOT/w2_wall_oracle"
  cp -f "${GOAL_SCRATCH}/w2_wall_oracle/"* "$OUT_ROOT/w2_wall_oracle/" 2>/dev/null || true
else
  echo "=== W2 wall multi-panel ===" | tee -a "$LOG"
  python3 "$ROOT/scripts/wall_coercive_oracle.py" \
    --repo "$ROOT" \
    --probe-cf "$PROBE" \
    --binary "$FLEXAIDDS_BINARY" \
    --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
    --panel "${WALL_PANEL[@]}" \
    --out-dir "$OUT_ROOT/w2_wall_oracle" \
    2>&1 | tee "$OUT_ROOT/w2_wall_oracle.log" || true
fi

if [[ -n "${GOAL_SCRATCH:-}" && -f "${GOAL_SCRATCH}/w1_elec/elec_oracle.json" ]]; then
  echo "=== W1.2 elec oracle SKIP (reuse GOAL_SCRATCH) ===" | tee -a "$LOG"
  mkdir -p "$OUT_ROOT/w1_elec"
  cp -f "${GOAL_SCRATCH}/w1_elec/"* "$OUT_ROOT/w1_elec/" 2>/dev/null || true
else
  echo "=== W1.2 elec oracle ===" | tee -a "$LOG"
  python3 "$ROOT/scripts/elec_native_cf_oracle.py" \
    --repo "$ROOT" \
    --probe-cf "$PROBE" \
    --binary "$FLEXAIDDS_BINARY" \
    --data-dir "$(dirname "$FLEXAIDDS_BINARY")" \
    --panel "${WALL_PANEL[@]}" \
    --out-dir "$OUT_ROOT/w1_elec" \
    2>&1 | tee "$OUT_ROOT/w1_elec.log" || true
fi

echo "=== W1.1 ACF_STRICT dock pilots ===" | tee -a "$LOG"
if [[ -x "$FLEXAIDDS_RUNNER" ]]; then
  # launch script refuses if v_autonomous still in ps
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" off 2>&1 | tee "$OUT_ROOT/acf_off.log" || true
  bash "$ROOT/ops/launch_acf_strict_pilot.sh" on 2>&1 | tee "$OUT_ROOT/acf_on.log" || true
  # Resolve pilot OUT dirs from launcher logs / latest pilot_acf_strict_*
  python3 - <<'PY' || true
import re, json
from pathlib import Path
out_root = Path(__import__("os").environ["OUT_ROOT"])
home = Path.home() / "flexaidds_results"
goods = {"1HNN", "1HP0", "1HQ2"}
gap_set = {"1G9V", "1M2Z", "1N1M", "1J3J", "1K3U", "1L7F"}

def find_out(log_name, arm):
    log = out_root / log_name
    if log.is_file():
        m = re.search(r"out=(\S+)", log.read_text(errors="replace"))
        if m:
            return Path(m.group(1))
    cands = sorted(home.glob(f"pilot_acf_strict_{arm}_*"), key=lambda p: p.stat().st_mtime)
    return cands[-1] if cands else None

def genuine_map(camp: Path):
    res = {}
    if not camp or not camp.is_dir():
        return res
    for csvp in camp.glob("*/result.csv"):
        pdb = csvp.parent.name
        try:
            import csv
            with csvp.open() as fh:
                rows = list(csv.DictReader(fh))
            if not rows:
                continue
            r0 = rows[0]
            rms = None
            for k in ("rmsd_hungarian", "rmsd_h", "rmsd", "RMSD"):
                if k in r0 and r0[k] not in (None, ""):
                    try:
                        rms = float(r0[k])
                        break
                    except ValueError:
                        pass
            se = r0.get("seed_echo", "0")
            try:
                se_v = float(se) if se not in (None, "") else 0.0
            except ValueError:
                se_v = 0.0
            bcr = None
            for k in ("bcr", "BCR", "best_cluster_rmsd"):
                if k in r0 and r0[k] not in (None, ""):
                    try:
                        bcr = float(r0[k])
                        break
                    except ValueError:
                        pass
            genuine = rms is not None and rms <= 2.0 and se_v == 0.0
            res[pdb] = {"rmsd": rms, "bcr": bcr, "seed_echo": se_v, "genuine": genuine}
        except Exception as e:
            res[pdb] = {"error": str(e)}
    return res

off_out = find_out("acf_off.log", "off")
on_out = find_out("acf_on.log", "on")
off_m = genuine_map(off_out)
on_m = genuine_map(on_out)
goods_flip = []
for g in sorted(goods):
    go = off_m.get(g, {}).get("genuine")
    gn = on_m.get(g, {}).get("genuine")
    if go is True and gn is False:
        goods_flip.append(g)
gap_move = []
for g in sorted(gap_set):
    o, n = off_m.get(g, {}), on_m.get(g, {})
    if o and n and o.get("rmsd") is not None and n.get("rmsd") is not None:
        gap_move.append({
            "pdb": g,
            "rmsd_off": o["rmsd"], "rmsd_on": n["rmsd"],
            "bcr_off": o.get("bcr"), "bcr_on": n.get("bcr"),
            "genuine_off": o.get("genuine"), "genuine_on": n.get("genuine"),
        })
summary = {
    "acf_off_out": str(off_out) if off_out else None,
    "acf_on_out": str(on_out) if on_out else None,
    "goods_success_to_fail": goods_flip,
    "goods_non_regression": len(goods_flip) == 0 and bool(off_m or on_m),
    "gap_targets": gap_move,
    "n_off": len(off_m),
    "n_on": len(on_m),
}
(out_root / "w1_acf_ab_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
md = ["# W1.1 ACF_STRICT dock A/B summary\n",
      f"- off OUT: `{off_out}`\n", f"- on OUT: `{on_out}`\n",
      f"- goods success→fail flips: **{goods_flip or 'none'}**\n",
      f"- goods_non_regression: **{summary['goods_non_regression']}**\n",
      f"- targets off/on: {summary['n_off']}/{summary['n_on']}\n"]
if gap_move:
    md.append("\n| pdb | rmsd_off | rmsd_on | genuine_off | genuine_on |\n|-----|----------|---------|-------------|------------|\n")
    for r in gap_move:
        md.append(f"| {r['pdb']} | {r['rmsd_off']} | {r['rmsd_on']} | {r['genuine_off']} | {r['genuine_on']} |\n")
(out_root / "w1_acf_ab_summary.md").write_text("".join(md))
print("".join(md))
PY
else
  echo "WARN: no runner $FLEXAIDDS_RUNNER" | tee -a "$LOG"
fi

echo "=== W3 BCR one-variable pilot (1J3J/1K3U) ===" | tee -a "$LOG"
export W3_BCR_OUT="$OUT_ROOT/w3_bcr_pilot"
bash "$ROOT/ops/launch_w3_bcr_pilot.sh" 2>&1 | tee "$OUT_ROOT/w3_bcr_pilot.log" || true
if [[ -f "$OUT_ROOT/w3_bcr_pilot/w3_bcr_summary.md" ]]; then
  cp -f "$OUT_ROOT/w3_bcr_pilot/w3_bcr_summary."* "$OUT_ROOT/" 2>/dev/null || true
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
  cp -f "$OUT_ROOT/w1_acf_ab_summary."* "$GOAL_SCRATCH/w1_acf_strict_pilot/" 2>/dev/null || true
  # Evaluator / flip gate also looks at GOAL_SCRATCH root for these names
  cp -f "$OUT_ROOT/w1_acf_ab_summary."* "$GOAL_SCRATCH/" 2>/dev/null || true
  mkdir -p "$GOAL_SCRATCH/w3_sampling"
  cp -f "$OUT_ROOT/w3_bcr_pilot/w3_bcr_summary."* "$GOAL_SCRATCH/w3_sampling/" 2>/dev/null || true
  cp -f "$OUT_ROOT/w3_bcr_pilot/w3_bcr_summary."* "$GOAL_SCRATCH/" 2>/dev/null || true
  cp -f "$OUT_ROOT/w3_bcr_summary."* "$GOAL_SCRATCH/" 2>/dev/null || true
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
aj = root / "w1_acf_ab_summary.json"
if aj.is_file():
    d = json.loads(aj.read_text())
    lines.append(f"- acf goods_non_regression: {d.get('goods_non_regression')} flips={d.get('goods_success_to_fail')}\n")
    lines.append(f"- acf off/on targets: {d.get('n_off')}/{d.get('n_on')}\n")
bj = root / "w3_bcr_pilot" / "w3_bcr_summary.json"
if bj.is_file():
    d = json.loads(bj.read_text())
    lines.append(f"- w3_bcr pairs: {len(d.get('pairs') or [])}\n")
(sc / "PILOTS_DONE.md").write_text("".join(lines))
print("".join(lines))
PY
fi

echo "DONE $(date -u +%Y-%m-%dT%H:%MZ) under $OUT_ROOT" | tee -a "$LOG"
