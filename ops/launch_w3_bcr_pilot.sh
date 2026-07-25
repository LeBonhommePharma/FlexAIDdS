#!/usr/bin/env bash
# W3 one-variable sampling pilot: default-OFF vs BOOM_INTERVAL×SIGMA_SCALE on 1J3J+1K3U.
# Run ONLY after v_autonomous baseline is free. Workers≤2.
#
# Usage:
#   export FLEXAIDDS_BINARY=.../build_wave0/FlexAIDdS
#   export FLEXAIDDS_RUNNER=.../build_wave0/benchmark_datasets
#   bash ops/launch_w3_bcr_pilot.sh
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Refuse while autonomous full85 still holds the box (path-based on macOS).
# Match only real FlexAIDdS/benchmark_datasets; monitors that mention the path are OK.
BASELINE="${FLEXAIDDS_BASELINE_OUT:-$HOME/flexaidds_results/v_autonomous_20260724_160919}"
if ps -axo command= 2>/dev/null | grep -F -- "$BASELINE" \
    | grep -v grep | grep -v wait_baseline | grep -v wave_pilots \
    | grep -v manual_pilot | grep -v force_pilots | grep -v pilot_ \
    | grep -E 'FlexAIDdS|benchmark_datasets' >/dev/null 2>&1; then
  echo "REFUSE: baseline still running at $BASELINE" >&2
  exit 92
fi

BIN="${FLEXAIDDS_BINARY:-}"
RUNNER="${FLEXAIDDS_RUNNER:-$ROOT/build_wave0/benchmark_datasets}"
if [[ -z "$BIN" || ! -x "$BIN" ]]; then
  echo "error: set FLEXAIDDS_BINARY" >&2
  exit 2
fi
if [[ ! -x "$RUNNER" ]]; then
  echo "error: runner not found: $RUNNER" >&2
  exit 2
fi

STAMP=$(date +%Y%m%d_%H%M%S)
OUT_ROOT="${W3_BCR_OUT:-$HOME/flexaidds_results/pilot_w3_bcr_${STAMP}}"
mkdir -p "$OUT_ROOT"
ONLY="1J3J,1K3U"
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
export FLEXAIDDS_BINARY="$BIN"
# Never arm memetic (wall pilot failed)
unset FLEXAIDDS_MEMETIC FLEXAIDDS_WALL_PILOT_PASS || true

run_arm() {
  local arm="$1"
  local out="$OUT_ROOT/$arm"
  mkdir -p "$out"
  if [[ "$arm" == "default" ]]; then
    unset FLEXAIDDS_BOOM_INTERVAL FLEXAIDDS_BOOM_FRAC FLEXAIDDS_SIGMA_SCALE || true
  else
    # One-variable anti-collapse: shorter BOOM interval + tighter niche scale
    export FLEXAIDDS_BOOM_INTERVAL=50
    export FLEXAIDDS_SIGMA_SCALE=0.5
  fi
  python3 - <<PY
import json, os
from pathlib import Path
out = Path("$out")
env = {k: v for k, v in os.environ.items() if k.startswith("FLEXAIDDS_") or k == "OMP_NUM_THREADS"}
(out / "SCORING_PROVENANCE.json").write_text(json.dumps({
    "arm": "$arm",
    "binary": "$BIN",
    "panel": "$ONLY",
    "scoring_env": env,
    "workers": 2,
}, indent=2) + "\n")
print("wrote", out / "SCORING_PROVENANCE.json")
PY
  echo "=== W3 BCR arm=$arm out=$out ==="
  caffeinate -i "$RUNNER" \
    --benchmark astex \
    --mode autonomous \
    --output "$out" \
    --threads 2 \
    --omp-threads 1 \
    --only-codes "$ONLY" \
    --job-timeout-seconds 3600 \
    --force
}

run_arm default
run_arm boom_sigma

# Summarize BCR / elected rmsd / freq proxies from result.csv
python3 - "$OUT_ROOT" <<'PY'
import csv, json, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = []
for arm in ("default", "boom_sigma"):
    camp = root / arm
    for csvp in sorted(camp.glob("*/result.csv")):
        pdb = csvp.parent.name
        with csvp.open() as fh:
            r = next(csv.DictReader(fh), None)
        if not r:
            continue
        def fget(*keys):
            for k in keys:
                if k in r and r[k] not in (None, ""):
                    try:
                        return float(r[k])
                    except ValueError:
                        pass
            return None
        rows.append({
            "arm": arm,
            "pdb": pdb,
            "rmsd": fget("rmsd_hungarian", "rmsd_h", "rmsd"),
            "bcr": fget("bcr", "BCR", "best_cluster_rmsd"),
            "freq": fget("freq", "cluster_freq", "elected_freq"),
            "seed_echo": fget("seed_echo"),
        })

by = {}
for r in rows:
    by.setdefault(r["pdb"], {})[r["arm"]] = r

summary = {"out_root": str(root), "rows": rows, "pairs": []}
md = ["# W3 BCR one-variable pilot (default vs BOOM_INTERVAL=50 × SIGMA_SCALE=0.5)\n",
      f"OUT: `{root}`\n\n",
      "| pdb | bcr_default | bcr_boom_sigma | rmsd_default | rmsd_boom_sigma | freq_default | freq_boom |\n",
      "|-----|-------------|----------------|--------------|-----------------|--------------|----------|\n"]
for pdb in sorted(by):
    d = by[pdb].get("default", {})
    b = by[pdb].get("boom_sigma", {})
    pair = {
        "pdb": pdb,
        "bcr_default": d.get("bcr"), "bcr_boom_sigma": b.get("bcr"),
        "rmsd_default": d.get("rmsd"), "rmsd_boom_sigma": b.get("rmsd"),
        "freq_default": d.get("freq"), "freq_boom_sigma": b.get("freq"),
    }
    summary["pairs"].append(pair)
    md.append(
        f"| {pdb} | {d.get('bcr')} | {b.get('bcr')} | {d.get('rmsd')} | {b.get('rmsd')} | {d.get('freq')} | {b.get('freq')} |\n"
    )
(root / "w3_bcr_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
(root / "w3_bcr_summary.md").write_text("".join(md))
print("".join(md))
print("OUT_ROOT", root)
PY

echo "=== W3 BCR pilot done under $OUT_ROOT ==="
echo "$OUT_ROOT" > "$OUT_ROOT/OUT_ROOT.txt"
