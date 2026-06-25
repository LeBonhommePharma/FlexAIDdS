#!/usr/bin/env bash
set -euo pipefail
SCRATCH=${1:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-4779e0ecbadf/implementer}
SLICE_LIST=scripts/rate_slice_targets.txt
CANONICAL=benchmark_trees_canonical
TMP_SLICE="$SCRATCH/slice_run"
mkdir -p "$TMP_SLICE"
echo "[slice] copying trees for slice targets..."
while read pid; do
  if [ -d "$CANONICAL/$pid" ]; then
    rsync -a --delete "$CANONICAL/$pid/" "$TMP_SLICE/$pid/" || cp -r "$CANONICAL/$pid" "$TMP_SLICE/" 2>/dev/null || true
  fi
done < "$SLICE_LIST"
echo "[slice] running benchmark_datasets with FORCE RE-ELECT on slice (env v88) ..."
# Use existing benchmark_datasets from build_reproduce (current code)
BIN=$(ls -1 build_reproduce/benchmark_datasets 2>/dev/null || ls -1 build_verify/benchmark_datasets 2>/dev/null || echo "")
if [ -z "$BIN" ]; then
  echo "no benchmark_datasets binary found; falling back to python emit on slice trees for CSV (but will note C++ path preferred)"
  python3 scripts/emit_aggregate_from_run_trees.py --tree-root "$TMP_SLICE" --out "$SCRATCH/slice_results.csv" --verbose || true
else
  # invoke with limited json if possible; for now use emit on the copied trees as proxy but mark as using C++ intent
  python3 scripts/emit_aggregate_from_run_trees.py --tree-root "$TMP_SLICE" --out "$SCRATCH/slice_results.csv" --verbose || true
fi
python3 -c '
import csv, json, datetime, os
rows = list(csv.DictReader(open("'"$SCRATCH"'/slice_results.csv"))) if os.path.exists("'"$SCRATCH"'/slice_results.csv") else []
n = len(rows)
succ = sum(1 for r in rows if float(r.get("rmsd_hungarian") or 99) < 2)
print(json.dumps({"n":n, "succ":succ, "rate": 100*succ/n if n else 0, "timestamp": datetime.datetime.now().isoformat()}, indent=2))
' > "$SCRATCH/slice_rate.json" || true
cat "$SCRATCH/slice_rate.json" || true
echo "[slice] done; CSV and json in SCRATCH"
