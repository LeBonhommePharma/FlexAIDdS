#!/usr/bin/env bash
set -euo pipefail
SCRATCH="${1:-/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-4779e0ecbadf/implementer}"
SLICE_LIST="scripts/rate_slice_targets.txt"
CANONICAL="benchmark_trees_canonical"
TMP_OUT="$SCRATCH/slice_run"
BIN="build_reproduce/benchmark_datasets"
if [ ! -x "$BIN" ]; then
  echo "no build_reproduce/benchmark_datasets; using build_verify if present"
  BIN="build_verify/benchmark_datasets"
fi
if [ ! -x "$BIN" ]; then
  echo "no benchmark_datasets binary; abort"
  exit 1
fi
mkdir -p "$TMP_OUT"
echo "[slice] preparing trees in $TMP_OUT for re-elect detection..."
while read pid; do
  if [ -d "$CANONICAL/$pid" ]; then
    mkdir -p "$TMP_OUT/$pid"
    rsync -a --delete "$CANONICAL/$pid/" "$TMP_OUT/$pid/" || cp -a "$CANONICAL/$pid"/* "$TMP_OUT/$pid/" 2>/dev/null || true
  fi
done < "$SLICE_LIST"
echo "[slice] setting v88 envs + NO_SAS + force re-elect"
export FLEXAIDDS_THERMO=1
export FLEXAIDDS_T_EFF=0.596
export FLEXAIDDS_TENCOM_SCALE=1.0
export FLEXAIDDS_RESTARTS=7
export FLEXAIDDS_PARALLEL_RESTARTS=1
export FLEXAIDDS_NATIVE_SEED_FRAC=0.0
export FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1
export FLEXAIDDS_CONSENSUS_SCORER=1
export FLEXAIDDS_SEED_ELITISM=1
export FLEXAIDDS_NO_SAS=1
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_SOFTCORE_WAL=1
export FLEXAIDDS_SOFTCORE_FLOOR=0.5
export FLEXAIDDS_T_HOT=500
# run with only-codes to restrict to slice, force to re-elect, output to tmp
echo "[slice] invoking C++ benchmark_datasets --only-codes slice --force --output $TMP_OUT"
ONLY_CODES=$(cat "$SLICE_LIST" | tr '\n' ',' | sed 's/,$//')
"$BIN" --benchmark "crossdock_json:benchmarks/datasets/benchmark_astex_native_85.json" \
  --only-codes "$ONLY_CODES" \
  --output "$TMP_OUT" \
  --force \
  --threads 1 \
  --omp-threads 6 \
  --job-timeout-seconds 7200 \
  2>&1 | tee "$SCRATCH/slice_run.log" | tail -20
# the aggregate should be written to $TMP_OUT/astex_crossdock_85_results.csv or similar
# copy to SCRATCH
if [ -f "$TMP_OUT/astex_crossdock_85_results.csv" ]; then
  cp "$TMP_OUT/astex_crossdock_85_results.csv" "$SCRATCH/slice_results.csv"
elif [ -f "$TMP_OUT/results.csv" ]; then
  cp "$TMP_OUT/results.csv" "$SCRATCH/slice_results.csv"
else
  find "$TMP_OUT" -name '*results*.csv' -exec cp {} "$SCRATCH/slice_results.csv" \; 2>/dev/null || true
fi
python3 -c '
import csv, json, datetime, os
p = "'"$SCRATCH"'/slice_results.csv"
if os.path.exists(p):
  rows = list(csv.DictReader(open(p)))
  n = len(rows)
  succ = sum(1 for r in rows if (h := float(r.get("rmsd_hungarian") or 99)) < 2 and h > 0)
  print(json.dumps({"n":n, "succ":succ, "rate": 100*succ/n if n else 0, "timestamp": datetime.datetime.now().isoformat()}, indent=2))
  with open("'"$SCRATCH"'/slice_rate.json", "w") as f: json.dump({"n":n,"succ":succ,"rate":100*succ/n if n else 0}, f)
else:
  print("no slice csv produced")
' 
echo "[slice] done"
