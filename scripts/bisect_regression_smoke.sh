#!/usr/bin/env bash
# bisect_regression_smoke.sh — git bisect helper for v127→v130 regression (6 targets)
# Gate: ≥5/6 Hungarian success = good (exit 0); else bad (exit 1). Build fail = skip (125).
set -euo pipefail

REPO="${REPO:-$HOME/Projects/FlexAIDdS}"
cd "$REPO"
NCPU="$(sysctl -n hw.logicalcpu 2>/dev/null || echo 4)"
SHA="$(git rev-parse --short HEAD)"
OUT="/tmp/bisect_smoke_${SHA}"
BENCH_JSON="$REPO/benchmarks/datasets/benchmark_bisect_regression_6.json"
TARGETS=(1HQ2 1OF1 1OPK 1P2Y 1T40 1XM6)

echo "[bisect] checkout $SHA" >&2

if ! cmake --build build_lto --target FlexAIDdS benchmark_datasets -j"$NCPU" >/tmp/bisect_build_"${SHA}".log 2>&1; then
  echo "[bisect] build failed — skip" >&2
  tail -20 /tmp/bisect_build_"${SHA}".log >&2 || true
  exit 125
fi

rm -rf "$OUT"
mkdir -p "$OUT"

export FLEXAIDDS_BINARY="$REPO/build_lto/FlexAIDdS"
export FLEXAIDDS_DATA_DIR="$REPO/build_lto"
export FLEXAIDDS_ORACLE_SITE_DIR="$REPO/benchmarks/astex_diverse/astex_diverse"
export FLEXAIDDS_RESTARTS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
export FLEXAIDDS_CONSENSUS_SCORER=1
export FLEXAIDDS_SEED_ELITISM=1
export FLEXAIDDS_N_ELITE=1
export FLEXAIDDS_BUDGET_SCALE=1
export FLEXAIDDS_SOFTCORE_WAL=1
export FLEXAIDDS_SOFTCORE_FLOOR=0.5
export FLEXAIDDS_T_HOT=500
export FLEXAIDDS_NATIVE_SEED_FRAC=0.90
export FLEXAIDDS_VCT_R0=4
export FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0
export OMP_NUM_THREADS=1

"$REPO/build_lto/benchmark_datasets" \
  --benchmark "crossdock_json:${BENCH_JSON}" \
  --output "$OUT" \
  --threads 1 \
  --omp-threads 1 \
  --temperature 298 \
  --job-timeout-seconds 1800 \
  --mode oracle-ceiling \
  --force \
  --only-codes "1HQ2,1OF1,1OPK,1P2Y,1T40,1XM6" \
  >/tmp/bisect_run_"${SHA}".log 2>&1

PASS=0
for t in "${TARGETS[@]}"; do
  csv="$OUT/$t/result.csv"
  if [[ -f "$csv" ]]; then
    s="$(python3 -c "import csv; print(next(csv.DictReader(open('$csv')))['success'])")"
    rmsd="$(python3 -c "import csv; print(next(csv.DictReader(open('$csv')))['rmsd_hungarian'])")"
    echo "[bisect] $t success=$s rmsd=$rmsd" >&2
    if [[ "$s" == "1" ]]; then
      PASS=$((PASS + 1))
    fi
  else
    echo "[bisect] $t MISSING result.csv" >&2
  fi
done

echo "[bisect] $SHA: $PASS/6 pass (gate ≥5 good)" >&2
if [[ "$PASS" -ge 5 ]]; then
  exit 0
fi
exit 1