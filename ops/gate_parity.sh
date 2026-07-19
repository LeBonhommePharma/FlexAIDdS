#!/bin/bash
# METHODOLOGY.md §1 — parity gate. Usage: gate_parity.sh <engineA> <engineB> [target]
set -u
REPO="/Users/lp.more/Projects/FlexAIDdS"; cd "$REPO"
A="$1"; B="$2"; T="${3:-1G9V}"
CACHE=/tmp/ab_bench_cache/astex_diverse
[ -d "$CACHE/$T" ] || CACHE="$REPO/benchmarks/astex_diverse/astex_diverse"
rec="$CACHE/$T/${T}_apo.pdb"; lig="$CACHE/$T/${T}_ligand.sdf"
export FLEXAID_SEED=12345 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_RESTARTS=1 FLEXAIDDS_DATA_DIR="$REPO/build"
cfg=/tmp/parity.json
for eng in "$A" "$B"; do
  tag=$(basename "$eng"); out=/tmp/gate_parity/$tag/$T; mkdir -p "$out"
  OMP_NUM_THREADS=1 "$eng" "$rec" "$lig" -c "$cfg" -o "$out/d" > "$out/run.log" 2>&1
done
# compare 10 elected poses byte-identical
ta=$(basename "$A"); tb=$(basename "$B"); ok=1
for i in $(seq 0 9); do
  fa=/tmp/gate_parity/$ta/$T/d_${i}.pdb; fb=/tmp/gate_parity/$tb/$T/d_${i}.pdb
  [ -f "$fa" ] && [ -f "$fb" ] || { echo "MISSING pose $i"; ok=0; break; }
  cmp -s "$fa" "$fb" || { echo "POSE $i DIFFERS"; ok=0; }
done
[ "$ok" = 1 ] && echo "PARITY PASS ($T, 10/10 poses byte-identical)" || echo "PARITY FAIL ($T)"
