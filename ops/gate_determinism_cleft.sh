#!/bin/bash
# METHODOLOGY.md §2 — cleft-grid determinism. Usage: gate_determinism_cleft.sh <engine> [targets...]
set -u
REPO="/Users/lp.more/Projects/FlexAIDdS"; cd "$REPO"
ENG="$1"; shift; TARGETS="${*:-1G9V 1GM8 1GPK 1HP0 1HNN}"
CACHE=/tmp/ab_bench_cache/astex_diverse
[ -d "$CACHE/1G9V" ] || CACHE="$REPO/benchmarks/astex_diverse/astex_diverse"
export FLEXAID_SEED=12345 FLEXAIDDS_NO_SEC=1 FLEXAIDDS_RESTARTS=1 FLEXAIDDS_DATA_DIR="$REPO/build"
python3 -c "import json;c=json.load(open('/tmp/parity.json'));c['ga']['num_generations']=20;json.dump(c,open('/tmp/cleft_gate.json','w'))"
pass=0; total=0
for T in $TARGETS; do
  rec="$CACHE/$T/${T}_apo.pdb"; lig="$CACHE/$T/${T}_ligand.sdf"
  [ -f "$rec" ] || continue
  o=/tmp/gate_cleft/$T; mkdir -p "$o/a" "$o/b"
  OMP_NUM_THREADS=4 "$ENG" "$rec" "$lig" -c /tmp/cleft_gate.json -o "$o/a/x" >/dev/null 2>&1
  OMP_NUM_THREADS=4 "$ENG" "$rec" "$lig" -c /tmp/cleft_gate.json -o "$o/b/x" >/dev/null 2>&1
  total=$((total+1))
  if cmp -s "$o/a/x.rrg" "$o/b/x.rrg"; then echo "$T: rrg IDENTICAL"; pass=$((pass+1)); else echo "$T: rrg DIFFERS"; fi
done
echo "DETERMINISM $pass/$total grids reproducible"
