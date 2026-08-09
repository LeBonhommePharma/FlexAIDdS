#!/bin/bash
# Reproduce the cleft-grid nondeterminism that drives Astex campaign variance.
#
# The SURFNET probe merge in CleftDetector::generate_probes runs under
# `omp critical` in thread-arrival order. That order permutes cleftgrid index
# assignment, and cleftgrid index IS GA gene 0 — so the same seed addresses
# different anchors. Background CPU load is REQUIRED to expose it: on an idle
# box the arrival order repeats and the defect hides.
#
# Usage: ./repro_cleft_determinism.sh [PDBID]   (default 1GPK; 1G9V diverges harder)
set -u
set +m          # no job-control chatter when the load generators are killed
TARGET=${1:-1GPK}
REPO=${FLEXAIDDS_REPO:-$(cd "$(dirname "$0")/.." && pwd)}
BIN=${FLEXAIDDS_BINARY:-$REPO/build_tests/FlexAIDdS}
DATA=$REPO/benchmarks/astex_diverse/astex_diverse/$TARGET
WORK=$(mktemp -d)
REPS=${REPS:-6}
LOAD=${LOAD:-8}

[ -x "$BIN" ] || { echo "no engine at $BIN (set FLEXAIDDS_BINARY)"; exit 1; }
[ -d "$DATA" ] || { echo "no target data at $DATA"; exit 1; }

# Small budget: cleft detection and the .rrg write happen long before this matters.
cat > "$WORK/fast.json" <<'EOF'
{ "ga": { "num_chromosomes": 200, "num_generations": 30 } }
EOF

run () {  # $1=tag  $2=threads  $3=cleft_sort
  env FLEXAIDDS_CLEFT_SORT="$3" \
      FLEXAIDDS_ORACLE_SITE="$DATA/${TARGET}_binding_site.pdb" \
      OMP_NUM_THREADS="$2" OMP_PLACES=cores OMP_PROC_BIND=spread \
      FLEXAID_SEED=12345 FLEXAIDDS_PARALLEL_RESTARTS=0 FLEXAIDDS_SEED_ELITISM=0 \
      timeout 300 "$BIN" "$DATA/${TARGET}_apo.pdb" "$DATA/${TARGET}_ligand.sdf" \
      -c "$WORK/fast.json" -o "$WORK/$1" > "$WORK/$1.log" 2>&1
}

distinct () { md5 -q "$WORK/$1"*.rrg 2>/dev/null | sort -u | wc -l | tr -d ' '; }

hogs () {  # background CPU load perturbs OpenMP thread arrival order
  HOGPIDS=""
  for _ in $(seq "$LOAD"); do yes > /dev/null & HOGPIDS="$HOGPIDS $!"; done
}
unhogs () { kill $HOGPIDS 2>/dev/null; }

echo "target=$TARGET  engine=$BIN  reps=$REPS  load=$LOAD"
echo

hogs
for r in $(seq $REPS); do run "off_r$r" 4 0; done
unhogs
echo "gate OFF, 4 threads, under load : $(distinct off_r) distinct grids of $REPS"
md5 -q "$WORK"/off_r*.rrg | sort | uniq -c

hogs
for r in $(seq $REPS); do run "on_r$r" 4 1; done
unhogs
echo "gate ON,  4 threads, under load : $(distinct on_r) distinct grids of $REPS"

run t1_a 1 1; run t1_b 1 1
echo "gate ON,  1 thread              : $(distinct t1_) distinct grids of 2"
echo "gate ON,  1 vs 4 threads        : $(md5 -q "$WORK"/on_r1.rrg "$WORK"/t1_a.rrg | sort -u | wc -l | tr -d ' ') distinct (1 == thread-count invariant)"

echo
echo "artifacts: $WORK"
