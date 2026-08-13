#!/bin/bash
set -u
REPO=/Users/lp.more/Projects/FlexAIDdS
WT=$REPO/ab_mac_20260806T133329/wt_post_cpu
OUT=$REPO/determinism_e1
cd $REPO
export FLEXAIDDS_BINARY=$WT/build/FlexAID
export FLEXAIDDS_BENCHMARK_DATA=$REPO/benchmarks/astex_diverse/data
export FLEXAID_SEED=12345
export FLEXAIDDS_PARALLEL_REPRODUCE=0
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
export PYTHONPATH=$REPO/python
for rep in 1 2 3 4 5; do
  for T in 1 4; do
    D=$OUT/t${T}_rep${rep}
    [ -f "$D/DONE" ] && continue
    mkdir -p "$D"
    echo "[$(date +%H:%M:%S)] START t=$T rep=$rep free=$(df -g / | awk 'NR==2{print $4}')Gi" >> $OUT/driver.log
    /opt/homebrew/bin/python3 $REPO/benchmarks/run.py --dataset astex_diverse --tier 1 \
      --results-dir "$D" --datasets-dir $REPO/benchmarks/datasets \
      --workers 1 --omp-threads $T --verbose > "$D/run.log" 2>&1
    echo "rc=$?" > "$D/DONE"
    echo "[$(date +%H:%M:%S)] END   t=$T rep=$rep" >> $OUT/driver.log
  done
done
echo "ALL_DONE $(date)" >> $OUT/driver.log
