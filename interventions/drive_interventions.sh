#!/bin/bash
set -u
REPO=/Users/lp.more/Projects/FlexAIDdS
WT=$REPO/ab_mac_20260806T133329/wt_post_cpu
OUT=$REPO/interventions
cd $WT
export FLEXAIDDS_BENCHMARK_DATA=$WT/benchmarks/astex_diverse
export FLEXAID_SEED=12345
export FLEXAIDDS_PARALLEL_REPRODUCE=0
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
export PYTHONPATH=$WT/python
run_cell () {  # $1=arm  $2=rep
  local D=$OUT/$1_rep$2
  [ -f "$D/DONE" ] && return
  mkdir -p "$D"
  echo "[$(date +%H:%M:%S)] START $1 rep=$2 avail=$(df -g / | awk 'NR==2{print $4}')Gi" >> $OUT/driver.log
  /opt/homebrew/bin/python3 $WT/benchmarks/run.py --dataset astex_diverse --tier 1 \
    --results-dir "$D" --datasets-dir $WT/benchmarks/datasets \
    --workers 1 --omp-threads 4 --verbose > "$D/run.log" 2>&1
  echo "rc=$?" > "$D/DONE"
  echo "[$(date +%H:%M:%S)] END   $1 rep=$2" >> $OUT/driver.log
}
for rep in 1 2 3; do
  ( export FLEXAIDDS_BINARY=$WT/build/FlexAID;        export FLEXAIDDS_CLEFT_SORT=1; run_cell I1_cleftsort   $rep )
  ( export FLEXAIDDS_BINARY=$OUT/FlexAID_mifoff.sh;                                  run_cell I2_mifoff      $rep )
  ( export FLEXAIDDS_BINARY=$WT/build/FlexAID;        export FLEXAID_DETERMINISTIC=1; run_cell I3_determ     $rep )
done
echo "ALL_DONE $(date)" >> $OUT/driver.log
