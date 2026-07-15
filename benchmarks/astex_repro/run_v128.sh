#!/bin/bash
# v128 baseline — 2BYS fixed (dimer A+E, 3336 atoms, commit 8f7828872)
# Grand partition: 39/39 tests passing
export PATH=/opt/homebrew/bin:$PATH
export FLEXAIDDS_BINARY="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine/FlexAIDdS"
export FLEXAIDDS_DATA_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine"
export FLEXAIDDS_CLEFT_SPHERE_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/spheres"
export FLEXAIDDS_ORACLE_SITE_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse"
export FLEXAIDDS_RESTARTS=10
export FLEXAIDDS_PARALLEL_RESTARTS=0  # 18GB contract: one FlexAIDdS child at a time
export FLEXAIDDS_NO_SEC=1
# NO --force => skip_completed resumes the 17 targets already copied from v127
"/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine/benchmark_datasets" \
  --benchmark astex \
  --cache "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse" \
  --output "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/full_v128" \
  --mode autonomous \
  --threads 1 \
  --omp-threads 1 \
  --ga-population 1000 \
  --ga-generations 2000 \
  --job-timeout-seconds 10800
