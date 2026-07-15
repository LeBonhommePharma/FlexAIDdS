#!/bin/bash
export PATH=/opt/homebrew/bin:$PATH
export FLEXAIDDS_BINARY="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine/FlexAIDdS"
# ENERGY MATRIX (2026-07-10): pin FLEXAIDDS_DATA_DIR to the engine dir so BOTH
# the per-job dock AND the provenance recorder (DatasetRunner.cpp:4971) resolve
# the SAME MC_st0r5.2_6.dat (md5 9dc93717...) deterministically. Without this,
# provenance.json recorded matrix_path="" / matrix_md5="" (the <bin>/../WRK
# fallback does not exist) — a reproducibility hole. The interaction/energy
# matrix is the complementarity function's core input and must always be loaded
# and provenance-anchored, never silently defaulted.
export FLEXAIDDS_DATA_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine"
export FLEXAIDDS_CLEFT_SPHERE_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/spheres"
export FLEXAIDDS_ORACLE_SITE_DIR="/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse"
export FLEXAIDDS_RESTARTS=10
# CRITICAL for benchmark parity (fix 35af1b3f8): disable the entropy/stagnation
# early-exit so the GA uses the FULL generation budget with exploration boost,
# matching methods that always run to completion. Without this the GA terminated
# at ~gen 400 of 3500 (the 22.6% under-budget result).
export FLEXAIDDS_NO_SEC=1
#
# ── POSE-EMISSION FIX (2026-07-10) ───────────────────────────────────────────
# ROOT CAUSE of the 11/16 "zero-pose" targets: all 10 restarts fork CONCURRENTLY
# (FLEXAIDDS_PARALLEL_RESTARTS default ON) and SHARE ONE per_job_timeout budget
# (DatasetRunner.cpp:6044, per_job_timeout_s - elapsed). On an 18 GiB / ~0-free
# machine, 10 concurrent GAs swap and crawl; flexible ligands (n_gen_effective
# =3500) only reach gen ~2000-2750 when the 3600 s wall fires -> all restarts
# SIGKILLed mid-GA, BEFORE the clustering/emission stage -> num_poses=0. Rigid
# ligands (n_gen~2000) just finish in time, which is why 5/16 emitted poses.
# The cleft/grid definition is NOT the cause (ligand coverage=1.0, center-to-
# ligand <2.7 A, grid size uncorrelated with failure).
#
# ENV-ONLY FIX (no rebuild): raise the per-complex wall-clock budget from 3600 s
# to 10800 s (3 h) so even under contention every restart reaches its full
# generation budget AND runs the clustering/emission stage before SIGTERM.
# --job-timeout-seconds already exists (benchmark_datasets.cpp:592) and is
# compiled into the pinned engine binary (verified via strings).
#
# NOTE (RAM): the 18 GiB cap is still pressured by 10 concurrent GAs.
# FLEXAIDDS_PARALLEL_RESTARTS is a BOOLEAN (on/off), NOT a concurrency count, so
# there is no env-only way to cap concurrency at N<10. omp-threads lowered to 1
# to remove CPU oversubscription (10 restarts x 1 thread = 10 <= 11 cores),
# which lets each GA run faster and reach budget sooner. A true concurrency cap
# (batched forking) needs a code change + rebuild — proposed as a follow-up.
# NO --force  => skip_completed resumes targets that already have result.csv
"/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine/benchmark_datasets" --benchmark astex --cache "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse" --output "/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/full" \
  --mode autonomous --threads 1 --omp-threads 1 \
  --ga-population 1000 --ga-generations 2000 \
  --job-timeout-seconds 10800
