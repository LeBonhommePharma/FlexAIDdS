#!/usr/bin/env python3
"""
v110 — Multi-cleft GA: independent GA per octree subregion.
Autonomous mode (blinding=ON, seed_elitism=OFF, no oracle assistance).
FLEXAIDDS_MULTI_CLEFT=8 → --parallel-dock --parallel-dock-regions 8
Commit: d5aac417  Binary: fb41fa836f718ba6 (build_lto LTO)

Fixes in this commit:
  1. Global-index seeding zeroed for subregion GAs (heap corruption fix)
  2. chrom allocated before get_best_chromosome() (SIGSEGV fix)
  3. chrom_snapshot calloc'd for downstream pipeline (SIGSEGV fix)
  4. gene_lim allocated + set_gene_lim() for downstream pipeline (crash fix)

Comparison table targets:
  v50b  : 69/85 = 81.2% (single GA, oracle site, 5-restart consensus)
  v108  : 42/85 = 49.4% (oracle-ceiling mode, seed_elitism ON)
  v109  : TBD   (autonomous, seed_elitism OFF — killed by rate-limit)
  v110  : TBD   (multi-cleft, autonomous — this run)
"""
import os, sys, subprocess, datetime

REPO = os.path.expanduser("~/Projects/FlexAIDdS")
BIN  = os.path.join(REPO, "build_lto/FlexAIDdS")
DS   = os.path.join(REPO, "build_lto/benchmark_datasets")
TAG  = datetime.datetime.now().strftime("%Y%m%d")
OUT  = os.path.expanduser(
    f"~/Documents/PhD/Programs/FlexAIDdS/results/v110_{TAG}_multicleft_d5aac417"
)

env = os.environ.copy()
env["FLEXAIDDS_ORACLE_SITE_DIR"] = os.path.join(
    REPO, "benchmarks/astex_diverse/astex_diverse"
)
env["FLEXAIDDS_BINARY"]   = "build_lto/FlexAIDdS"
env["FLEXAIDDS_BUILD"]    = "build_lto"
env["FLEXAIDDS_MULTI_CLEFT"] = "8"   # → --parallel-dock --parallel-dock-regions 8

cmd = [
    DS,
    "--benchmark",   "astex",
    "--mode",        "autonomous",
    "--threads",     "1",          # 1 target at a time (multi-cleft uses OMP internally)
    "--omp-threads", "8",          # 8 OMP threads for region parallelism per target
    "--output",      OUT,
]

os.makedirs(OUT, exist_ok=True)
log = open(os.path.join(OUT, "benchmark_runner.log"), "w", buffering=1)

# os.setsid() double-fork — survives SIGHUP from terminal/session close
pid = os.fork()
if pid != 0:
    print(f"[v110] Launched PID {pid}")
    print(f"[v110] Output: {OUT}")
    print(f"[v110] FLEXAIDDS_MULTI_CLEFT=8 → 24 octree subregions per target")
    sys.exit(0)

os.setsid()
pid2 = os.fork()
if pid2 != 0:
    sys.exit(0)

# grandchild — detached
os.chdir(REPO)
proc = subprocess.Popen(
    cmd, env=env,
    stdout=log, stderr=log,
    start_new_session=True,
)
print(f"[v110] Worker PID {proc.pid}", flush=True)
proc.wait()
log.close()
