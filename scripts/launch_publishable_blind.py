#!/usr/bin/env python3
"""
Publishable blind-docking benchmark — autonomous mode.

Protocol:
  - Mode:    autonomous  (pose_seed_enabled=false, seed_fraction=0.0, blinding=ON)
  - Site:    oracle (FLEXAIDDS_ORACLE_SITE_DIR = crystal binding-site PDB)
  - Binary:  FLEXAIDDS_BINARY (or local pin under FLEXAIDDS_LOCAL_ROOT)
  - Dataset: astex (Astex Diverse; historically used for cross-dock 85 claims)
  - Restarts: 5
  - Threads:  2 concurrent targets × 5 OMP each

Comparison baseline (historical notes from v123 session):
  v50b autonomous : 69/85 = 81.2%  (same mode, older binary without recent fixes)
  v127 oracle-ceiling : 78/85 = 91.8%  (NOT publishable — 90% native seeds)
  rDock (literature) : 88.2%

This run establishes the honest publishable rate for a fixed binary.
Paths resolve from env vars / repo root — no machine-absolute hardcoding.
"""
from __future__ import annotations

import datetime
import os
import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    # scripts/ → repo root
    return Path(__file__).resolve().parents[1]


REPO = _repo_root()
LOCAL_ROOT = Path(
    os.environ.get("FLEXAIDDS_LOCAL_ROOT", str(Path.home() / "flexaidds_results"))
).expanduser()

# Prefer explicit binary; fall back to local-first pin used by claim staging.
BIN = os.environ.get("FLEXAIDDS_BINARY") or str(
    LOCAL_ROOT / "three_engine_entropy_q1" / "bin" / "C"
)
DATA_DIR = os.environ.get(
    "FLEXAIDDS_DATA_DIR",
    str(REPO / "build"),
)
DS = os.environ.get(
    "FLEXAIDDS_BENCHMARK_DATASETS",
    str(Path(DATA_DIR) / "benchmark_datasets"),
)

TAG = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUT = os.environ.get(
    "FLEXAIDDS_PUBLISHABLE_OUT",
    str(LOCAL_ROOT / "campaigns" / f"publishable_{TAG}_blind"),
)

env = os.environ.copy()
env["FLEXAIDDS_BINARY"] = BIN
env["FLEXAIDDS_DATA_DIR"] = DATA_DIR
env.setdefault(
    "FLEXAIDDS_ORACLE_SITE_DIR",
    str(REPO / "benchmarks" / "astex_diverse" / "astex_diverse"),
)
env.setdefault("FLEXAIDDS_RESTARTS", "5")
env.setdefault("FLEXAIDDS_PARALLEL_RESTARTS", "1")
env.setdefault("FLEXAIDDS_EVAL_SCALE_DIHEDRAL", "1")
env.setdefault("FLEXAIDDS_SEED_ELITISM", "1")
env.setdefault("FLEXAIDDS_N_ELITE", "1")
env.setdefault("FLEXAIDDS_BUDGET_SCALE", "1")
env.setdefault("FLEXAIDDS_SOFTCORE_WAL", "1")
env.setdefault("FLEXAIDDS_SOFTCORE_FLOOR", "0.5")
env.setdefault("FLEXAIDDS_T_HOT", "500")
env.setdefault("FLEXAIDDS_VCT_R0", "4")
env.setdefault("FLEXAIDDS_RECEPTOR_ROTAMER_PREP", "0")
env.setdefault("FLEXAIDDS_ALLOW_CONCURRENT", "1")
env.setdefault("FLEXAIDDS_CONSENSUS_SCORER", "1")
# NO FLEXAIDDS_NATIVE_SEED_FRAC — autonomous mode sets seed_fraction=0.0

cmd = [
    DS,
    "--benchmark",
    "astex",
    "--mode",
    "autonomous",  # blind poses — publishable
    "--restarts",
    "5",
    "--threads",
    "2",
    "--omp-threads",
    "5",
    "--output",
    OUT,
]

if not Path(DS).is_file() and not Path(DS).exists():
    print(f"[publishable_blind] ERROR: benchmark runner not found: {DS}", file=sys.stderr)
    print(
        "  Set FLEXAIDDS_BENCHMARK_DATASETS or build with ENABLE_BENCHMARK_DATASETS.",
        file=sys.stderr,
    )
    sys.exit(2)

os.makedirs(OUT, exist_ok=True)
log_path = Path(OUT) / "stdout.log"
log = open(log_path, "w", buffering=1)

pid = os.fork()
if pid != 0:
    print(f"[publishable_blind] Launched PID {pid}")
    print(f"[publishable_blind] Binary: {BIN}")
    print("[publishable_blind] Mode: autonomous (blind poses + oracle site)")
    print(f"[publishable_blind] Output: {OUT}")
    sys.exit(0)

os.setsid()
pid2 = os.fork()
if pid2 != 0:
    sys.exit(0)

os.chdir(REPO)
proc = subprocess.Popen(
    cmd,
    env=env,
    stdout=log,
    stderr=log,
    start_new_session=True,
)
log.write(f"[publishable_blind] Worker PID {proc.pid}\n")
log.flush()
proc.wait()
log.close()
