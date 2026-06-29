#!/usr/bin/env python3
"""Resume a stalled v111_science campaign (skip-completed targets)."""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BUILD = os.path.join(REPO, "build")
OUTPUT = os.path.expanduser("~/flexaidds_results/v111_science_20260626_0613")
RUNNER = os.path.join(BUILD, "benchmark_datasets")
JSON_PAIRS = os.path.join(REPO, "benchmarks/datasets/benchmark_astex_native_85.json")
MATRIX_V2 = os.path.join(REPO, "MC_st0r5.2_6_v2_science.dat")
ORACLE_DIR = os.path.join(REPO, "benchmarks/astex_diverse/astex_diverse")
DATA_DIR = BUILD if os.path.isfile(os.path.join(BUILD, "MC_st0r5.2_6.dat")) else REPO

env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY": os.path.join(BUILD, "FlexAIDdS"),
    "FLEXAIDDS_BUILD": BUILD,
    "FLEXAIDDS_REPO": REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
    "FLEXAIDDS_DATA_DIR": DATA_DIR,
    "FLEXAIDDS_RESTARTS": "7",
    "FLEXAIDDS_PARALLEL_RESTARTS": "1",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
    "FLEXAIDDS_CONSENSUS_SCORER": "1",
    "FLEXAIDDS_SEED_ELITISM": "1",
    "FLEXAIDDS_N_ELITE": "1",
    "FLEXAIDDS_BUDGET_SCALE": "1",
    "FLEXAIDDS_SOFTCORE_WAL": "1",
    "FLEXAIDDS_SOFTCORE_FLOOR": "0.5",
    "FLEXAIDDS_NATIVE_SEED_FRAC": "0.90",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
    "FLEXAIDDS_SCIENCE_FIXES": "1",
    "FLEXAIDDS_ENERGY_MATRIX": MATRIX_V2,
    "FLEXAIDDS_NEARMISS_SHARPEN": "1",
    "FLEXAIDDS_HBOND_ANGLE_GATE": "1",
    "FLEXAIDDS_T_HOT": "350",
    "FLEXAIDDS_SHARING_ALPHA": "6",
    "OMP_WAIT_POLICY": "passive",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "spread",
})

cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark", f"crossdock_json:{JSON_PAIRS}",
    "--mode", "oracle-ceiling",
    "--output", OUTPUT,
    "--threads", "4",
    "--temperature", "298",
    "--job-timeout-seconds", "1800",
]

resume_log = os.path.join(OUTPUT, "stdout_resume.log")
stderr_log = os.path.join(OUTPUT, "stderr_resume.log")
pid = launch_session_isolated(cmd, env, OUTPUT, stdout_log=resume_log, stderr_log=stderr_log, cwd=REPO)
print(f"Resumed v111_science → pid={pid}")
print(f"  output: {OUTPUT}")
print(f"  logs:   {resume_log}")