#!/usr/bin/env python3
"""
launch_v111_science.py — Oracle-ceiling Astex 85 with v111 science-fix bundle.

Interventions (all independently disableable via env):
  FLEXAIDDS_SCIENCE_FIXES=1       matrix v2 + hbond recalibration
  FLEXAIDDS_ENERGY_MATRIX=...     MC_st0r5.2_6_v2_science.dat
  FLEXAIDDS_NEARMISS_SHARPEN=1    r0=4.5, entropy=0.15, T_HOT=350, sharing_alpha=6
  FLEXAIDDS_HBOND_ANGLE_GATE=1    120° hard cutoff

Post-run: failure_classify.py + cf_ground_truth_audit.py

Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = os.environ.get(
    "FLEXAIDDS_REPO",
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
)
def _pick_build_dir() -> str:
    for name in ("build", "build_lto"):
        d = os.path.join(REPO, name)
        if os.path.isfile(os.path.join(d, "benchmark_datasets")):
            return d
    return os.path.join(REPO, "build_lto")


BUILD = _pick_build_dir()
BINARY_SRC = os.path.join(BUILD, "FlexAIDdS")
RUNNER = os.path.join(BUILD, "benchmark_datasets")
DATA_DIR = BUILD if os.path.isfile(os.path.join(BUILD, "MC_st0r5.2_6.dat")) else REPO
ORACLE_DIR = os.path.join(REPO, "benchmarks/astex_diverse/astex_diverse")
JSON_PAIRS = os.path.join(REPO, "benchmarks/datasets/benchmark_astex_native_85.json")
MATRIX_V2 = os.path.join(REPO, "MC_st0r5.2_6_v2_science.dat")
TAG = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT = os.path.expanduser(f"~/flexaidds_results/v111_science_{TAG}")
PROV_FILE = os.path.join(OUTPUT, "provenance.json")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    for p in (BINARY_SRC, RUNNER, ORACLE_DIR, JSON_PAIRS, MATRIX_V2):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    os.makedirs(OUTPUT, exist_ok=True)

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY": BINARY_SRC,
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
        # v111 science bundle
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

    print(f"v111 science launch")
    print(f"  git:    {git_commit}")
    print(f"  matrix: {MATRIX_V2}")
    print(f"  output: {OUTPUT}")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version": "v111_science",
        "launched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_commit": git_commit,
        "description": (
            "Oracle-ceiling Astex 85 with v111 science fixes: matrix v2 Priority-1, "
            "hbond angular recalibration (weight=-3.5, sigma_angle=20, angle_gate=120), "
            "near-miss sharpen (r0=4.5, entropy=0.15, T_HOT=350, sharing_alpha=6)."
        ),
        "binary": BINARY_SRC,
        "binary_sha256": sha256(BINARY_SRC),
        "energy_matrix": MATRIX_V2,
        "output_dir": OUTPUT,
        "pid": child_pid,
        "env_flags": [
            "FLEXAIDDS_SCIENCE_FIXES",
            "FLEXAIDDS_ENERGY_MATRIX",
            "FLEXAIDDS_NEARMISS_SHARPEN",
            "FLEXAIDDS_HBOND_ANGLE_GATE",
            "FLEXAIDDS_T_HOT",
            "FLEXAIDDS_SHARING_ALPHA",
        ],
        "post_run": [
            f"python3 scripts/failure_classify.py {OUTPUT}",
            f"python3 scripts/cf_ground_truth_audit.py {OUTPUT}",
        ],
    }
    with open(PROV_FILE, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"  pid: {child_pid}")
    print(f"  prov: {PROV_FILE}")


if __name__ == "__main__":
    main()