#!/usr/bin/env python3
# launch_holo_smoke_single.py — single-target holo pipeline verification smoke
#
# Usage:
#   python3 scripts/launch_holo_smoke_single.py --target 1G9V
#   python3 scripts/launch_holo_smoke_single.py --target 1TW6
#
# Data fix @ bf8cf1d2. Protocol: v127 (r0=4, rotamer off, consensus, oracle-ceiling).
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = "/Users/lp.more/Projects/FlexAIDdS"
BUILD = f"{REPO}/build_lto"
BINARY_SRC = f"{BUILD}/FlexAIDdS"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
DATA_COMMIT = "bf8cf1d2"

TARGET_CONFIG = {
    "1G9V": {
        "json": f"{REPO}/benchmarks/datasets/benchmark_astex_smoke_1_1g9v_holo.json",
        "binary": "/tmp/FlexAIDdS_holo_smoke_1g9v",
        "cache": "cache_holo_smoke_1g9v",
        "required_paths": (
            f"{ORACLE_DIR}/1G9V/1G9V_apo.pdb",
            f"{ORACLE_DIR}/1G9V/1G9V_holo.pdb",
        ),
        "gate": "CF_native ~+50 (not +209); BCR < 2 Å",
    },
    "1TW6": {
        "json": f"{REPO}/benchmarks/datasets/benchmark_astex_smoke_1_1tw6_holo.json",
        "binary": "/tmp/FlexAIDdS_holo_smoke_1tw6",
        "cache": "cache_holo_smoke_1tw6",
        "required_paths": (f"{ORACLE_DIR}/1TW6/1TW6_holo.pdb",),
        "gate": "no intra-clash sentinel in stderr; RMSD improved vs v109/v127",
    },
}

REF_RUNS = (
    ("v109", "v109_20260626_tier1_consensus5r"),
    ("v127", "v127_20260629_0139_optB_full85"),
    ("v130", "v130_20260629_0548_sulfo_expB_full85"),
)

ENV_SNAPSHOT_KEYS = (
    "FLEXAIDDS_BINARY",
    "FLEXAIDDS_ORACLE_SITE_DIR",
    "FLEXAIDDS_RESTARTS",
    "FLEXAIDDS_PARALLEL_RESTARTS",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL",
    "FLEXAIDDS_CONSENSUS_SCORER",
    "FLEXAIDDS_SEED_ELITISM",
    "FLEXAIDDS_N_ELITE",
    "FLEXAIDDS_BUDGET_SCALE",
    "FLEXAIDDS_SOFTCORE_WAL",
    "FLEXAIDDS_SOFTCORE_FLOOR",
    "FLEXAIDDS_T_HOT",
    "FLEXAIDDS_NATIVE_SEED_FRAC",
    "FLEXAIDDS_VCT_R0",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
    "FLEXAIDDS_DATA_DIR",
    "FLEXAIDDS_BENCH_CACHE",
    "FLEXAIDDS_ALLOW_CONCURRENT",
)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def ref_snapshot(pdb: str) -> dict:
    import csv

    snap = {}
    for label, run_dir in REF_RUNS:
        p = RESULTS_DIR / run_dir / pdb / "result.csv"
        if p.exists():
            r = list(csv.DictReader(open(p)))[0]
            snap[label] = {
                "success": r.get("success"),
                "rmsd_hungarian": r.get("rmsd_hungarian"),
                "rmsd_to_crystal": r.get("rmsd_to_crystal"),
                "best_cluster_rmsd": r.get("best_cluster_rmsd"),
                "cf_native": r.get("cf_native"),
                "pose_source": r.get("pose_source"),
            }
    return snap


def launch(target: str) -> int:
    cfg = TARGET_CONFIG[target]
    json_pairs = cfg["json"]
    binary = cfg["binary"]

    for p in (BINARY_SRC, RUNNER, ORACLE_DIR, json_pairs, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")
    for p in cfg["required_paths"]:
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing holo data path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    native = json.load(open(json_pairs))
    assert len(native["pairs"]) == 1
    assert native["pairs"][0]["receptor_id"] == target

    shutil.copy2(BINARY_SRC, binary)
    os.chmod(binary, 0o755)

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"holo_smoke_{tag}_{target}")
    cache = str(RESULTS_DIR / cfg["cache"])
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "1")

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":                binary,
        "FLEXAIDDS_BUILD":                 BUILD,
        "FLEXAIDDS_REPO":                  REPO,
        "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
        "FLEXAIDDS_RESTARTS":              "5",
        "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
        "FLEXAIDDS_CONSENSUS_SCORER":      "1",
        "FLEXAIDDS_SEED_ELITISM":          "1",
        "FLEXAIDDS_N_ELITE":               "1",
        "FLEXAIDDS_BUDGET_SCALE":          "1",
        "FLEXAIDDS_SOFTCORE_WAL":          "1",
        "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
        "FLEXAIDDS_T_HOT":                 "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
        "FLEXAIDDS_VCT_R0":                "4",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_DATA_DIR":              DATA_DIR,
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           cache,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    })
    for k in (
        "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
        "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
        "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX",
        "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
        "FLEXAIDDS_PRIORITY_TARGETS",
    ):
        env.pop(k, None)

    cmd = [
        "caffeinate", "-i",
        RUNNER,
        "--benchmark",           f"crossdock_json:{json_pairs}",
        "--output",              output,
        "--threads",             bench_threads,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               cache,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print(f"\nLaunching holo smoke — {target} @ HEAD build_lto")
    print(f"  commit      : {git_commit}")
    print(f"  data_commit : {DATA_COMMIT}")
    print(f"  gate        : {cfg['gate']}")
    print(f"  output      : {output}")
    print(f"  threads     : {bench_threads}")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov_doc = {
        "version":       f"holo_smoke_{target}",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    git_commit,
        "data_commit":   DATA_COMMIT,
        "target":        target,
        "gate":          cfg["gate"],
        "binary":        binary,
        "binary_sha256": sha256(binary),
        "runner_sha256": sha256(RUNNER),
        "matrix_md5":    md5(f"{DATA_DIR}/MC_st0r5.2_6.dat"),
        "json_pairs":    json_pairs,
        "output_dir":    output,
        "cache_dir":     cache,
        "pid":           child_pid,
        "reference_runs": {
            f"{label}_dir": run_dir for label, run_dir in REF_RUNS
        },
        "reference_per_target": ref_snapshot(target),
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nHolo smoke {target} launched pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--target", required=True, choices=sorted(TARGET_CONFIG))
    args = parser.parse_args()
    launch(args.target)


if __name__ == "__main__":
    main()