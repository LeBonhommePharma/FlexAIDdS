#!/usr/bin/env python3
# launch_v129_sulfo_smoke.py — sulfonamide SybylTyper fix smoke (1HNN + 1T9B)
#
# Binary: HEAD build_lto @ 04ff1735 (–SO2NH– → N.am/S.O live VCT rows).
# Protocol: identical to v127 (r0=4, rotamer off, consensus ON, oracle-ceiling).
#
# Compare against:
#   v124 @ 15b536f8 — v124_full85_20260626_0413_consensus_guard
#   v127 @ ba5364d3 — v127_20260629_0139_optB_full85
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
BINARY = "/tmp/FlexAIDdS_v129_sulfo"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_smoke_2_sulfo.json"
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")

SMOKE_TARGETS = "1HNN,1T9B"
REF_V124_DIR = "v124_full85_20260626_0413_consensus_guard"
REF_V127_DIR = "v127_20260629_0139_optB_full85"

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
    for label, run_dir in [
        ("v124", REF_V124_DIR),
        ("v127", REF_V127_DIR),
    ]:
        p = RESULTS_DIR / run_dir / pdb / "result.csv"
        if p.exists():
            r = list(csv.DictReader(open(p)))[0]
            snap[label] = {
                "success": r.get("success"),
                "rmsd_hungarian": r.get("rmsd_hungarian"),
                "best_cluster_rmsd": r.get("best_cluster_rmsd"),
                "cf_native": r.get("cf_native"),
                "shannon_entropy": r.get("shannon_entropy"),
                "pose_source": r.get("pose_source"),
            }
    return snap


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-build-check", action="store_true")
    args = parser.parse_args()

    for p in (BINARY_SRC, RUNNER, ORACLE_DIR, JSON_PAIRS, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 2

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"v129_{tag}_smoke2_sulfo")
    cache = str(RESULTS_DIR / "cache_v129_sulfo_smoke2")
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":                BINARY,
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
        "--benchmark",           f"crossdock_json:{JSON_PAIRS}",
        "--output",              output,
        "--threads",             bench_threads,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               cache,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print("\nLaunching v129 sulfonamide smoke — 1HNN + 1T9B @ HEAD build_lto")
    print(f"  commit   : {git_commit}")
    print(f"  targets  : {SMOKE_TARGETS}")
    print(f"  output   : {output}")
    print(f"  threads  : {bench_threads} (concurrent OK)")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    ref = {pdb: ref_snapshot(pdb) for pdb in SMOKE_TARGETS.split(",")}

    prov_doc = {
        "version":       "v129_smoke2_sulfo",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    git_commit,
        "description": (
            "Sulfonamide SybylTyper fix smoke: 1HNN, 1T9B with HEAD binary "
            "(04ff1735 –SO2NH– → N.am/S.O live VCT rows). "
            "v127 protocol (r0=4, rotamer off, consensus, oracle-ceiling). "
            "Success criterion: top-1 RMSD < 2 Å (BCR) vs v127 failures "
            "(1HNN 11.06 Å, 1T9B 11.02 Å)."
        ),
        "smoke_targets": SMOKE_TARGETS,
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "output_dir":     output,
        "cache_dir":      cache,
        "pid":            child_pid,
        "reference_runs": {
            "v124_dir": REF_V124_DIR,
            "v127_dir": REF_V127_DIR,
            "per_target": ref,
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nv129 sulfonamide smoke launched pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


if __name__ == "__main__":
    main()