#!/usr/bin/env python3
# launch_v129b_1hnn_expB.py — Experiment B: 1HNN with SKF/ligand-centered oracle site
#
# Fixes SAH-offset oracle centroid (8.4 A from SKF) that guided GA into ~11 A false
# minimum. Uses ligand-centered site (DatasetRunner write_ligand_centered_site logic).
# Binary: build_lto @ 04ff1735 (sulfonamide SybylTyper). No rebuild.
#
# Compare: v129_20260629_0358_smoke2_sulfo/1HNN (FAIL 11.31 A)
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

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
BINARY = "/tmp/FlexAIDdS_v129b_1hnn"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_smoke_1_hnn_expB.json"
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
REF_V129_DIR = "v129_20260629_0358_smoke2_sulfo"

GIT_COMMIT = "04ff1735"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    for p in (BINARY_SRC, RUNNER, JSON_PAIRS, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    native = json.load(open(JSON_PAIRS))
    site = native["pairs"][0]["oracle_site_pdb"]
    if not os.path.exists(site):
        sys.exit(f"ERROR: ligand-centered site missing: {site}")

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"v129b_{tag}_1hnn_expB_ligand_site")
    cache = str(RESULTS_DIR / "cache_v129b_1hnn_expB")
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "1")

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":                BINARY,
        "FLEXAIDDS_BUILD":                 BUILD,
        "FLEXAIDDS_REPO":                  REPO,
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
        "FLEXAIDDS_VCT_NORM", "FLEXAIDDS_ORACLE_SITE_DIR",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX", "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
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

    print("\nLaunching v129b Experiment B — 1HNN ligand-centered oracle site")
    print(f"  commit   : {GIT_COMMIT}")
    print(f"  site     : {site}")
    print(f"  output   : {output}")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov_doc = {
        "version":       "v129b_1hnn_expB_ligand_site",
        "experiment":    "B",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    GIT_COMMIT,
        "description": (
            "Experiment B: 1HNN with ligand-centered oracle site (SKF COM) instead of "
            "SAH-offset Get_Cleft site. Addresses ~8.4 A oracle/SKF mismatch that "
            "guided GA into SAH-channel false minimum (v129 FAIL 11.31 A). "
            "Sulfonamide typing @ 04ff1735 unchanged."
        ),
        "oracle_site_pdb": site,
        "oracle_site_sah_backup": (
            f"{REPO}/benchmarks/astex_diverse/astex_diverse/1HNN/"
            "1HNN_binding_site_SAH_backup.pdb"
        ),
        "binary":         BINARY,
        "binary_sha256":  sha256(BINARY),
        "json_pairs":     JSON_PAIRS,
        "output_dir":     output,
        "cache_dir":      cache,
        "pid":            child_pid,
        "reference_run":  REF_V129_DIR,
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nv129b Experiment B launched pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


if __name__ == "__main__":
    main()