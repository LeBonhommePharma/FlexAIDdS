#!/usr/bin/env python3
# launch_v130_full85.py — sulfonamide remap + 1HNN expB oracle site (full Astex 85)
#
# Bundles v129 (04ff1735 SybylTyper –SO2NH– → N.am/S.O) with v129b expB fix:
# 1HNN uses ligand-centered oracle site (SKF COM) via benchmark_astex_native_85_v130.json.
# Protocol: v127 (Option B logsumexp, consensus ON, oracle-ceiling, r0=4, rotamer off).
#
# Success gate vs v127 (78/85): expect +1 from 1HNN; 1T9B may flip if elitism protects
# correctly under sulfo remap (v129 smoke: ini_elitism, GA collapsed).
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = "/Users/lp.more/Projects/FlexAIDdS"
BUILD = f"{REPO}/build_lto"
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY = "/tmp/FlexAIDdS_v130"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85_v130.json"
HNN_LIGAND_SITE = (
    f"{ORACLE_DIR}/1HNN/1HNN_ligand_centered_site.pdb"
)
RESULTS_DIR = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"

TAG = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT = f"{RESULTS_DIR}/v130_{TAG}_sulfo_expB_full85"
CACHE = f"{RESULTS_DIR}/cache_v130_sulfo_expB"
PROV = f"{OUTPUT}/launch_provenance.json"

REF_V127_DIR = "v127_20260629_0139_optB_full85"
REF_V129_DIR = "v129_20260629_0358_smoke2_sulfo"
REF_V129B_DIR = "v129b_20260629_0447_1hnn_expB_ligand_site"
MIN_COMMIT = "04ff1735"

BENCH_THREADS = os.environ.get("FLEXAIDDS_BENCH_THREADS", "4")

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


def main():
    for p in (
        BINARY_SRC,
        RUNNER,
        ORACLE_DIR,
        JSON_PAIRS,
        HNN_LIGAND_SITE,
        f"{DATA_DIR}/MC_st0r5.2_6.dat",
    ):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", MIN_COMMIT, "HEAD"],
        cwd=REPO,
    )
    if ancestor.returncode != 0:
        sys.exit(
            f"ERROR: HEAD ({git_commit}) does not include sulfonamide fix {MIN_COMMIT}; "
            "rebuild build_lto first"
        )

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 85
    hnn = next(p for p in native["pairs"] if p["receptor_id"] == "1HNN")
    assert hnn["oracle_site_pdb"] == HNN_LIGAND_SITE

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
        "FLEXAIDDS_BENCH_CACHE":           CACHE,
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
        "--output",              OUTPUT,
        "--threads",             BENCH_THREADS,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               CACHE,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(OUTPUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)

    print("\nLaunching v130 full-85 — sulfonamide remap + 1HNN expB oracle site")
    print(f"  commit   : {git_commit}")
    print(f"  output   : {OUTPUT}")
    print(f"  cache    : {CACHE}")
    print(f"  json     : {JSON_PAIRS}")
    print(f"  1HNN site: {HNN_LIGAND_SITE}")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v130_sulfo_expB_full85",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    git_commit,
        "description": (
            "Full Astex 85: v127 protocol + 04ff1735 sulfonamide SybylTyper remap + "
            "1HNN ligand-centered oracle site (expB). Validates bundled v129/v129b "
            "fixes on all 85 targets. Compare vs v127 (78/85)."
        ),
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "oracle_site_overrides": {
            "1HNN": {
                "site": HNN_LIGAND_SITE,
                "replaces": f"{ORACLE_DIR}/1HNN/1HNN_binding_site.pdb",
                "reason": "SAH-offset Get_Cleft site ~8.4 A from SKF COM",
            }
        },
        "output_dir":     OUTPUT,
        "cache_dir":      CACHE,
        "pid":            child_pid,
        "reference_runs": {
            "v127_dir": REF_V127_DIR,
            "v129_smoke": REF_V129_DIR,
            "v129b_1hnn": REF_V129B_DIR,
        },
        "benchmark": {
            "threads":             BENCH_THREADS,
            "temperature_k":       298,
            "job_timeout_seconds": 7200,
            "mode":                "oracle-ceiling",
            "n_pairs":             85,
        },
        "audit_notes": {
            "v129_1T9B": (
                "Smoke pass via ini_elitism (seed_echo=1); GA collapsed (positive CF). "
                "Valid oracle-ceiling count; not autonomous search win."
            ),
            "v129_1HNN": "Sulfo remap alone FAIL 11.53 A; expB site required.",
            "expected_delta_vs_v127": "+1 (1HNN) to +2 (+1T9B if elitism holds)",
            "count_metric": "per-target result.csv; report pose_source + seed_echo",
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(PROV, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\nv130 launched (commit={git_commit}) pid={child_pid}")
    print(f"  prov: {PROV}")
    return child_pid


if __name__ == "__main__":
    main()