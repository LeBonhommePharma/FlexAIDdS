#!/usr/bin/env python3
# launch_v127_full85.py — HEAD full-85 after v126 logsumexp + H-bond/VCT scoring
#
# Queued by queue_after_v124_v126.py when v124 resume and v126 complete.
# Protocol matches v126 (Option B logsumexp, consensus ON, oracle-ceiling, 5-restart)
# but uses current build_lto binary (ba5364d3 H-bond cosine gate + C.ar stacking).
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import os
import sys
import subprocess
import hashlib
import json
import datetime
import shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO        = "/Users/lp.more/Projects/FlexAIDdS"
BUILD       = f"{REPO}/build_lto"
BINARY_SRC  = f"{BUILD}/FlexAIDdS"
BINARY      = "/tmp/FlexAIDdS_v127"
RUNNER      = f"{BUILD}/benchmark_datasets"
DATA_DIR    = BUILD
ORACLE_DIR  = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS  = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85.json"
RESULTS_DIR = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"

TAG    = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT = f"{RESULTS_DIR}/v127_{TAG}_optB_full85"
CACHE  = f"{RESULTS_DIR}/cache_v127_optB"
PROV   = f"{OUTPUT}/launch_provenance.json"

BENCH_THREADS = os.environ.get("FLEXAIDDS_BENCH_THREADS", "4")

# Keys recorded in launch_provenance.json for thesis / reproducibility audit.
# Parent benchmark_datasets inherits this env for all 85 targets (not canary-only).
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
    for p in (BINARY_SRC, RUNNER, ORACLE_DIR, JSON_PAIRS,
              f"{DATA_DIR}/MC_st0r5.2_6.dat"):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    git_commit = subprocess.check_output(
        ["git", "log", "--oneline", "-1"], cwd=REPO, text=True
    ).strip().split()[0]

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 85

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
        "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX", "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
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

    print("\nLaunching v127 full-85 — logsumexp + H-bond/VCT scoring (HEAD)")
    print(f"  output   : {OUTPUT}")
    print(f"  cache    : {CACHE}")
    print(f"  threads  : {BENCH_THREADS}")
    print(f"  consensus: ON")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v127_optB_full85",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    git_commit,
        "description": (
            "Full Astex 85 queued after v124 resume + v126 complete. "
            "Option B logsumexp composite, consensus ON, oracle-ceiling, "
            "5-restart SMFREE. HEAD binary includes H-bond cosine gate + "
            "C.ar stacking discriminator (ba5364d3). "
            "First-principles free-energy selector benchmark (log Z - alpha*H); "
            "not a v50b degenerate-fitness reproduction."
        ),
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "output_dir":     OUTPUT,
        "cache_dir":      CACHE,
        "pid":            child_pid,
        "benchmark": {
            "threads":             BENCH_THREADS,
            "temperature_k":       298,
            "job_timeout_seconds": 7200,
            "mode":                "oracle-ceiling",
            "n_pairs":             85,
        },
        "audit_notes": {
            "code_default_FLEXAIDDS_CONSENSUS_SCORER": "0 (since ce8f3368 v125 Option A)",
            "launch_override_FLEXAIDDS_CONSENSUS_SCORER": "1",
            "consensus_applies_to": "all 85 targets via parent benchmark_datasets env",
            "verify_post_run": "grep -c '\\[CONSENSUS\\]' stderr.log should equal completed target count",
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(PROV, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\nv127 launched (commit={git_commit}) pid={child_pid}")
    print(f"  prov: {PROV}")
    return child_pid


if __name__ == "__main__":
    main()