#!/usr/bin/env python3
# launch_v128_v50b_repro.py — pinned efc4f5d v50b reproduction (v128)
#
# Queued by queue_after_v127_v128.py after v127 completes.
# Replays the v50b consensus protocol exactly: SMFREE overflow selector,
# no oracle-ceiling CLI flag, matrix MD5 72d7c739, binary SHA dbfaca09…
#
# Gate: >=71/85 successes on per-target result.csv (v50b reference: 71/85).
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = "/Users/lp.more/Projects/FlexAIDdS"
WORKTREE = f"{REPO}/../FlexAIDdS_v128_repro"
BUILD = f"{WORKTREE}/build_lto"
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY = "/tmp/FlexAIDdS_v128"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85.json"
RESULTS_DIR = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"
BUILD_SCRIPT = f"{REPO}/scripts/build_v128_repro.sh"

GIT_COMMIT = "efc4f5d"
EXP_ENGINE_SHA = "dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA = "53fa471cfe3a55b2b071bf87e2181caba889ee92124553199df436275d714781"
EXP_MATRIX_MD5 = "72d7c7396702331d96ff12d18f831796"
SUCCESS_GATE = 71

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


def ensure_build(skip_build: bool) -> None:
    need = [BINARY_SRC, RUNNER, f"{DATA_DIR}/MC_st0r5.2_6.dat"]
    if all(os.path.exists(p) for p in need):
        engine_sha = sha256(BINARY_SRC)
        runner_sha = sha256(RUNNER)
        matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")
        if (
            engine_sha == EXP_ENGINE_SHA
            and runner_sha == EXP_RUNNER_SHA
            and matrix_md5 == EXP_MATRIX_MD5
        ):
            return
        if skip_build:
            sys.exit(
                "ERROR: v128 artifacts present but fingerprint mismatch "
                "(re-run without --skip-build)"
            )

    if skip_build:
        sys.exit(f"ERROR: missing v128 build artifacts under {BUILD}")

    print(f"Building v128 repro via {BUILD_SCRIPT}")
    subprocess.check_call(["bash", BUILD_SCRIPT], cwd=REPO)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="require pre-built worktree artifacts; do not invoke build script",
    )
    args = parser.parse_args()

    ensure_build(args.skip_build)

    for p in (BINARY_SRC, RUNNER, ORACLE_DIR, JSON_PAIRS, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    engine_sha = sha256(BINARY_SRC)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

    if engine_sha != EXP_ENGINE_SHA:
        sys.exit(f"ERROR: engine SHA mismatch\n  got  {engine_sha}\n  want {EXP_ENGINE_SHA}")
    if runner_sha != EXP_RUNNER_SHA:
        sys.exit(f"ERROR: runner SHA mismatch\n  got  {runner_sha}\n  want {EXP_RUNNER_SHA}")
    if matrix_md5 != EXP_MATRIX_MD5:
        sys.exit(f"ERROR: matrix MD5 mismatch\n  got  {matrix_md5}\n  want {EXP_MATRIX_MD5}")

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 85
    for pair in native["pairs"]:
        assert pair["receptor_id"] == pair["ligand_id"], (
            f"Non-native pair: {pair['receptor_id']} != {pair['ligand_id']}"
        )

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = f"{RESULTS_DIR}/v128_{tag}_v50b_repro"
    cache = f"{RESULTS_DIR}/cache_v128_v50b_repro"
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "4")

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

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
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           cache,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    })
    for k in (
        "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
        "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
        "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX", "FLEXAIDDS_PRIORITY_TARGETS",
        "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
    ):
        env.pop(k, None)

    cmd = [
        "caffeinate", "-i",
        RUNNER,
        "--benchmark",           f"crossdock_json:{JSON_PAIRS}",
        "--output",              output,
        "--threads",             bench_threads,
        "--temperature",         "298",
        "--job-timeout-seconds", "5400",
        "--cache",               cache,
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print("\nLaunching v128 v50b reproduction — pinned efc4f5d")
    print(f"  worktree : {WORKTREE}")
    print(f"  output   : {output}")
    print(f"  cache    : {cache}")
    print(f"  gate     : >={SUCCESS_GATE}/85 per-target result.csv")
    print(f"  mode     : v50b (no oracle-ceiling CLI)")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov_doc = {
        "version":       "v128_v50b_repro",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit":    GIT_COMMIT,
        "worktree":        WORKTREE,
        "description": (
            "Pinned efc4f5d v50b reproduction baseline. Identical env to "
            "launch_v50b.py: consensus ON, seed elitism, 5-restart SMFREE, "
            "RECEPTOR_ROTAMER_PREP=0, 298 K, no oracle-ceiling CLI. "
            "Isolates HEAD drift vs v50b 71/85 (83.5%). Not a thesis headline "
            "number — ablation anchor only."
        ),
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "output_dir":     output,
        "cache_dir":      cache,
        "pid":            child_pid,
        "success_gate":   {"min_success": SUCCESS_GATE, "total": 85},
        "benchmark": {
            "threads":             bench_threads,
            "temperature_k":       298,
            "job_timeout_seconds": 5400,
            "mode":                "v50b_autonomous_oracle_site",
            "oracle_ceiling_cli":  False,
            "n_pairs":             85,
        },
        "audit_notes": {
            "reference_run": "v50b_20260614_consensus5r",
            "reference_success": "71/85 (83.5%)",
            "matrix_head_md5": "9dc93717 (HEAD build_lto — not used here)",
            "count_metric": "per-target result.csv only; ignore stale aggregate CSV",
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nv128 launched (commit={GIT_COMMIT}) pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


if __name__ == "__main__":
    main()