#!/usr/bin/env python3
# launch_v46.py — daemonized launcher for the v46 native self-docking benchmark.
#
# v46 = v43 config (native Astex Diverse 85, oracle mode) + current binary
#       (Fixes A+B+C committed + log1p(N) cluster bonus uncommitted in build).
#
# Differences from v43:
#   FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0   (explicit; v43 relied on default=off,
#                                        v44 set it ON and hurt results → always
#                                        be explicit from now on)
#   run_name  = v46_20260613_native_logpop
#   binary    = build_lto/benchmark_datasets (SHA 12d6d9bd…, same as v45 build)
#
# All other params IDENTICAL to v43:
#   --benchmark astex          (native self-docking, 85 targets)
#   RESTARTS=3, N_ELITE=1, SEED_ELITISM=1, BUDGET_SCALE=1
#   SOFTCORE_WAL=1, SOFTCORE_FLOOR=0.5, T_HOT=500
#   NATIVE_SEED_FRAC=0.90, PRIORITY_TARGETS=1Q4G
#   --threads 10, --temperature 298, --job-timeout-seconds 5400
#
# Runs in parallel with v45 cross-docking benchmark (PID 7484, 10 workers).
# 11 logical cores total; OMP_WAIT_POLICY=passive → OS time-shares gracefully.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import datetime
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

# ── Paths ─────────────────────────────────────────────────────────────────────
REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY     = f"{BUILD}/FlexAIDdS"
RUNNER     = f"{BUILD}/benchmark_datasets"
DATA_DIR   = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT     = os.path.expanduser("~/flexaidds_results/v46_20260613_native_logpop")
PROV_FILE  = f"{OUTPUT}/launch_provenance.json"

# ── Provenance anchors (v45 build binary, current HEAD 9b3fec9) ───────────────
EXP_ENGINE_SHA = "dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA = "12d6d9bd19470eccc7e33dd828e13e8a9db50c1516a4ac28153a1b38ce1bd7f7"
EXP_MATRIX_MD5 = "72d7c7396702331d96ff12d18f831796"
GIT_COMMIT     = "9b3fec9"


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


# ── Pre-flight checks ─────────────────────────────────────────────────────────
for p in (BINARY, RUNNER, ORACLE_DIR, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing required path: {p}")

engine_sha = sha256(BINARY)
runner_sha = sha256(RUNNER)
matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

print(f"engine SHA256     : {engine_sha}")
print(f"runner SHA256     : {runner_sha}")
print(f"matrix MD5        : {matrix_md5}")

if engine_sha != EXP_ENGINE_SHA:
    sys.exit(f"ERROR: engine SHA mismatch\n  got  {engine_sha}\n  want {EXP_ENGINE_SHA}")
if runner_sha != EXP_RUNNER_SHA:
    sys.exit(f"ERROR: runner SHA mismatch\n  got  {runner_sha}\n  want {EXP_RUNNER_SHA}")
if matrix_md5 != EXP_MATRIX_MD5:
    sys.exit(f"ERROR: matrix MD5 mismatch\n  got  {matrix_md5}\n  want {EXP_MATRIX_MD5}")

# NOTE: v45 cross-docking benchmark is intentionally running in parallel (PID 7484).
# We do NOT abort on finding a running benchmark_datasets — parallel execution
# of native (v46) + cross-docking (v45) is the explicit goal.

# ── Environment ───────────────────────────────────────────────────────────────
env = dict(os.environ)
env.update({
    # Identical to v43 ──────────────────────────────────────────────────────
    "FLEXAIDDS_BINARY":                BINARY,
    "FLEXAIDDS_BUILD":                 BUILD,
    "FLEXAIDDS_REPO":                  REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
    "FLEXAIDDS_RESTARTS":              "3",
    "FLEXAIDDS_SEED_ELITISM":          "1",
    "FLEXAIDDS_N_ELITE":               "1",
    "FLEXAIDDS_BUDGET_SCALE":          "1",
    "FLEXAIDDS_SOFTCORE_WAL":          "1",
    "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
    "FLEXAIDDS_T_HOT":                 "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
    "FLEXAIDDS_DATA_DIR":              DATA_DIR,
    "FLEXAIDDS_PRIORITY_TARGETS":      "1Q4G",   # near-miss from v43
    # v46 addition: explicit OFF (v43 defaulted off, v44 set ON → hurt results)
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
    # OMP tuning (same as v45) ──────────────────────────────────────────────
    "OMP_WAIT_POLICY":                 "passive",
    "OMP_PLACES":                      "cores",
    "OMP_PROC_BIND":                   "spread",
})

# Clear ablation/tuning vars that might bleed in from the calling shell
for k in (
    "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
):
    env.pop(k, None)

# ── Command ───────────────────────────────────────────────────────────────────
# 10 threads: same as v43/v45. OS time-shares 11 logical cores across both
# benchmarks; OMP_WAIT_POLICY=passive prevents spin-wait CPU burn when idle.
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark",           "astex",   # native self-docking (same as v43)
    "--output",              OUTPUT,
    "--threads",             "10",
    "--temperature",         "298",
    "--job-timeout-seconds", "5400",
]

# ── Launch (Fix B: double-fork daemon, start_new_session=True) ────────────────
if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)

    print(f"\nLaunching v46 native Astex Diverse benchmark …")
    print(f"  benchmark:  astex (native self-docking, 85 targets)")
    print(f"  oracle dir: {ORACLE_DIR}")
    print(f"  output:     {OUTPUT}")
    print(f"  parallel with v45 cross-docking (PID 7484)")

    child_pid = launch_session_isolated(
        cmd,
        env,
        OUTPUT,
        cwd=REPO,
    )

    # ── Write provenance JSON ─────────────────────────────────────────────────
    prov = {
        "run_name":                   "v46_20260613_native_logpop",
        "version":                    "v46",
        "launched_at":                datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":                 GIT_COMMIT,
        "description": (
            "Native self-docking Astex Diverse 85 (oracle mode). "
            "Replicates v43 config exactly with FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0 "
            "made explicit (v44 set it ON and degraded results). "
            "Binary includes Fixes A+B+C (committed) + log1p(N) cluster bonus "
            "(uncommitted in native_score.cpp). "
            "Runs in parallel with v45 cross-docking benchmark."
        ),
        "oracle_mode":                True,
        "receptor_rotamer_prep":      False,
        "soft_wall_cutoff_angstrom":  0.40,
        "benchmark":                  "astex",
        "binary":                     BINARY,
        "binary_sha256":              engine_sha,
        "runner":                     RUNNER,
        "runner_sha256":              runner_sha,
        "matrix":                     f"{DATA_DIR}/MC_st0r5.2_6.dat",
        "matrix_md5":                 matrix_md5,
        "oracle_dir":                 ORACLE_DIR,
        "output_dir":                 OUTPUT,
        "pid":                        child_pid,
        "parallel_with": {
            "v45_pid":    7484,
            "v45_output": "~/flexaidds_results/v45_20260613_crossdock85",
        },
        "env_snapshot": {
            k: env[k] for k in (
                "FLEXAIDDS_ORACLE_SITE_DIR",
                "FLEXAIDDS_RESTARTS", "FLEXAIDDS_SEED_ELITISM",
                "FLEXAIDDS_N_ELITE", "FLEXAIDDS_BUDGET_SCALE",
                "FLEXAIDDS_SOFTCORE_WAL", "FLEXAIDDS_SOFTCORE_FLOOR",
                "FLEXAIDDS_T_HOT", "FLEXAIDDS_NATIVE_SEED_FRAC",
                "FLEXAIDDS_DATA_DIR",
                "FLEXAIDDS_PRIORITY_TARGETS",
                "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
            )
        },
    }
    with open(PROV_FILE, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\n✓ v46 native Astex Diverse launched (commit={GIT_COMMIT})")
    print(f"  pid:     {child_pid}  →  {OUTPUT}/benchmark.pid")
    print(f"  prov:    {PROV_FILE}")
    print(f"  monitor: tail -f {OUTPUT}/stdout.log")
    print(f"  errors:  tail -f {OUTPUT}/stderr.log")
