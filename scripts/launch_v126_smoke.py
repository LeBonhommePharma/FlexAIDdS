#!/usr/bin/env python3
# launch_v126_smoke.py — Option B (logsumexp fix): smoke 6 key targets
#
# Identical to v124 (consensus-guard, oracle-ceiling, 5-restart, SMFREE,
# r0=4, hbond_search ON, no rotamer prep, sas_weight=1.0) but with the
# logsumexp-stable boltzmann_composite (Option B fix).
# Consensus scorer is back ON (now real signal, not emergency fallback).
# Priority targets: the 6 canary targets from v124/v125 smoke.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import os, sys, subprocess, hashlib, json, datetime, shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO        = "/Users/lp.more/Projects/FlexAIDdS"
BUILD       = f"{REPO}/build_lto"
BINARY_SRC  = f"{BUILD}/FlexAIDdS"
BINARY      = "/tmp/FlexAIDdS_v126"
RUNNER      = f"{BUILD}/benchmark_datasets"
DATA_DIR    = BUILD
ORACLE_DIR  = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS  = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85.json"
RESULTS_DIR = "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results"

TAG    = datetime.datetime.now().strftime("%Y%m%d_%H%M")
OUTPUT = f"{RESULTS_DIR}/v126_{TAG}_optB_smoke"
CACHE  = f"{RESULTS_DIR}/cache_v126_optB"
PROV   = f"{OUTPUT}/launch_provenance.json"

BENCH_THREADS = os.environ.get("FLEXAIDDS_BENCH_THREADS", "4")
SMOKE_TARGETS = "1R55,1G9V,1OF6,1T46,1XOZ,1Y6R"

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

print(f"git commit        : {git_commit}")
print(f"engine SHA256     : {engine_sha}")
print(f"runner SHA256     : {runner_sha}")
print(f"matrix MD5        : {matrix_md5}")
print(f"binary stamped at : {BINARY}")

native = json.load(open(JSON_PAIRS))
assert len(native["pairs"]) == 85, "Expected 85 pairs in native JSON"
print("Native JSON verified: 85 self-docking pairs")

env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":                BINARY,
    "FLEXAIDDS_BUILD":                 BUILD,
    "FLEXAIDDS_REPO":                  REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
    "FLEXAIDDS_RESTARTS":              "5",
    "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
    "FLEXAIDDS_CONSENSUS_SCORER":      "1",   # ON — now real signal (Option B fixes overflow)
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
    "FLEXAIDDS_PRIORITY_TARGETS":      SMOKE_TARGETS,
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
]

if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)
    os.makedirs(CACHE, exist_ok=True)

    print("\nLaunching v126 smoke — Option B: logsumexp stable composite")
    print(f"  output   : {OUTPUT}")
    print(f"  cache    : {CACHE}")
    print(f"  targets  : {SMOKE_TARGETS}")
    print(f"  consensus: ON  (real signal, overflow fixed by logsumexp)")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":      "v126_optB_smoke",
        "launched_at":  datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":   git_commit,
        "description": (
            "Option B smoke: logsumexp-stable boltzmann_composite replaces "
            "exp(-CF/kT) overflow. Consensus scorer ON (now real thermodynamic "
            "signal, not emergency fallback). 6 canary targets: 1R55,1G9V,1OF6,1T46,1XOZ,1Y6R. "
            "Expected: all 6 success where v125 (consensus OFF) got -1.0 on 3."
        ),
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "smoke_targets":  SMOKE_TARGETS,
        "pid":            child_pid,
    }
    with open(PROV, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\nv126 smoke launched (commit={git_commit})")
    print(f"  pid    : {child_pid}")
    print(f"  prov   : {PROV}")
    print(f"  log    : tail -f {OUTPUT}/stdout.log")
    print(f"  errors : tail -f {OUTPUT}/stderr.log")
