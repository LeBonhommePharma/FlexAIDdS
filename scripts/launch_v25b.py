#!/usr/bin/env python3
# launch_v25b.py — Lever 2 scoring-fix mini-benchmark over the 35 v24 failures.
#
# v25b = v24 (PSHARE alpha=4, CF clustering, default matrix, oracle sites,
# Fix A/B selection) with ONE scoring change: the VCT exp(-r/r0) decay length
# is raised 4.0 -> 7.0 Å via FLEXAIDDS_VCT_R0.
#
# Rationale (Step 1-4 audit): r0=4.0 (P9) crushed CF.com to the noise floor
# (1LPZ com -1.3, total CF -3.8). The orientation-independent SAS baseline
# (~-104) then dominates the total CF, the GA fitness landscape goes flat, and
# the population collapses to the seed clone (1LPZ: 8447 clashed, only 2 unique
# survive). A gentler decay restores the com gradient relative to the SAS floor
# without the 58x compression. r0=7.0 -> exp(-3.5/7)=0.61 vs exp(-3.5/4)=0.42
# for a 3.5 A contact; distal (~8 A) contacts keep 0.32 vs 0.13.
#
# The change is RUNNER-ONLY: the engine already reads vct_dist_weight_r0 from
# dock_config.json (vcfunction.cpp:388), so docking physics is otherwise the
# v24 engine. Output is a SEPARATE dir; no competing benchmark may run (shared
# grid cache). Same double-fork + setsid + caffeinate detachment as launch_v24.
import os, sys, signal, subprocess, hashlib, shutil

REPO   = "/Users/lp.more/Projects/FlexAIDdS"
BUILD  = f"{REPO}/build"
BINARY = f"{BUILD}/FlexAIDdS"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT = os.path.expanduser("~/flexaidds_results/v25b_scoringfix")
CODES_FILE = f"{OUTPUT}/v25b_codes_35.txt"
PIDFILE = f"{OUTPUT}/v25b.pid"
STDOUT_LOG = f"{OUTPUT}/v25b_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v25b_stderr.log"

VCT_R0 = os.environ.get("FLEXAIDDS_VCT_R0", "7.0")

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

for p in (BINARY, RUNNER, ORACLE_DIR, CODES_FILE):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing {p}")
print(f"engine SHA: {sha256(BINARY)}")
print(f"runner SHA: {sha256(RUNNER)}")
print(f"FLEXAIDDS_VCT_R0 = {VCT_R0}")

# never run a competing benchmark (corrupts the shared grid cache)
ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    sys.exit(f"ERROR: a benchmark is already running (pids {ps.stdout.split()}) — abort")

env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY": BINARY,
    "FLEXAIDDS_BUILD": BUILD,
    "FLEXAIDDS_REPO": REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
    "FLEXAIDDS_VCT_R0": VCT_R0,           # Lever 2: gentler exp(-r/r0) decay
    "OMP_WAIT_POLICY": "passive",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "spread",
})
# clear unrelated opt-in diagnostics so v25b differs from v24 ONLY by r0
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON", "FLEXAIDDS_VCT_NORM"):
    env.pop(k, None)

cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark", "astex",
    "--only-codes", CODES_FILE,
    "--output", OUTPUT,
    "--threads", "5",
    "--omp-threads", "2",
    "--temperature", "298",
    "--job-timeout-seconds", "5400",
]

def double_fork_launch():
    if os.fork() > 0:
        return
    os.setsid()
    if os.fork() > 0:
        os._exit(0)
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    os.chdir(REPO)
    with open(STDOUT_LOG, "w") as out, open(STDERR_LOG, "w") as err:
        p = subprocess.Popen(cmd, stdout=out, stderr=err, env=env,
                              start_new_session=True)
    with open(PIDFILE, "w") as pf:
        pf.write(str(p.pid) + "\n")
    os._exit(0)

if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)
    double_fork_launch()
    print(f"launched v25b (r0={VCT_R0}) -> {OUTPUT}")
