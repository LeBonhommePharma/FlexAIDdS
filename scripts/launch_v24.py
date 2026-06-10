#!/usr/bin/env python3
# launch_v24.py — daemonized launcher for the v24 "selection-fix" Astex 85 run.
#
# v24 = v23 (Arm A: PSHARE alpha=4, VCT exp(-r/4), 2000-gen P6 budget, CF
# clustering, default matrix, oracle binding sites) + two DatasetRunner-only
# post-processing fixes:
#   Fix A  success criterion uses the symmetry-corrected RMSD
#          min(rmsd_to_crystal, rmsd_hungarian) < 2.0  (was serial-only).
#   Fix B  frequency-gated cluster selection: drop degenerate (CF~0) poses,
#          prefer populated clusters (Frequency>1), pick min-CF within that pool
#          (was plain global min-CF over _0..19.pdb).
#
# The ENGINE (build/FlexAIDdS) is byte-identical to v23 — DatasetRunner.cpp links
# only into build/benchmark_datasets, so docking physics, grids and emitted poses
# are unchanged; only result selection/scoring/success differ. Engine SHA is
# therefore expected UNCHANGED; only the runner SHA moves.
#
# Same detachment contract as launch_v23.py / launch_v23_matrix_v1.py:
# double-fork + setsid, SIGHUP-immune, caffeinate -i, 5 workers x 2 OMP threads.
import os, sys, signal, subprocess, hashlib

REPO   = "/Users/lp.more/Projects/FlexAIDdS"
BUILD  = f"{REPO}/build"
BINARY = f"{BUILD}/FlexAIDdS"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT = os.path.expanduser("~/flexaidds_results/v24_20260610_selfix")
CODES_SRC  = os.path.expanduser(
    "~/flexaidds_results/v23_20260610_pshare_vct/v23_codes_85.txt")
CODES_FILE = f"{OUTPUT}/v24_codes_85.txt"
PIDFILE = f"{OUTPUT}/v24.pid"
STDOUT_LOG = f"{OUTPUT}/v24_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v24_stderr.log"

# Engine unchanged from v23; runner rebuilt with Fix A + Fix B.
EXP_ENGINE = "c1281359e43d5a6455d773d216784acb1c2a709debef8c212ec492926acbfb9f"
EXP_RUNNER = "5088c58b9cc29adab1f3ba0d6fa140c27455bcf6250244573b70207d05475784"

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ── pre-launch guards ────────────────────────────────────────────────
for p in (BINARY, RUNNER, ORACLE_DIR, CODES_SRC):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing {p}")
es, rs = sha256(BINARY), sha256(RUNNER)
print(f"engine SHA: {es}")
print(f"runner SHA: {rs}")
if es != EXP_ENGINE:
    sys.exit(f"ERROR: engine SHA mismatch (got {es}, want {EXP_ENGINE})")
if rs != EXP_RUNNER:
    sys.exit(f"ERROR: runner SHA mismatch (got {rs}, want {EXP_RUNNER}) "
             "— rebuild benchmark_datasets")
# never run a competing benchmark (corrupts the shared grid cache)
ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark astex"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    sys.exit(f"ERROR: a benchmark is already running (pids {ps.stdout.split()}) — abort")

env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY": BINARY,
    "FLEXAIDDS_BUILD": BUILD,
    "FLEXAIDDS_REPO": REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,   # enables oracle (cache != source tree)
    "OMP_WAIT_POLICY": "passive",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "spread",
})
# CF clustering + default matrix are the production defaults — clear all opt-in
# diagnostics so v24 differs from v23 Arm A ONLY by the runner fixes.
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON"):
    env.pop(k, None)

# 11 logical CPUs: 5 workers × 2 OMP threads = 10 ≤ 11.
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark", "astex",
    "--only-codes", CODES_FILE,        # 85 codes, full Astex Diverse (one code/line)
    "--output", OUTPUT,
    "--threads", "5",
    "--omp-threads", "2",
    "--temperature", "298",
    "--job-timeout-seconds", "5400",   # 90 min/job headroom for the 2000-gen P6 budget
]

def double_fork_launch():
    if os.fork() > 0:
        return
    os.setsid()                        # new session, detach controlling tty
    if os.fork() > 0:
        os._exit(0)
    signal.signal(signal.SIGHUP, signal.SIG_IGN)   # immune to terminal hangup
    os.chdir(REPO)
    with open(STDOUT_LOG, "w") as out, open(STDERR_LOG, "w") as err:
        p = subprocess.Popen(cmd, stdout=out, stderr=err, env=env,
                              start_new_session=True)
    with open(PIDFILE, "w") as pf:
        pf.write(str(p.pid) + "\n")
    os._exit(0)

if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)
    import shutil
    if not os.path.exists(CODES_FILE):
        shutil.copy(CODES_SRC, CODES_FILE)
    double_fork_launch()
    print(f"launched v24 -> {OUTPUT}")
