#!/usr/bin/env python3
# launch_v23.py — daemonized launcher for the v23 full Astex 85 oracle benchmark.
#
# v23 = v22 + P4..P10 optimization pass:
#   P4  best-of-N oracle cluster RMSD recorded per target
#   P5  PSHARE niche sharing alpha=4 + periodic BOOM random injection (every 100 gens)
#   P6  GA generation budget 500->2000 base (x ceil(n_genes/4)) + ADAPTVGA
#   P8  dual SEC termination (energy + gene-space joint convergence)
#   P9  VCT distance-weighted contacts exp(-r/r0), r0=4.0 A
# (P7 fine-grid and P10 matrix-swap are opt-in diagnostics, NOT enabled here.)
#
# Binary: build/FlexAIDdS (Release + IPO/LTO, Metal OFF), git faede51.
# Double-fork + setsid detachment, SIGHUP-immune, caffeinate-wrapped, CF clustering,
# oracle binding sites via FLEXAIDDS_ORACLE_SITE_DIR.
import os, sys, signal, subprocess

REPO   = "/Users/lp.more/Projects/FlexAIDdS"
BUILD  = f"{REPO}/build"
BINARY = f"{BUILD}/FlexAIDdS"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT = os.path.expanduser("~/flexaidds_results/v23_20260610_pshare_vct")
CODES_FILE = f"{OUTPUT}/v23_codes_85.txt"
PIDFILE = f"{OUTPUT}/v23.pid"
STDOUT_LOG = f"{OUTPUT}/v23_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v23_stderr.log"

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
# CF clustering is the default — deliberately NOT setting FLEXAIDDS_USE_DP.
# Diagnostics OFF for the production run.
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON"):
    env.pop(k, None)

# 11 logical CPUs: 5 workers × 2 OMP threads = 10 ≤ 11.
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark", "astex",
    "--only-codes", CODES_FILE,        # 85 codes, full Astex Diverse (file: one code/line)
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
    # Stage the codes file inside the durable output dir (not /tmp).
    import shutil
    src = "/tmp/v19_codes_85.txt"
    if os.path.exists(src) and not os.path.exists(CODES_FILE):
        shutil.copy(src, CODES_FILE)
    double_fork_launch()
    print("parent returning")
