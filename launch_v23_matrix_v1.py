#!/usr/bin/env python3
# launch_v23_matrix_v1.py — Arm B of the VCT matrix A/B experiment.
#
# Identical to launch_v23.py (same binary, oracle sites, P4..P10 pass, 2000-gen
# budget, CF clustering, temp 298) EXCEPT the scoring matrix: FLEXAIDDS_DATA_DIR
# points the engine at v23_matrix_v1_datadir/, whose MC_st0r5.2_6.dat IS the
# conservatively rescaled FA_matrix_v1.dat. AMINO.def / NUCLEOTIDES.def in that
# dir are byte-identical to canonical, so the matrix is the single variable.
#
# Arm A (default matrix) = the already-complete v23_20260610_pshare_vct run.
#
# Binary: build/FlexAIDdS  SHA256 c1281359e43d5a6455d773d216784acb1c2a709debef8c212ec492926acbfb9f
# Runner: build/benchmark_datasets SHA256 f1e5357257848aefd7f08a49d015f0881608568c80757d4f93b8285d6bca1cf9
import os, sys, signal, subprocess, hashlib

REPO   = "/Users/lp.more/Projects/FlexAIDdS"
BUILD  = f"{REPO}/build"
BINARY = f"{BUILD}/FlexAIDdS"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
DATA_DIR   = os.path.expanduser("~/flexaidds_results/v23_matrix_v1_datadir")
OUTPUT = os.path.expanduser("~/flexaidds_results/v23_matrix_v1_20260610")
CODES_SRC  = os.path.expanduser("~/flexaidds_results/v23_20260610_pshare_vct/v23_codes_85.txt")
CODES_FILE = f"{OUTPUT}/codes_85.txt"
PIDFILE = f"{OUTPUT}/run.pid"
STDOUT_LOG = f"{OUTPUT}/benchmark.log"
STDERR_LOG = f"{OUTPUT}/stderr.log"

EXP_ENGINE = "c1281359e43d5a6455d773d216784acb1c2a709debef8c212ec492926acbfb9f"
EXP_RUNNER = "f1e5357257848aefd7f08a49d015f0881608568c80757d4f93b8285d6bca1cf9"

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ── pre-launch guards ────────────────────────────────────────────────
for p in (BINARY, RUNNER, ORACLE_DIR, CODES_SRC,
          f"{DATA_DIR}/MC_st0r5.2_6.dat", f"{DATA_DIR}/AMINO.def"):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing {p}")
es, rs = sha256(BINARY), sha256(RUNNER)
print(f"engine SHA: {es}")
print(f"runner SHA: {rs}")
if es != EXP_ENGINE: sys.exit("ERROR: engine SHA mismatch")
if rs != EXP_RUNNER: sys.exit("ERROR: runner SHA mismatch")
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
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
    "FLEXAIDDS_DATA_DIR": DATA_DIR,            # <-- Arm B: swap in FA_matrix_v1
    "OMP_WAIT_POLICY": "passive",
    "OMP_PLACES": "cores",
    "OMP_PROC_BIND": "spread",
})
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON"):
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
    import shutil
    if not os.path.exists(CODES_FILE):
        shutil.copy(CODES_SRC, CODES_FILE)
    double_fork_launch()
    print(f"launched Arm B -> {OUTPUT}")
