#!/usr/bin/env python3
# launch_v34_ctrl.py — daemonized launcher for the v34 CONTROL run (old matrix).
#
# v34_ctrl = identical binary + identical env as v34 (softcore_floor compiled in,
#            SOFTCORE_WAL, T_HOT=500, RESTARTS=3, SEED_ELITISM=1, N_ELITE=1)
#            BUT using the ORIGINAL (pre-v34) MC_st0r5.2_6.dat via FLEXAIDDS_DATA_DIR.
#
# Purpose: apples-to-apples matrix comparison.  The ONLY difference from v34 is:
#   FLEXAIDDS_DATA_DIR=/tmp/data_ctrl_v34   (contains original matrix, backed up
#   from /tmp/MC_original_pre_v34.dat before any patching was applied)
#
# Entries that differ between v34 and v34_ctrl:
#   [10-13]  −86.64 (old)  vs  −15.0  (v34)
#   [10-35]    0.0  (old)  vs  −175.0 (v34)
#   [12-35]    0.0  (old)  vs  −110.0 (v34)
#   [13-40]  +33.99 (old)  vs  +90.0  (v34)
#   [14-40]  +43.24 (old)  vs  +90.0  (v34)
#   [15-40]  +29.56 (old)  vs  +90.0  (v34)
#
# Binary: bf6770a2 — same as v34 (softcore_floor always compiled in;
#   the floor only activates inside the SOFTCORE_WAL block which is ON in both)
#
import os, sys, signal, subprocess, hashlib

CTRL_DATA_DIR = "/tmp/data_ctrl_v34"

# Sanity: ctrl data dir must have the original matrix
_ctrl_mat = f"{CTRL_DATA_DIR}/MC_st0r5.2_6.dat"
if not os.path.exists(_ctrl_mat):
    sys.exit(
        f"ERROR: control data dir missing: {_ctrl_mat}\n"
        "       Run: cp /tmp/MC_original_pre_v34.dat /tmp/data_ctrl_v34/MC_st0r5.2_6.dat\n"
        "       and: cp build_lto/AMINO.def build_lto/NUCLEOTIDES.def /tmp/data_ctrl_v34/"
    )
import hashlib as _hlib
with open(_ctrl_mat, "rb") as _f:
    _ctrl_md5 = __import__('hashlib').md5(_f.read()).hexdigest()
EXPECTED_ORIG_MD5 = "204b75ef31b69e4a14deecf8a48c3f71"
if _ctrl_md5 != EXPECTED_ORIG_MD5:
    sys.exit(
        f"ERROR: ctrl matrix MD5 mismatch (got {_ctrl_md5}, want {EXPECTED_ORIG_MD5})\n"
        "       The control data dir must contain the PRE-PATCH matrix."
    )

REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY     = f"{BUILD}/FlexAIDdS"
RUNNER     = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT     = os.path.expanduser("~/flexaidds_results/v34_ctrl_20260612_oldmatrix")
CODES_CSV  = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"
CODES_FILE = f"{OUTPUT}/v34_ctrl_codes_85.txt"
PIDFILE    = f"{OUTPUT}/v34_ctrl.pid"
STDOUT_LOG = f"{OUTPUT}/v34_ctrl_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v34_ctrl_stderr.log"

EXP_ENGINE = "bf6770a2c14f8a8182fb73629de9782180334e8247935881079d03b67dde0eb7"
EXP_RUNNER = "3f059681dc5a779345040b3759345d2469624b68f6cb42961de51249ea98031c"

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

for p in (BINARY, RUNNER, ORACLE_DIR, CODES_CSV):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing {p}")

engine_sha = sha256(BINARY)
runner_sha = sha256(RUNNER)
print(f"engine SHA: {engine_sha}")
print(f"runner SHA: {runner_sha}")

if engine_sha != EXP_ENGINE:
    sys.exit(
        f"ERROR: engine SHA mismatch (got {engine_sha}, want {EXP_ENGINE})\n"
        "       Rebuild: cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build_lto --target FlexAIDdS -j8"
    )
if runner_sha != EXP_RUNNER:
    sys.exit(
        f"ERROR: runner SHA mismatch (got {runner_sha}, want {EXP_RUNNER})\n"
        "       Rebuild: cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build_lto --target benchmark_datasets -j8"
    )

ps = subprocess.run(["pgrep", "-x", "benchmark_datasets"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    running = ps.stdout.split()
    if len(running) >= 2:
        sys.exit(
            f"ERROR: 2+ benchmarks already running (pids {running}) — abort"
        )

def make_codes_file(csv_path, out_path):
    codes = []
    with open(csv_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            pdb_id = line.split(",")[0].strip()
            if pdb_id:
                codes.append(pdb_id)
    with open(out_path, "w") as f:
        f.write("\n".join(codes) + "\n")
    return len(codes)

env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":          BINARY,
    "FLEXAIDDS_BUILD":           BUILD,
    "FLEXAIDDS_REPO":            REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,
    "FLEXAIDDS_RESTARTS":        "3",
    "FLEXAIDDS_SEED_ELITISM":    "1",
    "FLEXAIDDS_N_ELITE":         "1",
    "FLEXAIDDS_BUDGET_SCALE":    "1",
    "FLEXAIDDS_SOFTCORE_WAL":    "1",
    "FLEXAIDDS_T_HOT":           "500",
    "FLEXAIDDS_DATA_DIR":        CTRL_DATA_DIR,   # ← ONLY difference from v34
    "OMP_WAIT_POLICY":           "passive",
    "OMP_PLACES":                "cores",
    "OMP_PROC_BIND":             "spread",
})
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
          "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
          "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
          "FLEXAIDDS_RING_FLEX"):
    env.pop(k, None)

cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark",           "astex",
    "--only-codes",          CODES_FILE,
    "--output",              OUTPUT,
    "--threads",             "5",
    "--omp-threads",         "2",
    "--temperature",         "298",
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
    n_codes = make_codes_file(CODES_CSV, CODES_FILE)
    print(f"codes: {n_codes} targets → {CODES_FILE}")
    double_fork_launch()
    print(f"launched v34_ctrl (OLD matrix via FLEXAIDDS_DATA_DIR={CTRL_DATA_DIR}) → {OUTPUT}")
    print(f"  monitor: tail -f {STDOUT_LOG}")
    print(f"  pid:     {PIDFILE}")
    print(f"  est. wall time: ~12.5 h")
