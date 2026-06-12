#!/usr/bin/env python3
# launch_v34.py — daemonized launcher for the v34 matrix-fix + softcore-floor Astex 85 run.
#
# v34 = v33-base (RESTARTS=3, N_ELITE=1, SEED_ELITISM=1, BUDGET_SCALE=1,
#                 SOFTCORE_WAL, T_HOT=500)
#       + SOFTCORE_FLOOR (default 0.5 × cr; hard r^-12 below that, parabola above)
#       + MATRIX PATCHES in MC_st0r5.2_6.dat:
#           [10-13]  −86.64 → −15.0   (N.AR vs backbone C=O: removed spurious H-bond)
#           [10-35]    0.0  → −175.0  (N.AR vs Zn: primary coordination)
#           [12-35]    0.0  → −110.0  (N.PL3 vs Zn: moderate coordination)
#           [13-40]  +33.99 → +90.0   (backbone C=O burial: penalty parity with N)
#           [14-40]  +43.24 → +90.0   (O.3 burial: penalty parity)
#           [15-40]  +29.56 → +90.0   (carboxylate O burial: penalty parity)
#       + metal_coord_enabled: already unconditional in DatasetRunner (no-op change)
#
# Hypothesis: [10-13] removal clears the dominant false-minimum driver in N-rich
#   pockets (His, Trp, aromatic ligands) — the parabola was letting buried N.AR
#   pairs score "free" below 50% cr. Zn coordination fixes rescue metalloprotein
#   targets. Burial penalty symmetry fixes carboxylate/hydroxyl scoring asymmetry.
#
# v33 result: 65/85 = 76.5%  (T_hot=500K anneal)
# v34 hypothesis: net +2..+5 from matrix fixes; softcore floor closes the
#   parabola loophole for catastrophically buried large-atom pairs.
#
# Control: launch_v34_ctrl.py — identical binary + env, but OLD matrix via
#   FLEXAIDDS_DATA_DIR=/tmp/data_ctrl_v34  (apples-to-apples comparison)
#
# Binary: bf6770a2 — rebuilt Jun 12 18:31 with vcfunction.cpp softcore_floor patch
#
import os, sys, signal, subprocess, hashlib

REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY     = f"{BUILD}/FlexAIDdS"
RUNNER     = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT     = os.path.expanduser("~/flexaidds_results/v34_20260612_metalfix_softfloor_n10fix")
CODES_CSV  = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"
CODES_FILE = f"{OUTPUT}/v34_codes_85.txt"
PIDFILE    = f"{OUTPUT}/v34.pid"
STDOUT_LOG = f"{OUTPUT}/v34_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v34_stderr.log"

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

ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark astex"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    # In the parallel launch scenario both v34 + v34_ctrl are launched in sequence;
    # allow up to 1 already-running benchmark (the first one we just started).
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
    "FLEXAIDDS_SOFTCORE_WAL":    "1",   # activates the parabola softcore
    "FLEXAIDDS_T_HOT":           "500",
    "OMP_WAIT_POLICY":           "passive",
    "OMP_PLACES":                "cores",
    "OMP_PROC_BIND":             "spread",
    # SOFTCORE_FLOOR uses compiled-in default 0.5; no env var needed unless overriding
})
# Ensure the patched matrix in build_lto is used (default data dir = binary dir)
# No FLEXAIDDS_DATA_DIR → uses build_lto/MC_st0r5.2_6.dat (already patched)
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
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
    print(f"launched v34 (SOFTCORE_FLOOR + MATRIX_FIX + T_HOT=500) → {OUTPUT}")
    print(f"  monitor: tail -f {STDOUT_LOG}")
    print(f"  pid:     {PIDFILE}")
    print(f"  est. wall time: ~12.5 h")
