#!/usr/bin/env python3
# launch_v25.py — daemonized launcher for the v25 "multi-restart pooling" Astex 85 run.
#
# v25 = v24 (PSHARE alpha=4, VCT exp(-r/4), 2000-gen P6 budget, CF clustering,
# default matrix, oracle binding sites, Fix A RMSD symmetry-correction, Fix B
# frequency-gated cluster selection) + one new DatasetRunner-only lever:
#
#   Multi-restart pooling (FLEXAIDDS_RESTARTS=3):
#     The GA engine is run N=3 independent times per target, each with a
#     different random seed (restart 0 → clock, 1 → seed 7919, 2 → seed 15838).
#     All emitted PDB poses from all restarts are pooled before Fix B cluster
#     selection, directly expanding the oracle ceiling from 53→~70+ by exploiting
#     run-to-run variance in the GA's conformational sampling.  This is a pure
#     DatasetRunner change — the FlexAIDdS engine binary is UNCHANGED.
#
#   Oracle ceiling context (v24 baseline):
#     50/85 = 58.8% sub-2Å.  Oracle ceiling = 53/85 (62.4%).
#     32 targets collapsed to 1-2 clusters and missed near-native geometry entirely.
#     Three independent 1000-gen (or scaled) runs × the full GA diversity give
#     3 independent random explorations — oracle ceiling push target ≥ 70/85.
#
# WALL-TIME ESTIMATE:
#   v24: 85 targets ÷ 5 workers × ~15 min/target ≈ 4.25 h
#   v25: same but 3× per target ≈ 12.5 h  (restarts run sequentially per slot)
#   Per-restart timeout: 5400 s (90 min/restart). Total per-target budget: ~270 min.
#
# Same detachment contract as launch_v24.py:
#   double-fork + setsid, SIGHUP-immune, caffeinate -i, 5 workers × 2 OMP threads.
#
import os, sys, signal, subprocess, hashlib, shutil

REPO      = "/Users/lp.more/Projects/FlexAIDdS"
BUILD     = f"{REPO}/build"
BINARY    = f"{BUILD}/FlexAIDdS"
RUNNER    = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT    = os.path.expanduser("~/flexaidds_results/v25_20260610_diversity")
# Authoritative code list lives in the repo — no dependency on prior run dirs.
CODES_CSV = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"
CODES_FILE = f"{OUTPUT}/v25_codes_85.txt"
PIDFILE   = f"{OUTPUT}/v25.pid"
STDOUT_LOG = f"{OUTPUT}/v25_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v25_stderr.log"

# v25: runner SHA is intentionally NOT pinned here — benchmark_datasets was just
# rebuilt with multi-restart pooling code.  The SHA is printed for provenance.
# Pin it once the first build is verified by checking v25_benchmark.log.
EXP_ENGINE = "c1281359e43d5a6455d773d216784acb1c2a709debef8c212ec492926acbfb9f"  # unchanged from v24

def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

# ── pre-launch guards ────────────────────────────────────────────────────────
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
        "       FlexAIDdS engine should be byte-identical to v24 — only "
        "benchmark_datasets was rebuilt."
    )

# Warn if runner SHA matches v24 (means rebuild didn't happen)
V24_RUNNER_SHA = "5088c58b9cc29adab1f3ba0d6fa140c27455bcf6250244573b70207d05475784"
if runner_sha == V24_RUNNER_SHA:
    sys.exit(
        f"ERROR: runner SHA matches v24 ({runner_sha}) — benchmark_datasets "
        "was not rebuilt with multi-restart pooling code.\n"
        "       Run: cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build --target benchmark_datasets -j4"
    )

# Never run a competing benchmark (corrupts the shared grid/pose cache).
ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark astex"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    sys.exit(
        f"ERROR: a benchmark is already running (pids {ps.stdout.split()}) — "
        "abort to avoid cache corruption"
    )

# ── Generate v25 codes file from repo CSV (if not already present) ──────────
def make_codes_file(csv_path, out_path):
    """Extract PDB-ID column (col 0, skip header) from astex_diverse_set.csv."""
    codes = []
    with open(csv_path) as f:
        for i, line in enumerate(f):
            if i == 0:
                continue  # skip header
            pdb_id = line.split(",")[0].strip()
            if pdb_id:
                codes.append(pdb_id)
    with open(out_path, "w") as f:
        f.write("\n".join(codes) + "\n")
    return len(codes)

# ── Environment ──────────────────────────────────────────────────────────────
env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":        BINARY,
    "FLEXAIDDS_BUILD":         BUILD,
    "FLEXAIDDS_REPO":          REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR": ORACLE_DIR,   # oracle pocket mode
    "FLEXAIDDS_RESTARTS":      "3",            # ← multi-restart pooling (v25 key lever)
    "OMP_WAIT_POLICY":         "passive",
    "OMP_PLACES":              "cores",
    "OMP_PROC_BIND":           "spread",
})
# v25 differs from v24 ONLY by FLEXAIDDS_RESTARTS — clear all other opt-in flags
# so the run is a clean oracle-mode multi-restart benchmark.
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
          "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM"):
    env.pop(k, None)

# ── Command ──────────────────────────────────────────────────────────────────
# 5 workers × 2 OMP threads = 10 cores.
# Per-restart timeout 5400 s (90 min); 3 restarts → up to 270 min/target.
# skip-completed=true: if v25 is interrupted, re-running resumes from where it stopped.
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark",       "astex",
    "--only-codes",      CODES_FILE,   # 85 Astex Diverse targets
    "--output",          OUTPUT,
    "--threads",         "5",
    "--omp-threads",     "2",
    "--temperature",     "298",
    "--job-timeout-seconds", "5400",   # per-restart timeout (90 min); total ≈ 3× for 3 restarts
]

# ── Double-fork daemoniser ───────────────────────────────────────────────────
def double_fork_launch():
    if os.fork() > 0:
        return                    # parent returns immediately
    os.setsid()                   # detach from controlling terminal
    if os.fork() > 0:
        os._exit(0)               # intermediate exits
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
    n_codes = make_codes_file(CODES_CSV, CODES_FILE)
    print(f"codes: {n_codes} targets → {CODES_FILE}")
    double_fork_launch()
    print(f"launched v25 (RESTARTS=3) -> {OUTPUT}")
    print(f"  monitor: tail -f {STDOUT_LOG}")
    print(f"  pid:     {PIDFILE}")
    print(f"  est. wall time: ~12.5 h (3x v24)")
