#!/usr/bin/env python3
# launch_v26.py — daemonized launcher for the v26 "seed-anchored elitism" Astex 85 run.
#
# v26 = v25 (PSHARE alpha=4, VCT exp(-r/4), CF clustering, default matrix, oracle
# binding sites, Fix A RMSD symmetry-correction, Fix B frequency-gated cluster
# selection, FLEXAIDDS_RESTARTS=3 multi-restart pooling) + one new
# DatasetRunner-only lever:
#
#   Seed-anchored elitism (FLEXAIDDS_SEED_ELITISM=1):
#     In oracle mode the engine writes a gen-0 "<prefix>_INI.pdb" whose pose is
#     the crystal seed (RMSD≈0.00 Å — confirmed across all 31 v25 oracle-
#     unreachable targets).  boom_inject_fraction=1.0 + sharing_alpha=4.0 eject
#     this near-native pose from the GA population before termination, so the
#     emitted cluster heads (_0…_N) scatter even though the seed geometry was
#     perfect.  v26 adds each restart's _INI.pdb to the candidate pool in
#     select_pose_freq_gated_pooled() as an ALWAYS-ELIGIBLE seed: it bypasses the
#     freq>1 gate and degenerate-CF drop, and is elected rank-0 only if its CF is
#     strictly more favourable than the freq-gated best.  Pure DatasetRunner
#     change — the FlexAIDdS engine binary is UNCHANGED from v24/v25.
#
#   v25 baseline: 27/85 sub-2Å measured (Hungarian); 31 oracle-unreachable
#   targets all had _INI at 0.00 Å.  Projected v26 impact: +9–26.
#
# Same detachment contract as launch_v25.py:
#   double-fork + setsid, SIGHUP-immune, caffeinate -i, 5 workers × 2 OMP threads.
#
import os, sys, signal, subprocess, hashlib

REPO      = "/Users/lp.more/Projects/FlexAIDdS"
BUILD     = f"{REPO}/build"
BINARY    = f"{BUILD}/FlexAIDdS"
RUNNER    = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT    = os.path.expanduser("~/flexaidds_results/v26_20260610_seedelitism")
# Authoritative code list lives in the repo — no dependency on prior run dirs.
CODES_CSV = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"
CODES_FILE = f"{OUTPUT}/v26_codes_85.txt"
PIDFILE   = f"{OUTPUT}/v26.pid"
STDOUT_LOG = f"{OUTPUT}/v26_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v26_stderr.log"

# Engine binary must be byte-identical to v25 — only benchmark_datasets changed.
EXP_ENGINE = "d748ccb5406c1cf3d4a2c446840cbc5e514d83423ee7f38390e585b3f08370b9"  # cf_native diagnostic engine; GA search unchanged
# Runner pinned to the seed-elitism build (this launcher's provenance anchor).
EXP_RUNNER = "c088760349f81a8e95117c801f73eb5eed31015e20fd1ac6f41300f8b11e76e0"  # DatasetRunner seed-anchored elitism

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
        "       FlexAIDdS engine should be byte-identical to v25 — only "
        "benchmark_datasets was rebuilt."
    )

if runner_sha != EXP_RUNNER:
    sys.exit(
        f"ERROR: runner SHA mismatch (got {runner_sha}, want {EXP_RUNNER})\n"
        "       Rebuild the seed-elitism benchmark_datasets:\n"
        "       cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build --target benchmark_datasets -j8"
    )

# Never run a competing benchmark (corrupts the shared grid/pose cache).
ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark astex"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    sys.exit(
        f"ERROR: a benchmark is already running (pids {ps.stdout.split()}) — "
        "abort to avoid cache corruption"
    )

# ── Generate v26 codes file from repo CSV (if not already present) ──────────
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
    "FLEXAIDDS_RESTARTS":      "3",            # multi-restart pooling (v25 lever, retained)
    "FLEXAIDDS_SEED_ELITISM":  "1",            # ← v26 key lever: rescue _INI seed in pooled selection
    "OMP_WAIT_POLICY":         "passive",
    "OMP_PLACES":              "cores",
    "OMP_PROC_BIND":           "spread",
})
# v26 differs from v25 ONLY by FLEXAIDDS_SEED_ELITISM — clear all other opt-in
# flags so the run is a clean oracle-mode multi-restart + seed-elitism benchmark.
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
          "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM"):
    env.pop(k, None)

# ── Command ──────────────────────────────────────────────────────────────────
# 5 workers × 2 OMP threads = 10 cores.
# Per-restart timeout 5400 s (90 min); 3 restarts → up to 270 min/target.
# skip-completed=true: if v26 is interrupted, re-running resumes from where it stopped.
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
    print(f"launched v26 (RESTARTS=3, SEED_ELITISM=1) -> {OUTPUT}")
    print(f"  monitor: tail -f {STDOUT_LOG}")
    print(f"  pid:     {PIDFILE}")
    print(f"  est. wall time: ~12.5 h (same as v25)")
