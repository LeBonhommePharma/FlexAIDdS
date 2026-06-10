#!/usr/bin/env python3
# launch_v27.py — daemonized launcher for the v27 "true GA-internal elitism" Astex 85 run.
#
# v27 = v26 (PSHARE alpha=4, VCT exp(-r/4), CF clustering, default matrix, oracle
# binding sites, Fix A RMSD symmetry-correction, Fix B frequency-gated cluster
# selection, FLEXAIDDS_RESTARTS=3 multi-restart pooling, FLEXAIDDS_SEED_ELITISM=1
# post-hoc _INI rescue) + two NEW levers that move the fix from DatasetRunner
# post-processing INTO the GA itself:
#
#   1. True GA-internal elitism (FLEXAIDDS_N_ELITE=1):
#        gaboom.cpp now snapshots the n_elite lowest-CF individuals EVERY
#        generation *before* boom injection and reproduce()/niche-sharing, then
#        restores them over the worst of the freshly reproduced population.  The
#        elites are additionally exempt from the PSHARE/SMFREE sharing fitness
#        reduction.  This guarantees the running best (the near-native seed in
#        oracle mode) can never be ejected by boom_inject_fraction=1.0 +
#        sharing_alpha=4.0 — the exact mechanism that scattered v25's cluster
#        heads despite a perfect 0.00 Å seed.  Unlike v26 (which rescued the seed
#        in select_pose_freq_gated_pooled AFTER the GA discarded it), v27 keeps
#        it in the population the whole run, so the emitted cluster heads are
#        built around it.  Engine binary CHANGED vs v25/v26.
#
#   2. High-DoF eval-budget scaling (FLEXAIDDS_BUDGET_SCALE=1):
#        DatasetRunner widens n_gen by max(1.0, n_genes/7.0) for high-DoF ligands
#        (n_genes >= 14, i.e. >=10 rotatable bonds) so the extra torsional
#        dimensions get enough GA generations for the elite seed to anchor a
#        converged basin.  Stacks on the existing ceil(n_genes/4) scaling.
#
#   v26 baseline: seed-anchored elitism rescued v25 oracle-unreachables in 3/3
#   smoke (1W2G/2D3U/1Q1G → 0.00 Å).  v27 projected impact: lift the high-DoF
#   unreachables (e.g. 1HQ2) that the post-hoc rescue still scattered, by keeping
#   the elite in-population across the full search.
#
# Same detachment contract as launch_v26.py:
#   double-fork + setsid, SIGHUP-immune, caffeinate -i, 5 workers × 2 OMP threads.
#
import os, sys, signal, subprocess, hashlib

REPO      = "/Users/lp.more/Projects/FlexAIDdS"
BUILD     = f"{REPO}/build_lto"               # ← v27 builds in build_lto (build/ is the running v26)
BINARY    = f"{BUILD}/FlexAIDdS"
RUNNER    = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
OUTPUT    = os.path.expanduser("~/flexaidds_results/v27_20260610_gaelitism")
# Authoritative code list lives in the repo — no dependency on prior run dirs.
CODES_CSV = f"{REPO}/benchmarks/astex_diverse/astex_diverse_set.csv"
CODES_FILE = f"{OUTPUT}/v27_codes_85.txt"
PIDFILE   = f"{OUTPUT}/v27.pid"
STDOUT_LOG = f"{OUTPUT}/v27_benchmark.log"
STDERR_LOG = f"{OUTPUT}/v27_stderr.log"

# v27 rebuilt BOTH binaries (GA-internal elitism touches the engine + the runner
# echoes n_elite / budget_scale).  Pin both to the build_lto SHAs.
EXP_ENGINE = "40a7b3378625a538ba900ac6fe17ad78bec95f2fa0fe82393d87c2d52cdac8e5"  # FlexAIDdS — GA-internal elitism
EXP_RUNNER = "31160ef6c8776407b406176b790140deac74f67f84b3d34dbb3ffa1e90292cae"  # benchmark_datasets — n_elite echo + budget scale

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
        "       Rebuild the GA-internal-elitism engine:\n"
        "       cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build_lto --target FlexAIDdS -j8"
    )

if runner_sha != EXP_RUNNER:
    sys.exit(
        f"ERROR: runner SHA mismatch (got {runner_sha}, want {EXP_RUNNER})\n"
        "       Rebuild the v27 benchmark_datasets:\n"
        "       cd /Users/lp.more/Projects/FlexAIDdS && "
        "/opt/homebrew/bin/cmake --build build_lto --target benchmark_datasets -j8"
    )

# Never run a competing benchmark (corrupts the shared grid/pose cache).
# This guard intentionally blocks v27 from launching while v26 is still running.
ps = subprocess.run(["pgrep", "-f", "benchmark_datasets --benchmark astex"],
                    capture_output=True, text=True)
if ps.stdout.strip():
    sys.exit(
        f"ERROR: a benchmark is already running (pids {ps.stdout.split()}) — "
        "abort to avoid cache corruption (wait for v26 to finish before launching v27)"
    )

# ── Generate v27 codes file from repo CSV (if not already present) ──────────
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
    "FLEXAIDDS_SEED_ELITISM":  "1",            # v26 lever, retained: rescue _INI in pooled selection
    "FLEXAIDDS_N_ELITE":       "1",            # ← v27 key lever: true GA-internal elitism
    "FLEXAIDDS_BUDGET_SCALE":  "1",            # ← v27 lever: widen n_gen for high-DoF ligands
    "OMP_WAIT_POLICY":         "passive",
    "OMP_PLACES":              "cores",
    "OMP_PROC_BIND":           "spread",
})
# Clear all other opt-in flags so the run is a clean oracle-mode multi-restart +
# seed-elitism + GA-internal-elitism benchmark.
for k in ("FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID", "FLEXAIDDS_DATA_DIR",
          "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
          "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
          "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC"):
    env.pop(k, None)

# ── Command ──────────────────────────────────────────────────────────────────
# 5 workers × 2 OMP threads = 10 cores.
# Per-restart timeout 5400 s (90 min); 3 restarts → up to 270 min/target.
# skip-completed=true: if v27 is interrupted, re-running resumes from where it stopped.
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
    print(f"launched v27 (RESTARTS=3, SEED_ELITISM=1, N_ELITE=1, BUDGET_SCALE=1) -> {OUTPUT}")
    print(f"  monitor: tail -f {STDOUT_LOG}")
    print(f"  pid:     {PIDFILE}")
    print(f"  est. wall time: ~12.5 h (same as v26)")
