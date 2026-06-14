#!/usr/bin/env python3
# launch_v50a.py — daemonized launcher for the v50a consensus benchmark.
#
# v50 = CLEAN AUTONOMOUS validation run on Astex re-docking 85 (self-docking / FLRP,
#       benchmark_astex_native_85.json — same dataset as v43..v49).
#
# Three levers (all committed in DatasetRunner.cpp, commit efc4f5d):
#   Change 1 — BCR gate is now DIAGNOSTIC-ONLY (permanent, no env toggle):
#     The v49 gate substituted the oracle-best cluster as the reported result
#     when the freq-gated selector missed a near-native cluster. That required
#     the crystal ligand (oracle-assisted), inflating the headline. v50 demotes
#     it to a pure observer: it logs [BCR-DIAGNOSTIC] when it WOULD have fired
#     but no longer overwrites rmsd_to_crystal / rmsd_hungarian / best_score.
#     The best_cluster_rmsd diagnostic column is still written to the CSV.
#     => the reported number is now a TRUE autonomous docking result.
#
#   Change 2 — Lever 1: dihedral eval-budget scaling (FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1):
#     n_gen_effective = n_gen_base * max(1.0, n_flex_bonds / 4.0). High-DoF
#     ligands get an eval budget matched to their search-space dimensionality
#     (v18 ablation: 4->11 dihedral genes needs >=2.75x more evaluations).
#
#   Change 3 — Lever 3: cross-restart cluster consensus re-ranking
#     (FLEXAIDDS_CONSENSUS_SCORER=1): a blind (no crystal / no oracle) signal.
#     Pool every emitted cluster pose across restarts; a candidate's consensus =
#     number of OTHER restart prefixes with a pose within 1.5A (Hungarian heavy-
#     atom RMSD). Re-rank consensus desc, CF asc; the winner is reported. If N
#     independent thermodynamic trajectories converge to the same basin, it is
#     more likely the true free-energy minimum than a basin only one restart
#     finds. Needs >=2 restart prefixes -> FLEXAIDDS_RESTARTS=3 here.
#
# v50a = 3 parallel restarts. v50b = same with 5 restarts.
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
#
import os
import sys
import subprocess
import hashlib
import json
import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

# ── Paths ────────────────────────────────────────────────────────────────────
REPO        = "/Users/lp.more/Projects/FlexAIDdS"
BUILD       = f"{REPO}/build_lto"
BINARY      = f"{BUILD}/FlexAIDdS"
RUNNER      = f"{BUILD}/benchmark_datasets"
DATA_DIR    = BUILD
ORACLE_DIR  = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS  = f"{REPO}/benchmarks/datasets/benchmark_astex_native_85.json"
OUTPUT      = os.path.expanduser("~/flexaidds_results/v50a_20260614_consensus")
PROV_FILE   = f"{OUTPUT}/launch_provenance.json"

# ── Concurrency knobs (env-overridable; defaults preserve safe single-instance) ─
# FLEXAIDDS_BENCH_THREADS  : concurrent FlexAIDdS workers (default 10).
# FLEXAIDDS_BENCH_CACHE    : isolated cache dir ("" → runner default, shared).
# FLEXAIDDS_ALLOW_CONCURRENT=1 : skip the single-instance pgrep guard so this run
#   may share the machine with another benchmark_datasets process.  Only safe
#   when paired with a DISTINCT --cache + --output (else the shared cache is
#   corrupted — see project_benchmark_harness_singleinstance).
BENCH_THREADS    = os.environ.get("FLEXAIDDS_BENCH_THREADS", "10")
BENCH_CACHE      = os.environ.get("FLEXAIDDS_BENCH_CACHE", "")
ALLOW_CONCURRENT = os.environ.get("FLEXAIDDS_ALLOW_CONCURRENT") == "1"

# ── Provenance anchors ────────────────────────────────────────────────────────
EXP_ENGINE_SHA  = "dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA  = "53fa471cfe3a55b2b071bf87e2181caba889ee92124553199df436275d714781"
EXP_MATRIX_MD5  = "72d7c7396702331d96ff12d18f831796"
GIT_COMMIT      = "efc4f5d"


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


# ── Pre-flight checks ─────────────────────────────────────────────────────────
for p in (BINARY, RUNNER, ORACLE_DIR, JSON_PAIRS, f"{DATA_DIR}/MC_st0r5.2_6.dat"):
    if not os.path.exists(p):
        sys.exit(f"ERROR: missing required path: {p}")

engine_sha = sha256(BINARY)
runner_sha = sha256(RUNNER)
matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")

print(f"engine SHA256     : {engine_sha}")
print(f"runner SHA256     : {runner_sha}")
print(f"matrix MD5        : {matrix_md5}")

if engine_sha != EXP_ENGINE_SHA:
    sys.exit(f"ERROR: engine SHA mismatch\n  got  {engine_sha}\n  want {EXP_ENGINE_SHA}")
if runner_sha != EXP_RUNNER_SHA:
    sys.exit(f"ERROR: runner SHA mismatch\n  got  {runner_sha}\n  want {EXP_RUNNER_SHA}")
if matrix_md5 != EXP_MATRIX_MD5:
    sys.exit(f"ERROR: matrix MD5 mismatch\n  got  {matrix_md5}\n  want {EXP_MATRIX_MD5}")

ps = subprocess.run(
    ["pgrep", "-x", "benchmark_datasets"],
    capture_output=True, text=True,
)
if ps.stdout.strip() and not ALLOW_CONCURRENT:
    sys.exit(
        f"ERROR: benchmark_datasets already running "
        f"(pids {ps.stdout.split()}) — abort to avoid collision"
    )
if ps.stdout.strip() and ALLOW_CONCURRENT:
    if not BENCH_CACHE:
        sys.exit(
            "ERROR: FLEXAIDDS_ALLOW_CONCURRENT=1 requires FLEXAIDDS_BENCH_CACHE "
            "(a distinct cache dir) to avoid corrupting the shared cache of the "
            f"already-running run (pids {ps.stdout.split()})"
        )
    print(
        f"NOTE: concurrent run permitted (other pids {ps.stdout.split()}); "
        f"isolated cache={BENCH_CACHE}, threads={BENCH_THREADS}"
    )

# Verify native self-docking JSON (receptor == ligand for every pair).
native = json.load(open(JSON_PAIRS))
assert len(native["pairs"]) == 85, "Expected 85 pairs in native JSON"
for pair in native["pairs"]:
    assert pair["receptor_id"] == pair["ligand_id"], \
        f"Non-native pair found: {pair['receptor_id']} != {pair['ligand_id']}"
print(f"Native JSON verified: 85 self-docking pairs ✓")

# ── Environment ───────────────────────────────────────────────────────────────
env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":                BINARY,
    "FLEXAIDDS_BUILD":                 BUILD,
    "FLEXAIDDS_REPO":                  REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
    "FLEXAIDDS_RESTARTS":              "3",
    "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
    # v50 levers:
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",   # Lever 1 (Change 2)
    "FLEXAIDDS_CONSENSUS_SCORER":      "1",   # Lever 3 (Change 3)
    # carried over from v49:
    "FLEXAIDDS_SEED_ELITISM":          "1",
    "FLEXAIDDS_N_ELITE":               "1",
    "FLEXAIDDS_BUDGET_SCALE":          "1",
    "FLEXAIDDS_SOFTCORE_WAL":          "1",
    "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
    "FLEXAIDDS_T_HOT":                 "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
    "FLEXAIDDS_DATA_DIR":              DATA_DIR,
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
    "OMP_WAIT_POLICY":                 "passive",
    "OMP_PLACES":                      "cores",
    "OMP_PROC_BIND":                   "spread",
})

for k in (
    "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
    "FLEXAIDDS_PRIORITY_TARGETS",
):
    env.pop(k, None)

# ── Command ───────────────────────────────────────────────────────────────────
cmd = [
    "caffeinate", "-i",
    RUNNER,
    "--benchmark",            f"crossdock_json:{JSON_PAIRS}",
    "--output",               OUTPUT,
    "--threads",              BENCH_THREADS,
    "--temperature",          "298",
    "--job-timeout-seconds",  "5400",
]
if BENCH_CACHE:
    cmd += ["--cache", BENCH_CACHE]

# ── Launch ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)

    print(f"\nLaunching v50a clean-autonomous consensus benchmark ...")
    print(f"  native JSON:    {JSON_PAIRS}")
    print(f"  oracle dir:     {ORACLE_DIR}")
    print(f"  output:         {OUTPUT}")
    print(f"  levers:         BCR=diagnostic-only | EVAL_SCALE_DIHEDRAL=1 | CONSENSUS_SCORER=1")
    print(f"  restarts:       3 (parallel)")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v50a_consensus",
        "launched_at":   datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":    GIT_COMMIT,
        "description":   (
            "Clean autonomous Astex re-docking 85 (self-docking / FLRP) validation. BCR gate is "
            "diagnostic-only (no oracle substitution of the reported result). "
            "Lever 1 (FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1): n_gen scaled by "
            "max(1.0, n_flex_bonds/4.0). Lever 3 (FLEXAIDDS_CONSENSUS_SCORER=1): "
            "cross-restart cluster consensus re-ranking (blind, delta=1.5A). "
            "3 parallel restarts."
        ),
        "oracle_mode":               True,
        "receptor_rotamer_prep":     False,
        "bcr_gate_enabled":          False,
        "bcr_gate_mode":             "diagnostic_only",
        "eval_scale_dihedral":       True,
        "consensus_scorer":          True,
        "consensus_delta_angstrom":  1.5,
        "n_restarts":                3,
        "parallel_restarts":         True,
        "benchmark":     f"crossdock_json:{JSON_PAIRS}",
        "binary":        BINARY,
        "binary_sha256": engine_sha,
        "runner":        RUNNER,
        "runner_sha256": runner_sha,
        "matrix":        f"{DATA_DIR}/MC_st0r5.2_6.dat",
        "matrix_md5":    matrix_md5,
        "oracle_dir":    ORACLE_DIR,
        "json_pairs":    JSON_PAIRS,
        "output_dir":    OUTPUT,
        "pid":           child_pid,
        "env_snapshot": {
            k: env[k] for k in (
                "FLEXAIDDS_ORACLE_SITE_DIR",
                "FLEXAIDDS_RESTARTS", "FLEXAIDDS_PARALLEL_RESTARTS",
                "FLEXAIDDS_EVAL_SCALE_DIHEDRAL", "FLEXAIDDS_CONSENSUS_SCORER",
                "FLEXAIDDS_SEED_ELITISM",
                "FLEXAIDDS_N_ELITE", "FLEXAIDDS_BUDGET_SCALE",
                "FLEXAIDDS_SOFTCORE_WAL", "FLEXAIDDS_SOFTCORE_FLOOR",
                "FLEXAIDDS_T_HOT", "FLEXAIDDS_NATIVE_SEED_FRAC",
                "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
                "FLEXAIDDS_DATA_DIR",
            )
        },
    }
    with open(PROV_FILE, "w") as f:
        json.dump(prov, f, indent=2)
        f.write("\n")

    print(f"\n✓ v50a_consensus launched (commit={GIT_COMMIT})")
    print(f"  pid:     {child_pid}  ->  {OUTPUT}/benchmark.pid")
    print(f"  prov:    {PROV_FILE}")
    print(f"  monitor: tail -f {OUTPUT}/stdout.log")
    print(f"  errors:  tail -f {OUTPUT}/stderr.log")
