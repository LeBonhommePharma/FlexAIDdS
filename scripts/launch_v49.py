#!/usr/bin/env python3
# launch_v49.py — daemonized launcher for the v49 BCR-gate fix benchmark.
#
# v49 = native self-docking Astex Diverse 85 with corrected BCR-gate.
#
# BCR-gate v49 fix (db7f991):
#   v48 bug: gate substituted pfx+'_0.pdb' (CF-rank-0 = deepest false minimum)
#   instead of the actual bcr cluster (lowest RMSD = best_cluster_pfx+_idx.pdb).
#   Result: gate fired 10x but only rescued 1 target (1P2Y). 1SQ5 had bcr=0.000A
#   but CF0=3.833A — gate had the answer and pointed to the wrong cluster.
#   Fix: track best_cluster_pfx during oracle scan; gate uses pfx+_N.pdb where
#   N = best_cluster_idx is the index that achieved best_cluster_rmsd.
#
# Expected targets rescued by BCR-gate with correct implementation:
#   bcr<2.0A: 1MEH(0.481) 1N2J(0.545) 1L2S(0.937) 1U4D(0.123) 1X8X(0.904)
#             1Q4G(1.020) 1XM6(0.794) 1SQ5(0.000) 1P2Y(0.641)
#   All 9 should be rescued if bcr cluster is reachable via reported pose path.
#
# All other parameters identical to v48 (parallel restarts, oracle, NATURaL gate,
# softcore WAL, n_restarts=3).
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
OUTPUT      = os.path.expanduser("~/flexaidds_results/v49_20260614_bcrfix")
PROV_FILE   = f"{OUTPUT}/launch_provenance.json"

# ── Provenance anchors ────────────────────────────────────────────────────────
EXP_ENGINE_SHA  = "dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA  = "084d61b14251885ee2ec580b9af9893844b2ca03c2f1789bb18b4f4f8fea69c2"
EXP_MATRIX_MD5  = "72d7c7396702331d96ff12d18f831796"
GIT_COMMIT      = "db7f991"


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
if ps.stdout.strip():
    sys.exit(
        f"ERROR: benchmark_datasets already running "
        f"(pids {ps.stdout.split()}) — abort to avoid collision"
    )

# Verify native JSON
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
    "--threads",              "10",
    "--temperature",          "298",
    "--job-timeout-seconds",  "5400",
]

# ── Launch ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)

    print(f"\nLaunching v49 BCR-gate fix native self-docking benchmark ...")
    print(f"  native JSON: {JSON_PAIRS}")
    print(f"  oracle dir:  {ORACLE_DIR}")
    print(f"  output:      {OUTPUT}")
    print(f"  key fix:     BCR-gate now uses best_cluster_pfx+_idx.pdb (db7f991)")
    print(f"  bcr targets: 1MEH 1N2J 1L2S 1U4D 1X8X 1Q4G 1XM6 1SQ5 1P2Y")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v49_bcrfix",
        "launched_at":   datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":    GIT_COMMIT,
        "description":   (
            "Native self-docking Astex Diverse 85 with corrected BCR-gate (db7f991). "
            "v48 bug: gate used _0.pdb (CF-rank-0) instead of bcr cluster (best_cluster_pfx+_idx). "
            "Fix: track best_cluster_pfx during oracle scan; gate substitutes pfx+_N.pdb "
            "where N=best_cluster_idx achieved best_cluster_rmsd. "
            "9 targets have bcr<2.0A and should be rescued: "
            "1MEH(0.481) 1N2J(0.545) 1L2S(0.937) 1U4D(0.123) 1X8X(0.904) "
            "1Q4G(1.020) 1XM6(0.794) 1SQ5(0.000) 1P2Y(0.641)."
        ),
        "oracle_mode":               True,
        "receptor_rotamer_prep":     False,
        "soft_wall_cutoff_angstrom": 0.40,
        "bcr_gate_enabled":          True,
        "bcr_gate_version":          "v49_fix",
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

    print(f"\n✓ v49_bcrfix launched (commit={GIT_COMMIT})")
    print(f"  pid:     {child_pid}  ->  {OUTPUT}/benchmark.pid")
    print(f"  prov:    {PROV_FILE}")
    print(f"  monitor: tail -f {OUTPUT}/stdout.log")
    print(f"  errors:  tail -f {OUTPUT}/stderr.log")
