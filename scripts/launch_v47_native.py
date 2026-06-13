#!/usr/bin/env python3
# launch_v47_native.py — daemonized launcher for the v47 NATIVE self-docking benchmark.
#
# v47_native = first corrected native self-docking run after v46 DatasetRunner cache bug.
#
# Root cause of v46 failure (2/85 = 2.4%):
#   DatasetRunner.load_astex_diverse() used ~/.flexaidds/benchmarks/ cache for ligand_path.
#   73/85 targets: cache SDF missing → no docking. 10/85: wrong ligand (HEM cofactor)
#   → count_delta > 2 → RMSD=999.
#
# Fix (no DatasetRunner.cpp changes needed):
#   Use crossdock_json: loader with benchmark_astex_native_85.json — receptor_id == ligand_id
#   for each pair. Explicit project SDF paths bypass cache entirely.
#   All 85 project SDFs verified present. 1TW6 uses 1TW6.pdb (no _apo.pdb).
#
# Config mirrors v43 (current record 68/85) + v44 settings:
#   --benchmark crossdock_json:<native_json>
#   --threads 10, --temperature 298, --job-timeout-seconds 5400
#   FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0  (OFF — net -1 in v44, keep off)
#   SOFTCORE_WAL=1, SOFTCORE_FLOOR=0.5 (v43 setting = current record)
#   T_HOT=500, NATIVE_SEED_FRAC=0.90, RESTARTS=3, N_ELITE=1, SEED_ELITISM=1
#   caffeinate -i, start_new_session=True (SIGHUP-immune double-fork)
#
# Binary: same LTO build as v43/v44/v45/v46. Do NOT rebuild.
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
OUTPUT      = os.path.expanduser("~/flexaidds_results/v47_native_20260613")
PROV_FILE   = f"{OUTPUT}/launch_provenance.json"

# ── Provenance anchors ────────────────────────────────────────────────────────
EXP_ENGINE_SHA  = "dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA  = "12d6d9bd19470eccc7e33dd828e13e8a9db50c1516a4ac28153a1b38ce1bd7f7"
EXP_MATRIX_MD5  = "72d7c7396702331d96ff12d18f831796"
GIT_COMMIT      = "0f8f160"


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

# Verify native JSON — quick sanity check
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
    "FLEXAIDDS_SEED_ELITISM":          "1",
    "FLEXAIDDS_N_ELITE":               "1",
    "FLEXAIDDS_BUDGET_SCALE":          "1",
    "FLEXAIDDS_SOFTCORE_WAL":          "1",
    "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
    "FLEXAIDDS_T_HOT":                 "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
    "FLEXAIDDS_DATA_DIR":              DATA_DIR,
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",      # OFF — net -1 in v44
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

    print(f"\nLaunching v47 NATIVE self-docking benchmark …")
    print(f"  native JSON: {JSON_PAIRS}")
    print(f"  oracle dir:  {ORACLE_DIR}")
    print(f"  output:      {OUTPUT}")

    child_pid = launch_session_isolated(cmd, env, OUTPUT, cwd=REPO)

    prov = {
        "version":       "v47_native",
        "launched_at":   datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":    GIT_COMMIT,
        "description":   (
            "Native self-docking Astex Diverse 85. Fixes v46 DatasetRunner cache bug "
            "by using crossdock_json loader with benchmark_astex_native_85.json — "
            "explicit project SDF paths, receptor_id==ligand_id."
        ),
        "oracle_mode":               True,
        "receptor_rotamer_prep":     False,
        "soft_wall_cutoff_angstrom": 0.40,
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
                "FLEXAIDDS_RESTARTS", "FLEXAIDDS_SEED_ELITISM",
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

    print(f"\n✓ v47_native launched (commit={GIT_COMMIT})")
    print(f"  pid:     {child_pid}  →  {OUTPUT}/benchmark.pid")
    print(f"  prov:    {PROV_FILE}")
    print(f"  monitor: tail -f {OUTPUT}/stdout.log")
    print(f"  errors:  tail -f {OUTPUT}/stderr.log")
