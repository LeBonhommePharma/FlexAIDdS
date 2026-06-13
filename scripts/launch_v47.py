#!/usr/bin/env python3
# launch_v47.py — daemonized launcher for the v47 official cross-docking benchmark.
#
# v47 = official first cross-dock benchmark of Astex Diverse (85 pairs).
#       Queued as follow-up to v45 (initial trial run) and v46 (native self-dock).
#       Fire this once v45 + v46 results are collected and reviewed.
#
# Cross-docking semantics (identical to v45):
#   receptor[i] docks with ligand[(i+1) % 85]  (alphabetical cyclic shift)
#   No receptor is docked against its native ligand → purely cross-docking.
#
# Oracle mode ON — native binding site PDB provides seed placement geometry
# even though the docked ligand is non-native.  The oracle constrains the GA
# search to the known pocket; this isolates pose-generation capability from
# binding-site detection.
#
# receptor_rotamer_prep DISABLED — there is no native ligand geometry to guide
# sidechain prep against (ghost occupancy from a different compound is noise).
#
# Config mirrors v45 exactly:
#   --benchmark crossdock_json:<pairs_file>   (85-pair JSON bijection)
#   --threads 10, --temperature 298, --job-timeout-seconds 5400
#   FLEXAIDDS_RECEPTOR_ROTAMER_PREP=0         (OFF — explicit, cross-dock)
#   SOFTCORE_WAL=1, SOFTCORE_FLOOR=0.5, T_HOT=500
#   NATIVE_SEED_FRAC=0.90, RESTARTS=3, N_ELITE=1, SEED_ELITISM=1
#   BUDGET_SCALE=1
#   caffeinate -i wrapper (prevent macOS idle sleep during long run)
#   start_new_session=True (Fix B: POSIX session isolation)
#   Signal handler includes SIGHUP (Fix A, already committed in DatasetRunner.cpp)
#
# Binary: same LTO build as v45/v46 — do NOT rebuild before launching.
#
# Commits in this build:
#   91e5922  Fix: SIGHUP handler + graceful drain + launch session isolation
#   0f8f160  Add: FLEXAIDDS_RECEPTOR_ROTAMER_PREP env var override
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
#
import os
import sys
import subprocess
import hashlib
import json
import datetime

# Import shared launch helper (Fix B: POSIX double-fork + start_new_session=True)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

# ── Paths ────────────────────────────────────────────────────────────────────
REPO       = "/Users/lp.more/Projects/FlexAIDdS"
BUILD      = f"{REPO}/build_lto"
BINARY     = f"{BUILD}/FlexAIDdS"
RUNNER     = f"{BUILD}/benchmark_datasets"
DATA_DIR   = BUILD
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_crossdock_85.json"
OUTPUT     = os.path.expanduser("~/flexaidds_results/v47_20260613_crossdock_official")
PROV_FILE  = f"{OUTPUT}/launch_provenance.json"

# ── Provenance anchors ────────────────────────────────────────────────────────
# Same binary as v45/v46 — do NOT rebuild.  SHAs must match exactly.
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

# ── Environment ───────────────────────────────────────────────────────────────
env = dict(os.environ)
env.update({
    "FLEXAIDDS_BINARY":                BINARY,
    "FLEXAIDDS_BUILD":                 BUILD,
    "FLEXAIDDS_REPO":                  REPO,
    "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,   # oracle binding sites for seed placement
    "FLEXAIDDS_RESTARTS":              "3",
    "FLEXAIDDS_SEED_ELITISM":          "1",
    "FLEXAIDDS_N_ELITE":               "1",
    "FLEXAIDDS_BUDGET_SCALE":          "1",
    "FLEXAIDDS_SOFTCORE_WAL":          "1",           # enable distance-based softcore
    "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
    "FLEXAIDDS_T_HOT":                 "500",
    "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",        # 90% seeds from oracle site
    "FLEXAIDDS_DATA_DIR":              DATA_DIR,
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",           # OFF — cross-docking: no native lig geom
    "OMP_WAIT_POLICY":                 "passive",
    "OMP_PLACES":                      "cores",
    "OMP_PROC_BIND":                   "spread",
})

# Clear ablation/tuning vars that might bleed in from the calling shell
for k in (
    "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
    "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
    "FLEXAIDDS_VCT_R0", "FLEXAIDDS_VCT_NORM",
    "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
    "FLEXAIDDS_RING_FLEX",
    "FLEXAIDDS_PRIORITY_TARGETS",          # no priority targeting for cross-docking
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

# ── Launch (Fix B: double-fork daemon, start_new_session=True) ────────────────
if __name__ == "__main__":
    os.makedirs(OUTPUT, exist_ok=True)

    print(f"\nLaunching v47 official cross-docking benchmark …")
    print(f"  pairs file: {JSON_PAIRS}")
    print(f"  oracle dir: {ORACLE_DIR}")
    print(f"  output:     {OUTPUT}")

    child_pid = launch_session_isolated(
        cmd,
        env,
        OUTPUT,
        cwd=REPO,
    )

    # ── Write provenance JSON ─────────────────────────────────────────────────
    prov = {
        "version":       "v47",
        "launched_at":   datetime.datetime.utcnow().isoformat() + "Z",
        "git_commit":    GIT_COMMIT,
        "description":   (
            "Official first cross-docking benchmark — 85-pair bijection: "
            "receptor[i] → ligand[(i+1)%85] (cyclic shift, alphabetical order "
            "over Astex Diverse). Queued after v45 trial + v46 native results."
        ),
        "oracle_mode":              True,
        "receptor_rotamer_prep":    False,
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

    print(f"\n✓ v47 crossdock_official launched (commit={GIT_COMMIT})")
    print(f"  pid:     {child_pid}  →  {OUTPUT}/benchmark.pid")
    print(f"  prov:    {PROV_FILE}")
    print(f"  monitor: tail -f {OUTPUT}/stdout.log")
    print(f"  errors:  tail -f {OUTPUT}/stderr.log")
