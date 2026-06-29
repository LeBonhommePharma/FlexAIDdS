#!/usr/bin/env python3
# launch_v131_smoke12.py — v131 Lane A smoke (12 targets) before full-85
#
# Binary: FlexAIDdS_v131_safe worktree via build_v131_safe.sh
#   Base 82ad51f4 + cherry-picks 04ff1735 (sulfo) + bf8cf1d2 (holo data)
#   SoA OFF, pre-27e68e51 gaboom parallel regression
#
# Protocol: v127 (r0=4, consensus ON, native_seed_frac=0.90, rotamer off)
# JSON: benchmark_astex_smoke_12_v131.json (1TW6_holo, 1HNN expB site)
#
# Gate: >=10/12 Hungarian success AND 0/3 v130 regression guards fail
#       (1HQ2, 1S3V, 1T40 must stay PASS vs v127)
#
# Usage:
#   python3 scripts/launch_v131_smoke12.py
#   python3 scripts/launch_v131_smoke12.py --skip-build
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = "/Users/lp.more/Projects/FlexAIDdS"
GIT_ROOT = subprocess.check_output(
    ["git", "-C", SCRIPT_DIR, "rev-parse", "--show-toplevel"], text=True
).strip()


def resolve_worktree() -> str:
    candidates = [
        os.environ.get("FLEXAIDDS_V131_WORKTREE", ""),
        f"{REPO}/../FlexAIDdS_v131_safe",
        "/Users/lp.more/.grok/worktrees/projects-flexaidds/FlexAIDdS_v131_safe",
    ]
    for path in candidates:
        if path and os.path.isdir(os.path.join(path, "build_lto")):
            return path
    for path in candidates:
        if path and os.path.isdir(path):
            return path
    return f"{REPO}/../FlexAIDdS_v131_safe"


WORKTREE = resolve_worktree()
BUILD = f"{WORKTREE}/build_lto"
BUILD_SCRIPT = os.path.join(SCRIPT_DIR, "build_v131_safe.sh")
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY = "/tmp/FlexAIDdS_v131_safe"
RUNNER = f"{BUILD}/benchmark_datasets"
DATA_DIR = BUILD
JSON_PAIRS = f"{GIT_ROOT}/benchmarks/datasets/benchmark_astex_smoke_12_v131.json"


def resolve_oracle_dir() -> str:
    rel = "benchmarks/astex_diverse/astex_diverse"
    for root in (REPO, WORKTREE, GIT_ROOT):
        path = os.path.join(root, rel)
        if os.path.isdir(path):
            return path
    return os.path.join(REPO, rel)


def resolve_oracle_asset(*parts: str) -> str:
    for root in (REPO, WORKTREE, GIT_ROOT):
        path = os.path.join(root, "benchmarks", "astex_diverse", "astex_diverse", *parts)
        if os.path.exists(path):
            return path
    return os.path.join(REPO, "benchmarks", "astex_diverse", "astex_diverse", *parts)


ORACLE_DIR = resolve_oracle_dir()
TW6_HOLO = resolve_oracle_asset("1TW6", "1TW6_holo.pdb")
HNN_LIGAND_SITE = resolve_oracle_asset("1HNN", "1HNN_ligand_centered_site.pdb")
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")

BASE_COMMIT = "82ad51f4"
CHERRY_PICKS = ("04ff1735", "bf8cf1d2")
SMOKE_TARGETS = (
    "1G9V,1TW6,1GPK,1T9B,"
    "1HQ2,1S3V,1T40,"
    "1M2Z,1R58,1XM6,"
    "1HNN,1N2V"
)
REGRESSION_GUARD = ("1HQ2", "1S3V", "1T40")
V127_GUARD = ("1M2Z", "1R58", "1XM6")
SUCCESS_GATE_MIN = 10
REF_V127_DIR = "v127_20260629_0139_optB_full85"
REF_V130_DIR = "v130_20260629_0548_sulfo_expB_full85"

ENV_SNAPSHOT_KEYS = (
    "FLEXAIDDS_BINARY",
    "FLEXAIDDS_ORACLE_SITE_DIR",
    "FLEXAIDDS_RESTARTS",
    "FLEXAIDDS_PARALLEL_RESTARTS",
    "FLEXAIDDS_EVAL_SCALE_DIHEDRAL",
    "FLEXAIDDS_CONSENSUS_SCORER",
    "FLEXAIDDS_SEED_ELITISM",
    "FLEXAIDDS_N_ELITE",
    "FLEXAIDDS_BUDGET_SCALE",
    "FLEXAIDDS_SOFTCORE_WAL",
    "FLEXAIDDS_SOFTCORE_FLOOR",
    "FLEXAIDDS_T_HOT",
    "FLEXAIDDS_NATIVE_SEED_FRAC",
    "FLEXAIDDS_VCT_R0",
    "FLEXAIDDS_RECEPTOR_ROTAMER_PREP",
    "FLEXAIDDS_DATA_DIR",
    "FLEXAIDDS_BENCH_CACHE",
    "FLEXAIDDS_ALLOW_CONCURRENT",
)


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


def run_build(force: bool = False) -> None:
    if not os.path.isfile(BUILD_SCRIPT):
        sys.exit(f"ERROR: missing build script: {BUILD_SCRIPT}")
    cmd = ["bash", BUILD_SCRIPT]
    if force:
        cmd.append("--force")
    print(f"Running: {' '.join(cmd)}")
    subprocess.check_call(cmd, cwd=REPO)


def worktree_head() -> str:
    if not os.path.isdir(WORKTREE):
        return ""
    return subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], cwd=WORKTREE, text=True
    ).strip()


def ref_snapshot(pdb: str) -> dict:
    import csv

    snap = {}
    for label, run_dir in [
        ("v127", REF_V127_DIR),
        ("v130", REF_V130_DIR),
    ]:
        p = RESULTS_DIR / run_dir / pdb / "result.csv"
        if p.exists():
            r = list(csv.DictReader(open(p)))[0]
            snap[label] = {
                "success": r.get("success"),
                "rmsd_hungarian": r.get("rmsd_hungarian"),
                "rmsd_to_crystal": r.get("rmsd_to_crystal"),
                "best_cluster_rmsd": r.get("best_cluster_rmsd"),
                "cf_native": r.get("cf_native"),
                "pose_source": r.get("pose_source"),
            }
    return snap


def main():
    parser = argparse.ArgumentParser(description="v131 Lane A smoke-12 launcher")
    parser.add_argument("--skip-build", action="store_true",
                        help="Require existing v131_safe worktree build artifacts")
    parser.add_argument("--force-build", action="store_true",
                        help="Pass --force to build_v131_safe.sh")
    args = parser.parse_args()

    if not args.skip_build:
        run_build(force=args.force_build)
    else:
        for p in (BINARY_SRC, RUNNER):
            if not os.path.exists(p):
                sys.exit(
                    f"ERROR: --skip-build but missing {p}; run build_v131_safe.sh first"
                )

    for p in (
        ORACLE_DIR,
        JSON_PAIRS,
        TW6_HOLO,
        HNN_LIGAND_SITE,
        f"{DATA_DIR}/MC_st0r5.2_6.dat",
    ):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    native = json.load(open(JSON_PAIRS))
    assert len(native["pairs"]) == 12
    tw6 = next(p for p in native["pairs"] if p["receptor_id"] == "1TW6")
    hnn = next(p for p in native["pairs"] if p["receptor_id"] == "1HNN")
    assert tw6["receptor_pdb"].endswith("1TW6/1TW6_holo.pdb")
    assert os.path.isfile(tw6["receptor_pdb"])
    assert hnn["oracle_site_pdb"].endswith("1HNN/1HNN_ligand_centered_site.pdb")
    assert os.path.isfile(hnn["oracle_site_pdb"])

    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{DATA_DIR}/MC_st0r5.2_6.dat")
    wt_head = worktree_head()

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"v131_{tag}_smoke12_safe")
    cache = str(RESULTS_DIR / "cache_v131_smoke12_safe")
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

    env = dict(os.environ)
    env.update({
        "FLEXAIDDS_BINARY":                BINARY,
        "FLEXAIDDS_BUILD":                 BUILD,
        "FLEXAIDDS_REPO":                  REPO,
        "FLEXAIDDS_ORACLE_SITE_DIR":       ORACLE_DIR,
        "FLEXAIDDS_RESTARTS":              "5",
        "FLEXAIDDS_PARALLEL_RESTARTS":     "1",
        "FLEXAIDDS_EVAL_SCALE_DIHEDRAL":   "1",
        "FLEXAIDDS_CONSENSUS_SCORER":      "1",
        "FLEXAIDDS_SEED_ELITISM":          "1",
        "FLEXAIDDS_N_ELITE":               "1",
        "FLEXAIDDS_BUDGET_SCALE":          "1",
        "FLEXAIDDS_SOFTCORE_WAL":          "1",
        "FLEXAIDDS_SOFTCORE_FLOOR":        "0.5",
        "FLEXAIDDS_T_HOT":                 "500",
        "FLEXAIDDS_NATIVE_SEED_FRAC":      "0.90",
        "FLEXAIDDS_VCT_R0":                "4",
        "FLEXAIDDS_RECEPTOR_ROTAMER_PREP": "0",
        "FLEXAIDDS_DATA_DIR":              DATA_DIR,
        "FLEXAIDDS_ALLOW_CONCURRENT":      "1",
        "FLEXAIDDS_BENCH_CACHE":           cache,
        "OMP_WAIT_POLICY":                 "passive",
        "OMP_PLACES":                      "cores",
        "OMP_PROC_BIND":                   "spread",
    })
    for k in (
        "FLEXAIDDS_USE_DP", "FLEXAIDDS_FINE_GRID",
        "FLEXAIDDS_FORCE_RIGID", "FLEXAIDDS_USE_SHANNON",
        "FLEXAIDDS_VCT_NORM",
        "FLEXAIDDS_SHARING_ALPHA", "FLEXAIDDS_BOOM_FRAC",
        "FLEXAIDDS_RING_FLEX",
        "FLEXAIDDS_THERMO", "FLEXAIDDS_HVIB",
        "FLEXAIDDS_PRIORITY_TARGETS", "FLEXAIDDS_FREQSEL",
    ):
        env.pop(k, None)

    cmd = [
        "caffeinate", "-i",
        RUNNER,
        "--benchmark",           f"crossdock_json:{JSON_PAIRS}",
        "--output",              output,
        "--threads",             bench_threads,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               cache,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print("\nLaunching v131 Lane A smoke-12 — v127-safe binary + sulfo + holo data")
    print(f"  worktree  : {WORKTREE} @ {wt_head}")
    print(f"  base      : {BASE_COMMIT} + cherry-picks {' '.join(CHERRY_PICKS)}")
    print(f"  targets   : {SMOKE_TARGETS}")
    print(f"  output    : {output}")
    print(f"  gate      : >={SUCCESS_GATE_MIN}/12 pass, 0/3 regression guards fail")
    print(f"  threads   : {bench_threads}")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    ref = {pdb: ref_snapshot(pdb) for pdb in SMOKE_TARGETS.split(",")}

    prov_doc = {
        "version":       "v131_smoke12_safe",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "worktree":      WORKTREE,
        "worktree_head": wt_head,
        "base_commit":   BASE_COMMIT,
        "cherry_picks":  list(CHERRY_PICKS),
        "description": (
            "v131 Lane A smoke before full-85: v127-safe binary (82ad51f4 + sulfo + holo data), "
            "SoA OFF, v127 protocol (r0=4, consensus ON). Validates data/sulfo gains without "
            "v130 perf-regression binary (27e68e51 gaboom parallel / SoA churn)."
        ),
        "smoke_targets": SMOKE_TARGETS,
        "smoke_roles": {
            "data_fix": ["1G9V"],
            "holo_receptor": ["1TW6"],
            "sulfo_remap": ["1GPK", "1T9B"],
            "v130_regression_guard": list(REGRESSION_GUARD),
            "v127_only_guard": list(V127_GUARD),
            "informational": ["1HNN", "1N2V"],
        },
        "success_gate": {
            "min_pass": SUCCESS_GATE_MIN,
            "max_regression_guard_fail": 0,
            "regression_guard_targets": list(REGRESSION_GUARD),
            "full85_launch_if": (
                f">={SUCCESS_GATE_MIN}/12 Hungarian success AND "
                "0/3 regression guards (1HQ2,1S3V,1T40) fail"
            ),
        },
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs":     JSON_PAIRS,
        "output_dir":     output,
        "cache_dir":      cache,
        "pid":            child_pid,
        "reference_runs": {
            "v127_dir": REF_V127_DIR,
            "v130_dir": REF_V130_DIR,
            "per_target": ref,
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nv131 smoke-12 launched pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


if __name__ == "__main__":
    main()