#!/usr/bin/env python3
# launch_v131_safe_full85.py — v131 Lane A full-85 (v127-safe binary + data/sulfo fixes)
#
# Prerequisite: v131 smoke-12 gate PASS (>=10/12, 0/3 regression guards fail).
#   python3 scripts/launch_v131_smoke12.py
#
# Binary: FlexAIDdS_v131_safe worktree (82ad51f4 + 04ff1735 + bf8cf1d2, SoA OFF)
# Protocol: v127 (r0=4, consensus ON, native_seed_frac=0.90, rotamer off)
# JSON: benchmark_astex_native_85_v131.json with runtime path patch for holo/data targets
#
# Usage:
#   python3 scripts/launch_v131_safe_full85.py --skip-build
#   python3 scripts/launch_v131_safe_full85.py --skip-build --ignore-smoke-gate
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

import argparse
import csv
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated
from v131_safe_common import (
    REPO,
    git_root,
    patch_manifest,
    resolve_oracle_dir,
    resolve_worktree,
    scrub_env,
    validate_lane_a_assets,
    v127_protocol_env,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GIT_ROOT = git_root(SCRIPT_DIR)
WORKTREE = resolve_worktree()
BUILD = f"{WORKTREE}/build_lto"
BUILD_SCRIPT = os.path.join(SCRIPT_DIR, "build_v131_safe.sh")
BINARY_SRC = f"{BUILD}/FlexAIDdS"
BINARY = "/tmp/FlexAIDdS_v131_safe"
RUNNER = f"{BUILD}/benchmark_datasets"
JSON_SRC = os.path.join(GIT_ROOT, "benchmarks/datasets/benchmark_astex_native_85_v131.json")
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")

BASE_COMMIT = "82ad51f4"
CHERRY_PICKS = ("04ff1735", "bf8cf1d2")
REF_V127_DIR = "v127_20260629_0139_optB_full85"
REF_V109_DIR = "v109_20260626_tier1_consensus5r"
SMOKE_GATE_MIN = 10
REGRESSION_GUARD = ("1HQ2", "1S3V", "1T40")

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
    cmd = ["bash", BUILD_SCRIPT]
    if force:
        cmd.append("--force")
    subprocess.check_call(cmd, cwd=REPO)


def worktree_head() -> str:
    if not os.path.isdir(WORKTREE):
        return ""
    return subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"], cwd=WORKTREE, text=True
    ).strip()


def latest_smoke_summary() -> dict | None:
    dirs = sorted(RESULTS_DIR.glob("v131_*_smoke12_safe"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not dirs:
        return None
    run_dir = dirs[0]
    done = []
    pass_ = 0
    reg_fail = []
    for td in run_dir.iterdir():
        if not td.is_dir() or not (td / "result.csv").exists():
            continue
        done.append(td.name)
        r = list(csv.DictReader(open(td / "result.csv")))[0]
        ok = r.get("success") == "1"
        if ok:
            pass_ += 1
        if td.name in REGRESSION_GUARD and not ok:
            reg_fail.append(td.name)
    return {
        "run_dir": run_dir.name,
        "completed": len(done),
        "pass": pass_,
        "regression_guard_fail": reg_fail,
        "gate_pass": pass_ >= SMOKE_GATE_MIN and len(reg_fail) == 0,
    }


def main():
    parser = argparse.ArgumentParser(description="v131 Lane A full-85 launcher")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--force-build", action="store_true")
    parser.add_argument(
        "--ignore-smoke-gate",
        action="store_true",
        help="Launch even if smoke-12 gate has not passed (not recommended)",
    )
    args = parser.parse_args()

    if not args.skip_build:
        run_build(force=args.force_build)
    else:
        for p in (BINARY_SRC, RUNNER):
            if not os.path.exists(p):
                sys.exit(f"ERROR: --skip-build but missing {p}")

    smoke = latest_smoke_summary()
    if not args.ignore_smoke_gate:
        if smoke is None:
            sys.exit("ERROR: no v131 smoke-12 run found; run launch_v131_smoke12.py first")
        if not smoke["gate_pass"]:
            sys.exit(
                f"ERROR: smoke gate FAIL ({smoke['pass']}/12 pass, "
                f"regression guards failed: {smoke['regression_guard_fail']}). "
                f"Fix manifest/binary and re-smoke, or pass --ignore-smoke-gate."
            )

    if not os.path.isfile(JSON_SRC):
        sys.exit(f"ERROR: missing JSON manifest: {JSON_SRC}")

    validate_lane_a_assets(WORKTREE, GIT_ROOT)

    manifest = patch_manifest(json.load(open(JSON_SRC)), WORKTREE, GIT_ROOT)
    assert len(manifest["pairs"]) == 85

    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_v131_safe_full85.json", delete=False
    ) as tmp:
        json.dump(manifest, tmp, indent=2)
        tmp.write("\n")
        json_pairs = tmp.name

    oracle_dir = resolve_oracle_dir(WORKTREE, GIT_ROOT)
    shutil.copy2(BINARY_SRC, BINARY)
    os.chmod(BINARY, 0o755)

    engine_sha = sha256(BINARY)
    runner_sha = sha256(RUNNER)
    matrix_md5 = md5(f"{BUILD}/MC_st0r5.2_6.dat")
    wt_head = worktree_head()

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"v131_{tag}_safe_full85")
    cache = str(RESULTS_DIR / "cache_v131_safe_full85")
    prov = f"{output}/launch_provenance.json"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "4")

    env = scrub_env(dict(os.environ))
    env.update(v127_protocol_env(BINARY, BUILD, cache, oracle_dir))

    cmd = [
        "caffeinate", "-i",
        RUNNER,
        "--benchmark",           f"crossdock_json:{json_pairs}",
        "--output",              output,
        "--threads",             bench_threads,
        "--temperature",         "298",
        "--job-timeout-seconds", "7200",
        "--cache",               cache,
        "--mode",                "oracle-ceiling",
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print("\nLaunching v131 Lane A full-85 — v127-safe binary + sulfo + holo data")
    print(f"  worktree  : {WORKTREE} @ {wt_head}")
    print(f"  output    : {output}")
    print(f"  smoke ref : {smoke['run_dir'] if smoke else 'none'}")
    print(f"  threads   : {bench_threads}")

    child_pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov_doc = {
        "version":       "v131_safe_full85",
        "launched_at":   datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "worktree":      WORKTREE,
        "worktree_head": wt_head,
        "base_commit":   BASE_COMMIT,
        "cherry_picks":  list(CHERRY_PICKS),
        "description": (
            "v131 Lane A full-85: v127-safe binary (82ad51f4 + sulfo + holo data), "
            "SoA OFF, v127 protocol. Manifest paths patched for 1G9V_apo, 1TW6_holo, 1HNN expB."
        ),
        "success_gate":  ">=80/85 Hungarian (beat v109 on comparable native_seed protocol)",
        "binary":         BINARY,
        "binary_sha256":  engine_sha,
        "runner_sha256":  runner_sha,
        "matrix_md5":     matrix_md5,
        "json_pairs_src": JSON_SRC,
        "json_pairs":     json_pairs,
        "output_dir":     output,
        "cache_dir":      cache,
        "pid":            child_pid,
        "smoke_gate":     smoke,
        "reference_runs": {
            "v127_dir": REF_V127_DIR,
            "v109_record": REF_V109_DIR,
        },
        "manifest_overrides": {
            "1G9V": "fixed 1G9V_apo.pdb (chain D stripped)",
            "1TW6": "1TW6_holo.pdb (inhibitor chains stripped)",
            "1HNN": "1HNN_ligand_centered_site.pdb (expB)",
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
    }
    with open(prov, "w") as f:
        json.dump(prov_doc, f, indent=2)
        f.write("\n")

    print(f"\nv131 safe full-85 launched pid={child_pid}")
    print(f"  prov: {prov}")
    return child_pid


if __name__ == "__main__":
    main()