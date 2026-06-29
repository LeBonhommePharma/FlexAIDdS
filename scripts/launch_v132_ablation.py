#!/usr/bin/env python3
# launch_v132_ablation.py — single-step launcher for v132 oracle ablation ladder
#
# Steps (one isolated knob each vs consensus_on baseline):
#   consensus_on   — HEAD + v131 data + consensus ON + r0=4 (audit primary)
#   safe_binary    — v131_safe binary only
#   logsumexp_only — a4056163 binary only (no H-bond/VCT patch)
#   hbond_zero     — FLEXAIDDS_HBOND_WEIGHT=0 on HEAD
#
# Usage:
#   python3 scripts/launch_v132_ablation.py consensus_on
#   python3 scripts/launch_v132_ablation.py safe_binary --skip-build
#   python3 scripts/launch_v132_ablation.py --list
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from lib_launch import launch_session_isolated
from lib_worker_orders import assert_campaign_allowed
from lib_v132_ablation import (
    ENV_SNAPSHOT_KEYS,
    RESULTS,
    STRIP_ENV_KEYS,
    AblationStep,
    base_oracle_env,
    build_steps,
    step_by_id,
)

BENCH_THREADS = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_short(repo: str) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", repo, "rev-parse", "--short", "HEAD"], text=True
        ).strip()
    except Exception:
        return "unknown"


def ensure_build(step: AblationStep, skip_build: bool) -> None:
    need = [step.binary_src, step.runner_path, f"{step.data_dir}/MC_st0r5.2_6.dat"]
    if all(os.path.isfile(p) for p in need):
        return
    if skip_build:
        missing = [p for p in need if not os.path.isfile(p)]
        sys.exit(f"ERROR: --skip-build but missing: {missing}")

    if step.build_script and os.path.isfile(step.build_script):
        print(f"Building via {step.build_script}", flush=True)
        subprocess.check_call(["bash", step.build_script], cwd=str(SCRIPT_DIR.parent))
        return

    print(f"Building {step.build_dir} …", flush=True)
    subprocess.check_call(
        ["cmake", "--build", step.build_dir, "-j8", "--target", "FlexAIDdS", "benchmark_datasets"],
        cwd=str(SCRIPT_DIR.parent),
    )


def validate_paths(step: AblationStep) -> None:
    for p in (
        step.binary_src,
        step.runner_path,
        step.json_pairs,
        f"{step.data_dir}/MC_st0r5.2_6.dat",
    ):
        if not os.path.exists(p):
            sys.exit(f"ERROR: missing required path: {p}")

    if step.min_commit:
        rc = subprocess.run(
            ["git", "merge-base", "--is-ancestor", step.min_commit, "HEAD"],
            cwd=step.git_cwd,
        )
        if rc.returncode != 0:
            sys.exit(
                f"ERROR: {step.git_cwd} HEAD missing required ancestor {step.min_commit}"
            )

    native = json.loads(Path(step.json_pairs).read_text())
    assert len(native["pairs"]) == 85, "v132 ablation ladder requires full Astex-85 JSON"


def launch_step(step: AblationStep, *, skip_build: bool) -> int:
    ensure_build(step, skip_build)
    validate_paths(step)

    shutil.copy2(step.binary_src, step.binary_path)
    os.chmod(step.binary_path, 0o755)

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M")
    output = RESULTS / f"v132_{stamp}_{step.step_id}_full85"
    cache = RESULTS / step.cache_suffix
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    env = dict(os.environ)
    env.update(base_oracle_env(
        binary=step.binary_path,
        build=step.build_dir,
        data_dir=step.data_dir,
        cache=str(cache),
    ))
    env.update(step.env_overrides)
    for k in STRIP_ENV_KEYS:
        env.pop(k, None)

    cmd = [
        "caffeinate",
        "-i",
        str(step.runner_path),
        "--benchmark",
        f"crossdock_json:{step.json_pairs}",
        "--output",
        str(output),
        "--threads",
        BENCH_THREADS,
        "--temperature",
        "298",
        "--job-timeout-seconds",
        "7200",
        "--cache",
        str(cache),
        "--mode",
        "oracle-ceiling",
    ]

    print(f"\nLaunching v132 ablation step: {step.step_id}")
    print(f"  label     : {step.label}")
    print(f"  knob      : {step.ablation_knob}")
    print(f"  delta     : {step.ablation_delta_vs_consensus_on}")
    print(f"  git       : {git_short(step.git_cwd)} @ {step.git_cwd}")
    print(f"  output    : {output}")
    print(f"  threads   : {BENCH_THREADS}")

    child_pid = launch_session_isolated(cmd, env, str(output), cwd=str(SCRIPT_DIR.parent))

    prov = {
        "version": step.label,
        "ladder": "v132_oracle_ablation",
        "step_id": step.step_id,
        "launched_at": datetime.datetime.now(datetime.timezone.utc)
        .isoformat()
        .replace("+00:00", "Z"),
        "git_commit": git_short(step.git_cwd),
        "git_cwd": step.git_cwd,
        "description": step.description,
        "audit_rationale": step.audit_rationale,
        "ablation_knob": step.ablation_knob,
        "ablation_delta_vs_consensus_on": step.ablation_delta_vs_consensus_on,
        "success_gate": step.success_gate,
        "binary": step.binary_path,
        "binary_src": step.binary_src,
        "binary_sha256": sha256(Path(step.binary_path)),
        "runner_sha256": sha256(Path(step.runner_path)),
        "matrix_md5": md5(Path(f"{step.data_dir}/MC_st0r5.2_6.dat")),
        "json_pairs": step.json_pairs,
        "output_dir": str(output),
        "cache_dir": str(cache),
        "pid": child_pid,
        "reference_runs": step.reference_runs,
        "benchmark": {
            "threads": BENCH_THREADS,
            "temperature_k": 298,
            "job_timeout_seconds": 7200,
            "mode": "oracle-ceiling",
            "n_pairs": 85,
        },
        "env_snapshot": {k: env[k] for k in ENV_SNAPSHOT_KEYS if k in env},
        "audit_quote": (
            "If you want the highest honest oracle number next, the audits point toward "
            "either restoring consensus ON for oracle-ceiling campaigns, or running an "
            "ablation ladder (consensus, hbond, logsumexp, binary safe vs HEAD) rather "
            "than another single combined knob turn."
        ),
    }
    prov_path = output / "launch_provenance.json"
    prov_path.write_text(json.dumps(prov, indent=2) + "\n")

    print(f"\nv132/{step.step_id} launched pid={child_pid}")
    print(f"  prov: {prov_path}")
    return child_pid


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "step",
        nargs="?",
        help="Ablation step id (consensus_on, safe_binary, logsumexp_only, hbond_zero)",
    )
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--list", action="store_true", help="List ladder steps and exit")
    args = parser.parse_args()

    if args.list or not args.step:
        print("v132 oracle ablation ladder:")
        for step in build_steps():
            print(f"  {step.step_id:16}  {step.ablation_knob:22}  {step.label}")
        return 0

    assert_campaign_allowed("v132_ablation_step", script_name=__file__)
    launch_step(step_by_id(args.step), skip_build=args.skip_build)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())