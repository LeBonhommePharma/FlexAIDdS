#!/usr/bin/env python3
# launch_vcontacts_bisect_smoke.py — smoke-12 for one Vcontacts bisect variant
#
# Variants:
#   safe         — v131_safe worktree (82ad51f4 + sulfo + holo, SoA OFF) [control]
#   head_soa_off — HEAD build_bisect_soa_off
#   head_soa_on  — HEAD build_bisect_soa_on (suspect: default SoA ON)
#
# Usage:
#   python3 scripts/launch_vcontacts_bisect_smoke.py safe
#   python3 scripts/launch_vcontacts_bisect_smoke.py head_soa_off --skip-build
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

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

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
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

# Pin perfornance-swarm worktree — never inherit caller cwd / Projects/FlexAIDdS.
GIT_ROOT = Path(
    os.environ.get("FLEXAIDDS_GIT_ROOT", git_root(str(SCRIPT_DIR)))
)
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
SMOKE_JSON = GIT_ROOT / "benchmarks/datasets/benchmark_astex_smoke_12_v131.json"
REGRESSION_GUARD = ("1HQ2", "1S3V", "1T40")
BUILD_BISECT_SCRIPT = SCRIPT_DIR / "build_vcontacts_bisect.sh"
BUILD_SAFE_SCRIPT = SCRIPT_DIR / "build_v131_safe.sh"

VARIANTS = {
    "safe": {
        "label": "bisect_safe_control",
        "description": "v131_safe @ 82ad51f4+sulfo+holo, pre-27e68e51 Vcontacts/gaboom",
        "suspect": "control — expect nearest v127 78/85 behaviour",
    },
    "head_soa_off": {
        "label": "bisect_head_soa_off",
        "description": "HEAD binary, FLEXAIDS_USE_SOA_DISTANCES=OFF at compile time",
        "suspect": "isolates SoA compile flag; inv_d12 + parallel reproduce still active",
    },
    "head_soa_on": {
        "label": "bisect_head_soa_on",
        "description": "HEAD binary, FLEXAIDS_USE_SOA_DISTANCES=ON (production default)",
        "suspect": "Wave-2 suspect: Vcontacts SoA hot path + PR4 parity commits",
    },
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_short(cwd: str) -> str:
    return subprocess.check_output(
        ["git", "-C", cwd, "rev-parse", "--short", "HEAD"], text=True
    ).strip()


def resolve_paths(variant: str) -> tuple[str, str, str, str]:
    if variant == "safe":
        wt = resolve_worktree()
        build = f"{wt}/build_lto"
        binary_src = f"{build}/FlexAIDdS"
        runner = f"{build}/benchmark_datasets"
        stamped = "/tmp/FlexAIDdS_bisect_safe"
        git_cwd = wt
    else:
        build = str(GIT_ROOT / f"build_bisect_{'soa_off' if variant == 'head_soa_off' else 'soa_on'}")
        binary_src = f"{build}/FlexAIDdS"
        runner = f"{build}/benchmark_datasets"
        stamped = f"/tmp/FlexAIDdS_bisect_{variant}"
        git_cwd = str(GIT_ROOT)
    return binary_src, runner, stamped, git_cwd


def ensure_build(variant: str, skip_build: bool) -> None:
    binary_src, runner, _, _ = resolve_paths(variant)
    need = [binary_src, runner]
    if all(os.path.isfile(p) for p in need):
        return
    if skip_build:
        sys.exit(f"ERROR: --skip-build but missing {need}")

    build_env = {**os.environ, "FLEXAIDDS_GIT_ROOT": str(GIT_ROOT)}
    if variant == "safe":
        subprocess.check_call(
            ["bash", str(BUILD_SAFE_SCRIPT)], cwd=str(GIT_ROOT), env=build_env
        )
    else:
        subprocess.check_call(
            ["bash", str(BUILD_BISECT_SCRIPT)],
            cwd=str(GIT_ROOT),
            env=build_env,
        )


def patched_smoke_json(worktree: str) -> str:
    native = patch_manifest(json.loads(SMOKE_JSON.read_text()), worktree, GIT_ROOT)
    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_bisect_smoke12.json", delete=False
    ) as tmp:
        json.dump(native, tmp, indent=2)
        tmp.write("\n")
        return tmp.name


def launch(variant: str, *, skip_build: bool) -> int:
    meta = VARIANTS[variant]
    ensure_build(variant, skip_build)
    binary_src, runner, stamped, git_cwd = resolve_paths(variant)
    worktree = resolve_worktree() if variant == "safe" else str(GIT_ROOT)

    validate_lane_a_assets(worktree, GIT_ROOT)

    shutil.copy2(binary_src, stamped)
    os.chmod(stamped, 0o755)

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M")
    output = RESULTS / f"vcontacts_bisect_{stamp}_{variant}"
    cache = RESULTS / f"cache_vcontacts_bisect_{variant}"
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    json_path = patched_smoke_json(worktree)
    oracle_dir = resolve_oracle_dir(worktree, GIT_ROOT)
    data_dir = str(Path(binary_src).parent)
    env = scrub_env(dict(os.environ))
    env.update(v127_protocol_env(stamped, data_dir, str(cache), oracle_dir))

    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")
    cmd = [
        "caffeinate", "-i",
        runner,
        "--benchmark", f"crossdock_json:{json_path}",
        "--output", str(output),
        "--threads", bench_threads,
        "--temperature", "298",
        "--job-timeout-seconds", "7200",
        "--cache", str(cache),
        "--mode", "oracle-ceiling",
    ]

    print(f"\nLaunching Vcontacts bisect smoke: {variant}")
    print(f"  {meta['description']}")
    print(f"  suspect: {meta['suspect']}")
    print(f"  output:  {output}")

    pid = launch_session_isolated(cmd, env, str(output), cwd=str(GIT_ROOT))

    prov = {
        "bisect": "vcontacts_wave2",
        "variant": variant,
        "label": meta["label"],
        "description": meta["description"],
        "suspect": meta["suspect"],
        "launched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_commit": git_short(git_cwd),
        "binary": stamped,
        "binary_sha256": sha256(Path(stamped)),
        "output_dir": str(output),
        "cache_dir": str(cache),
        "pid": pid,
        "reference": {
            "v109_record": "80/85",
            "v127_baseline": "78/85",
            "v130_observed": "73/85",
            "priority": "bisect Vcontacts before v132 ablation or new knob turns",
        },
        "regression_guard": list(REGRESSION_GUARD),
    }
    (output / "launch_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"pid={pid} prov={output / 'launch_provenance.json'}")
    return pid


def harvest_report(output: Path) -> dict:
    rows = []
    for rf in sorted(output.glob("*/result.csv")):
        tid = rf.parent.name
        r = next(csv.DictReader(rf.open()))
        ok = r.get("success") == "1"
        rows.append({
            "target": tid,
            "success": ok,
            "rmsd_hungarian": r.get("rmsd_hungarian"),
            "pose_source": r.get("pose_source"),
        })
    guard_fail = [
        x["target"] for x in rows
        if x["target"] in REGRESSION_GUARD and not x["success"]
    ]
    return {
        "n_complete": len(rows),
        "n_success": sum(1 for x in rows if x["success"]),
        "regression_guard_fail": guard_fail,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=sorted(VARIANTS))
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--harvest", type=Path, help="Write bisect_smoke_report.json for output dir")
    args = parser.parse_args()

    if args.harvest:
        report = harvest_report(args.harvest)
        path = args.harvest / "bisect_smoke_report.json"
        path.write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2))
        return 0

    launch(args.variant, skip_build=args.skip_build)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())