#!/usr/bin/env python3
# launch_vcontacts_commit_bisect_smoke.py — smoke-12 for single-commit revert variants
#
# Usage:
#   python3 scripts/launch_vcontacts_commit_bisect_smoke.py revert_f9c80fe5
#   python3 scripts/launch_vcontacts_commit_bisect_smoke.py revert_d2295cf0 --skip-build
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
import tempfile
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from launch_vcontacts_bisect_smoke import (
    REGRESSION_GUARD,
    harvest_report,
    sha256,
)
from lib_launch import launch_session_isolated
from v131_safe_common import (
    git_root,
    patch_manifest,
    resolve_oracle_dir,
    resolve_worktree,
    scrub_env,
    validate_lane_a_assets,
    v127_protocol_env,
)

GIT_ROOT = Path(os.environ.get("FLEXAIDDS_GIT_ROOT", git_root(str(SCRIPT_DIR))))
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
SMOKE_JSON = GIT_ROOT / "benchmarks/datasets/benchmark_astex_smoke_12_v131.json"
BUILD_SCRIPT = SCRIPT_DIR / "build_vcontacts_commit_bisect.sh"

VARIANTS = {
    "revert_f9c80fe5": {
        "sha": "f9c80fe5",
        "label": "bisect_revert_f9c80fe5",
        "description": "HEAD (27e68e51 reverted) minus f9c80fe5 SoA sqrdist parity",
    },
    "revert_d2295cf0": {
        "sha": "d2295cf0",
        "label": "bisect_revert_d2295cf0",
        "description": "HEAD minus d2295cf0 thread-local coord cache atoms",
    },
    "revert_d4d68592": {
        "sha": "d4d68592",
        "label": "bisect_revert_d4d68592",
        "description": "HEAD minus d4d68592 PR4 scalar-identical Vcontacts loop",
    },
}


def paths(variant: str) -> tuple[str, str, str, str]:
    wt = str(GIT_ROOT.parent / f"FlexAIDdS_{variant}")
    build = f"{wt}/build_bisect"
    main_build = str(GIT_ROOT / "build_lto")
    return (
        f"{build}/FlexAIDdS",
        f"{main_build}/benchmark_datasets",
        f"/tmp/FlexAIDdS_{variant}",
        wt,
    )


def ensure_build(variant: str, skip_build: bool) -> None:
    binary, runner, _, _ = paths(variant)
    if os.path.isfile(binary) and os.path.isfile(runner):
        return
    if skip_build:
        sys.exit(f"ERROR: --skip-build but missing {binary} or {runner}")
    env = {**os.environ, "FLEXAIDDS_GIT_ROOT": str(GIT_ROOT)}
    subprocess.check_call(["bash", str(BUILD_SCRIPT), variant], env=env)


def launch(variant: str, *, skip_build: bool) -> int:
    meta = VARIANTS[variant]
    ensure_build(variant, skip_build)
    binary_src, runner, stamped, wt = paths(variant)
    validate_lane_a_assets(resolve_worktree(), GIT_ROOT)

    shutil.copy2(binary_src, stamped)
    os.chmod(stamped, 0o755)

    manifest = patch_manifest(
        json.loads(SMOKE_JSON.read_text()), resolve_worktree(), GIT_ROOT
    )
    with tempfile.NamedTemporaryFile(
        mode="w", suffix="_commit_bisect_smoke12.json", delete=False
    ) as tmp:
        json.dump(manifest, tmp, indent=2)
        tmp.write("\n")
        json_path = tmp.name

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M")
    output = RESULTS / f"vcontacts_commit_bisect_{stamp}_{variant}"
    cache = RESULTS / f"cache_vcontacts_commit_bisect_{variant}"
    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    oracle_dir = resolve_oracle_dir(resolve_worktree(), GIT_ROOT)
    data_dir = str(GIT_ROOT / "build_lto")
    env = scrub_env(dict(os.environ))
    env.update(v127_protocol_env(stamped, data_dir, str(cache), oracle_dir))

    cmd = [
        "caffeinate", "-i",
        runner,
        "--benchmark", f"crossdock_json:{json_path}",
        "--output", str(output),
        "--threads", os.environ.get("FLEXAIDDS_BENCH_THREADS", "1"),
        "--temperature", "298",
        "--job-timeout-seconds", "7200",
        "--cache", str(cache),
        "--mode", "oracle-ceiling",
    ]

    print(f"\nLaunching commit bisect smoke: {variant}")
    print(f"  revert: {meta['sha']}")
    print(f"  {meta['description']}")
    print(f"  output: {output}")

    pid = launch_session_isolated(cmd, env, str(output), cwd=str(GIT_ROOT))
    prov = {
        "bisect": "vcontacts_commit_revert",
        "variant": variant,
        "reverted_sha": meta["sha"],
        "label": meta["label"],
        "description": meta["description"],
        "launched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "worktree": wt,
        "binary": stamped,
        "binary_sha256": sha256(Path(stamped)),
        "output_dir": str(output),
        "regression_guard": list(REGRESSION_GUARD),
        "pid": pid,
    }
    (output / "launch_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    print(f"pid={pid}")
    return pid


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("variant", choices=sorted(VARIANTS))
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--harvest", type=Path)
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