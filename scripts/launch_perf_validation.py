#!/usr/bin/env python3
# launch_perf_validation.py — queue-safe tier-1 perf validation (paired timing + SoA gate)
#
# Runs scripts/run_perf_validation_smoke.py as a detached daemon:
#   post_p0_scalar (build/, SoA OFF) + post_p0_soa (build_soa/) + pre_p0_scalar
#
# Queue-safe: nice -n 19, OMP_NUM_THREADS=1, FLEXAIDDS_ALLOW_CONCURRENT=1,
# caffeinate -i, isolated session (lib_launch double-fork).
#
# Usage:
#   python3 scripts/launch_perf_validation.py [--skip-build]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = Path(__file__).resolve().parents[1]
BUILD_SCALAR = REPO / "build"
BUILD_SOA = REPO / "build_soa"
BUILD_PRE_P0 = Path("/Users/lp.more/.grok/worktrees/flexaidds-pre-p0/build")
RESULTS_ROOT = Path(
    os.environ.get(
        "FLEXAIDDS_PERF_VALIDATION_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/perf_swarm_validation",
    )
)
SMOKE_SCRIPT = REPO / "scripts" / "run_perf_validation_smoke.py"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_short(repo: Path) -> str:
    return subprocess.check_output(
        ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"], text=True
    ).strip()


def ensure_build(skip_build: bool) -> None:
    need = [
        BUILD_SCALAR / "FlexAIDdS",
        BUILD_SCALAR / "benchmark_datasets",
        BUILD_SOA / "FlexAIDdS",
        BUILD_SOA / "benchmark_datasets",
        BUILD_PRE_P0 / "FlexAIDdS",
        BUILD_PRE_P0 / "benchmark_datasets",
    ]
    if all(p.is_file() for p in need):
        return
    if skip_build:
        missing = [str(p) for p in need if not p.is_file()]
        sys.exit(f"ERROR: missing build artifacts: {missing}")

    for label, build_dir, targets in (
        ("scalar", BUILD_SCALAR, ("FlexAIDdS", "benchmark_datasets")),
        ("soa", BUILD_SOA, ("FlexAIDdS", "benchmark_datasets")),
    ):
        print(f"Building {label} ({build_dir}) …", flush=True)
        subprocess.check_call(
            ["cmake", "--build", str(build_dir), "-j8", "--target", *targets],
            cwd=str(REPO),
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--nice", type=int, default=19)
    args = parser.parse_args()

    ensure_build(args.skip_build)

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M")
    output = RESULTS_ROOT / f"tier1_paired_{stamp}"
    output.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env.update(
        {
            "FLEXAIDDS_ALLOW_CONCURRENT": "1",
            "FLEXAIDDS_RESTARTS": "1",
            "FLEXAIDDS_PARALLEL_RESTARTS": "0",
            "OMP_NUM_THREADS": "1",
            "OMP_PLACES": "cores",
            "OMP_PROC_BIND": "spread",
            "OMP_WAIT_POLICY": "passive",
        }
    )

    cmd = [
        "caffeinate",
        "-i",
        sys.executable,
        str(SMOKE_SCRIPT),
        "--results-root",
        str(RESULTS_ROOT),
        "--build-scalar",
        str(BUILD_SCALAR),
        "--build-soa",
        str(BUILD_SOA),
        "--build-pre-p0",
        str(BUILD_PRE_P0),
        f"--nice",
        str(args.nice),
    ]

    print(f"\nQueueing tier-1 perf validation → {output.name}", flush=True)
    print(f"  scalar:   {BUILD_SCALAR} @ {git_short(REPO)}", flush=True)
    print(f"  soa:      {BUILD_SOA}", flush=True)
    print(f"  pre_p0:   {BUILD_PRE_P0}", flush=True)
    print(f"  nice:     {args.nice}", flush=True)

    pid = launch_session_isolated(
        cmd,
        env,
        str(output),
        cwd=str(REPO),
        stdout_log=str(output / "launcher_stdout.log"),
        stderr_log=str(output / "launcher_stderr.log"),
    )

    provenance = {
        "launched_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "pid": pid,
        "campaign": output.name,
        "output_root": str(output),
        "repo": str(REPO),
        "git": {
            "commit": subprocess.check_output(
                ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
            ).strip(),
            "short": git_short(REPO),
        },
        "builds": {
            "scalar": {
                "dir": str(BUILD_SCALAR),
                "flexaidds_sha256": sha256(BUILD_SCALAR / "FlexAIDdS"),
            },
            "soa": {
                "dir": str(BUILD_SOA),
                "flexaidds_sha256": sha256(BUILD_SOA / "FlexAIDdS"),
            },
            "pre_p0": {
                "dir": str(BUILD_PRE_P0),
                "flexaidds_sha256": sha256(BUILD_PRE_P0 / "FlexAIDdS"),
            },
        },
        "queue_safe": {
            "nice": args.nice,
            "omp_threads": 1,
            "bench_threads": 1,
            "allow_concurrent": True,
            "caffeinate": True,
        },
        "benchmark_json": str(REPO / "benchmarks/perf_swarm/tier1_paired_5.json"),
        "expected_artifacts": [
            "post_p0_scalar/",
            "post_p0_soa/",
            "pre_p0_scalar/",
            "validation_report.json",
        ],
    }
    prov_path = output / "launch_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"Queued pid={pid}", flush=True)
    print(f"Provenance: {prov_path}", flush=True)
    print(f"Monitor: tail -f {output}/launcher_stderr.log", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())