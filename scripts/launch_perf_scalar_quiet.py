#!/usr/bin/env python3
# launch_perf_scalar_quiet.py — scalar-only tier-1 timing (quiet window)
#
# Single pass: post_p0_scalar on tier1_paired_5 (5 targets), SoA OFF.
# Intended to run when no other benchmark_datasets jobs are active.
#
# Usage:
#   python3 scripts/launch_perf_scalar_quiet.py [--skip-build]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated

REPO = Path(__file__).resolve().parents[1]
BUILD_SCALAR = REPO / "build"
BENCH_JSON = REPO / "benchmarks/perf_swarm/tier1_paired_5.json"
RESULTS_ROOT = Path(
    os.environ.get(
        "FLEXAIDDS_PERF_VALIDATION_ROOT",
        "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/perf_swarm_validation",
    )
)
TIMING_RE = re.compile(
    r"TIMING SUMMARY:\s+\d+\s+gens timed,\s+avg\s+([\d.]+)\s+ms/gen"
)
REF_REPORT = RESULTS_ROOT / "tier1_paired_20260629_0854" / "validation_report.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_info(repo: Path) -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
        short = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"], text=True
        ).strip()
        return {"commit": commit, "short": short}
    except Exception:
        return {"commit": "unknown", "short": "unknown"}


def ensure_build(skip_build: bool) -> None:
    need = [BUILD_SCALAR / "FlexAIDdS", BUILD_SCALAR / "benchmark_datasets"]
    if all(p.is_file() for p in need):
        return
    if skip_build:
        missing = [str(p) for p in need if not p.is_file()]
        sys.exit(f"ERROR: missing build artifacts: {missing}")
    subprocess.check_call(
        ["cmake", "--build", str(BUILD_SCALAR), "-j8", "--target", "FlexAIDdS", "benchmark_datasets"],
        cwd=str(REPO),
    )


def harvest_timings(out_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for log in sorted(out_dir.rglob("stderr.log")):
        try:
            text = log.read_text(errors="replace")
        except OSError:
            continue
        m = TIMING_RE.search(text)
        if not m:
            continue
        rows.append(
            {
                "target": log.parent.name,
                "avg_ms_per_gen": float(m.group(1)),
                "log_path": str(log),
            }
        )
    rows.sort(key=lambda r: r["target"])
    return rows


def compare_to_reference(timings: list[dict]) -> dict | None:
    if not REF_REPORT.is_file():
        return None
    ref = json.loads(REF_REPORT.read_text())
    ref_map = {
        r["target"]: r["avg_ms_per_gen"]
        for r in ref.get("timings", {}).get("post_scalar", [])
    }
    rows = []
    deltas = []
    for row in timings:
        tgt = row["target"]
        if tgt not in ref_map:
            continue
        ref_ms = ref_map[tgt]
        cur_ms = row["avg_ms_per_gen"]
        delta_pct = 100.0 * (cur_ms - ref_ms) / ref_ms if ref_ms > 0 else 0.0
        deltas.append(delta_pct)
        rows.append(
            {
                "target": tgt,
                "ref_0854_ms_per_gen": ref_ms,
                "current_ms_per_gen": cur_ms,
                "delta_pct": round(delta_pct, 2),
            }
        )
    return {
        "reference": str(REF_REPORT),
        "paired_targets": len(rows),
        "median_delta_pct": round(sorted(deltas)[len(deltas) // 2], 2) if deltas else None,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--nice", type=int, default=19)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Override output directory (default: tier1_scalar_quiet_<stamp>)",
    )
    args = parser.parse_args()

    ensure_build(args.skip_build)

    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%d_%H%M")
    output = args.output_root or (RESULTS_ROOT / f"tier1_scalar_quiet_{stamp}")
    label_dir = output / "post_p0_scalar"
    label_dir.mkdir(parents=True, exist_ok=True)

    runner = BUILD_SCALAR / "benchmark_datasets"
    flexaid_bin = BUILD_SCALAR / "FlexAIDdS"

    env = os.environ.copy()
    env.update(
        {
            "FLEXAIDDS_BUILD": str(BUILD_SCALAR),
            "FLEXAIDDS_BINARY": str(flexaid_bin),
            "FLEXAIDDS_REPO": str(REPO),
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
        "nice",
        f"-n{args.nice}",
        str(runner),
        "--benchmark",
        f"crossdock_json:{BENCH_JSON}",
        "--output",
        str(label_dir),
        "--threads",
        "1",
        "--omp-threads",
        "1",
        "--ga-generations",
        "88",
        "--ga-population",
        "100",
        "--job-timeout-seconds",
        "1800",
        "--mode",
        "oracle-ceiling",
        "--force",
    ]

    print(f"\nQueueing scalar-only quiet tier-1 → {output.name}", flush=True)
    print(f"  build: {BUILD_SCALAR} @ {git_info(REPO)['short']}", flush=True)
    print(f"  nice:  {args.nice}", flush=True)

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
        "mode": "scalar_only_quiet",
        "repo": str(REPO),
        "git": git_info(REPO),
        "build_scalar": {
            "dir": str(BUILD_SCALAR),
            "flexaidds_sha256": sha256(flexaid_bin),
        },
        "queue_safe": {
            "nice": args.nice,
            "omp_threads": 1,
            "bench_threads": 1,
            "allow_concurrent": False,
            "caffeinate": True,
        },
        "benchmark_json": str(BENCH_JSON),
        "reference_report": str(REF_REPORT),
    }
    prov_path = output / "launch_provenance.json"
    prov_path.write_text(json.dumps(provenance, indent=2) + "\n")

    print(f"Queued pid={pid}", flush=True)
    print(f"Provenance: {prov_path}", flush=True)
    print(f"Monitor: tail -f {output}/launcher_stderr.log", flush=True)
    return 0


def write_timing_report(output: Path) -> Path:
    """Post-run helper: harvest timings and write timing_report.json."""
    label_dir = output / "post_p0_scalar"
    timings = harvest_timings(label_dir)
    ms_vals = [r["avg_ms_per_gen"] for r in timings]
    median_ms = sorted(ms_vals)[len(ms_vals) // 2] if ms_vals else None
    report = {
        "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "output_root": str(output),
        "git": git_info(REPO),
        "timings": timings,
        "median_ms_per_gen": median_ms,
        "n_targets": len(timings),
        "vs_reference_0854": compare_to_reference(timings),
    }
    path = output / "timing_report.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    return path


if __name__ == "__main__":
    raise SystemExit(main())