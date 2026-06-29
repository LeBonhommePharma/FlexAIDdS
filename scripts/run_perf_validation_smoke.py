#!/usr/bin/env python3
"""Queue-safe tier-1 validation: SoA accuracy gate + paired pre/post-P0 timing."""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BENCH_JSON = REPO / "benchmarks/perf_swarm/tier1_paired_5.json"
TIMING_RE = re.compile(
    r"TIMING SUMMARY:\s+\d+\s+gens timed,\s+avg\s+([\d.]+)\s+ms/gen"
)
PRE_P0_SHA = "04ff1735"  # parent of Wave-1 P0 (27e68e51)


def git_info() -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
        ).strip()
        return {"commit": commit}
    except Exception:
        return {"commit": "unknown"}


def harvest_run(out_dir: Path) -> list[dict]:
    rows: list[dict] = []
    if not out_dir.is_dir():
        return rows
    for log in out_dir.rglob("stderr.log"):
        try:
            text = log.read_text(errors="replace")
        except OSError:
            continue
        m = TIMING_RE.search(text)
        if not m:
            continue
        target = log.parent.name
        rows.append(
            {
                "target": target,
                "avg_ms_per_gen": float(m.group(1)),
                "log_path": str(log),
            }
        )
    rows.sort(key=lambda r: r["target"])
    return rows


def load_results(out_dir: Path) -> dict[str, dict]:
    """Load per-target RMSD from benchmark_datasets CSV output."""
    out: dict[str, dict] = {}
    csv_candidates = sorted(out_dir.glob("astex_*_results.csv"))
    if not csv_candidates:
        return out
    import csv

    with csv_candidates[0].open(newline="") as fh:
        for row in csv.DictReader(fh):
            tid = str(row.get("pdb_id", "")).upper()
            if not tid:
                continue
            try:
                rmsd = float(row.get("rmsd_to_crystal", "nan"))
            except ValueError:
                rmsd = float("nan")
            out[tid] = {
                "rmsd_best": rmsd,
                "selection_method": row.get("pose_source", ""),
                "success": row.get("success", ""),
            }
    return out


def run_benchmark(
    *,
    label: str,
    build_dir: Path,
    out_root: Path,
    nice_level: int,
) -> Path:
    runner = build_dir / "benchmark_datasets"
    if not runner.is_file():
        raise FileNotFoundError(f"missing benchmark_datasets in {build_dir}")

    out_dir = out_root / label
    out_dir.mkdir(parents=True, exist_ok=True)
    flexaid_bin = build_dir / "FlexAIDdS"
    if not flexaid_bin.is_file():
        raise FileNotFoundError(f"missing FlexAIDdS in {build_dir}")

    env = os.environ.copy()
    env["FLEXAIDDS_BUILD"] = str(build_dir)
    env["FLEXAIDDS_BINARY"] = str(flexaid_bin)
    env["FLEXAIDDS_REPO"] = str(REPO)
    env["FLEXAIDDS_RESTARTS"] = "1"
    env["FLEXAIDDS_PARALLEL_RESTARTS"] = "0"
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_PLACES"] = "cores"
    env["OMP_PROC_BIND"] = "spread"
    env["OMP_WAIT_POLICY"] = "passive"

    cmd = [
        "nice",
        f"-n{nice_level}",
        str(runner),
        "--benchmark",
        f"crossdock_json:{BENCH_JSON}",
        "--output",
        str(out_dir),
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
    print(f"[run] {label}: {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, cwd=str(REPO))
    if proc.returncode != 0:
        raise RuntimeError(f"{label} benchmark_datasets exited {proc.returncode}")
    return out_dir


def compare_accuracy(scalar_dir: Path, soa_dir: Path) -> dict:
    scalar = load_results(scalar_dir)
    soa = load_results(soa_dir)
    common = sorted(set(scalar) & set(soa))
    rows = []
    max_rmsd_delta = 0.0
    for tid in common:
        s = scalar[tid]
        o = soa[tid]
        r_s = float(s.get("rmsd_best", float("nan")))
        r_o = float(o.get("rmsd_best", float("nan")))
        d = abs(r_o - r_s) if r_s == r_s and r_o == r_o else float("nan")
        if d == d:
            max_rmsd_delta = max(max_rmsd_delta, d)
        rows.append(
            {
                "target": tid,
                "scalar_rmsd": r_s,
                "soa_rmsd": r_o,
                "rmsd_delta": d,
                "scalar_method": s.get("selection_method"),
                "soa_method": o.get("selection_method"),
            }
        )
    gate_pass = max_rmsd_delta <= 0.05 if rows else False
    return {
        "paired_targets": len(common),
        "max_rmsd_delta": round(max_rmsd_delta, 4),
        "gate_pass": gate_pass,
        "tolerance_angstrom": 0.05,
        "rows": rows,
    }


def compare_paired_timing(pre_dir: Path, post_dir: Path) -> dict:
    pre = {r["target"]: r["avg_ms_per_gen"] for r in harvest_run(pre_dir)}
    post = {r["target"]: r["avg_ms_per_gen"] for r in harvest_run(post_dir)}
    common = sorted(set(pre) & set(post))
    rows = []
    speedups = []
    for tid in common:
        b, c = pre[tid], post[tid]
        sp = b / c if c > 0 else 1.0
        speedups.append(sp)
        rows.append(
            {
                "target": tid,
                "pre_p0_ms_per_gen": b,
                "post_p0_ms_per_gen": c,
                "speedup": round(sp, 3),
                "delta_pct": round(100.0 * (c - b) / b, 2),
            }
        )
    rows.sort(key=lambda r: r["speedup"], reverse=True)
    median_sp = sorted(speedups)[len(speedups) // 2] if speedups else 1.0
    return {
        "paired_targets": len(common),
        "median_speedup": round(median_sp, 3),
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/perf_swarm_validation"),
    )
    parser.add_argument(
        "--build-scalar",
        type=Path,
        default=Path("/Users/lp.more/Projects/FlexAIDdS/build_lto"),
        help="Post-P0 scalar binary dir (default: production build_lto)",
    )
    parser.add_argument("--build-soa", type=Path, default=REPO / "build_soa")
    parser.add_argument(
        "--build-pre-p0",
        type=Path,
        default=Path("/Users/lp.more/.grok/worktrees/flexaidds-pre-p0/build"),
    )
    parser.add_argument("--nice", type=int, default=19)
    parser.add_argument(
        "--skip-run",
        action="store_true",
        help="Only compare existing result directories",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=None,
        help="Existing run directory (required with --skip-run)",
    )
    args = parser.parse_args()

    if args.skip_run:
        if not args.output_root:
            parser.error("--output-root is required when --skip-run is set")
        out_root = args.output_root
    else:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M")
        out_root = args.results_root / f"tier1_paired_{stamp}"
        out_root.mkdir(parents=True, exist_ok=True)

    if not args.skip_run:
        scalar_dir = run_benchmark(
            label="post_p0_scalar",
            build_dir=args.build_scalar,
            out_root=out_root,
            nice_level=args.nice,
        )
        soa_dir = run_benchmark(
            label="post_p0_soa",
            build_dir=args.build_soa,
            out_root=out_root,
            nice_level=args.nice,
        )
        pre_dir = run_benchmark(
            label="pre_p0_scalar",
            build_dir=args.build_pre_p0,
            out_root=out_root,
            nice_level=args.nice,
        )
    else:
        scalar_dir = out_root / "post_p0_scalar"
        soa_dir = out_root / "post_p0_soa"
        pre_dir = out_root / "pre_p0_scalar"

    report = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "git": git_info(),
        "pre_p0_sha": PRE_P0_SHA,
        "benchmark_json": str(BENCH_JSON),
        "output_root": str(out_root),
        "soa_accuracy_gate": compare_accuracy(scalar_dir, soa_dir),
        "paired_timing": compare_paired_timing(pre_dir, scalar_dir),
        "timings": {
            "pre_p0": harvest_run(pre_dir),
            "post_scalar": harvest_run(scalar_dir),
            "post_soa": harvest_run(soa_dir),
        },
    }

    report_path = out_root / "validation_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    soa = report["soa_accuracy_gate"]
    timing = report["paired_timing"]
    print(f"\nSoA gate: {'PASS' if soa['gate_pass'] else 'FAIL'} "
          f"(max RMSD delta {soa['max_rmsd_delta']} Å, n={soa['paired_targets']})")
    print(f"Paired timing: median {timing['median_speedup']}x "
          f"(n={timing['paired_targets']})")
    print(f"Wrote {report_path}")
    return 0 if soa["gate_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())