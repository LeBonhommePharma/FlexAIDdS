#!/usr/bin/env python3
"""Deterministic SoA vs scalar tier-1 accuracy gate (paired RMSD delta <= 0.05 A)."""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_BENCH = REPO / "benchmarks/perf_swarm/tier1_paired_5.json"
GATE_TOL = 0.05


def git_info() -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"], text=True
        ).strip()
        return {"commit": commit}
    except Exception:
        return {"commit": "unknown"}


def load_results(out_dir: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    csv_candidates = sorted(out_dir.glob("astex_*_results.csv"))
    if not csv_candidates:
        return out
    with csv_candidates[0].open(newline="") as fh:
        for row in csv.DictReader(fh):
            tid = str(row.get("pdb_id", "")).upper()
            if not tid:
                continue
            try:
                rmsd = float(row.get("rmsd_to_crystal", "nan"))
            except ValueError:
                rmsd = float("nan")
            try:
                nposes = int(float(row.get("num_poses", "0")))
            except ValueError:
                nposes = 0
            out[tid] = {
                "rmsd": rmsd,
                "num_poses": nposes,
                "pose_source": row.get("pose_source", ""),
                "success": row.get("success", ""),
            }
    return out


def run_benchmark(
    *,
    label: str,
    build_dir: Path,
    out_root: Path,
    bench_json: Path,
    seed: int,
    nice_level: int,
) -> Path:
    runner = build_dir / "benchmark_datasets"
    flexaid = build_dir / "FlexAIDdS"
    if not runner.is_file():
        raise FileNotFoundError(f"missing benchmark_datasets in {build_dir}")
    if not flexaid.is_file():
        raise FileNotFoundError(f"missing FlexAIDdS in {build_dir}")

    out_dir = out_root / label
    out_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["FLEXAIDDS_BUILD"] = str(build_dir)
    env["FLEXAIDDS_BINARY"] = str(flexaid)
    env["FLEXAIDDS_REPO"] = str(REPO)
    env["FLEXAIDDS_RESTARTS"] = "1"
    env["FLEXAIDDS_PARALLEL_RESTARTS"] = "0"
    # gaboom.cpp: when ga.seed==0 in dock JSON, FLEXAID_SEED sets srand/mt19937
    env["FLEXAID_SEED"] = str(seed)
    env["OMP_NUM_THREADS"] = "1"
    env["OMP_PLACES"] = "cores"
    env["OMP_PROC_BIND"] = "spread"
    env["OMP_WAIT_POLICY"] = "passive"

    cmd = [
        "nice",
        f"-n{nice_level}",
        str(runner),
        "--benchmark",
        f"crossdock_json:{bench_json}",
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
    print(f"[gate] {label} seed={seed}: {' '.join(cmd)}", flush=True)
    proc = subprocess.run(cmd, env=env, cwd=str(REPO))
    if proc.returncode != 0:
        raise RuntimeError(f"{label} exited {proc.returncode}")
    return out_dir


def compare_pair(scalar_dir: Path, soa_dir: Path) -> dict:
    scalar = load_results(scalar_dir)
    soa = load_results(soa_dir)
    common = sorted(set(scalar) & set(soa))
    rows = []
    max_delta = 0.0
    complete_pairs = 0

    for tid in common:
        s = scalar[tid]
        o = soa[tid]
        rs, ro = s["rmsd"], o["rmsd"]
        complete = (
            rs == rs
            and ro == ro
            and rs >= 0.0
            and ro >= 0.0
            and s["num_poses"] > 0
            and o["num_poses"] > 0
        )
        if complete:
            complete_pairs += 1
            d = abs(ro - rs)
            max_delta = max(max_delta, d)
        else:
            d = float("nan")
        rows.append(
            {
                "target": tid,
                "scalar_rmsd": rs,
                "soa_rmsd": ro,
                "rmsd_delta": d if d == d else None,
                "complete": complete,
                "scalar_poses": s["num_poses"],
                "soa_poses": o["num_poses"],
            }
        )

    gate_pass = complete_pairs > 0 and max_delta <= GATE_TOL
    return {
        "paired_targets": len(common),
        "complete_pairs": complete_pairs,
        "max_rmsd_delta": round(max_delta, 4),
        "gate_pass": gate_pass,
        "tolerance_angstrom": GATE_TOL,
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build-scalar", type=Path, default=REPO / "build")
    parser.add_argument("--build-soa", type=Path, default=REPO / "build_soa")
    parser.add_argument("--benchmark", type=Path, default=DEFAULT_BENCH)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--nice", type=int, default=19)
    args = parser.parse_args()

    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M")
    out_root = args.out_root or (REPO / "results/perf_swarm" / f"soa_gate_{stamp}")
    out_root.mkdir(parents=True, exist_ok=True)

    repeat_reports = []
    overall_max = 0.0
    all_complete = True

    for rep in range(args.repeats):
        seed = args.seed + rep
        rep_dir = out_root / f"repeat_{rep}_seed{seed}"
        scalar_dir = run_benchmark(
            label="scalar",
            build_dir=args.build_scalar,
            out_root=rep_dir,
            bench_json=args.benchmark,
            seed=seed,
            nice_level=args.nice,
        )
        soa_dir = run_benchmark(
            label="soa",
            build_dir=args.build_soa,
            out_root=rep_dir,
            bench_json=args.benchmark,
            seed=seed,
            nice_level=args.nice,
        )
        gate = compare_pair(scalar_dir, soa_dir)
        repeat_reports.append({"seed": seed, "output": str(rep_dir), "gate": gate})
        overall_max = max(overall_max, gate["max_rmsd_delta"])
        if gate["complete_pairs"] < gate["paired_targets"]:
            all_complete = False

    report = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "git": git_info(),
        "benchmark_json": str(args.benchmark),
        "build_scalar": str(args.build_scalar),
        "build_soa": str(args.build_soa),
        "base_seed": args.seed,
        "repeats": args.repeats,
        "output_root": str(out_root),
        "overall_max_rmsd_delta": round(overall_max, 4),
        "gate_pass": all_complete and overall_max <= GATE_TOL,
        "tolerance_angstrom": GATE_TOL,
        "repeats_detail": repeat_reports,
    }

    report_path = out_root / "soa_parity_gate_report.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")

    print(
        f"\nSoA parity gate: {'PASS' if report['gate_pass'] else 'FAIL'} "
        f"(overall max Δ={report['overall_max_rmsd_delta']} Å, repeats={args.repeats})"
    )
    print(f"Wrote {report_path}")
    return 0 if report["gate_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())