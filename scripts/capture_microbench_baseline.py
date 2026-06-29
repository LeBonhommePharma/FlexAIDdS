#!/usr/bin/env python3
"""Capture microbenchmark text output into a perf baseline JSON artifact."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/perf_swarm"

TENCOM_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
)
VCFBATCH_SPEEDUP_RE = re.compile(r"Speedup\s*:\s*([\d.]+)", re.I)


def git_info(repo: Path) -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"], text=True
        ).strip()
        dirty = bool(
            subprocess.check_output(
                ["git", "-C", str(repo), "status", "--porcelain"], text=True
            ).strip()
        )
        branch = subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "--abbrev-ref", "HEAD"], text=True
        ).strip()
        return {"commit": commit, "dirty": dirty, "branch": branch}
    except Exception as exc:  # noqa: BLE001
        return {"commit": "unknown", "dirty": None, "branch": "unknown", "error": str(exc)}


def parse_tencom(text: str, *, reference_n_res: int = 200) -> dict[str, float]:
    rows: dict[int, tuple[float, float]] = {}
    for line in text.splitlines():
        match = TENCOM_ROW_RE.match(line)
        if not match:
            continue
        n_res = int(match.group(1))
        rows[n_res] = (float(match.group(2)), float(match.group(4)))

    if reference_n_res in rows:
        build_ms, sample_ms = rows[reference_n_res]
    elif rows:
        n_res = sorted(rows)[len(rows) // 2]
        build_ms, sample_ms = rows[n_res]
        reference_n_res = n_res
    else:
        return {}

    return {
        "reference_n_res": reference_n_res,
        "build_ms_full": build_ms,
        "sample_ms_full": sample_ms,
    }


def parse_vcfbatch(text: str) -> dict[str, float]:
    match = VCFBATCH_SPEEDUP_RE.search(text)
    if not match:
        return {}
    return {"speedup_vs_scalar": float(match.group(1))}


def load_existing(path: Path) -> dict:
    if path.is_file():
        return json.loads(path.read_text())
    return {}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tencom", type=Path, help="benchmark_tencom stdout capture")
    parser.add_argument("--vcfbatch", type=Path, help="benchmark_vcfbatch stdout capture")
    parser.add_argument("--label", type=str, default="linux_cpu")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--merge",
        type=Path,
        default=None,
        help="Existing baseline JSON to merge (preserves dock_timings_harvested)",
    )
    parser.add_argument("--reference-n-res", type=int, default=200)
    args = parser.parse_args()

    merge_path = args.merge or (args.out / f"baseline_{args.label}.json")
    payload = load_existing(merge_path)

    benchmarks: list[dict] = []
    if args.tencom and args.tencom.is_file():
        tencom = parse_tencom(
            args.tencom.read_text(errors="replace"),
            reference_n_res=args.reference_n_res,
        )
        if tencom:
            benchmarks.append(
                {
                    "name": "tencom",
                    "reference_n_res": tencom.pop("reference_n_res"),
                    "metrics": {
                        "build_ms_full": tencom["build_ms_full"],
                        "sample_ms_full": tencom["sample_ms_full"],
                    },
                    "status": "captured",
                }
            )

    if args.vcfbatch and args.vcfbatch.is_file():
        vcf = parse_vcfbatch(args.vcfbatch.read_text(errors="replace"))
        if vcf:
            benchmarks.append(
                {
                    "name": "vcfbatch",
                    "args": [200, 20],
                    "metrics": vcf,
                    "status": "captured",
                }
            )

    if not benchmarks:
        print("No benchmark metrics parsed; nothing written.")
        return 1

    payload.update(
        {
            "schema_version": "1.0.0",
            "recorded_at": datetime.now(timezone.utc).isoformat(),
            "git": git_info(REPO),
            "phase0_mode": "ci_microbench",
            "platform": args.label,
            "benchmarks": benchmarks,
        }
    )
    payload.setdefault("dock_timings_harvested", {"count": 0, "records": []})
    payload.setdefault("notes", [])
    payload["notes"] = list(payload["notes"]) + [
        "Microbench captured via scripts/capture_microbench_baseline.py",
    ]

    args.out.mkdir(parents=True, exist_ok=True)
    out_path = args.out / f"baseline_{args.label}.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {out_path} ({len(benchmarks)} benchmark entries)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())