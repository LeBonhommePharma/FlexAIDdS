#!/usr/bin/env python3
"""Compare microbenchmark text output against a perf baseline JSON artifact."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

THRESHOLD_PCT = 5.0


TENCOM_ROW_RE = re.compile(
    r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)"
)
VCFBATCH_SPEEDUP_RE = re.compile(r"Speedup\s*:\s*([\d.]+)", re.I)
DEFAULT_TENCOM_N_RES = 200


def parse_vcfbatch(text: str) -> dict[str, float]:
    m = VCFBATCH_SPEEDUP_RE.search(text)
    return {"speedup_vs_scalar": float(m.group(1))} if m else {}


def parse_tencom(text: str, *, reference_n_res: int = DEFAULT_TENCOM_N_RES) -> dict[str, float]:
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
    else:
        return {}

    return {"build_ms_full": build_ms, "sample_ms_full": sample_ms}


def delta_pct(current: float, baseline: float) -> float:
    if baseline <= 0:
        return 0.0
    return 100.0 * (current - baseline) / baseline


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--tencom", type=Path)
    parser.add_argument("--vcfbatch", type=Path)
    parser.add_argument("--threshold", type=float, default=THRESHOLD_PCT)
    args = parser.parse_args()

    baseline = json.loads(args.baseline.read_text())
    failures: list[str] = []

    checks: list[tuple[str, dict[str, float]]] = []
    if args.tencom and args.tencom.is_file():
        checks.append(("tencom", parse_tencom(args.tencom.read_text())))
    if args.vcfbatch and args.vcfbatch.is_file():
        checks.append(("vcfbatch", parse_vcfbatch(args.vcfbatch.read_text())))

    bench_map = {b["name"]: b for b in baseline.get("benchmarks", []) if isinstance(b, dict)}

    for name, metrics in checks:
        ref = bench_map.get(name, {})
        ref_metrics = ref.get("metrics", {})
        for key, val in metrics.items():
            base_val = ref_metrics.get(key)
            if base_val is None:
                print(f"[skip] {name}.{key}: no baseline")
                continue
            dp = delta_pct(val, base_val)
            status = "pass" if abs(dp) <= args.threshold else "FAIL"
            print(f"[{status}] {name}.{key}: current={val} baseline={base_val} delta={dp:+.1f}%")
            if status == "FAIL":
                failures.append(f"{name}.{key}")

    if failures:
        print(f"Regression: {', '.join(failures)}", file=sys.stderr)
        return 1
    print("All compared metrics within threshold.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())