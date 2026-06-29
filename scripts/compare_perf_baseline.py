#!/usr/bin/env python3
"""Compare microbenchmark text output against a perf baseline JSON artifact."""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

THRESHOLD_PCT = 5.0


def parse_vcfbatch(text: str) -> dict[str, float]:
    m = re.search(r"speedup:\s*([\d.]+)", text, re.I)
    return {"speedup_vs_scalar": float(m.group(1))} if m else {}


def parse_tencom(text: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for key, pat in [
        ("build_ms_full", r"build.*?full.*?([\d.]+)\s*ms"),
        ("sample_ms_full", r"sample.*?full.*?([\d.]+)\s*ms"),
    ]:
        m = re.search(pat, text, re.I)
        if m:
            out[key] = float(m.group(1))
    return out


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