#!/usr/bin/env python3
"""Harvest read-only performance baselines without running new benchmarks."""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "results/perf_swarm"
DEFAULT_RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
TIMING_RE = re.compile(
    r"TIMING SUMMARY:\s+\d+\s+gens timed,\s+avg\s+([\d.]+)\s+ms/gen"
)


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


def harvest_timings(
    results_root: Path,
    *,
    log_name: str = "stderr.log",
    campaign_filter: str | None = None,
) -> list[dict]:
    harvested: list[dict] = []
    if not results_root.is_dir():
        return harvested

    for log in results_root.rglob(log_name):
        if campaign_filter and campaign_filter not in str(log):
            continue
        try:
            text = log.read_text(errors="replace")
        except OSError:
            continue
        match = TIMING_RE.search(text)
        if not match:
            continue

        rel = log.relative_to(results_root)
        parts = rel.parts
        campaign = parts[0] if len(parts) > 2 else results_root.name
        target = "/".join(parts[1:-1]) if len(parts) > 2 else parts[0]

        job_key = "/".join(parts[:-1])
        entry_key = job_key
        for marker in ("/results_resume_missing_198/", "/results/"):
            idx = job_key.find(marker)
            if idx >= 0:
                entry_key = job_key[idx + len(marker) :]
                break

        harvested.append(
            {
                "campaign": campaign,
                "target": target,
                "job_key": job_key,
                "entry_key": entry_key,
                "log_path": str(log),
                "avg_ms_per_gen": float(match.group(1)),
            }
        )

    harvested.sort(key=lambda row: row["avg_ms_per_gen"], reverse=True)
    return harvested


def active_queue() -> list[str]:
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return []
    lines = []
    for line in ps.strip().splitlines():
        if "caffeinate" in line:
            continue
        lines.append(line.strip())
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--campaign-filter", type=str, default=None)
    parser.add_argument("--label", type=str, default="macos_metal")
    args = parser.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    harvested = harvest_timings(
        args.results_root, campaign_filter=args.campaign_filter
    )
    queue = active_queue()
    git = git_info(REPO)

    payload = {
        "schema_version": "1.0.0",
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "git": git,
        "phase0_mode": "harvested_readonly",
        "queue_interference": "avoided",
        "active_benchmark_processes": queue,
        "platform": args.label,
        "results_root": str(args.results_root),
        "dock_timings_harvested": {
            "count": len(harvested),
            "slowest_5": harvested[:5],
            "fastest_5": harvested[-5:] if len(harvested) >= 5 else harvested,
            "records": harvested,
        },
        "notes": [
            "Dock timings harvested from stderr.log TIMING SUMMARY lines.",
            "Read-only harvest; does not launch new dock jobs.",
        ],
    }

    out_path = args.out / f"baseline_{args.label}.json"
    out_path.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {out_path}")
    print(f"Harvested {len(harvested)} timing records")


if __name__ == "__main__":
    main()