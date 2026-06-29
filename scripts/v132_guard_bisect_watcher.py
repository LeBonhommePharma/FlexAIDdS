#!/usr/bin/env python3
"""Poll v132 guard-bisect ladder; write report when all arms complete."""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
MONITOR = RESULTS / "campaign_monitor.log"
POLL_SEC = 90


def log(msg: str) -> None:
    ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}\n"
    with open(MONITOR, "a") as fh:
        fh.write(line)
    print(line, end="")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("parent_dir", help="Guard bisect ladder parent directory")
    args = parser.parse_args()
    parent = Path(args.parent_dir).resolve()
    manifest_path = parent / "ladder_manifest.json"
    report_script = SCRIPTS / "launch_v132_guard_bisect.py"

    log(f"v132 GUARD BISECT WATCHER started — {parent.name}")

    while True:
        if not manifest_path.is_file():
            log("waiting for ladder_manifest.json")
            time.sleep(POLL_SEC)
            continue

        manifest = json.loads(manifest_path.read_text())
        arms = manifest.get("arms", [])
        if not arms:
            time.sleep(POLL_SEC)
            continue

        statuses = []
        all_done = True
        for entry in arms:
            run_dir = Path(entry["output_dir"])
            n = len(list(run_dir.glob("*/result.csv")))
            statuses.append(f"{entry['arm']}={n}/2")
            if n < 2:
                all_done = False

        expected = len(manifest.get("expected_arms", ("fixb_crg", "no_fixb", "no_crg", "no_fixb_no_crg")))

        if not all_done or len(arms) < expected:
            suffix = ""
            if all_done and len(arms) < expected:
                suffix = f" (ladder {len(arms)}/{expected} arms queued)"
            log(f"v132 guard_bisect progress: {', '.join(statuses)}{suffix}")
            time.sleep(POLL_SEC)
            continue

        status = manifest.get("status")
        if status == "failed":
            log(f"v132 guard_bisect ladder failed ({len(arms)}/{expected} arms)")
        else:
            log(f"v132 guard_bisect ladder complete ({len(arms)}/{expected} arms)")

        import subprocess

        r = subprocess.run(
            [sys.executable, str(report_script), "--report", str(parent)],
            capture_output=True,
            text=True,
            cwd=str(REPO),
        )
        report = (r.stdout or r.stderr).strip()
        log(f"v132 guard_bisect REPORT: {parent}/v132_guard_bisect_report.txt")
        for line in report.splitlines():
            if line.startswith(("##", "BEST_ARM", "VERDICT", "   1")):
                log(f"  {line}")
        return 0


if __name__ == "__main__":
    sys.exit(main())