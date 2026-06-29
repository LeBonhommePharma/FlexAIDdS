#!/usr/bin/env python3
"""Poll v132 isolation-4 run; write report when 4/4 result.csv land."""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import os
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
MONITOR = RESULTS / "campaign_monitor.log"
TARGETS = ("1HQ2", "1OF1", "1T40", "1HNN")
GUARD_PAIR = ("1HQ2", "1T40")
POLL_SEC = 60


def log(msg: str) -> None:
    ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}\n"
    with open(MONITOR, "a") as fh:
        fh.write(line)
    print(line, end="")


def count_done(run_dir: Path) -> int:
    return len(list(run_dir.glob("*/result.csv")))


def pid_alive(run_dir: Path) -> bool:
    pid_file = run_dir / "benchmark.pid"
    if not pid_file.exists():
        return False
    try:
        pid = int(pid_file.read_text().strip())
        os.kill(pid, 0)
        return True
    except (ValueError, OSError):
        return False


def load_rows(run_dir: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for rf in run_dir.glob("*/result.csv"):
        row = dict(next(csv.DictReader(open(rf))))
        rows[row["pdb_id"]] = row
    return rows


def write_report(run_dir: Path, report_path: Path) -> str:
    rows = load_rows(run_dir)
    pass_n = sum(1 for t in TARGETS if rows.get(t, {}).get("success") == "1")
    guard_ok = all(rows.get(t, {}).get("success") == "1" for t in GUARD_PAIR)
    verdict = (
        "VARIANCE — rerun full smoke-12 solo (1HQ2+1T40 recovered)"
        if guard_ok
        else "SELECTOR REGRESSION — bisect Fix-B / elitism / BCR-gate on 1HQ2+1T40"
    )

    lines = [
        f"v132 isolation-4 report — {run_dir.name}",
        f"generated: {datetime.datetime.now(datetime.UTC).isoformat().replace('+00:00', 'Z')}",
        "",
        f"PASS: {pass_n}/4",
        f"Guard pair (1HQ2+1T40): {'PASS' if guard_ok else 'FAIL'}",
        f"VERDICT: {verdict}",
        "",
        "Per-target:",
    ]
    for tid in TARGETS:
        r = rows.get(tid, {})
        ok = r.get("success") == "1"
        lines.append(
            f"  {tid}: {'PASS' if ok else 'FAIL':4} "
            f"rmsd_h={r.get('rmsd_hungarian','?')} "
            f"rmsd_c={r.get('rmsd_to_crystal','?')} "
            f"bcr={r.get('best_cluster_rmsd','?')} "
            f"src={r.get('pose_source','')}"
        )
    lines.append("")
    lines.append("1HNN: consistent science blocker — site/sulfo variants regardless of guard verdict.")

    text = "\n".join(lines) + "\n"
    report_path.write_text(text)
    return text


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", help="Isolation-4 output directory")
    args = parser.parse_args()
    run_dir = Path(args.run_dir).resolve()
    check_script = SCRIPTS / "launch_v132_isolation4.py"
    report_path = run_dir / "v132_isolation4_report.txt"

    log(f"v132 ISOLATION4 WATCHER started — {run_dir.name}")

    while True:
        done = count_done(run_dir)
        alive = pid_alive(run_dir)
        if done < 4 and alive:
            log(f"v132 isolation4 {done}/4 result.csv — pid alive")
            time.sleep(POLL_SEC)
            continue
        if done < 4:
            log(f"v132 isolation4 STOPPED incomplete {done}/4")
            return 1

        if check_script.is_file():
            r = subprocess.run(
                [sys.executable, str(check_script), "--check", str(run_dir)],
                capture_output=True,
                text=True,
                cwd=str(REPO),
            )
            check_out = (r.stdout or r.stderr).strip()
        else:
            check_out = "check script missing"

        report = write_report(run_dir, report_path)
        log(f"v132 isolation4 DONE: {check_out.replace(chr(10), ' | ')}")
        log(f"v132 isolation4 REPORT: {report_path}")
        for line in report.splitlines():
            if line.startswith(("PASS:", "Guard", "VERDICT:", "  1")):
                log(f"  {line}")
        return 0


if __name__ == "__main__":
    sys.exit(main())