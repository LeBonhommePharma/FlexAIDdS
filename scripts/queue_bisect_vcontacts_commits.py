#!/usr/bin/env python3
# queue_bisect_vcontacts_commits.py — serial smoke-12 for commit-revert variants
#
# Order: revert_f9c80fe5 → revert_d2295cf0 → revert_d4d68592
# Writes vcontacts_commit_bisect_summary.json
#
# Usage:
#   python3 scripts/queue_bisect_vcontacts_commits.py
#   python3 scripts/queue_bisect_vcontacts_commits.py --daemon
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from launch_vcontacts_bisect_smoke import harvest_report
from launch_vcontacts_commit_bisect_smoke import VARIANTS

REPO = SCRIPT_DIR.parent
GIT_ROOT = Path(os.environ.get("FLEXAIDDS_GIT_ROOT", str(REPO)))
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
LAUNCH = SCRIPT_DIR / "launch_vcontacts_commit_bisect_smoke.py"
LOG = RESULTS / "queue_bisect_vcontacts_commits.log"
STATE = RESULTS / "queue_bisect_vcontacts_commits.state.json"
SUMMARY = RESULTS / "vcontacts_commit_bisect_summary.json"
# d2295cf0 skipped: file-level + full git revert both conflict; coord-cache only (not Vcontacts).
ORDER = ("revert_f9c80fe5", "revert_d4d68592")


def log(msg: str) -> None:
    line = f"[{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}] {msg}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def save_state(**kwargs) -> None:
    st = json.loads(STATE.read_text()) if STATE.exists() else {}
    st.update(kwargs)
    st["updated_at"] = datetime.now(timezone.utc).isoformat()
    STATE.write_text(json.dumps(st, indent=2) + "\n")


def active_foreign_benchmarks() -> list[str]:
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return []
    lines = []
    for line in ps.strip().splitlines():
        if "caffeinate" in line:
            continue
        if "vcontacts_commit_bisect" in line or "queue_bisect_vcontacts_commits" in line:
            continue
        lines.append(line.strip())
    return lines


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def smoke_complete(out_dir: Path, expected: int = 12) -> bool:
    if not out_dir.is_dir():
        return False
    done = sum(
        1 for d in out_dir.iterdir()
        if d.is_dir() and (d / "result.csv").is_file()
    )
    return done >= expected


def launch_one(variant: str) -> tuple[int, Path | None]:
    env = {**os.environ, "FLEXAIDDS_GIT_ROOT": str(GIT_ROOT)}
    proc = subprocess.run(
        [sys.executable, str(LAUNCH), variant],
        cwd=str(REPO),
        capture_output=True,
        text=True,
        env=env,
    )
    if proc.stdout.strip():
        log(proc.stdout.strip())
    if proc.stderr.strip():
        log(f"{variant} stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        return -1, None
    dirs = sorted(
        RESULTS.glob(f"vcontacts_commit_bisect_*_{variant}"),
        key=lambda p: p.stat().st_mtime,
    )
    out = dirs[-1] if dirs else None
    pid = 0
    if out and (out / "launch_provenance.json").is_file():
        pid = int(json.loads((out / "launch_provenance.json").read_text()).get("pid") or 0)
    return pid, out


def run(poll_s: int, resume_from: str | None = None, settle_s: int = 60) -> int:
    save_state(status="watching", order=list(ORDER))
    while active_foreign_benchmarks():
        n = len(active_foreign_benchmarks())
        log(f"waiting: foreign_benchmarks={n}")
        time.sleep(poll_s)
    log(f"Queue quiet — settling {settle_s}s before commit bisect")
    time.sleep(settle_s)
    save_state(status="running", order=list(ORDER))
    results: list[dict] = []
    variants = ORDER
    if resume_from:
        if resume_from not in ORDER:
            log(f"ERROR: unknown --resume-from {resume_from}")
            return 1
        variants = ORDER[ORDER.index(resume_from):]
        log(f"Resuming commit bisect from {resume_from}")
    for variant in variants:
        save_state(status="launching", current_variant=variant)
        pid, out = launch_one(variant)
        if pid <= 0 or not out:
            save_state(status="launch_failed", variant=variant)
            return 1
        while not smoke_complete(out):
            if pid > 0 and not pid_alive(pid):
                # benchmark_datasets parent may exit while workers continue.
                done = sum(
                    1 for d in out.iterdir()
                    if d.is_dir() and (d / "result.csv").is_file()
                )
                log(f"{variant} parent exited; progress {done}/12")
            time.sleep(poll_s)
        report = harvest_report(out)
        (out / "bisect_smoke_report.json").write_text(json.dumps(report, indent=2) + "\n")
        entry = {
            "variant": variant,
            "reverted_sha": VARIANTS[variant]["sha"],
            **report,
            "output_dir": str(out),
        }
        results.append(entry)
        log(f"done {variant}: {report['n_success']}/12 guard_fail={report['regression_guard_fail']}")
        save_state(status="variant_done", last=entry)

    best = max(results, key=lambda r: (r["n_success"], -len(r["regression_guard_fail"])))
    summary = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "context": {
            "27e68e51": "already reverted on master (150402ac)",
            "v131_safe_full85": "80/85 @ v131_20260629_1401_safe_full85",
            "vcontacts_soA_bisect_best": "safe 7/12",
        },
        "variants": results,
        "best_variant": best,
        "recommendation": (
            f"cherry_pick_or_keep_revert_{best['reverted_sha']}"
            if best["n_success"] >= 8 and not best["regression_guard_fail"]
            else "keep_v131_safe_hold_further_reverts"
        ),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2) + "\n")
    log(f"Wrote {SUMMARY}")
    save_state(status="done", summary_path=str(SUMMARY), best=best)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poll", type=int, default=120)
    parser.add_argument("--settle", type=int, default=120)
    parser.add_argument("--resume-from", choices=ORDER)
    parser.add_argument("--daemon", action="store_true")
    args = parser.parse_args()
    if args.daemon:
        from lib_launch import launch_session_isolated

        watcher = RESULTS / "queue_bisect_vcontacts_commits_watcher"
        watcher.mkdir(parents=True, exist_ok=True)
        cmd = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--poll",
            str(args.poll),
            "--settle",
            str(args.settle),
        ]
        if args.resume_from:
            cmd.extend(["--resume-from", args.resume_from])
        pid = launch_session_isolated(
            ["caffeinate", "-i", *cmd],
            os.environ.copy(),
            str(watcher),
            cwd=str(REPO),
        )
        save_state(status="daemon_started", watcher_pid=pid)
        print(f"commit bisect watcher pid={pid}")
        return 0
    return run(args.poll, resume_from=args.resume_from, settle_s=args.settle)


if __name__ == "__main__":
    raise SystemExit(main())