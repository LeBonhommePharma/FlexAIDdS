#!/usr/bin/env python3
# queue_bisect_vcontacts.py — Priority #1: bisect Vcontacts/SoA before any new full-85
#
# Pauses v132 ablation ladder. After queue clears, runs smoke-12 bisect variants:
#   safe → head_soa_off → head_soa_on
# Then writes vcontacts_bisect_summary.json and optionally queues v131_safe_full85
# if the best variant beats head_soa_on and passes regression guards.
#
# Usage:
#   python3 scripts/queue_bisect_vcontacts.py --daemon
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
from lib_launch import launch_session_isolated
from launch_vcontacts_bisect_smoke import VARIANTS, harvest_report

REPO = SCRIPT_DIR.parent
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
LAUNCH_SMOKE = SCRIPT_DIR / "launch_vcontacts_bisect_smoke.py"
LAUNCH_SAFE_FULL85 = SCRIPT_DIR / "launch_v131_safe_full85.py"
QUEUE_LOG = RESULTS / "queue_bisect_vcontacts.log"
STATE_FILE = RESULTS / "queue_bisect_vcontacts.state.json"
V132_STATE = RESULTS / "queue_v132_ablation_ladder.state.json"

BISECT_ORDER = ("safe", "head_soa_off", "head_soa_on")
SMOKE_GATE_MIN = 8  # conservative vs 10/12 until Vcontacts fixed


def log(msg: str) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    with open(QUEUE_LOG, "a") as f:
        f.write(line + "\n")


def save_state(**kwargs) -> None:
    state = {}
    if STATE_FILE.exists():
        state = json.loads(STATE_FILE.read_text())
    state.update(kwargs)
    state["updated_at"] = datetime.now(timezone.utc).isoformat()
    STATE_FILE.write_text(json.dumps(state, indent=2) + "\n")


def pid_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def pause_v132_watcher() -> None:
    if not V132_STATE.exists():
        return
    st = json.loads(V132_STATE.read_text())
    wpid = int(st.get("watcher_pid") or 0)
    if pid_alive(wpid):
        log(f"Pausing v132 ablation watcher pid={wpid} (Vcontacts bisect takes priority)")
        try:
            os.kill(wpid, 15)
        except OSError:
            pass
    st["status"] = "paused_for_vcontacts_bisect"
    st["paused_at"] = datetime.now(timezone.utc).isoformat()
    V132_STATE.write_text(json.dumps(st, indent=2) + "\n")


def active_foreign_benchmarks() -> list[str]:
    try:
        ps = subprocess.check_output(["pgrep", "-fl", "benchmark_datasets"], text=True)
    except subprocess.CalledProcessError:
        return []
    lines = []
    for line in ps.strip().splitlines():
        if "caffeinate" in line:
            continue
        if "vcontacts_bisect" in line or "queue_bisect" in line:
            continue
        lines.append(line.strip())
    return lines


def launch_variant(variant: str) -> tuple[int, Path | None]:
    proc = subprocess.run(
        [sys.executable, str(LAUNCH_SMOKE), variant],
        cwd=str(REPO),
        capture_output=True,
        text=True,
    )
    if proc.stdout.strip():
        log(proc.stdout.strip())
    if proc.stderr.strip():
        log(f"{variant} stderr: {proc.stderr.strip()}")
    if proc.returncode != 0:
        return -1, None

    dirs = sorted(
        RESULTS.glob(f"vcontacts_bisect_*_{variant}"),
        key=lambda p: p.stat().st_mtime,
    )
    out = dirs[-1] if dirs else None
    pid = 0
    if out and (out / "launch_provenance.json").is_file():
        pid = int(json.loads((out / "launch_provenance.json").read_text()).get("pid") or 0)
    return pid, out


def wait_pid(pid: int, poll_s: int) -> None:
    while pid_alive(pid):
        time.sleep(poll_s)


def run_bisect(poll_s: int, settle_s: int, queue_full85: bool) -> int:
    pause_v132_watcher()
    save_state(status="watching", bisect_order=list(BISECT_ORDER))

    while active_foreign_benchmarks():
        n = len(active_foreign_benchmarks())
        log(f"waiting: foreign_benchmarks={n}")
        time.sleep(poll_s)

    log(f"Queue quiet — settling {settle_s}s before Vcontacts bisect")
    time.sleep(settle_s)

    results: list[dict] = []
    for variant in BISECT_ORDER:
        save_state(status="launching", current_variant=variant)
        pid, out_dir = launch_variant(variant)
        if pid <= 0 or not out_dir:
            save_state(status="launch_failed", variant=variant)
            return 1
        wait_pid(pid, poll_s)
        report = harvest_report(out_dir)
        (out_dir / "bisect_smoke_report.json").write_text(
            json.dumps(report, indent=2) + "\n"
        )
        entry = {
            "variant": variant,
            "label": VARIANTS[variant]["label"],
            "suspect": VARIANTS[variant]["suspect"],
            "output_dir": str(out_dir),
            **report,
        }
        results.append(entry)
        log(
            f"done {variant}: {report['n_success']}/12 success, "
            f"guard_fail={report['regression_guard_fail']}"
        )
        save_state(status="variant_done", last=entry)

    best = max(results, key=lambda r: (r["n_success"], -len(r["regression_guard_fail"])))
    head_on = next(r for r in results if r["variant"] == "head_soa_on")
    recommendation = "revert_to_safe_binary"
    if best["variant"] == "head_soa_on":
        recommendation = "soa_on_ok_unexpected_investigate"
    elif best["variant"] == "head_soa_off":
        recommendation = "ship_head_with_Soa_OFF_and_revert_inv_d12_if_needed"
    elif best["variant"] == "safe":
        recommendation = "use_v131_safe_full85_hold_HEAD_Vcontacts_changes"

    summary = {
        "recorded_at": datetime.now(timezone.utc).isoformat(),
        "priority": "bisect Vcontacts.cpp before new knob turns or v132 ablation",
        "context": {
            "v109_record": "80/85 (94.1%)",
            "v127_baseline": "78/85",
            "v130_observed": "73/85",
            "v131_observed": "72/84 incomplete",
            "vcontacts_commits_since_82ad51f4": [
                "27e68e51 inv_d12 + gaboom parallel",
                "f9c80fe5 SoA double sqrdist parity",
                "d4d68592 PR4 scalar-identical loop",
            ],
        },
        "variants": results,
        "best_variant": best,
        "recommendation": recommendation,
        "next_action": (
            "launch_v131_safe_full85"
            if best["n_success"] >= SMOKE_GATE_MIN and not best["regression_guard_fail"]
            else "fix_Vcontacts_then_rerun_bisect"
        ),
    }
    summary_path = RESULTS / "vcontacts_bisect_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    log(f"Wrote {summary_path}")
    log(f"Best: {best['variant']} ({best['n_success']}/12) → {recommendation}")

    if (
        queue_full85
        and best["n_success"] >= SMOKE_GATE_MIN
        and not best["regression_guard_fail"]
        and best["variant"] in ("safe", "head_soa_off")
    ):
        log("Queueing v131_safe_full85 (--ignore-smoke-gate) after bisect pass")
        proc = subprocess.run(
            [
                sys.executable,
                str(LAUNCH_SAFE_FULL85),
                "--skip-build",
                "--ignore-smoke-gate",
            ],
            cwd=str(REPO),
            capture_output=True,
            text=True,
        )
        log(proc.stdout.strip() or "(no stdout)")
        if proc.stderr.strip():
            log(proc.stderr.strip())

    save_state(status="bisect_done", summary_path=str(summary_path), best=best)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poll", type=int, default=120)
    parser.add_argument("--settle", type=int, default=120)
    parser.add_argument("--no-full85", action="store_true")
    parser.add_argument("--daemon", action="store_true")
    args = parser.parse_args()

    if args.daemon:
        watcher_log = RESULTS / "queue_bisect_vcontacts_watcher"
        watcher_log.mkdir(parents=True, exist_ok=True)
        cmd = [
            "caffeinate",
            "-i",
            sys.executable,
            str(Path(__file__).resolve()),
            "--poll",
            str(args.poll),
            "--settle",
            str(args.settle),
        ]
        if args.no_full85:
            cmd.append("--no-full85")
        pid = launch_session_isolated(
            cmd,
            os.environ.copy(),
            str(watcher_log),
            cwd=str(REPO),
            stdout_log=str(watcher_log / "watcher_stdout.log"),
            stderr_log=str(watcher_log / "watcher_stderr.log"),
        )
        save_state(status="watcher_daemon_started", watcher_pid=pid)
        print(f"Vcontacts bisect watcher pid={pid}")
        print(f"Log: {QUEUE_LOG}")
        return 0

    return run_bisect(args.poll, args.settle, not args.no_full85)


if __name__ == "__main__":
    raise SystemExit(main())