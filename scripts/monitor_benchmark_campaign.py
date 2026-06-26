#!/usr/bin/env python3
"""
monitor_benchmark_campaign.py — read-only polling monitor for benchmark runs.

Watches result directories, logs progress, and runs post-hoc audits when a
campaign completes. Never modifies or signals running benchmark processes.

Usage:
    python3 scripts/monitor_benchmark_campaign.py --config <monitor.json>
    python3 scripts/monitor_benchmark_campaign.py --v111 <dir> --baseline <dir>

Copyright 2026 Le Bonhomme Pharma. Apache-2.0.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(SCRIPT_DIR)
POLL_INTERVAL_S = 120
TOTAL_TARGETS = 85


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def count_results(root: str) -> tuple[int, int, int]:
    """Return (done, sub2, fail) from per-target result.csv files."""
    done = sub2 = fail = 0
    for rp in glob.glob(os.path.join(root, "*/result.csv")):
        try:
            with open(rp, newline="") as f:
                row = next(csv.DictReader(f))
            r = float(row.get("rmsd_hungarian") or row.get("rmsd_to_crystal") or -1)
            if not (r >= 0 and r < 998):
                continue
            done += 1
            if r < 2.0:
                sub2 += 1
            elif r >= 2.0:
                fail += 1
        except (OSError, StopIteration, ValueError, TypeError):
            pass
    return done, sub2, fail


def active_targets(root: str) -> list[dict]:
    """Targets with dirs but no result.csv yet."""
    out = []
    if not os.path.isdir(root):
        return out
    for d in glob.glob(os.path.join(root, "*")):
        if not os.path.isdir(d):
            continue
        pid = os.path.basename(d)
        if pid.endswith(".json") or pid.endswith(".log"):
            continue
        if os.path.isfile(os.path.join(d, "result.csv")):
            continue
        slog = os.path.join(d, "stdout.log")
        cf = None
        if os.path.isfile(slog):
            try:
                import re
                cfs = re.findall(r"cf=\s*([-\d.]+)", open(slog).read())
                if cfs:
                    cf = float(cfs[-1])
            except OSError:
                pass
        out.append({"pdb_id": pid, "last_cf": cf})
    return sorted(out, key=lambda x: x["pdb_id"])


def runner_alive(pid: int | None) -> bool | None:
    if not pid:
        return None
    try:
        os.kill(pid, 0)
        return True
    except OSError:
        return False


def read_pid_file(root: str) -> int | None:
    path = os.path.join(root, "benchmark.pid")
    if not os.path.isfile(path):
        return None
    try:
        return int(open(path).read().strip())
    except (OSError, ValueError):
        return None


def run_post_audit(root: str, label: str, log_path: str) -> bool:
    """Run failure_classify + cf_ground_truth_audit once per campaign."""
    marker = os.path.join(root, ".monitor_audit_done")
    if os.path.isfile(marker):
        return False
    if count_results(root)[0] < TOTAL_TARGETS:
        return False

    science = REPO
    for script in ("failure_classify.py", "cf_ground_truth_audit.py"):
        cmd = [sys.executable, os.path.join(science, "scripts", script), root]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=300)
        except subprocess.CalledProcessError as e:
            _log(log_path, f"[{label}] WARN audit {script} failed: {e.stderr[:500]}")
            return False
        except subprocess.TimeoutExpired:
            _log(log_path, f"[{label}] WARN audit {script} timed out")
            return False

    with open(marker, "w") as f:
        f.write(utc_now() + "\n")
    _log(log_path, f"[{label}] Post-run audits complete.")
    return True


def run_ab_compare(v111: str, baseline: str, log_path: str) -> None:
    marker = os.path.join(v111, ".monitor_ab_compare_done")
    if os.path.isfile(marker):
        return
    if count_results(v111)[0] < TOTAL_TARGETS:
        return
    proxy = os.path.expanduser("~/flexaidds_results/v68_20260615_3fixes_oracle")
    if not os.path.isdir(proxy):
        proxy = baseline
    script = os.path.join(REPO, "scripts", "ab_compare_report.py")
    if not os.path.isfile(script):
        return
    out = os.path.join(v111, "ab_compare_v111_vs_proxy.txt")
    try:
        r = subprocess.run(
            [sys.executable, script, "v68_proxy", proxy, "v111", v111],
            capture_output=True, text=True, timeout=120,
        )
        with open(out, "w") as f:
            f.write(r.stdout or "")
            if r.stderr:
                f.write("\n--- stderr ---\n" + r.stderr)
        _log(log_path, f"[v111] A/B report written: {out}")
        with open(marker, "w") as f:
            f.write(utc_now() + "\n")
    except (subprocess.TimeoutExpired, OSError) as e:
        _log(log_path, f"[v111] WARN ab_compare failed: {e}")


def _log(log_path: str, msg: str) -> None:
    line = f"{utc_now()} {msg}"
    print(line, flush=True)
    try:
        with open(log_path, "a") as f:
            f.write(line + "\n")
    except OSError:
        pass


def snapshot_campaign(root: str, label: str, pid: int | None) -> dict:
    done, sub2, fail = count_results(root)
    active = active_targets(root)
    alive = runner_alive(pid)
    status = "running"
    if done >= TOTAL_TARGETS:
        status = "complete"
    elif alive is False and done < TOTAL_TARGETS:
        status = "stalled"
    return {
        "label": label,
        "root": root,
        "status": status,
        "done": done,
        "total": TOTAL_TARGETS,
        "sub2": sub2,
        "fail": fail,
        "pct": round(100.0 * done / TOTAL_TARGETS, 1),
        "runner_pid": pid,
        "runner_alive": alive,
        "active_count": len(active),
        "active_sample": active[:6],
        "updated_at": utc_now(),
    }


def format_status(s: dict) -> str:
    lines = [
        f"{s['label']}: {s['done']}/{s['total']} ({s['pct']}%) "
        f"sub2={s['sub2']} fail={s['fail']} status={s['status']}"
    ]
    if s.get("runner_alive") is not None:
        lines[0] += f" pid={s.get('runner_pid')} alive={s['runner_alive']}"
    for a in s.get("active_sample") or []:
        cf = a.get("last_cf")
        cf_s = f"{cf:+.1f}" if cf is not None else "?"
        flag = " [clash]" if cf is not None and cf > 0 else ""
        lines.append(f"  → {a['pdb_id']} cf={cf_s}{flag}")
    return "\n".join(lines)


def monitor_loop(campaigns: list[dict], log_path: str, interval: int) -> None:
    state_path = os.path.join(os.path.dirname(log_path), "monitor_state.json")
    _log(log_path, f"Monitor started — polling every {interval}s")
    _log(log_path, f"Campaigns: {[c['label'] for c in campaigns]}")

    last_done = {c["label"]: -1 for c in campaigns}
    all_complete = False

    while not all_complete:
        snapshots = []
        for c in campaigns:
            root = c["root"]
            label = c["label"]
            pid = c.get("pid") or read_pid_file(root)
            snap = snapshot_campaign(root, label, pid)
            snapshots.append(snap)

            if snap["done"] != last_done[label]:
                _log(log_path, format_status(snap))
                last_done[label] = snap["done"]

            if snap["status"] == "complete":
                run_post_audit(root, label, log_path)

        v111 = next((s for s in snapshots if s["label"] == "v111_science"), None)
        baseline = next((s for s in snapshots if s["label"] == "baseline_8196829"), None)
        if v111 and v111["status"] == "complete":
            base_dir = baseline["root"] if baseline else ""
            run_ab_compare(v111["root"], base_dir, log_path)

        try:
            with open(state_path, "w") as f:
                json.dump({"campaigns": snapshots, "updated_at": utc_now()}, f, indent=2)
                f.write("\n")
        except OSError:
            pass

        all_complete = all(s["status"] == "complete" for s in snapshots)
        if all_complete:
            _log(log_path, "All campaigns complete. Monitor exiting.")
            summary_path = os.path.join(os.path.dirname(log_path), "monitor_final_summary.md")
            _write_final_summary(snapshots, summary_path, log_path)
            break

        time.sleep(interval)


def _write_final_summary(snapshots: list[dict], path: str, log_path: str) -> None:
    lines = ["# Benchmark Campaign Monitor — Final Summary", "", f"Generated: {utc_now()}", ""]
    for s in snapshots:
        lines.append(f"## {s['label']}")
        lines.append(f"- **Dir:** `{s['root']}`")
        lines.append(f"- **Result:** {s['sub2']}/{s['total']} sub-2 Å ({100*s['sub2']/s['total']:.1f}%)")
        lines.append(f"- **Failures:** {s['fail']}")
        audit = os.path.join(s["root"], "cf_audit_report.json")
        if os.path.isfile(audit):
            lines.append(f"- **CF audit:** `{audit}`")
        lines.append("")
    try:
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        _log(log_path, f"Final summary: {path}")
    except OSError:
        pass


def main() -> None:
    parser = argparse.ArgumentParser(description="Monitor FlexAIDdS benchmark campaigns")
    parser.add_argument("--v111", default=os.path.expanduser(
        "~/flexaidds_results/v111_science_20260626_0613"))
    parser.add_argument("--baseline", default=os.path.expanduser(
        "~/flexaidds_results/baseline_8196829_audit"))
    parser.add_argument("--interval", type=int, default=POLL_INTERVAL_S)
    parser.add_argument("--log", default=os.path.expanduser(
        "~/flexaidds_results/benchmark_monitor.log"))
    args = parser.parse_args()

    campaigns = [
        {"label": "v111_science", "root": os.path.expanduser(args.v111)},
        {"label": "baseline_8196829", "root": os.path.expanduser(args.baseline)},
    ]
    monitor_loop(campaigns, os.path.expanduser(args.log), args.interval)


if __name__ == "__main__":
    main()