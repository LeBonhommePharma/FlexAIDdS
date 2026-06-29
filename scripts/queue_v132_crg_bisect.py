#!/usr/bin/env python3
"""Queue v132 CRG-parameter bisect ladder (solo) with watcher auto-report."""
from __future__ import annotations

import datetime
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
MONITOR = RESULTS / "campaign_monitor.log"

ARM_ORDER = ("crg_default", "crg_rmsd15", "crg_win8", "crg_softoff")


def log(msg: str) -> None:
    ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}\n"
    with open(MONITOR, "a") as fh:
        fh.write(line)
    print(line, end="")


def main() -> int:
    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    parent = RESULTS / f"v132_{tag}_crg_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)

    manifest = {
        "parent": str(parent),
        "arms": [],
        "targets": ["1HQ2", "1T40"],
        "expected_arms": list(ARM_ORDER),
        "reference": "v132_20260629_1657_guard_bisect_ladder",
        "status": "queued",
    }
    (parent / "ladder_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    log(
        "v132 CRG BISECT QUEUED: targets=1HQ2,1T40 "
        f"arms={','.join(ARM_ORDER)} | ladder={parent.name} | solo 2-thread"
    )

    ladder_log = parent / "ladder.log"
    launch = SCRIPTS / "launch_v132_crg_bisect.py"
    watcher = SCRIPTS / "v132_crg_bisect_watcher.py"

    with open(ladder_log, "a") as fh:
        ladder_proc = subprocess.Popen(
            ["caffeinate", "-i", sys.executable, str(launch), "--arm", "all", "--parent", str(parent)],
            cwd=str(REPO),
            stdout=fh,
            stderr=subprocess.STDOUT,
        )

    watcher_proc = subprocess.Popen(
        [sys.executable, str(watcher), str(parent)],
        cwd=str(REPO),
    )

    (parent / "queue_provenance.json").write_text(
        json.dumps(
            {
                "ladder_dir": str(parent),
                "ladder_pid": ladder_proc.pid,
                "watcher_pid": watcher_proc.pid,
                "queued_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
            },
            indent=2,
        )
        + "\n"
    )

    print(f"ladder_parent: {parent}")
    print(f"ladder_pid:    {ladder_proc.pid}")
    print(f"watcher_pid:   {watcher_proc.pid}")
    print(f"log:           {ladder_log}")
    return 0


if __name__ == "__main__":
    sys.exit(main())