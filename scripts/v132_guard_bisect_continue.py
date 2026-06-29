#!/usr/bin/env python3
"""Continue guard bisect ladder after arm1 (fixb_crg) completes."""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")

ARM1_DIR = Path(
    "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v132_20260629_1657_guard_bisect_fixb_crg"
)
REMAINING = ("no_fixb", "no_crg", "no_fixb_no_crg")


def wait_done(run_dir: Path, n: int = 2, poll: int = 60) -> bool:
    import os

    pid_file = run_dir / "benchmark.pid"
    while True:
        done = len(list(run_dir.glob("*/result.csv")))
        if done >= n:
            return True
        alive = False
        if pid_file.exists():
            try:
                os.kill(int(pid_file.read_text().strip()), 0)
                alive = True
            except (ValueError, OSError):
                pass
        if not alive and done < n:
            return False
        time.sleep(poll)


def main() -> int:
    tag = "20260629_1657"
    parent = RESULTS / f"v132_{tag}_guard_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)
    manifest_path = parent / "ladder_manifest.json"

    print(f"waiting for arm1 fixb_crg: {ARM1_DIR}")
    if not wait_done(ARM1_DIR):
        print("arm1 incomplete — abort")
        return 1

    manifest = {
        "parent": str(parent),
        "arms": [{"arm": "fixb_crg", "output_dir": str(ARM1_DIR)}],
        "targets": ["1HQ2", "1T40"],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    launch = SCRIPTS / "launch_v132_guard_bisect.py"
    for arm in REMAINING:
        print(f"\n=== launching arm {arm} ===")
        r = subprocess.run(
            [sys.executable, str(launch), "--arm", arm],
            cwd=str(REPO),
        )
        if r.returncode != 0:
            print(f"launch failed for {arm}")
            return r.returncode
        # find newest guard_bisect dir for this arm
        candidates = sorted(RESULTS.glob(f"v132_*_guard_bisect_{arm}"))
        if not candidates:
            print(f"no output dir for {arm}")
            return 1
        out_dir = candidates[-1]
        if not wait_done(out_dir):
            print(f"arm {arm} incomplete")
            break
        manifest["arms"].append({"arm": arm, "output_dir": str(out_dir)})
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    r = subprocess.run(
        [sys.executable, str(launch), "--report", str(parent)],
        cwd=str(REPO),
    )
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())