#!/usr/bin/env python3
"""Continue guard bisect ladder after arm1 (fixb_crg) completes."""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")

sys.path.insert(0, str(SCRIPTS))
from v132_common import wait_for_benchmark_done

ARM1_DIR = Path(
    "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v132_20260629_1657_guard_bisect_fixb_crg"
)
REMAINING = ("no_fixb", "no_crg", "no_fixb_no_crg")


def main() -> int:
    tag = "20260629_1657"
    parent = RESULTS / f"v132_{tag}_guard_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)
    manifest_path = parent / "ladder_manifest.json"

    print(f"waiting for arm1 fixb_crg: {ARM1_DIR}")
    if not wait_for_benchmark_done(ARM1_DIR):
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
        if not wait_for_benchmark_done(out_dir):
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