#!/usr/bin/env python3
"""Continue guard bisect ladder after arm1 (fixb_crg) completes."""
from __future__ import annotations

import argparse
import datetime
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPTS = REPO / "scripts"
RESULTS = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
MONITOR = RESULTS / "campaign_monitor.log"

sys.path.insert(0, str(SCRIPTS))
from v132_common import run_dir_has_active_docking, wait_for_benchmark_done

DEFAULT_ARM1 = Path(
    "/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v132_20260629_1657_guard_bisect_fixb_crg"
)
REMAINING = ("no_fixb", "no_crg", "no_fixb_no_crg")


def log(msg: str) -> None:
    ts = datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%d %H:%M:%S UTC")
    line = f"[{ts}] {msg}\n"
    with open(MONITOR, "a") as fh:
        fh.write(line)
    print(line, end="")


def arm1_complete(run_dir: Path) -> bool:
    return len(list(run_dir.glob("*/result.csv"))) >= 2


def newest_fixb_crg_dir() -> Path | None:
    candidates = sorted(RESULTS.glob("v132_*_guard_bisect_fixb_crg"))
    return candidates[-1] if candidates else None


def relaunch_fixb_crg() -> Path:
    launch = SCRIPTS / "launch_v132_guard_bisect.py"
    log("v132 guard_bisect: relaunching arm1 fixb_crg (prior arm incomplete)")
    r = subprocess.run(
        [sys.executable, str(launch), "--arm", "fixb_crg"],
        cwd=str(REPO),
    )
    if r.returncode != 0:
        raise RuntimeError(f"fixb_crg relaunch failed (exit={r.returncode})")
    out = newest_fixb_crg_dir()
    if out is None:
        raise RuntimeError("fixb_crg relaunch produced no output dir")
    log(f"v132 guard_bisect arm1 relaunched → {out.name}")
    return out


def resolve_arm1(arm1_dir: Path | None) -> Path:
    if arm1_dir is not None:
        return arm1_dir.resolve()
    if DEFAULT_ARM1.is_dir() and (arm1_complete(DEFAULT_ARM1) or run_dir_has_active_docking(DEFAULT_ARM1)):
        return DEFAULT_ARM1.resolve()
    newer = newest_fixb_crg_dir()
    if newer is not None and (arm1_complete(newer) or run_dir_has_active_docking(newer)):
        return newer.resolve()
    return DEFAULT_ARM1.resolve()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm1-dir", type=Path, default=None, help="Arm1 output directory")
    parser.add_argument("--tag", default="20260629_1657", help="Ladder parent tag suffix")
    parser.add_argument("--relaunch-arm1", action="store_true", help="Force fresh fixb_crg launch")
    args = parser.parse_args()

    tag = args.tag
    parent = RESULTS / f"v132_{tag}_guard_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)
    manifest_path = parent / "ladder_manifest.json"

    arm1_dir = resolve_arm1(args.arm1_dir)
    if args.relaunch_arm1 or (
        not arm1_complete(arm1_dir) and not run_dir_has_active_docking(arm1_dir)
    ):
        arm1_dir = relaunch_fixb_crg()

    print(f"waiting for arm1 fixb_crg: {arm1_dir}")
    log(f"v132 guard_bisect continue: waiting arm1 {arm1_dir.name}")
    if not wait_for_benchmark_done(arm1_dir, stall_grace=900):
        log(f"v132 guard_bisect arm1 INCOMPLETE — {arm1_dir.name}")
        print("arm1 incomplete — abort")
        return 1

    manifest = {
        "parent": str(parent),
        "arms": [{"arm": "fixb_crg", "output_dir": str(arm1_dir)}],
        "targets": ["1HQ2", "1T40"],
        "expected_arms": ["fixb_crg", "no_fixb", "no_crg", "no_fixb_no_crg"],
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
        if not wait_for_benchmark_done(out_dir, stall_grace=900):
            log(f"v132 guard_bisect arm {arm} INCOMPLETE — {out_dir.name}")
            print(f"arm {arm} incomplete")
            break
        manifest["arms"].append({"arm": arm, "output_dir": str(out_dir)})
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    manifest["status"] = "complete" if len(manifest["arms"]) >= 4 else "failed"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    r = subprocess.run(
        [sys.executable, str(launch), "--report", str(parent)],
        cwd=str(REPO),
    )
    log(f"v132 guard_bisect ladder finished status={manifest['status']}")
    return r.returncode


if __name__ == "__main__":
    sys.exit(main())