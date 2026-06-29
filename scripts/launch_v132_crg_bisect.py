#!/usr/bin/env python3
# launch_v132_crg_bisect.py — CRG parameter bisect on guard pair (1HQ2, 1T40)
#
# Hold CONSENSUS=1, CRG=1; sweep CRG_RMSD_MAX / CRG_CF_WINDOW on fixb_crg baseline.
#
# Arms:
#   crg_default  — RMSD_MAX=2.5  CF_WINDOW=15  (v132 default)
#   crg_rmsd15   — RMSD_MAX=1.5  CF_WINDOW=15
#   crg_win8     — RMSD_MAX=2.5  CF_WINDOW=8
#   crg_softoff  — RMSD_MAX=0    CF_WINDOW=15  (CRG on, gate disabled)
#
# Usage:
#   python3 scripts/launch_v132_crg_bisect.py --arm crg_rmsd15
#   python3 scripts/launch_v132_crg_bisect.py --arm all
#   python3 scripts/launch_v132_crg_bisect.py --report <parent_dir>
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated
from v132_common import REPO, validate_manifest, v132_protocol_env, wait_for_benchmark_done

BUILD = f"{REPO}/build_lto"
BINARY = "/tmp/FlexAIDdS_v132_f7a0708f"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_guard_bisect_2_v132.json"
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
TARGETS = ("1HQ2", "1T40")
REFERENCE = "v132_20260629_1657_guard_bisect_ladder"
PINNED_COMMIT = "f7a0708f"

ARMS: dict[str, dict[str, str]] = {
    "crg_default": {
        "label": "CRG default (RMSD_MAX=2.5, CF_WINDOW=15)",
        "FLEXAIDDS_CRG_RMSD_MAX": "2.5",
        "FLEXAIDDS_CRG_CF_WINDOW": "15",
    },
    "crg_rmsd15": {
        "label": "CRG tighter gate (RMSD_MAX=1.5, CF_WINDOW=15)",
        "FLEXAIDDS_CRG_RMSD_MAX": "1.5",
        "FLEXAIDDS_CRG_CF_WINDOW": "15",
    },
    "crg_win8": {
        "label": "CRG narrow CF window (RMSD_MAX=2.5, CF_WINDOW=8)",
        "FLEXAIDDS_CRG_RMSD_MAX": "2.5",
        "FLEXAIDDS_CRG_CF_WINDOW": "8",
    },
    "crg_softoff": {
        "label": "CRG soft-off (RMSD_MAX=0, CF_WINDOW=15)",
        "FLEXAIDDS_CRG_RMSD_MAX": "0",
        "FLEXAIDDS_CRG_CF_WINDOW": "15",
    },
}
ARM_ORDER = ("crg_default", "crg_rmsd15", "crg_win8", "crg_softoff")


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_rows(run_dir: Path) -> dict[str, dict]:
    rows: dict[str, dict] = {}
    for rf in run_dir.glob("*/result.csv"):
        row = dict(next(csv.DictReader(open(rf))))
        rows[row["pdb_id"]] = row
    return rows


def arm_summary(run_dir: Path) -> dict:
    rows = load_rows(run_dir)
    out = {"dir": str(run_dir), "pass": 0, "targets": {}}
    for tid in TARGETS:
        r = rows.get(tid, {})
        ok = r.get("success") == "1"
        out["targets"][tid] = {
            "pass": ok,
            "rmsd_h": r.get("rmsd_hungarian", ""),
            "rmsd_c": r.get("rmsd_to_crystal", ""),
            "bcr": r.get("best_cluster_rmsd", ""),
            "source": r.get("pose_source", ""),
        }
        if ok:
            out["pass"] += 1
    out["guard_pass"] = all(out["targets"].get(t, {}).get("pass") for t in TARGETS)
    return out


def launch_arm(arm: str, parent: Path | None = None) -> tuple[Path, int]:
    if arm not in ARMS:
        raise ValueError(f"unknown arm: {arm}")

    validate_manifest(JSON_PAIRS)
    for path in (BINARY, RUNNER, JSON_PAIRS, f"{BUILD}/MC_st0r5.2_6.dat"):
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = RESULTS_DIR / f"v132_{tag}_crg_bisect_{arm}"
    cache = RESULTS_DIR / f"cache_v132_crg_bisect_{arm}"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

    env = dict(os.environ)
    env.update(v132_protocol_env(BINARY, BUILD, str(cache), ORACLE_DIR))
    env.update({k: v for k, v in ARMS[arm].items() if k.startswith("FLEXAIDDS_")})
    env["FLEXAIDDS_BINARY"] = BINARY
    env["FLEXAIDDS_CONSENSUS_SCORER"] = "1"
    env["FLEXAIDDS_CRG"] = "1"
    env["FLEXAIDDS_PRIORITY_TARGETS"] = ",".join(TARGETS)

    cmd = [
        "caffeinate", "-i", RUNNER,
        "--benchmark", f"crossdock_json:{JSON_PAIRS}",
        "--output", str(output),
        "--threads", bench_threads,
        "--omp-threads", "1",
        "--temperature", "298",
        "--job-timeout-seconds", "7200",
        "--cache", str(cache),
        "--mode", "oracle-ceiling",
        "--force",
    ]

    output.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    print(f"\nLaunch arm={arm} — {ARMS[arm]['label']}")
    print(f"  output : {output}")
    pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov = {
        "version": f"v132_crg_bisect_{arm}",
        "arm": arm,
        "label": ARMS[arm]["label"],
        "launched_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit": PINNED_COMMIT,
        "binary_path": BINARY,
        "binary_sha256": sha256(BINARY),
        "targets": list(TARGETS),
        "reference": REFERENCE,
        "output_dir": str(output),
        "pid": pid,
        "env_overrides": {
            "FLEXAIDDS_CONSENSUS_SCORER": "1",
            "FLEXAIDDS_CRG": "1",
            **{k: v for k, v in ARMS[arm].items() if k.startswith("FLEXAIDDS_")},
        },
    }
    if parent:
        prov["parent_ladder"] = str(parent)
    with open(output / "launch_provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2)
        fh.write("\n")

    print(f"  pid    : {pid}")
    return output, pid


def write_report(parent: Path) -> str:
    manifest = json.loads((parent / "ladder_manifest.json").read_text())
    lines = [
        f"v132 CRG parameter bisect — {parent.name}",
        f"generated: {datetime.datetime.now(datetime.UTC).isoformat().replace('+00:00', 'Z')}",
        f"targets: {','.join(TARGETS)}",
        f"reference: {REFERENCE}",
        f"baseline: CONSENSUS=1 CRG=1",
        "",
    ]
    best_arm = None
    best_pass = -1
    guard_arms: list[str] = []
    hq2_arms: list[str] = []
    t40_arms: list[str] = []

    for entry in manifest["arms"]:
        run_dir = Path(entry["output_dir"])
        if not run_dir.is_dir():
            lines.append(f"## {entry['arm']}: MISSING")
            continue
        s = arm_summary(run_dir)
        arm_cfg = ARMS[entry["arm"]]
        lines.append(f"## {entry['arm']} — {arm_cfg['label']}")
        lines.append(
            f"   knobs: RMSD_MAX={arm_cfg.get('FLEXAIDDS_CRG_RMSD_MAX', '?')} "
            f"CF_WINDOW={arm_cfg.get('FLEXAIDDS_CRG_CF_WINDOW', '?')}"
        )
        lines.append(f"   guard_pass={s['guard_pass']} ({s['pass']}/2)")
        for tid in TARGETS:
            t = s["targets"].get(tid, {})
            lines.append(
                f"   {tid}: {'PASS' if t.get('pass') else 'FAIL'} "
                f"rmsd_h={t.get('rmsd_h')} src={t.get('source')}"
            )
            if t.get("pass"):
                (hq2_arms if tid == "1HQ2" else t40_arms).append(entry["arm"])
        if s["guard_pass"]:
            guard_arms.append(entry["arm"])
        if s["pass"] > best_pass or (s["pass"] == best_pass and s.get("guard_pass")):
            best_pass = s["pass"]
            best_arm = entry["arm"]
        lines.append("")

    if best_arm:
        lines.append(f"BEST_ARM: {best_arm} ({best_pass}/2)")
        if guard_arms:
            lines.append(
                f"VERDICT: guard pair recovered — promote {guard_arms[0]} CRG settings to smoke-12"
            )
        elif hq2_arms and t40_arms:
            lines.append(
                f"VERDICT: CRG tradeoff persists (1HQ2 arms={','.join(hq2_arms)}; "
                f"1T40 arms={','.join(t40_arms)}) — try intermediate RMSD/window grid"
            )
        else:
            lines.append("VERDICT: no CRG knob recovered guards — escalate selector bisect")
    text = "\n".join(lines) + "\n"
    (parent / "v132_crg_bisect_report.txt").write_text(text)
    return text


def run_ladder(parent: Path | None = None) -> int:
    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    if parent is None:
        parent = RESULTS_DIR / f"v132_{tag}_crg_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "parent": str(parent),
        "arms": [],
        "targets": list(TARGETS),
        "expected_arms": list(ARM_ORDER),
        "reference": REFERENCE,
    }

    for arm in ARM_ORDER:
        print(f"\n=== ladder step: {arm} ===")
        out_dir, pid = launch_arm(arm, parent=parent)
        manifest["arms"].append({"arm": arm, "output_dir": str(out_dir), "pid": pid})
        with open(parent / "ladder_manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2)
            fh.write("\n")
        ok = wait_for_benchmark_done(out_dir, n=2, stall_grace=900)
        if not ok:
            manifest["status"] = "failed"
            with open(parent / "ladder_manifest.json", "w") as fh:
                json.dump(manifest, fh, indent=2)
                fh.write("\n")
            print(f"arm {arm} incomplete — stopping ladder")
            break
        s = arm_summary(out_dir)
        print(f"arm {arm} done: {s['pass']}/2 guard_pass={s['guard_pass']}")
    else:
        manifest["status"] = "complete"

    with open(parent / "ladder_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")

    report = write_report(parent)
    print(report)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=[*ARMS.keys(), "all"], default="all")
    parser.add_argument("--report", metavar="PARENT_DIR", help="Write bisect report from ladder manifest")
    parser.add_argument("--parent", metavar="PARENT_DIR", help="Ladder parent dir for --arm all")
    args = parser.parse_args()

    if args.report:
        print(write_report(Path(args.report)))
        return 0

    if args.arm == "all":
        parent = Path(args.parent).resolve() if args.parent else None
        return run_ladder(parent)

    out_dir, _ = launch_arm(args.arm)
    print(f"Launched single arm → {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())