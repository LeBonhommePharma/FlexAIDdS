#!/usr/bin/env python3
# launch_v132_guard_bisect.py — selector bisect on guard regressions (1HQ2, 1T40)
#
# Arms (one knob each vs v132 baseline):
#   fixb_crg        — CONSENSUS=1 CRG=1  (reproduce isolation4 failure)
#   no_fixb         — CONSENSUS=0 CRG=1  (test Fix-B suppressing ini_elitism)
#   no_crg          — CONSENSUS=1 CRG=0  (CRG interference control)
#   no_fixb_no_crg  — CONSENSUS=0 CRG=0  (maximal seed-election path)
#
# Usage:
#   python3 scripts/launch_v132_guard_bisect.py --arm no_fixb
#   python3 scripts/launch_v132_guard_bisect.py --arm all   # sequential ladder
#   python3 scripts/launch_v132_guard_bisect.py --report <parent_dir>
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

from __future__ import annotations

import argparse
import csv
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time
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
REFERENCE = "v132_20260629_1230_isolation4_crg"
PINNED_COMMIT = "f7a0708f"

ARMS: dict[str, dict[str, str]] = {
    "fixb_crg": {
        "label": "v132 baseline Fix-B+CRG",
        "FLEXAIDDS_CONSENSUS_SCORER": "1",
        "FLEXAIDDS_CRG": "1",
    },
    "no_fixb": {
        "label": "Fix-B OFF (consensus=0)",
        "FLEXAIDDS_CONSENSUS_SCORER": "0",
        "FLEXAIDDS_CRG": "1",
    },
    "no_crg": {
        "label": "CRG OFF",
        "FLEXAIDDS_CONSENSUS_SCORER": "1",
        "FLEXAIDDS_CRG": "0",
    },
    "no_fixb_no_crg": {
        "label": "Fix-B OFF + CRG OFF",
        "FLEXAIDDS_CONSENSUS_SCORER": "0",
        "FLEXAIDDS_CRG": "0",
    },
}
ARM_ORDER = ("fixb_crg", "no_fixb", "no_crg", "no_fixb_no_crg")


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
    output = RESULTS_DIR / f"v132_{tag}_guard_bisect_{arm}"
    cache = RESULTS_DIR / f"cache_v132_guard_bisect_{arm}"
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

    env = dict(os.environ)
    env.update(v132_protocol_env(BINARY, BUILD, str(cache), ORACLE_DIR))
    env.update({k: v for k, v in ARMS[arm].items() if k.startswith("FLEXAIDDS_")})
    env["FLEXAIDDS_BINARY"] = BINARY
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
        "version": f"v132_guard_bisect_{arm}",
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
        "env_overrides": {k: v for k, v in ARMS[arm].items() if k.startswith("FLEXAIDDS_")},
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
        f"v132 guard selector bisect — {parent.name}",
        f"generated: {datetime.datetime.now(datetime.UTC).isoformat().replace('+00:00', 'Z')}",
        f"targets: {','.join(TARGETS)}",
        f"reference: {REFERENCE}",
        "",
    ]
    best_arm = None
    best_pass = -1
    for entry in manifest["arms"]:
        run_dir = Path(entry["output_dir"])
        if not run_dir.is_dir():
            lines.append(f"## {entry['arm']}: MISSING")
            continue
        s = arm_summary(run_dir)
        lines.append(f"## {entry['arm']} — {ARMS[entry['arm']]['label']}")
        lines.append(f"   guard_pass={s['guard_pass']} ({s['pass']}/2)")
        for tid in TARGETS:
            t = s["targets"].get(tid, {})
            lines.append(
                f"   {tid}: {'PASS' if t.get('pass') else 'FAIL'} "
                f"rmsd_h={t.get('rmsd_h')} src={t.get('source')}"
            )
        if s["pass"] > best_pass:
            best_pass = s["pass"]
            best_arm = entry["arm"]
        lines.append("")

    if best_arm:
        lines.append(f"BEST_ARM: {best_arm} ({best_pass}/2)")
        if manifest["arms"][-1] and arm_summary(Path(manifest["arms"][-1]["output_dir"])).get("guard_pass"):
            lines.append("VERDICT: selector fix identified — promote winning arm to smoke-12")
        elif any(
            arm_summary(Path(e["output_dir"]))["guard_pass"]
            for e in manifest["arms"]
            if e["arm"] == "no_fixb"
        ):
            lines.append("VERDICT: Fix-B (consensus scorer) suppresses ini_elitism on guards")
        else:
            lines.append("VERDICT: no arm recovered guards — deeper selector bisect needed")
    text = "\n".join(lines) + "\n"
    (parent / "v132_guard_bisect_report.txt").write_text(text)
    return text


def run_ladder() -> int:
    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    parent = RESULTS_DIR / f"v132_{tag}_guard_bisect_ladder"
    parent.mkdir(parents=True, exist_ok=True)
    manifest = {"parent": str(parent), "arms": [], "targets": list(TARGETS)}

    for arm in ARM_ORDER:
        print(f"\n=== ladder step: {arm} ===")
        out_dir, pid = launch_arm(arm, parent=parent)
        manifest["arms"].append({"arm": arm, "output_dir": str(out_dir), "pid": pid})
        with open(parent / "ladder_manifest.json", "w") as fh:
            json.dump(manifest, fh, indent=2)
            fh.write("\n")
        ok = wait_for_benchmark_done(out_dir, n=2)
        if not ok:
            print(f"arm {arm} incomplete — stopping ladder")
            break
        s = arm_summary(out_dir)
        print(f"arm {arm} done: {s['pass']}/2 guard_pass={s['guard_pass']}")

    report = write_report(parent)
    print(report)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arm", choices=[*ARMS.keys(), "all"], default="all")
    parser.add_argument("--report", metavar="PARENT_DIR", help="Write bisect report from ladder manifest")
    args = parser.parse_args()

    if args.report:
        print(write_report(Path(args.report)))
        return 0

    if args.arm == "all":
        return run_ladder()

    out_dir, _ = launch_arm(args.arm)
    print(f"Launched single arm → {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())