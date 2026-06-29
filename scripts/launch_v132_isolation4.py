#!/usr/bin/env python3
# launch_v132_isolation4.py — 4-target isolation batch (smoke-1148 failures)
#
# Targets: 1HQ2, 1OF1, 1T40, 1HNN
# Protocol: identical to v132_20260629_1148_smoke12_crg (f7a0708f binary)
#
# Decision rule (from v132 plan):
#   - If 1HQ2 + 1T40 pass → variance; schedule another full smoke-12 solo.
#   - If they fail again → bisect selector on guard regressions.
#   - 1HNN consistent miss → science target regardless of guard noise.
#
# Usage:
#   python3 scripts/launch_v132_isolation4.py
#   python3 scripts/launch_v132_isolation4.py --check <run_dir>
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
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_launch import launch_session_isolated
from v132_common import REPO, validate_manifest, v132_protocol_env

BUILD = f"{REPO}/build_lto"
# Pin to the binary that produced valid v132_20260629_1148_smoke12_crg results.
BINARY = "/tmp/FlexAIDdS_v132_f7a0708f"
RUNNER = f"{BUILD}/benchmark_datasets"
ORACLE_DIR = f"{REPO}/benchmarks/astex_diverse/astex_diverse"
JSON_PAIRS = f"{REPO}/benchmarks/datasets/benchmark_astex_isolation_4_v132.json"
RESULTS_DIR = Path("/Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results")
TARGETS = ("1HQ2", "1OF1", "1T40", "1HNN")
REFERENCE_RUN = "v132_20260629_1148_smoke12_crg"
PINNED_COMMIT = "f7a0708f"


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def summarize(run_dir: Path) -> tuple[int, list[dict]]:
    rows: list[dict] = []
    for rf in sorted(run_dir.glob("*/result.csv")):
        rows.append(dict(next(csv.DictReader(open(rf)))))
    return len(rows), rows


def check_run(run_dir: Path) -> bool:
    n, rows = summarize(run_dir)
    if n < 4:
        print(f"incomplete {n}/4")
        return False
    by_id = {r["pdb_id"]: r for r in rows}
    print(f"isolation4 {run_dir.name}: {sum(1 for r in rows if r.get('success')=='1')}/4 pass")
    for tid in TARGETS:
        r = by_id.get(tid, {})
        ok = r.get("success") == "1"
        rmsd = r.get("rmsd_hungarian", "NA")
        src = r.get("pose_source", "")
        print(f"  {tid}: {'PASS' if ok else 'FAIL'} rmsd={rmsd} src={src}")
    guard = ("1HQ2", "1T40")
    guard_ok = all(by_id.get(g, {}).get("success") == "1" for g in guard)
    print(f"guard_check 1HQ2+1T40: {'PASS' if guard_ok else 'FAIL'}")
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", metavar="RUN_DIR", help="Summarize completed isolation run")
    args = parser.parse_args()

    if args.check:
        check_run(Path(args.check))
        return 0

    validate_manifest(JSON_PAIRS)
    for path in (BINARY, RUNNER, JSON_PAIRS, f"{BUILD}/MC_st0r5.2_6.dat"):
        if not os.path.exists(path):
            sys.exit(f"ERROR: missing {path}")

    tag = datetime.datetime.now().strftime("%Y%m%d_%H%M")
    output = str(RESULTS_DIR / f"v132_{tag}_isolation4_crg")
    cache = str(RESULTS_DIR / "cache_v132_isolation4")
    bench_threads = os.environ.get("FLEXAIDDS_BENCH_THREADS", "2")

    env = dict(os.environ)
    env.update(v132_protocol_env(BINARY, BUILD, cache, ORACLE_DIR))
    env["FLEXAIDDS_BINARY"] = BINARY
    env["FLEXAIDDS_PRIORITY_TARGETS"] = ",".join(TARGETS)

    cmd = [
        "caffeinate", "-i", RUNNER,
        "--benchmark", f"crossdock_json:{JSON_PAIRS}",
        "--output", output,
        "--threads", bench_threads,
        "--omp-threads", "1",
        "--temperature", "298",
        "--job-timeout-seconds", "7200",
        "--cache", cache,
        "--mode", "oracle-ceiling",
    ]

    os.makedirs(output, exist_ok=True)
    os.makedirs(cache, exist_ok=True)

    print("\nLaunching v132 isolation-4 — smoke1148 failure panel")
    print(f"  targets : {','.join(TARGETS)}")
    print(f"  binary  : {BINARY} ({sha256(BINARY)[:12]}…)")
    print(f"  ref run : {REFERENCE_RUN} @ {PINNED_COMMIT}")
    print(f"  output  : {output}")
    print(f"  threads : {bench_threads}")

    pid = launch_session_isolated(cmd, env, output, cwd=REPO)

    prov = {
        "version": "v132_isolation4_crg",
        "launched_at": datetime.datetime.now(datetime.UTC).isoformat().replace("+00:00", "Z"),
        "git_commit": PINNED_COMMIT,
        "binary_path": BINARY,
        "binary_sha256": sha256(BINARY),
        "json_pairs": JSON_PAIRS,
        "targets": list(TARGETS),
        "reference_run": REFERENCE_RUN,
        "output_dir": output,
        "pid": pid,
        "decision_rule": {
            "guard_pair": ["1HQ2", "1T40"],
            "pass_guard_again": "variance — rerun full smoke-12 solo",
            "fail_guard_again": "bisect selector (Fix-B / elitism / BCR-gate)",
            "1HNN": "consistent science blocker — deep dive regardless",
        },
    }
    with open(f"{output}/launch_provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2)
        fh.write("\n")

    print(f"v132 isolation-4 launched pid={pid}")
    return pid


if __name__ == "__main__":
    main()