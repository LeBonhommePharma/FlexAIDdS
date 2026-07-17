#!/usr/bin/env python3
"""launch_v135.py — Crystal-blind basin recovery election (BCR-proxy).

v135 enables DatasetRunner election recovery behind FLEXAIDDS_ELECTION_V135:

  * Include Frequency=1 cluster heads (no freq>1 gate)
  * Score temperature τ in CF arbitrary units (default 25; not mixed kT=0.592)
  * FO dual-suffix pose enumeration

Default claim ranking is unchanged when the flag is off (AGENTS.md).

Does NOT dual-launch into the live C0 claim OUT. Uses a dedicated v135
iCloud campaign namespace.

Usage:
  source ~/.flexaidds_env
  source scripts/use_icloud_benchmark_storage.sh
  python3 scripts/launch_v135.py --dry-run
  python3 scripts/launch_v135.py

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def repo_root() -> Path:
    env = os.environ.get("FLEXAIDDS_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
        return Path(out)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--manifest",
        default=None,
        help="crossdock JSON (default: queue astex_native_85.json)",
    )
    ap.add_argument("--pop", type=int, default=1000)
    ap.add_argument("--gen", type=int, default=2000)
    ap.add_argument("--tau", type=float, default=25.0, help="CF a.u. score temperature")
    args = ap.parse_args()

    root = repo_root()
    os.chdir(root)

    # Prefer iCloud queue / results when configured
    q = Path(
        os.environ.get(
            "FLEXAIDDS_QUEUE_ROOT",
            str(
                Path.home()
                / "Library/Mobile Documents/com~apple~CloudDocs"
                / "FlexAIDdS_benchmarks/queues/three_engine_entropy_q1"
            ),
        )
    )
    results = Path(
        os.environ.get(
            "FLEXAIDDS_RESULTS",
            str(
                Path.home()
                / "Library/Mobile Documents/com~apple~CloudDocs"
                / "FlexAIDdS_benchmarks/results"
            ),
        )
    )

    stamp = dt.datetime.utcnow().strftime("%Y%m%d")
    out = results / "campaigns" / f"v135_bcr_proxy_election_{stamp}"
    manifest = Path(
        args.manifest
        or (q / "inputs" / "astex_native_85.json")
    )
    runner = q / "bin" / "C" / "benchmark_datasets"
    binary = q / "bin" / "C" / "FlexAIDdS"
    # Prefer freshly built claim engine if present and newer worktree build
    for cand in (
        root / "build_claim" / "FlexAID",
        root / "build_claim" / "FlexAIDdS",
        root / "build_lto" / "FlexAID",
        root / "build" / "FlexAID",
    ):
        if cand.is_file():
            binary = cand
            break
    for cand in (
        root / "build_claim" / "benchmark_datasets",
        root / "build_lto" / "benchmark_datasets",
        root / "build" / "benchmark_datasets",
    ):
        if cand.is_file():
            runner = cand
            break

    data_dir = Path(os.environ.get("FLEXAIDDS_DATA_DIR", str(q / "data")))

    print("=== v135 BCR-proxy election launch ===")
    print(f"OUT={out}")
    print(f"MANIFEST={manifest}")
    print(f"RUNNER={runner}")
    print(f"BINARY={binary}")
    print(f"pop={args.pop} gen={args.gen} tau={args.tau} (CF a.u.)")
    print("FLEXAIDDS_ELECTION_V135=1 INCLUDE_SINGLETONS=1")

    if not runner.is_file():
        print(f"FAIL: runner missing: {runner}", file=sys.stderr)
        return 1
    if not binary.is_file():
        print(f"FAIL: binary missing: {binary}", file=sys.stderr)
        return 1
    if not manifest.is_file():
        print(f"FAIL: manifest missing: {manifest}", file=sys.stderr)
        return 1

    # Refuse dual-launch into live claim OUT
    claim = results / "campaigns" / "C0_full85_claim_g2000_popmod_20260715"
    if out.resolve() == claim.resolve():
        print("REFUSE: v135 must not write claim OUT", file=sys.stderr)
        return 91

    if args.dry_run:
        print("DRY-RUN OK")
        return 0

    out.mkdir(parents=True, exist_ok=True)
    try:
        git = subprocess.check_output(
            ["git", "log", "--oneline", "-1"], cwd=root, text=True
        ).strip()
    except Exception:
        git = "unknown"

    env = dict(os.environ)
    env.update(
        {
            "FLEXAIDDS_BINARY": str(binary),
            "FLEXAIDDS_DATA_DIR": str(data_dir),
            "FLEXAIDDS_RESTARTS": "5",
            "FLEXAIDDS_PARALLEL_RESTARTS": "0",
            "FLEXAIDDS_EVAL_SCALE_DIHEDRAL": "1",
            "EVAL_SCALE_DIHEDRAL": "1",
            "FLEXAIDDS_BUDGET_SCALE": "1",
            "FLEXAIDDS_NATIVE_SEED_FRAC": "0",
            "FLEXAIDDS_SEED_ELITISM": "0",
            "FLEXAIDDS_ELECTION_V135": "1",
            "FLEXAIDDS_ELECTION_SCORE_TAU": str(args.tau),
            "FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS": "1",
            "OMP_NUM_THREADS": "1",
            "OMP_WAIT_POLICY": "passive",
        }
    )
    # Claim-style: no native seed / no DP override unless user sets
    for k in (
        "FLEXAIDDS_FORCE_SEED",
        "FLEXAIDDS_USE_DP",
        "FLEXAIDDS_FREQSEL",  # v70 regression; leave off unless user enables
    ):
        env.pop(k, None)

    prov = {
        "version": "v135_bcr_proxy_election",
        "launched_at": dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git": git,
        "description": (
            "Crystal-blind basin recovery election: include Frequency=1 heads; "
            "score temperature τ in CF a.u. (default 25). CF is a scoring proxy "
            "not kcal/mol. Default claim ranking unchanged when flag off."
        ),
        "protocol": {
            "FLEXAIDDS_ELECTION_V135": "1",
            "FLEXAIDDS_ELECTION_SCORE_TAU": args.tau,
            "FLEXAIDDS_ELECTION_INCLUDE_SINGLETONS": "1",
            "pop_base": args.pop,
            "gen": args.gen,
            "EVAL_SCALE_DIHEDRAL": 1,
            "mode": "defined-cleft-redock",
            "seed_elitism": 0,
        },
        "binary_sha256": sha256(binary),
        "runner_sha256": sha256(runner),
        "manifest": str(manifest),
        "output": str(out),
    }
    (out / "PROVENANCE.json").write_text(json.dumps(prov, indent=2) + "\n")

    cmd = [
        "caffeinate",
        "-i",
        "-s",
        str(runner),
        "--benchmark",
        f"crossdock_json:{manifest}",
        "--mode",
        "defined-cleft-redock",
        "--output",
        str(out) + "/",
        "--threads",
        "1",
        "--omp-threads",
        "1",
        "--ga-population",
        str(args.pop),
        "--ga-generations",
        str(args.gen),
        "--temperature",
        "298",
        "--job-timeout-seconds",
        "10800",
    ]

    log = out / "v135_launch.log"
    print(f"Launching: {' '.join(cmd)}")
    print(f"log={log}")
    with log.open("w") as lf:
        lf.write("CMD " + " ".join(cmd) + "\n")
        lf.flush()
        proc = subprocess.Popen(
            cmd,
            env=env,
            cwd=str(root),
            stdout=lf,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    (out / "v135.pid").write_text(str(proc.pid) + "\n")
    print(f"STARTED pid={proc.pid}")
    print("Does not dual-launch claim OUT. Monitor separately.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
