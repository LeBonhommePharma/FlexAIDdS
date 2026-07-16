#!/usr/bin/env python3
"""Stage FlexAID / FlexAIDdS binaries for three-engine red-pair with claim checks.

Claim science requires a **true** A vs master split:
  - bin/A/FlexAID  = historical FlexAID pin (JCIM-era lineage)
  - bin/B/FlexAID  = current master rebuild
  - bin/C/FlexAID + FlexAIDdS = master (arm C FO@298K only after oracle PASS)

If A and B are the same SHA256, B0 is a **deterministic twin** of A (same STRTSEED
+ same code) and is **not** a scientific control.

Usage:
  # Stage master build into B and C; leave A untouched
  python3 scripts/stage_three_engine_bins.py --master-bin build/FlexAID \\
    --flexaidds-bin build/FlexAIDdS --dest ~/flexaidds_results/three_engine_entropy_q1/bin

  # Claim-mode check (exit 1 if A==B)
  python3 scripts/stage_three_engine_bins.py --check-claim-split \\
    --dest ~/flexaidds_results/three_engine_entropy_q1/bin

  # Force-stage same binary to A and B (diagnostic only; refuse with --claim)
  python3 scripts/stage_three_engine_bins.py --master-bin build/FlexAID --also-a

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_sha(repo: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "-C", str(repo), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def write_pin(arm_dir: Path, *, binary: Path, git: str, note: str) -> str:
    arm_dir.mkdir(parents=True, exist_ok=True)
    dest = arm_dir / "FlexAID"
    shutil.copy2(binary, dest)
    dest.chmod(0o755)
    digest = sha256_file(dest)
    (arm_dir / "SHA256").write_text(digest + "\n")
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    (arm_dir / "BUILD_PIN.md").write_text(
        f"# FlexAID pin (arm {arm_dir.name})\n"
        f"- built/staged: {ts}\n"
        f"- git: {git}\n"
        f"- sha256: {digest}\n"
        f"- source: {binary.resolve()}\n"
        f"- note: {note}\n"
    )
    return digest


def main(argv: Optional[list] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--dest",
        type=Path,
        default=Path(
            os.environ.get(
                "FLEXAIDDS_LOCAL_QUEUE",
                str(Path.home() / "flexaidds_results/three_engine_entropy_q1"),
            )
        )
        / "bin",
    )
    ap.add_argument("--master-bin", type=Path, default=None, help="Master FlexAID binary")
    ap.add_argument(
        "--historical-a-bin",
        type=Path,
        default=None,
        help="Historical FlexAID for arm A (required for claim split)",
    )
    ap.add_argument("--flexaidds-bin", type=Path, default=None)
    ap.add_argument(
        "--also-a",
        action="store_true",
        help="Also copy master into A (diagnostic only — breaks claim split)",
    )
    ap.add_argument(
        "--check-claim-split",
        action="store_true",
        help="Exit 0 only if A and B FlexAID SHAs differ",
    )
    ap.add_argument(
        "--claim",
        action="store_true",
        help="Refuse staging if A and B would be identical",
    )
    ap.add_argument(
        "--repo",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    args = ap.parse_args(argv)
    dest: Path = args.dest.expanduser().resolve()
    dest.mkdir(parents=True, exist_ok=True)
    git = git_sha(args.repo.expanduser())

    if args.check_claim_split:
        a = dest / "A" / "FlexAID"
        b = dest / "B" / "FlexAID"
        if not a.is_file() or not b.is_file():
            print("FAIL: missing bin/A or bin/B FlexAID", file=sys.stderr)
            return 1
        sa, sb = sha256_file(a), sha256_file(b)
        print(f"A={sa}")
        print(f"B={sb}")
        if sa == sb:
            print(
                "FAIL: A and B are the same binary — B0 is a deterministic twin of A; "
                "not a claim control. Stage a historical A pin with --historical-a-bin.",
                file=sys.stderr,
            )
            return 1
        print("OK: claim binary split (A ≠ B)")
        return 0

    if args.master_bin is None and args.historical_a_bin is None:
        print("ERROR: provide --master-bin and/or --historical-a-bin", file=sys.stderr)
        return 2

    digests = {}
    if args.historical_a_bin:
        p = args.historical_a_bin.expanduser()
        if not p.is_file():
            print(f"ERROR: missing {p}", file=sys.stderr)
            return 2
        digests["A"] = write_pin(
            dest / "A",
            binary=p,
            git=git,
            note="historical FlexAID pin for arm A (TEMPER0 CF)",
        )
        print(f"staged A {digests['A']}")

    if args.master_bin:
        p = args.master_bin.expanduser()
        if not p.is_file():
            print(f"ERROR: missing {p}", file=sys.stderr)
            return 2
        digests["B"] = write_pin(
            dest / "B",
            binary=p,
            git=git,
            note="master FlexAID for B0/B (TEMPER0 CF / TEMPER21 FO)",
        )
        print(f"staged B {digests['B']}")
        digests["C"] = write_pin(
            dest / "C",
            binary=p,
            git=git,
            note="master FlexAID for arm C FO@298K (oracle gate required)",
        )
        print(f"staged C FlexAID {digests['C']}")
        if args.also_a:
            digests["A"] = write_pin(
                dest / "A",
                binary=p,
                git=git,
                note="DIAGNOSTIC: master copied to A (NOT claim split)",
            )
            print(f"staged A (also-a) {digests['A']}")

    if args.flexaidds_bin:
        fb = args.flexaidds_bin.expanduser()
        if not fb.is_file():
            print(f"ERROR: missing {fb}", file=sys.stderr)
            return 2
        cdir = dest / "C"
        cdir.mkdir(parents=True, exist_ok=True)
        out = cdir / "FlexAIDdS"
        shutil.copy2(fb, out)
        out.chmod(0o755)
        print(f"staged C FlexAIDdS {sha256_file(out)}")

    a_path, b_path = dest / "A" / "FlexAID", dest / "B" / "FlexAID"
    if a_path.is_file() and b_path.is_file():
        sa, sb = sha256_file(a_path), sha256_file(b_path)
        status = {
            "A": sa,
            "B": sb,
            "identical": sa == sb,
            "claim_binary_split_ok": sa != sb,
            "git": git,
            "ts_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        }
        (dest / "BINARY_SPLIT.json").write_text(json.dumps(status, indent=2) + "\n")
        if sa == sb:
            print(
                "WARN: A and B SHAs identical — B0 is not an independent control",
                file=sys.stderr,
            )
            if args.claim:
                print("FAIL: --claim refuses identical A/B", file=sys.stderr)
                return 1
        else:
            print("OK: A ≠ B claim binary split")
    return 0


if __name__ == "__main__":
    sys.exit(main())
