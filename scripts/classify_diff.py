#!/usr/bin/env python3
"""Classify a git diff into engine / pack / tests / hygiene buckets.

Use this before reviewing a PR or a `git diff A B` that looks like "mostly
science docs." It does not decide Fix vs Add (that is the commit prefix);
it decides *what kind of tree* moved so a mixed PR can be split.

Examples:
  python3 scripts/classify_diff.py
  python3 scripts/classify_diff.py origin/main HEAD
  python3 scripts/classify_diff.py pre-swarm-baseline-20260813 origin/main
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

import check_repo_hygiene as hygiene  # noqa: E402

REPO_ROOT = hygiene.REPO_ROOT


def _name_only(a: str, b: str) -> list[str]:
    proc = subprocess.run(
        ["git", "diff", "--name-only", a, b],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    if proc.returncode != 0:
        raise SystemExit(proc.stderr.strip() or f"git diff failed: {a} {b}")
    return [ln for ln in proc.stdout.splitlines() if ln]


def format_report(buckets: dict[str, list[str]]) -> str:
    lines = ["=== change class (path buckets) ===", ""]
    nonempty = [(k, v) for k, v in buckets.items() if v]
    if not nonempty:
        return "=== change class (path buckets) ===\n(no files)\n"
    for key, files in nonempty:
        lines.append(f"{key} ({len(files)})")
        for rel in files:
            lines.append(f"  {rel}")
        lines.append("")
    engine = buckets.get("engine_critical") or []
    pack = buckets.get("science_pack") or []
    if engine and pack:
        lines.append("VERDICT: MIXED — split engine/bugfix from swarm/audit pack")
    elif engine:
        lines.append("VERDICT: engine-critical (can move coordinates or CF ranking)")
    elif pack:
        lines.append("VERDICT: science pack (docs/swarm or docs/audit; not engine)")
    elif buckets.get("benchmark"):
        lines.append("VERDICT: benchmark harness")
    elif buckets.get("ci_hygiene"):
        lines.append("VERDICT: CI / hygiene")
    elif buckets.get("docs_other"):
        lines.append("VERDICT: docs")
    else:
        lines.append("VERDICT: non-engine")
    lines.append("Intent (Fix vs Add vs Docs) is the commit prefix, not this tool.")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "commits",
        nargs="*",
        help="git diff A B (default: merge-base origin/main HEAD .. HEAD)",
    )
    parser.add_argument(
        "--strict-split",
        action="store_true",
        help="exit 1 if engine_critical and science_pack both nonempty",
    )
    args = parser.parse_args(argv)

    if len(args.commits) == 2:
        files = _name_only(args.commits[0], args.commits[1])
    elif len(args.commits) == 0:
        files, skip = hygiene._changed_files("origin/main")
        if skip:
            print(f"NOTE: {skip}", file=sys.stderr)
            files = []
    else:
        parser.error("give zero args or exactly two commits")

    globs = hygiene._load_science_globs(REPO_ROOT)
    buckets = hygiene.classify_paths(files, globs)
    print(format_report(buckets), end="")
    if args.strict_split and buckets["engine_critical"] and buckets["science_pack"]:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
