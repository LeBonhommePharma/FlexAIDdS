#!/usr/bin/env python3
"""Decide whether a PR can move Tier-1 docking_power_top1.

PoseBust is post-election physical validity (METHODOLOGY.md). Changing
LIB/PoseBust/** cannot change RMSD ranking, so the 4-target Astex quality
gate is noise against an aspirational 0.70 top-1 baseline.

The workflow still *starts* on LIB/** (required checks stay green rather
than path-filter skipped). This script then skips the dock.

Fail closed: if the diff cannot be classified, run the dock.

Exit 0 always (except argparse errors). Prints dock=true|false to stdout
and appends the same key to $GITHUB_OUTPUT when that env var is set.

Examples:
  python3 scripts/tier1_pr_needs_dock.py origin/main HEAD
  python3 scripts/tier1_pr_needs_dock.py --files LIB/PoseBust/Engine.cpp
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Paths that cannot move CF ranking / RMSD / search budget. A PR whose
# *entire* file list matches these is post-election (or docs/CI for the
# skip gate itself) and must not fire the 4-target Astex quality gate.
_SKIP_EXACT = frozenset(
    {
        "README.md",
        "scripts/tier1_pr_needs_dock.py",
        "tests/test_tier1_pr_needs_dock.py",
    }
)
_SKIP_PREFIXES = (
    "LIB/PoseBust/",
    "tests/",
    "python/tests/",
    "docs/",
    ".github/",
)


def is_post_election_path(rel: str) -> bool:
    """True when this path cannot move docking_power_top1."""
    rel = rel.replace("\\", "/")
    if rel in _SKIP_EXACT:
        return True
    return any(rel.startswith(prefix) for prefix in _SKIP_PREFIXES)


def files_need_dock(files: list[str]) -> bool:
    """True when any changed file can move RMSD ranking or the quality gate."""
    if not files:
        return False
    return any(not is_post_election_path(path) for path in files)


def _name_only(base: str, head: str) -> tuple[list[str], str | None]:
    proc = subprocess.run(
        ["git", "diff", "--name-only", base, head],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        timeout=30,
        check=False,
    )
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip() or f"git diff failed: {base} {head}"
        return [], err
    return [ln for ln in proc.stdout.splitlines() if ln], None


def emit_dock(dock: bool, reason: str) -> None:
    value = "true" if dock else "false"
    print(f"dock={value}")
    print(reason)
    github_output = os.environ.get("GITHUB_OUTPUT")
    if not github_output:
        return
    with open(github_output, "a", encoding="utf-8") as fh:
        fh.write(f"dock={value}\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "commits",
        nargs="*",
        help="git diff A B (default: origin/main HEAD)",
    )
    parser.add_argument(
        "--files",
        nargs="+",
        help="classify these paths instead of running git diff",
    )
    args = parser.parse_args(argv)

    if args.files is not None:
        files = args.files
        git_err = None
    elif len(args.commits) == 2:
        files, git_err = _name_only(args.commits[0], args.commits[1])
    elif len(args.commits) == 0:
        files, git_err = _name_only("origin/main", "HEAD")
    else:
        parser.error("give zero args, exactly two commits, or --files")

    if git_err:
        emit_dock(True, f"cannot classify diff ({git_err}); fail-closed to running the dock")
        return 0

    dock = files_need_dock(files)
    if dock:
        movers = [p for p in files if not is_post_election_path(p)]
        preview = ", ".join(movers[:8])
        extra = " ..." if len(movers) > 8 else ""
        emit_dock(True, f"PR can move docking_power_top1: {preview}{extra}")
    else:
        emit_dock(
            False,
            "PR is post-election / docs / PoseBust-only; skipping 4-target Astex dock "
            "(RMSD ranking unchanged; FlexAID C++ CI covers NativePoseQC)",
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
