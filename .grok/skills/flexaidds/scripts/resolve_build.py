#!/usr/bin/env python3
"""
resolve_build.py

Resolve FlexAIDdS CMake build directories by profile and validate that
production binaries are present.

Profiles:
  default -> build/
  lto     -> build_lto/
  metal   -> build_metal/

Usage:
  python3 .grok/skills/flexaidds/scripts/resolve_build.py --profile lto
  python3 .grok/skills/flexaidds/scripts/resolve_build.py --check --profile lto
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parents[2]

PROFILE_DIRS: Dict[str, str] = {
    "default": "build",
    "lto": "build_lto",
    "metal": "build_metal",
}

REQUIRED_BINARIES = (
    "FlexAIDdS",
    "benchmark_datasets",
    "cavity_detect_cli",
)


def resolve_build_dir(profile: str, repo_root: Path) -> Path:
    if profile not in PROFILE_DIRS:
        raise ValueError(f"unknown profile: {profile}")
    return (repo_root / PROFILE_DIRS[profile]).resolve()


def read_cmake_cache_value(cache_path: Path, key: str) -> Optional[str]:
    if not cache_path.is_file():
        return None
    pattern = re.compile(rf"^{re.escape(key)}:(?:BOOL|STRING|FILEPATH|PATH|INTERNAL)=(?:(.*))$")
    try:
        text = cache_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for line in text.splitlines():
        match = pattern.match(line.strip())
        if match:
            return match.group(1).strip()
    return None


def check_build(profile: str, repo_root: Path) -> Tuple[bool, str]:
    build_dir = resolve_build_dir(profile, repo_root)
    cache_path = build_dir / "CMakeCache.txt"
    lines = [
        "=== FlexAIDdS build resolve (--check) ===",
        f"profile: {profile}",
        f"build_dir: {build_dir}",
    ]

    metal_value = read_cmake_cache_value(cache_path, "FLEXAIDS_USE_METAL")
    if metal_value is None:
        lines.append("FLEXAIDS_USE_METAL: <not found in CMakeCache.txt>")
    else:
        lines.append(f"FLEXAIDS_USE_METAL: {metal_value} (from CMakeCache.txt)")

    all_ok = True
    for name in REQUIRED_BINARIES:
        binary = build_dir / name
        if binary.exists() and binary.is_file():
            lines.append(f"  OK: {name}")
        else:
            lines.append(f"  FAIL: {name} (missing)")
            all_ok = False

    if all_ok:
        lines.append("CHECK: PASS")
    else:
        lines.append("CHECK: FAIL")

    return all_ok, "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="resolve_build",
        description=(
            "Resolve FlexAIDdS build directories by profile "
            "(default/build, lto/build_lto, metal/build_metal)."
        ),
    )
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILE_DIRS),
        default="default",
        help="Build profile directory selector (default: build/).",
    )
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=REPO_ROOT,
        help="Repository root containing build directories.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Validate FlexAIDdS, benchmark_datasets, and cavity_detect_cli exist; "
            "print FLEXAIDS_USE_METAL from CMakeCache.txt when present."
        ),
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    repo_root = args.repo_root.resolve()

    if args.check:
        ok, report = check_build(args.profile, repo_root)
        print(report)
        return 0 if ok else 1

    print(resolve_build_dir(args.profile, repo_root))
    return 0


if __name__ == "__main__":
    sys.exit(main())