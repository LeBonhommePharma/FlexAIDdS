#!/usr/bin/env python3
"""
FlexAIDdS Source File Guard — Python implementation

Ensures that every source file that looks like it should be compiled is
actually referenced somewhere in the build system (CMakeLists.txt, *.cmake,
setup.py, pybind11 bindings, etc.).

This is the primary "idiot-proof" mechanism for the repository. It catches
the most common cause of mysterious link failures when adding new .cpp/.cu/.mm
files.

Usage (standalone):
    python scripts/validate_sources.py --root . --strict
    python scripts/validate_sources.py --warn-only

Usage from CMake (via cmake/ValidateSources.cmake):
    python scripts/validate_sources.py --root /path/to/repo --cmake-mode --strict

Usage from Python packaging (setup.py / pip install -e .):
    from scripts.validate_sources import validate_sources
    validate_sources(root=".", strict=True)

Exit codes:
    0 = clean (or only warnings)
    1 = orphaned sources found (in --strict mode)

The script is deliberately heuristic but extremely effective in practice for
this codebase. It trades theoretical perfection for zero configuration pain
and immediate value.
"""

from __future__ import annotations

import argparse
import fnmatch
import json
import os
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Set, List, Dict, Optional

# =============================================================================
# Configuration — easy to extend
# =============================================================================

# File extensions that represent compilable / buildable units we care about.
# Order does not matter. We focus on things that must appear in a target.
SOURCE_EXTENSIONS: Set[str] = {
    ".cpp", ".cc", ".c", ".cu", ".cuh",
    ".mm", ".metal", ".hip", ".h", ".hpp", ".hxx", ".inl"
}

# Directories (relative to repo root) that we recursively scan for sources.
DEFAULT_SCAN_DIRS: List[str] = [
    "LIB",
    "src",
    "tests",
    "python/bindings",
    "tools",
]

# Build definition files we scan for references to source files.
# We look inside these for any mention of the filenames.
BUILD_DEFINITION_PATTERNS: List[str] = [
    "**/CMakeLists.txt",
    "**/*.cmake",
    "setup.py",
    "pyproject.toml",
    "python/setup.py",
    "python/pyproject.toml",
    "python/bindings/*.cpp",   # pybind11 glue often lists sources
]

# Files that are allowed to exist without being built (legacy, data, examples,
# intentionally disabled code, etc.). Use gitignore-style globs.
# One pattern per line. Lines starting with # are comments.
DEFAULT_IGNORE_FILE = "build_sources.ignore"

# Very common legacy / special files that have historically been problematic
# or are intentionally not part of the main build.
BUILT_IN_IGNORES: List[str] = [
    "LIB/vendor/**",
    "LIB/old/**",
    "LIB/wif083.cpp",           # historical internal tool
    "LIB/python_bindings.cpp",  # legacy pybind11 stub (intentionally not used)
    "**/__pycache__/**",
    "**/build/**",
    "**/CMakeFiles/**",
    "site/**",                  # gh-pages static assets
]


@dataclass
class ValidationResult:
    orphans: List[Path] = field(default_factory=list)
    ignored_orphans: List[Path] = field(default_factory=list)
    scanned_files: int = 0
    referenced_files: Set[str] = field(default_factory=set)
    errors: List[str] = field(default_factory=list)


# =============================================================================
# Core logic
# =============================================================================

def find_candidate_sources(root: Path, extra_dirs: Iterable[str] = ()) -> List[Path]:
    """Recursively find all files with source extensions under the scan directories."""
    candidates: List[Path] = []
    scan_roots = list(DEFAULT_SCAN_DIRS) + list(extra_dirs)

    for rel_dir in scan_roots:
        base = root / rel_dir
        if not base.exists():
            continue
        for ext in SOURCE_EXTENSIONS:
            # Use rglob for clarity; performance is irrelevant here
            for p in base.rglob(f"*{ext}"):
                if p.is_file():
                    candidates.append(p)

    # Deduplicate while preserving order
    seen = set()
    unique = []
    for p in candidates:
        if p not in seen:
            seen.add(p)
            unique.append(p)
    return unique


def load_ignore_patterns(root: Path, ignore_file: Optional[Path] = None) -> List[str]:
    """Load ignore patterns from file + built-in list."""
    patterns: List[str] = list(BUILT_IN_IGNORES)

    if ignore_file is None:
        ignore_file = root / DEFAULT_IGNORE_FILE

    if ignore_file.exists():
        for line in ignore_file.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                patterns.append(line)

    return patterns


def is_ignored(path: Path, root: Path, patterns: List[str]) -> bool:
    """Return True if the path (relative to root) matches any ignore pattern."""
    try:
        rel = path.relative_to(root).as_posix()
    except ValueError:
        rel = path.as_posix()

    for pat in patterns:
        # Support both forward and backslash, and ** globs
        if fnmatch.fnmatch(rel, pat) or fnmatch.fnmatch(rel, pat.replace("\\", "/")):
            return True
        # Also try matching against just the basename for convenience
        if fnmatch.fnmatch(path.name, pat):
            return True
    return False


def collect_referenced_names(root: Path) -> Set[str]:
    """
    Scan all build definition files and collect every token that looks like
    a source file we care about.

    This is deliberately simple and greedy: if "MyNewKernel.cpp" appears
    anywhere in a CMakeLists.txt or setup.py, we consider it referenced.
    This catches 99% of real mistakes with almost zero false positives for
    the "I forgot to add my file" scenario.
    """
    referenced: Set[str] = set()

    # Build a big list of files to search
    search_files: List[Path] = []
    for pattern in BUILD_DEFINITION_PATTERNS:
        search_files.extend(root.glob(pattern))
        search_files.extend(root.glob("**/" + pattern))  # be thorough

    # Remove duplicates
    search_files = list({p.resolve() for p in search_files if p.is_file()})

    # Regex that matches typical source file references
    # Examples matched:
    #   Foo.cpp   LIB/Bar.cu   "something/Widget.mm"   MyKernel.cuh
    source_regex = re.compile(
        r'["\']?([A-Za-z0-9_./\\-]+\.(?:' +
        "|".join(re.escape(ext.lstrip(".")) for ext in SOURCE_EXTENSIONS) +
        r'))["\']?'
    )

    for f in search_files:
        try:
            text = f.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue

        for match in source_regex.finditer(text):
            name = match.group(1)
            # Normalize to forward slashes and basename for matching
            norm = name.replace("\\", "/")
            referenced.add(Path(norm).name)           # basename
            referenced.add(norm)                      # relative path form
            # Also store without leading LIB/ etc. for robustness
            if "/" in norm:
                referenced.add(norm.split("/")[-1])

    return referenced


def validate_sources(
    root: Path | str = ".",
    *,
    strict: bool = False,
    warn_only: bool = False,
    ignore_file: Optional[Path | str] = None,
    extra_scan_dirs: Iterable[str] = (),
    cmake_mode: bool = False,
) -> ValidationResult:
    """
    Main validation routine. Returns a ValidationResult.
    """
    root = Path(root).resolve()
    result = ValidationResult()

    candidates = find_candidate_sources(root, extra_scan_dirs)
    result.scanned_files = len(candidates)

    ignore_patterns = load_ignore_patterns(root, Path(ignore_file) if ignore_file else None)

    referenced = collect_referenced_names(root)
    result.referenced_files = referenced

    for cand in candidates:
        name = cand.name
        rel_str = cand.relative_to(root).as_posix() if cand.is_relative_to(root) else cand.as_posix()

        is_ref = (name in referenced) or (rel_str in referenced)

        if is_ref:
            continue

        if is_ignored(cand, root, ignore_patterns):
            result.ignored_orphans.append(cand)
            continue

        result.orphans.append(cand)

    return result


def format_report(result: ValidationResult, root: Path, strict: bool, cmake_mode: bool) -> str:
    """Produce human-readable (and CI-friendly) output."""
    lines: List[str] = []

    if not result.orphans:
        if cmake_mode:
            lines.append("FlexAID source validator: clean (no orphaned source files).")
        else:
            lines.append("✅  FlexAIDdS source guard — clean")
            lines.append(f"   Scanned {result.scanned_files} candidate source files. No orphans found.")
        return "\n".join(lines)

    # We have orphans
    severity = "ERROR" if strict else "WARNING"
    lines.append(f"{severity}: FlexAIDdS source guard found orphaned source files!")
    lines.append("")
    lines.append(f"  Scanned : {result.scanned_files} files")
    lines.append(f"  Orphans : {len(result.orphans)}")
    if result.ignored_orphans:
        lines.append(f"  (Ignored: {len(result.ignored_orphans)} via build_sources.ignore)")
    lines.append("")
    lines.append("Orphaned files (these exist on disk but are not referenced in any build file):")
    lines.append("")

    for orphan in sorted(result.orphans):
        try:
            rel = orphan.relative_to(root)
        except ValueError:
            rel = orphan

        lines.append(f"  • {rel}")

        # Give a helpful hint
        parent = rel.parent
        if rel.suffix in {".h", ".hpp", ".hxx", ".cuh"}:
            # Headers are often included transitively; suggest the module dir
            suggested = f"LIB/{parent.name}/CMakeLists.txt (if this header belongs to that module)" if parent.name else "the module's CMakeLists.txt"
        elif "LIB" in str(parent):
            suggested = f"LIB/{parent.name}/CMakeLists.txt" if parent.name else "LIB/CMakeLists.txt (or create a new module dir)"
        else:
            suggested = "the appropriate CMakeLists.txt or setup.py"
        lines.append(f"    Suggested fix: Add '{rel.name}' to {suggested}")
        lines.append("")

    lines.append("How to fix:")
    lines.append("  1. If this is a new module, create LIB/YourModule/CMakeLists.txt and list the files there.")
    lines.append("  2. Add one line `add_subdirectory(YourModule)` in LIB/CMakeLists.txt.")
    lines.append("  3. Or add the file to an existing module's CMakeLists.txt.")
    lines.append("  4. Re-run cmake (or pip install -e .).")
    lines.append("")
    lines.append("To temporarily allow an orphan (migration / legacy):")
    lines.append(f"  echo '{rel}' >> build_sources.ignore")
    lines.append("")

    if not strict:
        lines.append("This is a non-fatal warning. Switch to --strict (or STRICT in CMake) once clean.")

    return "\n".join(lines)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="FlexAIDdS build source validator — prevents missing sources in CMake targets."
    )
    parser.add_argument("--root", default=".", help="Repository root (default: current directory)")
    parser.add_argument("--strict", action="store_true", help="Fail with non-zero exit if orphans exist")
    parser.add_argument("--warn-only", action="store_true", help="Never fail, only warn (useful during migration)")
    parser.add_argument("--ignore-file", default=None, help="Path to ignore file (default: build_sources.ignore)")
    parser.add_argument("--extra-scan-dir", action="append", default=[], dest="extra_scan_dirs",
                        help="Additional directory to scan (can be repeated)")
    parser.add_argument("--cmake-mode", action="store_true",
                        help="Produce shorter, CMake-friendly output and always exit 0 unless strict+orphans")
    parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON instead of text")
    args = parser.parse_args(argv)

    root = Path(args.root).resolve()
    strict = args.strict and not args.warn_only

    result = validate_sources(
        root=root,
        strict=strict,
        warn_only=args.warn_only,
        ignore_file=args.ignore_file,
        extra_scan_dirs=args.extra_scan_dirs,
        cmake_mode=args.cmake_mode,
    )

    if args.json:
        payload = {
            "orphans": [str(p.relative_to(root)) for p in result.orphans],
            "ignored": [str(p.relative_to(root)) for p in result.ignored_orphans],
            "scanned": result.scanned_files,
            "clean": len(result.orphans) == 0,
        }
        print(json.dumps(payload, indent=2))
        return 0 if not strict or not result.orphans else 1

    report = format_report(result, root, strict, args.cmake_mode)
    print(report)

    if result.orphans and strict:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
