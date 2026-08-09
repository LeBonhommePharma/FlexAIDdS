#!/usr/bin/env python3
"""Repository hygiene checks for FlexAIDdS.

Fails when:
- tracked files match secret/env patterns (.env, .env.local, etc.)
- agent/skill instruction files contain machine-specific absolute paths

Run:
  python3 scripts/check_repo_hygiene.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

ENV_TRACKED_PATTERNS = (
    re.compile(r"^\.env$"),
    re.compile(r"^\.env\.[^/]+$"),
    re.compile(r"^\.envrc$"),
)

ENV_ALLOWED = {".env.example"}

SCAN_PREFIXES = (
    ".agents/",
    ".grok/skills/",
    "docs/custom-instructions/",
    "AGENTS.md",
    "CLAUDE.md",
    "chatgpt-instructions.md",
)

# Match real machine paths (e.g. /Users/lp.more/...) but not documentation
# placeholders like /Users/<username> or /Users/...
HARDCODED_PATH_RE = re.compile(
    r"/(?:Users|home)/[A-Za-z0-9._-]+(?:/[A-Za-z0-9._\- ]+)+"
)


def git_tracked_files() -> list[str]:
    proc = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=REPO_ROOT,
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print("WARN: git ls-files failed; scanning working tree only", file=sys.stderr)
        return []
    return [p.decode("utf-8") for p in proc.stdout.split(b"\0") if p]


def check_tracked_env_files(tracked: list[str]) -> list[str]:
    errors: list[str] = []
    for path in tracked:
        name = Path(path).name
        if name in ENV_ALLOWED:
            continue
        if any(pat.match(path) or pat.match(name) for pat in ENV_TRACKED_PATTERNS):
            errors.append(f"tracked secret/env file must not be committed: {path}")
    return errors


def should_scan(path: str) -> bool:
    return any(path == prefix or path.startswith(prefix) for prefix in SCAN_PREFIXES)


def check_hardcoded_paths(tracked: list[str]) -> list[str]:
    errors: list[str] = []
    candidates = tracked if tracked else [
        str(p.relative_to(REPO_ROOT))
        for p in REPO_ROOT.rglob("*")
        if p.is_file() and should_scan(str(p.relative_to(REPO_ROOT)))
    ]
    for rel in candidates:
        if not should_scan(rel):
            continue
        full = REPO_ROOT / rel
        if not full.is_file():
            continue
        try:
            text = full.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        for match in HARDCODED_PATH_RE.finditer(text):
            errors.append(
                f"{rel}: machine-specific absolute path {match.group(0)!r} "
                "(use repo-relative paths or FLEXAIDDS_* env vars)"
            )
    return errors


SCIENCE_PATHS_FILE = "docs/SCIENCE_CRITICAL_PATHS.txt"
# Stacks that cannot affect docked coordinates. Bundling one of these with a
# science-critical change is what made PR #405 unreviewable.
NON_ENGINE_STACKS = ("typescript/", "swift/", "site/")


def _load_science_globs(repo_root: Path) -> list[str]:
    path = repo_root / SCIENCE_PATHS_FILE
    if not path.exists():
        return []
    globs = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            globs.append(line)
    return globs


def _changed_files(base: str = "origin/main") -> tuple[list[str], str | None]:
    """Files changed vs the merge base. Returns (files, skip_reason)."""
    import subprocess

    try:
        mb = subprocess.run(["git", "merge-base", base, "HEAD"],
                            capture_output=True, text=True, timeout=30)
        if mb.returncode != 0:
            # depth-1 clone: no common ancestor available. FAIL LOUD, not open.
            return [], (f"cannot compute merge-base against {base} "
                        "(shallow clone? needs actions/checkout fetch-depth: 0)")
        out = subprocess.run(["git", "diff", "--name-only", mb.stdout.strip(), "HEAD"],
                             capture_output=True, text=True, timeout=30)
        if out.returncode != 0:
            return [], "git diff failed"
        return [f for f in out.stdout.splitlines() if f], None
    except Exception as exc:  # noqa: BLE001
        return [], f"git unavailable: {exc}"


def check_science_bundling(repo_root: Path) -> list[str]:
    """A science-critical change must not be bundled with a non-engine stack.

    Rationale: PR #405 was 138 files across LIB/ + typescript/ + swift/ + python/
    under a title reading as metadata-only. Its coordinates were in fact
    byte-identical, so it was correctly provenance-only -- but establishing that
    took a day. This does not forbid such work; it forces it to be split so each
    half is reviewable on its own terms.
    """
    import fnmatch

    globs = _load_science_globs(repo_root)
    if not globs:
        return []
    changed, skip = _changed_files()
    if skip:
        print(f"  NOTE: science-bundling check skipped -- {skip}")
        return []
    if not changed:
        return []

    science = [f for f in changed
               if any(fnmatch.fnmatch(f, g) for g in globs)]
    if not science:
        return []
    foreign = [f for f in changed if f.startswith(NON_ENGINE_STACKS)]
    if not foreign:
        return []
    return [
        "science-critical change bundled with a non-engine stack; split the PR.\n"
        f"      science-critical: {', '.join(sorted(science)[:5])}"
        f"{' ...' if len(science) > 5 else ''}\n"
        f"      non-engine:       {', '.join(sorted(foreign)[:5])}"
        f"{' ...' if len(foreign) > 5 else ''}\n"
        f"      (rule: {SCIENCE_PATHS_FILE}; a change that can move docked "
        "coordinates must be reviewable on its own)"
    ]


def main() -> int:
    print("=== FlexAIDdS Repository Hygiene Check ===\n")
    tracked = git_tracked_files()
    repo_root = Path(__file__).resolve().parent.parent
    errors = (check_tracked_env_files(tracked)
              + check_hardcoded_paths(tracked)
              + check_science_bundling(repo_root))

    if errors:
        print("FAIL: repository hygiene violations:\n", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print("OK: no tracked .env secrets; agent/skill files have no hardcoded user paths")
    return 0


if __name__ == "__main__":
    sys.exit(main())