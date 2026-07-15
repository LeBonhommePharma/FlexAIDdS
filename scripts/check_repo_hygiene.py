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


def main() -> int:
    print("=== FlexAIDdS Repository Hygiene Check ===\n")
    tracked = git_tracked_files()
    errors = check_tracked_env_files(tracked) + check_hardcoded_paths(tracked)

    if errors:
        print("FAIL: repository hygiene violations:\n", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        return 1

    print("OK: no tracked .env secrets; agent/skill files have no hardcoded user paths")
    return 0


if __name__ == "__main__":
    sys.exit(main())