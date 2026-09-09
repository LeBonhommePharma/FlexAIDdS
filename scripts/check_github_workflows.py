#!/usr/bin/env python3
"""Reject GitHub workflow YAML that GitHub itself will not parse.

GitHub fails the whole run in 0s (no jobs, workflow name becomes the file
path) when a workflow has duplicate mapping keys — most commonly a
duplicated job id. PyYAML's ``safe_load`` silently keeps the last
duplicate, so a parse gate that uses it cannot see the class of breakage
it exists to catch.

Usage:
  python3 scripts/check_github_workflows.py
  python3 scripts/check_github_workflows.py --root /path/to/repo
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "PyYAML is required: python3 -m pip install pyyaml"
    ) from exc

REPO_ROOT = Path(__file__).resolve().parent.parent


class UniqueKeyLoader(yaml.SafeLoader):
    """SafeLoader that errors on duplicate mapping keys."""


def _no_duplicates_constructor(loader, node, deep=False):
    mapping = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in mapping:
            mark = key_node.start_mark
            raise yaml.constructor.ConstructorError(
                None,
                None,
                f"duplicate YAML key {key!r} at line {mark.line + 1} "
                f"column {mark.column + 1}",
                mark,
            )
        mapping[key] = loader.construct_object(value_node, deep=deep)
    return mapping


UniqueKeyLoader.add_constructor(
    yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
    _no_duplicates_constructor,
)


def load_unique(text: str):
    """Parse YAML and raise on duplicate keys."""
    return yaml.load(text, Loader=UniqueKeyLoader)


def iter_github_yaml(root: Path) -> list[Path]:
    files: list[Path] = []
    workflows = root / ".github" / "workflows"
    if workflows.is_dir():
        files.extend(sorted(workflows.glob("*.yml")))
        files.extend(sorted(workflows.glob("*.yaml")))
    actions = root / ".github" / "actions"
    if actions.is_dir():
        files.extend(sorted(actions.glob("**/action.yml")))
        files.extend(sorted(actions.glob("**/action.yaml")))
    return files


def check_file(path: Path, root: Path) -> list[str]:
    rel = path.relative_to(root).as_posix()
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as exc:
        return [f"{rel}: cannot read: {exc}"]
    try:
        data = load_unique(text)
    except yaml.YAMLError as exc:
        return [f"{rel}: {exc}"]
    if not isinstance(data, dict):
        return [f"{rel}: parsed but is not a mapping"]
    errors: list[str] = []
    if path.parent.name == "workflows" or "/workflows/" in rel:
        # GitHub stores `on:` as a boolean under YAML 1.1; only `jobs` is
        # required for a runnable workflow file.
        if "jobs" not in data:
            errors.append(f"{rel}: parsed but has no 'jobs' mapping")
        else:
            jobs = data["jobs"]
            if not isinstance(jobs, dict) or not jobs:
                errors.append(f"{rel}: 'jobs' must be a non-empty mapping")
    elif path.name in {"action.yml", "action.yaml"}:
        if "runs" not in data:
            errors.append(f"{rel}: composite/javascript action missing 'runs'")
    return errors


def collect_errors(root: Path) -> list[str]:
    files = iter_github_yaml(root)
    if not files:
        return [f"{root}: no .github/workflows/*.yml files found"]
    errors: list[str] = []
    for path in files:
        errors.extend(check_file(path, root))
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=REPO_ROOT,
        help="repository root (default: inferred from this script)",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()
    errors = collect_errors(root)
    if errors:
        print("UNPARSEABLE GITHUB YAML:")
        print("\n".join(errors))
        return 1
    n = len(iter_github_yaml(root))
    print(f"ok: {n} GitHub YAML files parse with unique mapping keys")
    return 0


if __name__ == "__main__":
    sys.exit(main())
