#!/usr/bin/env python3
"""
validate_dataset_semantics.py — Fail-closed self-docking vs cross-docking YAML gate.

Every dataset YAML under benchmarks/datasets/ (and the Python package mirror when
present) must declare an explicit docking_mode. Metrics and structural_states
must not contradict that mode.

Modes (normative):
  self_docking       — redock cognate ligand into native holo receptor
  cross_docking      — dock into non-native / apo / alternative receptor
  affinity_scoring   — score/rank vs experimental affinity (pose recovery optional)
  virtual_screening  — actives vs decoys enrichment
  specialized        — custom protocols (must still not claim crossdock metrics
                       unless structural evidence supports it)

Usage:
  python3 .grok/skills/flexaidds/scripts/validate_dataset_semantics.py
  python3 .grok/skills/flexaidds/scripts/validate_dataset_semantics.py --path benchmarks/datasets
  python3 .grok/skills/flexaidds/scripts/validate_dataset_semantics.py --strict  # default
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:  # pragma: no cover
    print("PyYAML required: pip install pyyaml", file=sys.stderr)
    sys.exit(2)

VALID_MODES = frozenset(
    {
        "self_docking",
        "cross_docking",
        "affinity_scoring",
        "virtual_screening",
        "specialized",
    }
)

CROSS_STATE_TOKENS = frozenset({"apo", "alternative", "crossdock", "non_native", "nonnative"})
CROSS_METRIC_PREFIXES = (
    "crossdock_",
    "cross_dock_",
    "nonnative_",
    "non_native_",
)
CROSS_SLUG_MARKERS = (
    "nonnative",
    "non_native",
    "crossdock",
    "cross_dock",
    "posex_cd",
    "hap2",
)
SELF_SLUG_MARKERS = (
    "astex_diverse",
)


def _is_cross_metric(name: str) -> bool:
    n = name.lower()
    return any(n.startswith(p) for p in CROSS_METRIC_PREFIXES) or "crossdock" in n


def discover_repo_root(start: Path) -> Path | None:
    cur = start.resolve()
    for _ in range(12):
        if (cur / ".git").exists() or (cur / "AGENTS.md").is_file():
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    return None


def load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"{path}: root must be a mapping")
    return data


def validate_config(raw: dict[str, Any], *, path: Path) -> list[str]:
    """Return list of error strings (empty = OK)."""
    errors: list[str] = []
    slug = str(raw.get("slug") or path.stem).strip().lower()
    mode = raw.get("docking_mode")
    if mode is None or str(mode).strip() == "":
        errors.append(f"{path.name}: missing required field docking_mode")
        return errors
    mode = str(mode).strip().lower()
    if mode not in VALID_MODES:
        errors.append(
            f"{path.name}: docking_mode={mode!r} invalid; expected one of {sorted(VALID_MODES)}"
        )
        return errors

    states = [str(s).strip().lower() for s in (raw.get("structural_states") or ["holo"])]
    metrics = [str(m).strip() for m in (raw.get("metrics") or [])]
    cross_metrics = [m for m in metrics if _is_cross_metric(m)]
    has_cross_state = any(s in CROSS_STATE_TOKENS for s in states)

    # Slug vs mode hard bans
    if any(m in slug for m in SELF_SLUG_MARKERS) and mode == "cross_docking":
        errors.append(
            f"{path.name}: slug {slug!r} is a self-docking benchmark but docking_mode=cross_docking"
        )
    if any(m in slug for m in CROSS_SLUG_MARKERS) and mode == "self_docking":
        errors.append(
            f"{path.name}: slug {slug!r} indicates cross-docking but docking_mode=self_docking"
        )

    # Metrics vs mode
    if cross_metrics and mode != "cross_docking":
        errors.append(
            f"{path.name}: cross-dock metrics {cross_metrics} require docking_mode=cross_docking "
            f"(got {mode})"
        )
    if mode == "cross_docking" and not cross_metrics and not has_cross_state:
        errors.append(
            f"{path.name}: docking_mode=cross_docking needs crossdock_* metrics "
            f"and/or structural_states in {sorted(CROSS_STATE_TOKENS)}"
        )

    # States vs mode
    if has_cross_state and mode == "self_docking":
        errors.append(
            f"{path.name}: structural_states {states} include non-native tokens but "
            f"docking_mode=self_docking (use cross_docking)"
        )
    if mode == "self_docking" and states and set(states) - {"holo"}:
        # Allow only holo for pure self-docking
        extra = sorted(set(states) - {"holo"})
        errors.append(
            f"{path.name}: self_docking allows structural_states=['holo'] only; found {extra}"
        )

    # Name/description soft signals elevated to errors when contradictory
    blob = " ".join(
        str(raw.get(k, "")) for k in ("name", "description")
    ).lower()
    if mode == "self_docking" and ("cross-dock" in blob or "cross dock" in blob or "non-native" in blob):
        if "self-dock" not in blob and "native holo" not in blob and "redock" not in blob:
            errors.append(
                f"{path.name}: description/name suggests cross-docking but docking_mode=self_docking"
            )
    if mode == "cross_docking" and "self-dock" in blob and "cross" not in blob:
        errors.append(
            f"{path.name}: description suggests self-docking but docking_mode=cross_docking"
        )

    return errors


def validate_path(directory: Path) -> tuple[int, list[str]]:
    errors: list[str] = []
    yamls = sorted(directory.glob("*.yaml"))
    if not yamls:
        return 1, [f"No *.yaml found under {directory}"]
    for path in yamls:
        try:
            raw = load_yaml(path)
        except Exception as exc:  # noqa: BLE001 — surface all load failures
            errors.append(f"{path.name}: YAML load failed: {exc}")
            continue
        errors.extend(validate_config(raw, path=path))
    return (1 if errors else 0), errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--path",
        type=Path,
        action="append",
        dest="paths",
        help="Directory of dataset YAMLs (repeatable). Default: repo benchmarks/datasets",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        default=True,
        help="Exit non-zero on any error (default)",
    )
    parser.add_argument(
        "--no-strict",
        dest="strict",
        action="store_false",
        help="Report errors but exit 0",
    )
    args = parser.parse_args()

    repo = discover_repo_root(Path(__file__))
    paths = args.paths
    if not paths:
        if not repo:
            print("Cannot find repo root; pass --path", file=sys.stderr)
            return 2
        paths = [repo / "benchmarks" / "datasets"]
        pkg = repo / "python" / "flexaidds" / "dataset_runner" / "datasets"
        if pkg.is_dir():
            paths.append(pkg)

    all_errors: list[str] = []
    worst = 0
    for directory in paths:
        directory = directory.expanduser().resolve()
        if not directory.is_dir():
            all_errors.append(f"Not a directory: {directory}")
            worst = 1
            continue
        code, errs = validate_path(directory)
        if code:
            worst = 1
        if errs:
            print(f"=== {directory} ===")
            for e in errs:
                print(f"  FAIL: {e}")
            all_errors.extend(errs)
        else:
            n = len(list(directory.glob("*.yaml")))
            print(f"OK: {directory} ({n} YAML files, docking_mode semantics consistent)")

    if all_errors:
        print(f"\n{len(all_errors)} semantic error(s).", file=sys.stderr)
        return 1 if args.strict else 0
    print("\nVALIDATION PASSED: dataset docking_mode semantics OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
