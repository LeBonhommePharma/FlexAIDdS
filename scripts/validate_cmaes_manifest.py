#!/usr/bin/env python3
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""Validate cmaes_ab_manifest.json schema for the locked-arch harness.

Exit codes:
  0 — schema and invariants OK
  1 — validation failure (printed to stderr)
  2 — usage / I/O error
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any, Iterable, List, Mapping, Sequence

SCHEMA_VERSIONS = {"1.0"}
ABS_PATH_RE = re.compile(r"^(/Users/|/home/[^/]+/|/opt/homebrew/|[A-Za-z]:\\)")
REQUIRED_TOP = (
    "schema_version",
    "name",
    "complex",
    "eval_budget",
    "arms",
    "environment",
    "container",
)
REQUIRED_COMPLEX = ("code", "receptor", "ligand", "reference_ligand")
REQUIRED_BUDGET = ("target_evals", "population", "generations")
REQUIRED_ARM = ("id", "search", "array_task_id", "env", "output_prefix")
REQUIRED_CONTAINER = ("arch", "def", "sif", "binary")


def _err(msg: str, errors: List[str]) -> None:
    errors.append(msg)


def _require_mapping(obj: Any, path: str, errors: List[str]) -> Mapping[str, Any] | None:
    if not isinstance(obj, Mapping):
        _err(f"{path}: expected object, got {type(obj).__name__}", errors)
        return None
    return obj


def _require_keys(
    obj: Mapping[str, Any], keys: Iterable[str], path: str, errors: List[str]
) -> None:
    for k in keys:
        if k not in obj:
            _err(f"{path}: missing required key '{k}'", errors)


def _check_no_host_abs(value: str, path: str, errors: List[str]) -> None:
    """Repo-relative paths only for inputs; allow container-internal /opt/flexaidds/*."""
    if value.startswith("/opt/flexaidds/"):
        return
    if ABS_PATH_RE.match(value) or value.startswith("/Users/"):
        _err(f"{path}: machine-specific absolute path not allowed: {value!r}", errors)


def _check_rel_or_work(value: str, path: str, errors: List[str]) -> None:
    if not value or value.strip() != value:
        _err(f"{path}: empty or padded path", errors)
        return
    if value.startswith("/work/"):
        return
    if value.startswith("/"):
        # only allow /opt/flexaidds container paths elsewhere
        if not value.startswith("/opt/flexaidds/"):
            _err(f"{path}: must be repo-relative (or /work/…, /opt/flexaidds/…): {value!r}", errors)
        return
    if ".." in Path(value).parts:
        _err(f"{path}: parent traversal not allowed: {value!r}", errors)


def validate_manifest(data: Any, *, check_files: bool = False, root: Path | None = None) -> List[str]:
    errors: List[str] = []
    obj = _require_mapping(data, "$", errors)
    if obj is None:
        return errors

    _require_keys(obj, REQUIRED_TOP, "$", errors)

    ver = obj.get("schema_version")
    if ver not in SCHEMA_VERSIONS:
        _err(f"$.schema_version: unsupported {ver!r}; allowed={sorted(SCHEMA_VERSIONS)}", errors)

    # --- complex ---
    cx = _require_mapping(obj.get("complex"), "$.complex", errors)
    if cx is not None:
        _require_keys(cx, REQUIRED_COMPLEX, "$.complex", errors)
        code = cx.get("code")
        if not isinstance(code, str) or not re.fullmatch(r"[0-9][A-Za-z0-9]{3}", code or ""):
            _err(f"$.complex.code: expected PDB-like code (e.g. 1G9V), got {code!r}", errors)
        for key in ("receptor", "ligand", "reference_ligand"):
            if key in cx:
                if not isinstance(cx[key], str):
                    _err(f"$.complex.{key}: expected string", errors)
                else:
                    _check_no_host_abs(cx[key], f"$.complex.{key}", errors)
                    _check_rel_or_work(cx[key], f"$.complex.{key}", errors)
                    if check_files and root is not None:
                        p = root / cx[key]
                        if not p.is_file():
                            _err(f"$.complex.{key}: file not found under root: {p}", errors)

    # --- eval budget ---
    eb = _require_mapping(obj.get("eval_budget"), "$.eval_budget", errors)
    if eb is not None:
        _require_keys(eb, REQUIRED_BUDGET, "$.eval_budget", errors)
        try:
            pop = int(eb["population"])  # type: ignore[index]
            gen = int(eb["generations"])  # type: ignore[index]
            target = int(eb["target_evals"])  # type: ignore[index]
        except (KeyError, TypeError, ValueError) as exc:
            _err(f"$.eval_budget: population/generations/target_evals must be ints ({exc})", errors)
        else:
            if pop <= 0 or gen <= 0:
                _err("$.eval_budget: population and generations must be > 0", errors)
            product = pop * gen
            if product != target:
                _err(
                    f"$.eval_budget: population*generations={product} != target_evals={target}",
                    errors,
                )
            # Claim budget for GA vs CMA-ES A/B
            if target != 2_000_000:
                _err(
                    f"$.eval_budget.target_evals: expected 2000000 for claim A/B, got {target}",
                    errors,
                )
            if pop != 1000 or gen != 2000:
                _err(
                    "$.eval_budget: expected population=1000 and generations=2000 "
                    f"(got pop={pop} gen={gen})",
                    errors,
                )

    # --- arms ---
    arms = obj.get("arms")
    if not isinstance(arms, list) or len(arms) < 2:
        _err("$.arms: expected a list with at least 2 arms (ga, cmaes)", errors)
    else:
        seen_ids: set[str] = set()
        seen_search: set[str] = set()
        seen_tasks: set[int] = set()
        for i, arm in enumerate(arms):
            path = f"$.arms[{i}]"
            am = _require_mapping(arm, path, errors)
            if am is None:
                continue
            _require_keys(am, REQUIRED_ARM, path, errors)
            aid = am.get("id")
            search = am.get("search")
            if not isinstance(aid, str) or not aid:
                _err(f"{path}.id: non-empty string required", errors)
            elif aid in seen_ids:
                _err(f"{path}.id: duplicate id {aid!r}", errors)
            else:
                seen_ids.add(aid)
            if search not in ("ga", "cmaes"):
                _err(f"{path}.search: must be 'ga' or 'cmaes', got {search!r}", errors)
            else:
                seen_search.add(str(search))
            try:
                tid = int(am.get("array_task_id"))  # type: ignore[arg-type]
            except (TypeError, ValueError):
                _err(f"{path}.array_task_id: must be int", errors)
            else:
                if tid in seen_tasks:
                    _err(f"{path}.array_task_id: duplicate task id {tid}", errors)
                seen_tasks.add(tid)
            env = am.get("env")
            em = _require_mapping(env, f"{path}.env", errors)
            if em is not None:
                if "FLEXAIDDS_SEARCH" not in em:
                    _err(f"{path}.env: missing FLEXAIDDS_SEARCH", errors)
                elif str(em["FLEXAIDDS_SEARCH"]) != str(search):
                    _err(
                        f"{path}.env.FLEXAIDDS_SEARCH={em['FLEXAIDDS_SEARCH']!r} "
                        f"!= search={search!r}",
                        errors,
                    )
            op = am.get("output_prefix")
            if isinstance(op, str):
                _check_no_host_abs(op, f"{path}.output_prefix", errors)
                _check_rel_or_work(op, f"{path}.output_prefix", errors)
            else:
                _err(f"{path}.output_prefix: expected string", errors)
        if "ga" not in seen_search or "cmaes" not in seen_search:
            _err("$.arms: must include both search='ga' and search='cmaes'", errors)

    # --- environment ---
    env = _require_mapping(obj.get("environment"), "$.environment", errors)
    if env is not None:
        for key in ("FLEXAID_SEED", "FLEXAIDDS_NO_SEC", "FLEXAIDDS_DATA_DIR"):
            if key not in env:
                _err(f"$.environment: missing '{key}'", errors)
        data_dir = env.get("FLEXAIDDS_DATA_DIR")
        if isinstance(data_dir, str):
            if not data_dir.startswith("/opt/flexaidds/"):
                _err(
                    "$.environment.FLEXAIDDS_DATA_DIR: expected container path "
                    f"under /opt/flexaidds/, got {data_dir!r}",
                    errors,
                )
        seed = env.get("FLEXAID_SEED")
        if seed is not None and str(seed) != "12345":
            _err(
                f"$.environment.FLEXAID_SEED: claim seed is 12345, got {seed!r}",
                errors,
            )

    # --- container ---
    ct = _require_mapping(obj.get("container"), "$.container", errors)
    if ct is not None:
        _require_keys(ct, REQUIRED_CONTAINER, "$.container", errors)
        arch = ct.get("arch")
        if arch not in ("x86_64", "amd64"):
            _err(f"$.container.arch: expected x86_64 (locked-arch), got {arch!r}", errors)
        for key in ("def", "sif"):
            val = ct.get(key)
            if isinstance(val, str):
                _check_no_host_abs(val, f"$.container.{key}", errors)
                _check_rel_or_work(val, f"$.container.{key}", errors)
            else:
                _err(f"$.container.{key}: expected string", errors)
        binary = ct.get("binary")
        if isinstance(binary, str):
            if not binary.startswith("/opt/flexaidds/"):
                _err(
                    f"$.container.binary: expected /opt/flexaidds/... path, got {binary!r}",
                    errors,
                )
        toolchain = ct.get("toolchain")
        if toolchain is not None and str(toolchain) != "locked":
            _err(f"$.container.toolchain: expected 'locked', got {toolchain!r}", errors)

    # --- optional ga_config consistency ---
    ga = obj.get("ga_config")
    if ga is not None:
        gm = _require_mapping(ga, "$.ga_config", errors)
        if gm is not None and eb is not None:
            try:
                if int(gm.get("num_chromosomes", -1)) != int(eb["population"]):  # type: ignore[index]
                    _err("$.ga_config.num_chromosomes must equal eval_budget.population", errors)
                if int(gm.get("num_generations", -1)) != int(eb["generations"]):  # type: ignore[index]
                    _err("$.ga_config.num_generations must equal eval_budget.generations", errors)
            except (TypeError, ValueError, KeyError):
                _err("$.ga_config: num_chromosomes/num_generations must be ints", errors)

    # --- optional slurm block ---
    slurm = obj.get("slurm")
    if slurm is not None:
        sm = _require_mapping(slurm, "$.slurm", errors)
        if sm is not None:
            if sm.get("account_env") != "CC_ACCOUNT":
                _err("$.slurm.account_env: expected 'CC_ACCOUNT'", errors)
            arr = sm.get("array")
            if arr not in ("0-1", "0-1%2"):
                # soft check — must cover at least two tasks
                if not (isinstance(arr, str) and arr.startswith("0-")):
                    _err(f"$.slurm.array: unexpected array spec {arr!r}", errors)

    return errors


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate CMA-ES A/B harness manifest (schema + 2e6 budget invariants)."
    )
    parser.add_argument(
        "manifest",
        nargs="?",
        default=str(Path(__file__).resolve().parent / "cmaes_ab_manifest.json"),
        help="Path to cmaes_ab_manifest.json (default: alongside this script)",
    )
    parser.add_argument(
        "--check-files",
        action="store_true",
        help="Also verify complex input files exist relative to --root",
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Repo root for --check-files (default: walk up for CMakeLists.txt)",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    path = Path(args.manifest)
    if not path.is_file():
        print(f"ERROR: manifest not found: {path}", file=sys.stderr)
        return 2

    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        print(f"ERROR: invalid JSON: {exc}", file=sys.stderr)
        return 2

    root = args.root
    if args.check_files and root is None:
        here = path.resolve().parent
        root = None
        for cand in [here, *here.parents]:
            if (cand / "CMakeLists.txt").is_file():
                root = cand
                break
        if root is None:
            print(
                "ERROR: --check-files set but repo root (CMakeLists.txt) not found; pass --root",
                file=sys.stderr,
            )
            return 2

    errors = validate_manifest(data, check_files=args.check_files, root=root)
    if errors:
        print(f"FAIL: {path} ({len(errors)} error(s))", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        return 1

    eb = data["eval_budget"]
    print(
        f"OK: {path} schema_version={data.get('schema_version')} "
        f"complex={data['complex'].get('code')} "
        f"budget={eb['population']}x{eb['generations']}={eb['target_evals']} "
        f"arms={len(data['arms'])}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
