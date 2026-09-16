# Copyright 2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
# SPDX-License-Identifier: Apache-2.0
"""Tier-1 Astex CI gate contract: sampling vs ranking, plus replay receipts.

This module is the executable form of ``docs/TIER1_CI_CONTRACT.md`` and
``METHODOLOGY.md`` §0.1.1. It does not retune the CF energy matrix and does
not change pose election. Ranking metrics stay measured; they do not fail
the PR until a registered election rule is recorded.

Hard (blocking) vs advisory is declared per dataset in ``ci_gate_class``.
Datasets that omit the map keep the historical behaviour: every
``expected_baselines`` key is a blocking gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# Production CF matrix pin (claim / comparative / CI). See METHODOLOGY.md
# and .grok/skills/flexaidds/SKILL.md. The 72d7 packing-sweetened fork is
# forbidden here. This module never writes the matrix file.
MATRIX_NAME = "MC_st0r5.2_6.dat"
MATRIX_PIN_MD5 = "9dc93717dfed0698006d88dd6a9627bc"

RECEIPT_SCHEMA = "flexaidds.tier1_receipt.v1"
RECEIPT_SCHEMA_VERSION = 1
PINNED_TIER1_SEED = 20260816
PINNED_ROSTER = ("1gpk", "1mq6", "1xm6", "2cet")
PINNED_ANCHORS = ("1gpk", "1mq6")
DEFAULT_ELECTION_OBJECTIVE = "cf_minus_ts"
VALID_GATE_CLASSES = frozenset({"hard", "advisory"})
VALID_RANKING_STATUS = frozenset({"unregistered", "registered"})
# Ranking / ΔS metrics that must stay advisory while the election rule is
# unregistered. Default-hard would reintroduce the 35016255184 false fail.
RANKING_QUALITY_METRICS = (
    "docking_power_top1",
    "docking_power_top3",
    "entropy_rescue_rate",
)

# Thesis direction for a *future* registered election rule. Not a CI weight
# and not a Cartesian ligand-entropy term on the CF proxy.
THESIS_ELECTION_DIRECTION = (
    "pose-density clustering for conformational entropy; "
    "torsional-consistent tENCoM for vibrational entropy. "
    "Cartesian ligand entropy is not a CF weight."
)

RECEIPT_REQUIRED_KEYS = (
    "schema",
    "schema_version",
    "seed",
    "roster",
    "anchors",
    "drawn",
    "matrix_name",
    "matrix_md5",
    "matrix_pin_md5",
    "matrix_pin_ok",
    "election_objective",
    "ci_ranking_status",
    "hard_gates",
    "advisory_ranking",
    "verdict",
    "exit_code",
)


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    # .../python/flexaidds/dataset_runner/tier1_contract.py → repo root
    return here.parents[3]


def md5_file(path: Path) -> str:
    """Identity hash for the energy matrix. Not a security construct."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_file(path: Path) -> Optional[str]:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def find_matrix_file(
    *,
    binary: Optional[os.PathLike] = None,
    repo_root: Optional[os.PathLike] = None,
) -> Optional[Path]:
    """Locate ``MC_st0r5.2_6.dat`` next to the binary or in the repo."""
    candidates: List[Path] = []
    if binary:
        bpath = Path(binary)
        candidates.append(bpath.parent / MATRIX_NAME)
    root = Path(repo_root) if repo_root is not None else _repo_root()
    candidates.extend(
        [
            root / MATRIX_NAME,
            root / "WRK" / MATRIX_NAME,
            root / ".grok" / "skills" / "flexaidds" / "data" / MATRIX_NAME,
        ]
    )
    seen = set()
    for cand in candidates:
        try:
            resolved = cand.resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_file():
            return resolved
    return None


def assert_matrix_pin(
    *,
    binary: Optional[os.PathLike] = None,
    repo_root: Optional[os.PathLike] = None,
    expected_md5: Optional[str] = None,
    dry_run: bool = False,
) -> str:
    """Return the observed matrix MD5. Fail closed on pin mismatch.

    ``dry_run`` skips the check (synthetic poses never load the matrix).
    """
    pin = (expected_md5 or MATRIX_PIN_MD5).strip().lower()
    if dry_run:
        return pin
    path = find_matrix_file(binary=binary, repo_root=repo_root)
    if path is None:
        raise RuntimeError(
            f"Tier-1 matrix pin: {MATRIX_NAME} not found next to binary "
            f"{binary!s} or under the repository. Refusing to dock without "
            f"the {pin} CF matrix."
        )
    observed = md5_file(path)
    if observed.lower() != pin:
        raise RuntimeError(
            f"Tier-1 matrix pin mismatch: {path} md5={observed} "
            f"expected {pin}. Mixing matrix forks is mixing physics. "
            "The CI gate does not retune the CF matrix."
        )
    return observed


def normalize_gate_class(value: Any, *, metric: str = "") -> str:
    klass = str(value or "hard").strip().lower() or "hard"
    if klass not in VALID_GATE_CLASSES:
        where = f" for {metric!r}" if metric else ""
        raise ValueError(
            f"ci_gate_class{where} must be 'hard' or 'advisory', got {value!r}"
        )
    return klass


def parse_ci_gate_class(raw: Any) -> Dict[str, str]:
    if not raw:
        return {}
    if not isinstance(raw, Mapping):
        raise ValueError("ci_gate_class must be a mapping of metric -> hard|advisory")
    return {
        str(key): normalize_gate_class(val, metric=str(key))
        for key, val in raw.items()
    }


def gate_class_for(config: Any, metric: str) -> str:
    classes = getattr(config, "ci_gate_class", None) or {}
    return normalize_gate_class(classes.get(metric, "hard"), metric=metric)


def split_regression_flags(
    config: Any,
    flags: Optional[Mapping[str, bool]],
) -> Tuple[Dict[str, bool], Dict[str, bool]]:
    """Split ``check_regressions`` flags into (blocking, advisory).

    Metrics with no ``ci_gate_class`` entry default to **hard** so existing
    datasets keep fail-closed quality gates. Advisory flags are reported
    but must not fail the process.
    """
    blocking: Dict[str, bool] = {}
    advisory: Dict[str, bool] = {}
    for metric, flagged in (flags or {}).items():
        if not flagged:
            continue
        if gate_class_for(config, metric) == "advisory":
            advisory[metric] = True
        else:
            blocking[metric] = True
    return blocking, advisory


def quality_exit_code(
    *,
    inconclusive: Sequence[str],
    blocking: Mapping[str, bool],
    dry_run: bool = False,
) -> int:
    """Map gate outcomes to the DatasetRunner CLI exit code.

    0 = pass (advisory ranking misses do not fail)
    1 = hard quality regression
    3 = inconclusive (#326 liveness / productivity / completeness)
    """
    if dry_run:
        return 0
    if inconclusive:
        return 3
    if any(blocking.values()):
        return 1
    return 0


def verdict_label(exit_code: int) -> str:
    return {0: "pass", 1: "fail", 3: "inconclusive"}.get(int(exit_code), "fail")


def normalize_roster(values: Iterable[Any]) -> List[str]:
    return [str(t).strip().lower() for t in values if str(t).strip()]


def assert_tier1_roster_lock(
    config: Any,
    scheduled_targets: Sequence[str],
    tier: int,
) -> None:
    """Refuse to dock a different draw than the YAML-locked roster."""
    if int(tier) != 1:
        return
    expected = normalize_roster(getattr(config, "tier1_expected_roster", None) or [])
    if not expected:
        return
    got = normalize_roster(scheduled_targets)
    if got == expected:
        return
    seed = "?"
    resolve = getattr(config, "resolve_tier1_seed", None)
    if callable(resolve):
        try:
            seed = str(resolve())
        except Exception:
            seed = "?"
    raise RuntimeError(
        f"Tier-1 roster replay failed (seed={seed}): expected {expected}, "
        f"got {got}. Update tier1_expected_roster together with tier1_seed "
        "and FLEXAIDDS_TIER1_SEED — never dock an unlocked draw against "
        "pinned sampling baselines."
    )


def assert_election_objective(config: Any, ranking_objective: str) -> None:
    pinned = str(getattr(config, "ci_election_objective", "") or "").strip()
    if not pinned:
        return
    if pinned != ranking_objective:
        raise RuntimeError(
            f"Tier-1 election objective lock: yaml ci_election_objective="
            f"{pinned!r} but runner ranking_objective={ranking_objective!r}. "
            "Refusing to elect under a different objective than the receipt."
        )


def assert_unregistered_ranking_is_advisory(config: Any) -> None:
    """Refuse to silently re-harden aspirational docking_power while unregistered.

    ``ci_ranking_status: unregistered`` is the contract that 0.70 / 0.85 are
    measured, not merge-blocking. A YAML edit that drops those metrics from
    ``ci_gate_class`` would default them to hard and recreate the CI miss.
    """
    status = str(getattr(config, "ci_ranking_status", "") or "").strip().lower()
    if status != "unregistered":
        return
    baselines = getattr(config, "expected_baselines", None) or {}
    hardened = [
        metric
        for metric in RANKING_QUALITY_METRICS
        if metric in baselines and gate_class_for(config, metric) != "advisory"
    ]
    if hardened:
        raise ValueError(
            "ci_ranking_status=unregistered but ranking metrics "
            f"{hardened} are not advisory. Mark them advisory in "
            "ci_gate_class, or set ci_ranking_status=registered after "
            "recording an election rule (docs/TIER1_CI_CONTRACT.md)."
        )


def _drawn_from_roster(roster: Sequence[str], anchors: Sequence[str]) -> List[str]:
    anchor_set = {a.lower() for a in anchors}
    return [t for t in roster if t.lower() not in anchor_set]


def build_receipt(
    *,
    seed: int,
    roster: Sequence[str],
    anchors: Sequence[str] = PINNED_ANCHORS,
    matrix_md5: str,
    matrix_path: Optional[str] = None,
    binary: str = "",
    binary_sha256: Optional[str] = None,
    election_objective: str = DEFAULT_ELECTION_OBJECTIVE,
    ci_ranking_status: str = "unregistered",
    git_sha: str = "",
    hard_gates: Optional[Mapping[str, Any]] = None,
    advisory_ranking: Optional[Mapping[str, Any]] = None,
    inconclusive: Optional[Sequence[str]] = None,
    blocking: Optional[Mapping[str, bool]] = None,
    exit_code: Optional[int] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    """Assemble a replayable Tier-1 receipt (JSON-serialisable)."""
    roster_n = normalize_roster(roster)
    anchors_n = normalize_roster(anchors)
    status = str(ci_ranking_status or "unregistered").strip().lower()
    if status not in VALID_RANKING_STATUS:
        raise ValueError(
            f"ci_ranking_status must be unregistered|registered, got {ci_ranking_status!r}"
        )
    blocking = dict(blocking or {})
    inconclusive = list(inconclusive or [])
    if exit_code is None:
        exit_code = quality_exit_code(inconclusive=inconclusive, blocking=blocking)
    payload: Dict[str, Any] = {
        "schema": RECEIPT_SCHEMA,
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "seed": int(seed),
        "seed_env": os.environ.get("FLEXAIDDS_TIER1_SEED"),
        "ga_seed": os.environ.get("FLEXAID_SEED") or None,
        "roster": roster_n,
        "anchors": anchors_n,
        "drawn": _drawn_from_roster(roster_n, anchors_n),
        "dataset": "astex_diverse",
        "tier": 1,
        "matrix_name": MATRIX_NAME,
        "matrix_path": matrix_path,
        "matrix_md5": str(matrix_md5).lower(),
        "matrix_pin_md5": MATRIX_PIN_MD5,
        "matrix_pin_ok": str(matrix_md5).lower() == MATRIX_PIN_MD5,
        "binary": binary,
        "binary_sha256": binary_sha256,
        "election_objective": election_objective,
        "ci_ranking_status": status,
        "git_sha": git_sha,
        "hard_gates": dict(hard_gates or {}),
        "advisory_ranking": dict(advisory_ranking or {}),
        "inconclusive": inconclusive,
        "blocking_regressions": [k for k, v in blocking.items() if v],
        "verdict": verdict_label(exit_code),
        "exit_code": int(exit_code),
        "thesis_election_direction": THESIS_ELECTION_DIRECTION,
        "notes": [
            "CF matrix was not retuned; CI pin is 9dc9 "
            f"({MATRIX_PIN_MD5}).",
            "docking_power_top1 0.70 is aspirational / advisory until a "
            "registered election rule is ready. Do not quote it as a "
            "receipted FlexAIDdS rate.",
            "A ranking miss with near-natives in the 10-pose library is "
            "sampling success + ranking collapse, not Softβ/entropy failure.",
            "Replay the draw: FLEXAIDDS_TIER1_SEED="
            f"{int(seed)} python -m flexaidds.dataset_runner.tier1_contract "
            "--replay-roster",
        ],
    }
    if extra:
        payload.update(dict(extra))
    return payload


def validate_receipt(
    data: Mapping[str, Any],
    *,
    require_matrix_pin: bool = True,
    require_roster: Sequence[str] = PINNED_ROSTER,
    require_seed: int = PINNED_TIER1_SEED,
) -> List[str]:
    """Return human-readable errors; empty means the receipt can be replayed."""
    errs: List[str] = []
    if not isinstance(data, Mapping):
        return ["receipt is not a JSON object"]
    if data.get("schema") != RECEIPT_SCHEMA:
        errs.append(
            f"schema {data.get('schema')!r} != {RECEIPT_SCHEMA!r}"
        )
    for key in RECEIPT_REQUIRED_KEYS:
        if key not in data:
            errs.append(f"missing key {key!r}")
    seed = data.get("seed")
    try:
        if int(seed) != int(require_seed):
            errs.append(f"seed {seed!r} != pinned {require_seed}")
    except (TypeError, ValueError):
        errs.append(f"seed {seed!r} is not an int")
    roster = normalize_roster(data.get("roster") or [])
    expected = normalize_roster(require_roster)
    if roster != expected:
        errs.append(f"roster {roster} != pinned {expected}")
    if require_matrix_pin:
        md5 = str(data.get("matrix_md5") or "").lower()
        if md5 != MATRIX_PIN_MD5:
            errs.append(f"matrix_md5 {md5!r} != pin {MATRIX_PIN_MD5}")
        if data.get("matrix_pin_ok") is not True:
            errs.append("matrix_pin_ok is not true")
    objective = data.get("election_objective")
    if objective != DEFAULT_ELECTION_OBJECTIVE:
        errs.append(
            f"election_objective {objective!r} != {DEFAULT_ELECTION_OBJECTIVE!r}"
        )
    return errs


def gate_snapshot(
    metrics: Mapping[str, Any],
    baselines: Mapping[str, Any],
    flags: Mapping[str, bool],
    names: Iterable[str],
) -> Dict[str, Any]:
    snap: Dict[str, Any] = {}
    for name in names:
        snap[name] = {
            "measured": metrics.get(name),
            "baseline": baselines.get(name),
            "regressed": bool(flags.get(name, False)),
        }
    return snap


def write_receipt(path: os.PathLike, payload: Mapping[str, Any]) -> Path:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n")
    return dest


def load_astex_config():
    """Load the package-copy Astex YAML (same draw algorithm as CI)."""
    from .runner import DatasetConfig

    pkg = Path(__file__).resolve().parent / "datasets" / "astex_diverse.yaml"
    canonical = _repo_root() / "benchmarks" / "datasets" / "astex_diverse.yaml"
    path = pkg if pkg.is_file() else canonical
    if not path.is_file():
        raise FileNotFoundError(f"astex_diverse.yaml not found at {pkg} or {canonical}")
    return DatasetConfig.from_yaml(path)


def replay_roster(seed: Optional[int] = None) -> List[str]:
    """Return the exact Tier-1 draw for ``seed`` (env wins inside DatasetConfig)."""
    if seed is not None:
        os.environ["FLEXAIDDS_TIER1_SEED"] = str(int(seed))
    cfg = load_astex_config()
    return list(cfg.tier1_targets())


def _main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-matrix",
        action="store_true",
        help="Fail unless the on-disk CF matrix matches the 9dc9 pin.",
    )
    parser.add_argument(
        "--binary",
        default=os.environ.get("FLEXAIDDS_BINARY", ""),
        help="FlexAID binary (matrix is resolved next to it, then repo root).",
    )
    parser.add_argument(
        "--replay-roster",
        action="store_true",
        help="Print the Tier-1 draw and fail if it is not the pinned roster.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=PINNED_TIER1_SEED,
        help=f"FLEXAIDDS_TIER1_SEED to replay (default {PINNED_TIER1_SEED}).",
    )
    parser.add_argument(
        "--validate-receipt",
        metavar="PATH",
        help="Validate a TIER1_RECEIPT.json against the pinned contract.",
    )
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.check_matrix:
        observed = assert_matrix_pin(binary=args.binary or None)
        print(f"matrix_pin_ok md5={observed} file pin={MATRIX_PIN_MD5}")

    if args.replay_roster:
        os.environ["FLEXAIDDS_TIER1_SEED"] = str(int(args.seed))
        got = replay_roster(args.seed)
        print(f"seed={int(args.seed)} roster={got}")
        if normalize_roster(got) != normalize_roster(PINNED_ROSTER):
            print(
                f"ERROR: draw {got} != pinned {list(PINNED_ROSTER)}",
                file=sys.stderr,
            )
            return 1
        cfg = load_astex_config()
        assert_tier1_roster_lock(cfg, got, tier=1)
        print("roster_lock_ok")

    if args.validate_receipt:
        path = Path(args.validate_receipt)
        data = json.loads(path.read_text(encoding="utf-8"))
        errs = validate_receipt(data)
        if errs:
            print("TIER1_RECEIPT invalid:", file=sys.stderr)
            for err in errs:
                print(f"  - {err}", file=sys.stderr)
            return 1
        print(f"receipt_ok path={path} verdict={data.get('verdict')}")

    if not (args.check_matrix or args.replay_roster or args.validate_receipt):
        parser.print_help()
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(_main())
