#!/usr/bin/env python3
"""10k bootstrap median of S_top10 success rates (3Dsig 2017 metric).

Per case: success if **any** of the top-10 ranked modes has RMSD ≤ 2.0 Å
(claim contract threshold; see docs/implementation/3dsig_red_pair_protocol.md).

Headline: **median** of B bootstrap resamples of the case set (with replacement).

Input formats (auto-detected):
  - Directory of per-PDB result.csv with **mode_rmsd_0..9** columns (required)
  - TSV/CSV with columns: pdb_id, rank, rmsd  (rank 0..9, emitted order)
  - JSON list of {pdb_id, rmsds: [..]} or {cases: {...}}

**Fail-closed:** arm-dir result.csv without mode_rmsd_* (or equivalent rank
columns) is rejected. Do **not** silently treat rmsd_top1 / rmsd_bcr as S_top10.

**Also fail-closed:** a case that never executed (``n_poses`` /
``restarts_finished`` present and zero, with no mode RMSDs) is a FAILED RUN,
not a miss. Scoring it as a miss fabricates a success rate — an arm that never
started would otherwise report a clean ``0.0000`` indistinguishable from a zero
the engine earned. Pass ``--allow-unrun-cases`` to drop them from N instead.

Usage:
  python3 scripts/bootstrap_3dsig_s_top10.py --cases cases.json --bootstraps 10000
  python3 scripts/bootstrap_3dsig_s_top10.py --arm-dir path/to/3dsig_r10
  python3 scripts/bootstrap_3dsig_s_top10.py --rank-table ranks.csv

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# Claim contract: success if RMSD ≤ 2.0 Å (not strict <)
DEFAULT_THRESH = 2.0
DEFAULT_BOOTSTRAPS = 10_000
TOP_N = 10

# Column aliases accepted for mode RMSDs in result.csv (prefer mode_rmsd_i)
MODE_RMSD_KEYS = (
    "mode_rmsd_{i}",
    "rmsd_{i}",
    "top{i}_rmsd",
    "rank{i}_rmsd",
)


class MissingModeRmsdError(ValueError):
    """Raised when result.csv lacks mode_rmsd_* (fail-closed for S_top10)."""


class UnrunCaseError(ValueError):
    """Raised when a case never executed (fail-closed for S_top10).

    A row with no poses and no finished restarts is a FAILED RUN, not a miss.
    Counting it as a miss silently deflates the success rate — a zero from an
    arm that never started is indistinguishable from a zero the engine earned.
    """


# Columns that witness whether the run actually executed. Absent columns are
# not evidence of anything, so an unrun verdict requires a present zero.
_EXECUTION_WITNESS_KEYS = ("n_poses", "restarts_finished")


def row_never_ran(row: dict, rmsds: Sequence[Optional[float]]) -> bool:
    """True when the row has no usable RMSD and a witness proves it never ran.

    Requires BOTH: every mode slot empty, and at least one execution witness
    column present and equal to zero. A row whose modes are merely empty is
    left alone — only a positive witness of non-execution disqualifies it.
    """
    if any(r is not None for r in rmsds):
        return False
    keys_lower = {k.lower(): k for k in row}
    for key in _EXECUTION_WITNESS_KEYS:
        col = keys_lower.get(key)
        if col is None:
            continue
        raw = (row.get(col) or "").strip()
        if not raw:
            continue
        try:
            if float(raw) == 0.0:
                return True
        except ValueError:
            continue
    return False


def _finite_rmsd(value: object) -> Optional[float]:
    if value is None or value == "" or value == "NA":
        return None
    try:
        v = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v) or v < 0.0:
        return None
    return v


def s_top10(
    rmsds: Sequence[Optional[float]], thresh: float = DEFAULT_THRESH
) -> bool:
    """True if any of the first 10 finite non-negative RMSDs is ≤ thresh."""
    for r in list(rmsds)[:TOP_N]:
        if r is None:
            continue
        if math.isfinite(r) and r >= 0.0 and r <= thresh:
            return True
    return False


def extract_mode_rmsds_from_row(
    row: dict, *, require_mode_columns: bool = True
) -> List[Optional[float]]:
    """Extract mode_rmsd_0..9 from a result.csv row.

    Prefer exact ``mode_rmsd_i`` keys. Accept documented aliases.
    If *require_mode_columns* and no mode/rank columns exist at all,
    raise MissingModeRmsdError (fail-closed — no BCR / top1 fallback).
    """
    keys_lower = {k.lower(): k for k in row}

    def _get_for_i(i: int) -> Tuple[bool, Optional[float]]:
        """Return (column_present, value)."""
        for template in MODE_RMSD_KEYS:
            key = template.format(i=i)
            if key in row:
                return True, _finite_rmsd(row[key])
            if key.lower() in keys_lower:
                return True, _finite_rmsd(row[keys_lower[key.lower()]])
        return False, None

    any_mode_col = False
    rmsds: List[Optional[float]] = []
    for i in range(TOP_N):
        present, val = _get_for_i(i)
        if present:
            any_mode_col = True
        rmsds.append(val)

    if not any_mode_col and require_mode_columns:
        # Explicit fail-closed: never treat rmsd_top1 / rmsd_bcr as S_top10
        sample_keys = sorted(row.keys())[:20]
        raise MissingModeRmsdError(
            "result.csv missing mode_rmsd_0..9 (or rmsd_i / top{i}_rmsd / "
            f"rank{i}_rmsd). Refusing BCR/top1 fallback. Columns sample: {sample_keys}"
        )
    return rmsds


def load_cases_json(path: Path) -> Dict[str, List[Optional[float]]]:
    data = json.loads(path.read_text())
    out: Dict[str, List[Optional[float]]] = {}
    if isinstance(data, dict) and "cases" in data:
        data = data["cases"]
    if isinstance(data, dict):
        for k, v in data.items():
            if isinstance(v, dict) and "rmsds" in v:
                out[str(k)] = [_finite_rmsd(x) for x in v["rmsds"]]
            elif isinstance(v, list):
                out[str(k)] = [_finite_rmsd(x) for x in v]
    elif isinstance(data, list):
        for row in data:
            pid = str(row.get("pdb_id") or row.get("pdb") or row["id"])
            out[pid] = [_finite_rmsd(x) for x in row["rmsds"]]
    return out


def load_rank_table(path: Path) -> Dict[str, List[Optional[float]]]:
    """pdb_id, rank, rmsd rows → top-10 rmsds sorted by **emitted rank** (not RMSD)."""
    by: Dict[str, List[Tuple[int, Optional[float]]]] = {}
    with path.open(newline="") as f:
        sample = f.read(2048)
        f.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        r = csv.DictReader(f, delimiter=delim)
        for row in r:
            keys = {k.lower(): k for k in row}
            pid = row[
                keys.get("pdb_id")
                or keys.get("pdb")
                or keys.get("case")
                or list(row)[0]
            ]
            rank = int(float(row[keys.get("rank") or list(row)[1]]))
            rmsd = _finite_rmsd(row[keys.get("rmsd") or list(row)[2]])
            by.setdefault(str(pid), []).append((rank, rmsd))
    out: Dict[str, List[Optional[float]]] = {}
    for pid, pairs in by.items():
        pairs.sort(key=lambda x: x[0])  # emitted rank order
        slots: List[Optional[float]] = [None] * TOP_N
        for rank, rmsd in pairs:
            if 0 <= rank < TOP_N:
                slots[rank] = rmsd
        # If ranks were 1-based sparse, also pack first 10 in rank order
        if all(v is None for v in slots) and pairs:
            packed = [rmsd for _, rmsd in pairs[:TOP_N]]
            slots = packed + [None] * (TOP_N - len(packed))
        out[pid] = slots
    return out


def load_arm_dir(
    arm_dir: Path,
    *,
    strict: bool = True,
    allow_unrun: bool = False,
) -> Dict[str, List[Optional[float]]]:
    """Scan result.csv under arm_dir for mode_rmsd_0..9.

    Fail-closed when *strict* (default): missing mode columns → error, not
    silent rmsd_top1 / rmsd_bcr substitution.

    Also fail-closed on cases that never executed (``n_poses``/
    ``restarts_finished`` present and zero, all mode slots empty). Pass
    *allow_unrun* to drop them with a warning instead — they are never
    scored as misses either way.
    """
    out: Dict[str, List[Optional[float]]] = {}
    unrun: List[str] = []
    unrun_detail: List[str] = []
    csv_paths = sorted(arm_dir.rglob("result.csv"))
    if not csv_paths:
        raise FileNotFoundError(f"no result.csv under {arm_dir}")

    errors: List[str] = []
    for csv_path in csv_paths:
        try:
            with csv_path.open(newline="") as f:
                rows = list(csv.DictReader(f))
        except OSError as e:
            errors.append(f"{csv_path}: {e}")
            continue
        if not rows:
            continue
        row = rows[0]
        pdb = (
            row.get("pdb_id") or row.get("receptor_id") or csv_path.parent.name
        ).upper()[:4]
        try:
            rmsds = extract_mode_rmsds_from_row(row, require_mode_columns=strict)
        except MissingModeRmsdError as e:
            errors.append(f"{csv_path}: {e}")
            continue
        if row_never_ran(row, rmsds):
            # NOT a miss: the run never executed. Never score it as failure.
            unrun.append(pdb)
            unrun_detail.append(f"{pdb} ({csv_path})")
            continue
        # Keep case even if all mode slots empty (counts as S_top10 fail)
        out[pdb] = rmsds

    if unrun and not allow_unrun:
        raise UnrunCaseError(
            f"S_top10 fail-closed: {len(unrun)} case(s) never executed "
            "(n_poses/restarts_finished = 0 with no mode RMSDs). Scoring them "
            "as misses would fabricate a success rate. Re-run those targets, "
            "or pass --allow-unrun-cases to exclude them from N. "
            f"First: {unrun_detail[0]}"
        )
    if unrun:
        print(
            f"warning: excluded {len(unrun)} never-executed case(s) from N: "
            + ", ".join(sorted(unrun)),
            file=sys.stderr,
        )

    if errors and strict:
        msg = (
            f"S_top10 fail-closed: {len(errors)} result.csv lack mode_rmsd_*; "
            "re-parse arms with scripts/parse_flexaid_arm_results.py. "
            f"First error: {errors[0]}"
        )
        raise MissingModeRmsdError(msg)
    if not out:
        raise MissingModeRmsdError(
            f"no usable mode_rmsd cases under {arm_dir} "
            f"({len(errors)} parse errors)"
        )
    return out


def bootstrap_median(
    case_success: Sequence[bool], n_boot: int, seed: int
) -> dict:
    """Bootstrap median success rate (and p05/p95) over cases with replacement."""
    rng = random.Random(seed)
    n = len(case_success)
    if n == 0:
        return {
            "n_cases": 0,
            "n_success": 0,
            "point": None,
            "median": None,
            "p05": None,
            "p95": None,
            "n_boot": n_boot,
        }
    successes = list(case_success)
    point = sum(1 for x in successes if x) / n
    samples: List[float] = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        samples.append(sum(1 for i in idx if successes[i]) / n)
    samples.sort()
    med = statistics.median(samples)
    # Inclusive empirical percentiles
    p05 = samples[max(0, int(math.floor(0.05 * (n_boot - 1))))]
    p95 = samples[min(n_boot - 1, int(math.ceil(0.95 * (n_boot - 1))))]
    return {
        "n_cases": n,
        "n_success": sum(1 for x in successes if x),
        "point": point,
        "median": med,
        "p05": p05,
        "p95": p95,
        "n_boot": n_boot,
    }


def compute_s_top10_stats(
    cases: Dict[str, List[Optional[float]]],
    *,
    thresh: float = DEFAULT_THRESH,
    n_boot: int = DEFAULT_BOOTSTRAPS,
    seed: int = 20170715,
) -> dict:
    """Unit-testable end-to-end: cases → S_top10 bootstrap stats."""
    ordered = sorted(cases.items(), key=lambda kv: kv[0])
    success = [s_top10(v, thresh) for _, v in ordered]
    stats = bootstrap_median(success, n_boot, seed)
    stats["metric"] = "S_top10"
    stats["thresh_A"] = thresh
    stats["thresh_op"] = "<="
    stats["n_pdbs_loaded"] = len(cases)
    stats["pdb_ids"] = [k for k, _ in ordered]
    stats["per_case_success"] = {
        k: bool(s_top10(v, thresh)) for k, v in ordered
    }
    return stats


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--cases", type=Path, help="JSON cases file")
    src.add_argument("--rank-table", type=Path, help="CSV/TSV pdb_id,rank,rmsd")
    src.add_argument(
        "--arm-dir",
        type=Path,
        help="Campaign arm directory with per-target result.csv (mode_rmsd_*)",
    )
    ap.add_argument(
        "--bootstraps",
        type=int,
        default=DEFAULT_BOOTSTRAPS,
        help=f"Bootstrap resamples (default {DEFAULT_BOOTSTRAPS})",
    )
    ap.add_argument("--seed", type=int, default=20170715)
    ap.add_argument(
        "--thresh",
        type=float,
        default=DEFAULT_THRESH,
        help=f"RMSD success threshold in Å (default {DEFAULT_THRESH}, inclusive ≤)",
    )
    ap.add_argument(
        "--allow-missing-modes",
        action="store_true",
        help="Do not fail when mode_rmsd_* missing (still never uses BCR as S_top10)",
    )
    ap.add_argument(
        "--allow-unrun-cases",
        action="store_true",
        help=(
            "Exclude never-executed cases (n_poses/restarts_finished = 0) from N "
            "with a warning, instead of failing. They are never scored as misses."
        ),
    )
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args(argv)

    try:
        if args.cases:
            cases = load_cases_json(args.cases)
        elif args.rank_table:
            cases = load_rank_table(args.rank_table)
        else:
            cases = load_arm_dir(
                args.arm_dir,
                strict=not args.allow_missing_modes,
                allow_unrun=args.allow_unrun_cases,
            )
    except (
        MissingModeRmsdError,
        UnrunCaseError,
        FileNotFoundError,
        OSError,
        json.JSONDecodeError,
    ) as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if not cases:
        print("error: no cases loaded", file=sys.stderr)
        return 2

    stats = compute_s_top10_stats(
        cases, thresh=args.thresh, n_boot=args.bootstraps, seed=args.seed
    )

    point = stats["point"]
    med = stats["median"]
    p05 = stats["p05"]
    p95 = stats["p95"]
    print(
        f"S_top10 point={point:.4f}  "
        f"bootstrap_median={med:.4f}  "
        f"p05–p95=[{p05:.4f},{p95:.4f}]  "
        f"n={stats['n_cases']} success={stats['n_success']}  "
        f"thresh=≤{args.thresh} Å  boot={args.bootstraps}"
    )
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        # per_case_success is useful; keep JSON serializable
        args.json_out.write_text(json.dumps(stats, indent=2) + "\n")
        print("wrote", args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
