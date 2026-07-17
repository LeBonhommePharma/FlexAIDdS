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
        for template in MODE_RMSD_KEYS:
            key = template.format(i=i)
            if key in row:
                return True, _finite_rmsd(row[key])
            if key.lower() in keys_lower:
                return True, _finite_rmsd(row[keys_lower[key.lower()]])
        return False, None

    any_mode_col = False
    out: List[Optional[float]] = []
    for i in range(TOP_N):
        present, val = _get_for_i(i)
        if present:
            any_mode_col = True
        out.append(val)

    if require_mode_columns and not any_mode_col:
        # Explicitly reject BCR / top1 as S_top10 stand-ins
        bcr_keys = ("rmsd_bcr", "bcr", "rmsd_top1", "top1_rmsd")
        has_bcr = any(k in row or k.lower() in keys_lower for k in bcr_keys)
        raise MissingModeRmsdError(
            "result.csv lacks mode_rmsd_0..9 (or aliases); refuse S_top10 "
            f"(bcr/top1 present={has_bcr} — must not be used as S_top10)"
        )
    return out


def load_arm_dir(arm_dir: Path, thresh: float = DEFAULT_THRESH) -> Dict[str, bool]:
    """Load per-pdb success flags from arm_dir/*/result.csv."""
    arm_dir = Path(arm_dir)
    if not arm_dir.is_dir():
        raise FileNotFoundError(f"arm-dir not found: {arm_dir}")
    cases: Dict[str, bool] = {}
    csvs = sorted(arm_dir.glob("*/result.csv"))
    if not csvs:
        raise FileNotFoundError(f"no result.csv under {arm_dir}")
    for p in csvs:
        pdb_id = p.parent.name
        with p.open(newline="") as f:
            rows = list(csv.DictReader(f))
        if not rows:
            raise ValueError(f"empty result.csv: {p}")
        row = rows[0]
        modes = extract_mode_rmsds_from_row(row, require_mode_columns=True)
        # Prefer explicit success_s_top10 if present and consistent; always recompute from modes
        ok = s_top10(modes, thresh=thresh)
        cases[pdb_id] = ok
    return cases


def load_rank_table(path: Path, thresh: float = DEFAULT_THRESH) -> Dict[str, bool]:
    """Load from table with pdb_id, rank, rmsd columns."""
    by_pdb: Dict[str, List[Optional[float]]] = {}
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f, delimiter="\t" if path.suffix.lower() == ".tsv" else ",")
        if not reader.fieldnames:
            raise ValueError("rank table has no header")
        fields = {h.lower(): h for h in reader.fieldnames}
        for req in ("pdb_id", "rank", "rmsd"):
            if req not in fields:
                raise ValueError(f"rank table missing column {req}")
        for row in reader:
            pdb = row[fields["pdb_id"]].strip()
            try:
                rank = int(float(row[fields["rank"]]))
            except (TypeError, ValueError):
                continue
            if rank < 0 or rank >= TOP_N:
                continue
            rms = _finite_rmsd(row[fields["rmsd"]])
            arr = by_pdb.setdefault(pdb, [None] * TOP_N)
            arr[rank] = rms
    if not by_pdb:
        raise ValueError("rank table produced zero cases")
    return {pdb: s_top10(rmsds, thresh=thresh) for pdb, rmsds in by_pdb.items()}


def load_cases_json(path: Path, thresh: float = DEFAULT_THRESH) -> Dict[str, bool]:
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict) and "cases" in data:
        data = data["cases"]
    cases: Dict[str, bool] = {}
    if isinstance(data, list):
        for item in data:
            pdb = str(item.get("pdb_id") or item.get("pdb") or item.get("id"))
            rmsds = item.get("rmsds") or item.get("mode_rmsds")
            if rmsds is None:
                raise ValueError(f"case {pdb} missing rmsds")
            cases[pdb] = s_top10(
                [_finite_rmsd(x) for x in rmsds], thresh=thresh
            )
    elif isinstance(data, dict):
        for pdb, rmsds in data.items():
            if isinstance(rmsds, dict) and "rmsds" in rmsds:
                rmsds = rmsds["rmsds"]
            cases[str(pdb)] = s_top10(
                [_finite_rmsd(x) for x in rmsds], thresh=thresh
            )
    else:
        raise ValueError("unsupported cases JSON shape")
    if not cases:
        raise ValueError("cases JSON empty")
    return cases


def bootstrap_median_rate(
    successes: Sequence[bool],
    bootstraps: int = DEFAULT_BOOTSTRAPS,
    seed: int = 0,
) -> Tuple[float, float, List[float]]:
    """Return (observed_rate, bootstrap_median, all bootstrap rates)."""
    n = len(successes)
    if n == 0:
        raise ValueError("empty success vector")
    obs = sum(1 for s in successes if s) / n
    if bootstraps <= 0 or n < 2:
        return obs, obs, [obs]
    rng = random.Random(seed)
    rates: List[float] = []
    idx = list(range(n))
    for _ in range(bootstraps):
        sample = [successes[rng.choice(idx)] for _ in range(n)]
        rates.append(sum(1 for s in sample if s) / n)
    return obs, statistics.median(rates), rates


def summarize(
    cases: Dict[str, bool],
    *,
    bootstraps: int,
    seed: int,
    thresh: float,
    source: str,
) -> dict:
    ids = sorted(cases.keys())
    succ = [cases[i] for i in ids]
    n = len(succ)
    hits = sum(1 for s in succ if s)
    obs, med, rates = bootstrap_median_rate(succ, bootstraps=bootstraps, seed=seed)
    rates_sorted = sorted(rates)
    def pct(p: float) -> Optional[float]:
        if not rates_sorted:
            return None
        k = min(len(rates_sorted) - 1, max(0, int(round(p * (len(rates_sorted) - 1)))))
        return rates_sorted[k]
    return {
        "metric": "S_top10",
        "threshold_A": thresh,
        "definition": "any of ranks 0..9 RMSD <= threshold",
        "source": source,
        "n_cases": n,
        "n_hits": hits,
        "observed_rate": obs,
        "bootstrap": {
            "n_resamples": bootstraps if n >= 2 else 0,
            "seed": seed,
            "median": med,
            "p05": pct(0.05),
            "p95": pct(0.95),
        },
        "fail_closed": True,
        "bcr_not_used_as_s_top10": True,
        "case_ids": ids,
        "case_success": {i: cases[i] for i in ids},
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--arm-dir", type=Path, help="Directory of pdb_id/result.csv")
    src.add_argument("--cases", type=Path, help="JSON cases file")
    src.add_argument("--rank-table", type=Path, help="CSV/TSV pdb_id,rank,rmsd")
    ap.add_argument("--bootstraps", type=int, default=DEFAULT_BOOTSTRAPS)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--thresh", type=float, default=DEFAULT_THRESH)
    ap.add_argument("--json-out", type=Path, default=None)
    args = ap.parse_args(argv)

    try:
        if args.arm_dir is not None:
            cases = load_arm_dir(args.arm_dir, thresh=args.thresh)
            source = f"arm-dir:{args.arm_dir}"
        elif args.cases is not None:
            cases = load_cases_json(args.cases, thresh=args.thresh)
            source = f"cases:{args.cases}"
        else:
            cases = load_rank_table(args.rank_table, thresh=args.thresh)
            source = f"rank-table:{args.rank_table}"
    except MissingModeRmsdError as e:
        print(f"FAIL closed: {e}", file=sys.stderr)
        return 2
    except (FileNotFoundError, ValueError) as e:
        print(f"FAIL: {e}", file=sys.stderr)
        return 1

    summary = summarize(
        cases,
        bootstraps=args.bootstraps,
        seed=args.seed,
        thresh=args.thresh,
        source=source,
    )
    # Human-readable headline
    print(
        f"S_top10 N={summary['n_cases']} hits={summary['n_hits']} "
        f"observed={summary['observed_rate']:.4f} "
        f"bootstrap_median={summary['bootstrap']['median']:.4f} "
        f"(B={summary['bootstrap']['n_resamples']}, seed={args.seed}, thresh={args.thresh}A)"
    )
    print(
        f"  p05={summary['bootstrap']['p05']} p95={summary['bootstrap']['p95']} "
        f"source={source}"
    )
    print("  fail_closed=True  (mode_rmsd_* required; BCR not used as S_top10)")
    if args.json_out:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        # drop bulky case map optional — keep for audit
        args.json_out.write_text(json.dumps(summary, indent=2) + "\n")
        print(f"  wrote {args.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
