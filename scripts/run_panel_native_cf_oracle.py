#!/usr/bin/env python3
"""Panel-wide native CF oracle gate for a classic-arm campaign OUT tree.

Scans one-level ``*/result.csv`` (CloudDocs-safe), pairs each target with work
and/or OUT poses, runs the same logic as ``native_cf_oracle_gate.py``.

Claim policy (full85 / Softβ / arm C):
  ranking and Softβ experiments are **forbidden** until this gate PASSes on
  a canary set (or full panel). See ``docs/implementation/softbeta_election_policy.md``.

Usage:
  python3 scripts/run_panel_native_cf_oracle.py \\
    --out-dir ~/flexaidds_results/campaigns/three_engine/A/3dsig_full85_scratch_3b2fa57cc \\
    --work-root ~/flexaidds_results/three_engine_entropy_q1/work_scratch_3b2fa57cc/A \\
    --json-out oracle_panel.json --status-out oracle_status.json

Exit:
  0  all evaluated targets PASS (or --min-pass rate met)
  1  any pathology FAIL
  2  usage
  3  no evaluable targets / all missing native CF
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from native_cf_oracle_gate import (  # noqa: E402
    EXIT_FAIL_PATHOLOGY,
    EXIT_MISSING_NATIVE,
    EXIT_PASS,
    EXIT_USAGE,
    run_gate,
)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Campaign OUT with <PDB>/result.csv",
    )
    ap.add_argument(
        "--work-root",
        type=Path,
        default=None,
        help="Optional work tree root with <PDB>/ subdirs",
    )
    ap.add_argument("--tolerance", type=float, default=0.0)
    ap.add_argument(
        "--min-pass-rate",
        type=float,
        default=1.0,
        help="Fraction of evaluable targets that must PASS (default 1.0)",
    )
    ap.add_argument("--json-out", type=Path, default=None)
    ap.add_argument(
        "--status-out",
        type=Path,
        default=None,
        help="Compact status JSON for claim-mode launchers (oracle_status.json)",
    )
    ap.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Max targets to evaluate (0 = all)",
    )
    args = ap.parse_args(argv)

    out_dir = args.out_dir.expanduser().resolve()
    if not out_dir.is_dir():
        print(f"ERROR: out-dir missing: {out_dir}", file=sys.stderr)
        return EXIT_USAGE

    csvs = sorted(
        p
        for p in out_dir.glob("*/result.csv")
        if "incomplete" not in p.parent.name
    )
    if args.limit and args.limit > 0:
        csvs = csvs[: args.limit]
    if not csvs:
        print(f"ERROR: no result.csv under {out_dir}", file=sys.stderr)
        return EXIT_MISSING_NATIVE

    rows: List[Dict[str, Any]] = []
    n_pass = n_fail = n_missing = 0
    for rc in csvs:
        pdb = rc.parent.name
        work = None
        if args.work_root is not None:
            cand = args.work_root.expanduser().resolve() / pdb
            if cand.is_dir():
                work = cand
        res = run_gate(
            work=work,
            results=rc.parent,
            poses_dir=rc.parent,
            pdb_id=pdb,
            tolerance=args.tolerance,
            require_poses=True,
        )
        d = res.as_dict()
        d["pdb_id"] = pdb
        rows.append(d)
        if res.exit_code == EXIT_PASS:
            n_pass += 1
        elif res.exit_code == EXIT_FAIL_PATHOLOGY:
            n_fail += 1
        else:
            n_missing += 1
        tag = "PASS" if res.ok else ("PATHOL" if res.exit_code == 1 else "MISS")
        print(
            f"{pdb:6s} {tag:6s} native={res.cf_native} best_ga={res.best_ga_cf} "
            f"gap={res.gap}  {res.source_native} / {res.source_ga}"
        )

    n_eval = n_pass + n_fail
    rate = (n_pass / n_eval) if n_eval else 0.0
    ranking_ok = n_eval > 0 and rate + 1e-12 >= args.min_pass_rate and n_fail == 0
    # Allow min-pass-rate soft mode: if min_pass_rate < 1, allow some fails
    if args.min_pass_rate < 1.0 - 1e-12:
        ranking_ok = n_eval > 0 and rate + 1e-12 >= args.min_pass_rate

    summary = {
        "out_dir": str(out_dir),
        "work_root": str(args.work_root) if args.work_root else None,
        "n_targets": len(rows),
        "n_pass": n_pass,
        "n_fail_pathology": n_fail,
        "n_missing_native": n_missing,
        "pass_rate": round(rate, 4),
        "min_pass_rate": args.min_pass_rate,
        "tolerance": args.tolerance,
        "ranking_allowed": bool(ranking_ok),
        "softbeta_allowed": bool(ranking_ok),
        "arm_c_fo298_allowed": bool(ranking_ok),
        "claim_eligible": False,  # claim also needs R=10 + binary split
        "message": (
            "PASS: native CF competitive — Softβ/arm-C may proceed when other gates OK"
            if ranking_ok
            else "FAIL: native CF not competitive — Softβ/arm-C/ranking FORBIDDEN"
        ),
    }
    payload = {"summary": summary, "targets": rows}

    if args.json_out:
        args.json_out.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.json_out.expanduser().write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.json_out}", file=sys.stderr)
    if args.status_out:
        args.status_out.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.status_out.expanduser().write_text(
            json.dumps(summary, indent=2) + "\n", encoding="utf-8"
        )
        print(f"wrote {args.status_out}", file=sys.stderr)

    print(
        f"\n# panel: pass={n_pass} pathol={n_fail} missing={n_missing} "
        f"rate={rate:.2%} ranking_allowed={ranking_ok}"
    )
    if n_eval == 0:
        return EXIT_MISSING_NATIVE
    if ranking_ok:
        return EXIT_PASS
    return EXIT_FAIL_PATHOLOGY


if __name__ == "__main__":
    sys.exit(main())
