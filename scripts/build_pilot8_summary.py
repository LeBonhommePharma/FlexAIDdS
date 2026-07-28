#!/usr/bin/env python3
"""Build a pilot8 results_summary JSON for evaluate_p3_pilot from real OUT trees.

Reads per-target result.csv under one or more arm OUT dirs and reports:
  n_targets, schema_ok (mode_rmsd_0..9 present), bcr_success, s_top10_success,
  per-arm and pooled counts. Never invents successes.

Usage:
  python3 scripts/build_pilot8_summary.py \\
    --arm-dir ~/flexaidds_results/campaigns/three_engine/A/comparative_pilot8_20260728 \\
    --arm-dir ~/flexaidds_results/campaigns/three_engine/B/comparative_pilot8_20260728 \\
    --json-out pilot_summary.json
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional


PILOT8 = ("1G9V", "1GPK", "1MEH", "1P62", "1Q4G", "1R9O", "1T40", "2BYS")
MODE_COLS = [f"mode_rmsd_{i}" for i in range(10)]


def _f(v: Any) -> Optional[float]:
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", "NA", "N/A", "nan", "None"):
        return None
    try:
        x = float(s)
    except (TypeError, ValueError):
        return None
    if x != x:
        return None
    return x


def _i(v: Any) -> Optional[int]:
    if v is None:
        return None
    s = str(v).strip()
    if s in ("", "NA", "N/A"):
        return None
    try:
        return int(float(s))
    except (TypeError, ValueError):
        return None


def read_result_csv(path: Path) -> Optional[Dict[str, str]]:
    if not path.is_file():
        return None
    try:
        with path.open(newline="", encoding="utf-8", errors="ignore") as fh:
            rows = list(csv.DictReader(fh))
    except OSError:
        return None
    if not rows:
        return None
    return {k: (v if v is not None else "") for k, v in rows[0].items()}


def schema_ok_row(row: Dict[str, str]) -> bool:
    present = 0
    for c in MODE_COLS:
        if c not in row:
            return False
        if _f(row.get(c)) is not None:
            present += 1
    # At least mode_rmsd_0 present and all columns exist
    return present >= 1 and all(c in row for c in MODE_COLS)


def bcr_success_row(row: Dict[str, str], thresh: float = 2.0) -> bool:
    for key in ("rmsd_bcr", "best_cluster_rmsd", "BCR"):
        v = _f(row.get(key))
        if v is not None:
            return v <= thresh
    return False


def s_top10_success_row(row: Dict[str, str], thresh: float = 2.0) -> bool:
    s = _i(row.get("success_s_top10"))
    if s is not None:
        return s == 1
    for c in MODE_COLS:
        v = _f(row.get(c))
        if v is not None and v <= thresh:
            return True
    return False


def scan_arm(arm_dir: Path, panel: List[str]) -> Dict[str, Any]:
    arm_dir = arm_dir.expanduser().resolve()
    targets: List[Dict[str, Any]] = []
    n_schema = 0
    n_bcr = 0
    n_s10 = 0
    for pdb in panel:
        rc = arm_dir / pdb / "result.csv"
        row = read_result_csv(rc)
        if row is None:
            targets.append({"pdb_id": pdb, "status": "missing"})
            continue
        ok = schema_ok_row(row)
        if ok:
            n_schema += 1
        bcr = bcr_success_row(row)
        s10 = s_top10_success_row(row)
        if bcr:
            n_bcr += 1
        if s10:
            n_s10 += 1
        targets.append(
            {
                "pdb_id": pdb,
                "status": "ok" if ok else "schema_incomplete",
                "schema_ok": ok,
                "bcr_success": bcr,
                "s_top10_success": s10,
                "rmsd_top1": _f(row.get("rmsd_top1")),
                "rmsd_bcr": _f(row.get("rmsd_bcr")),
                "matrix_md5": row.get("matrix_md5", ""),
                "seed_echo": _i(row.get("seed_echo")),
                "native_pose_seeded": _i(row.get("native_pose_seeded")),
                "result_csv": str(rc),
            }
        )
    n_done = sum(1 for t in targets if t.get("status") in ("ok", "schema_incomplete"))
    return {
        "arm_dir": str(arm_dir),
        "n_done": n_done,
        "n_schema_ok": n_schema,
        "bcr_success": n_bcr,
        "s_top10_success": n_s10,
        "targets": targets,
    }


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--arm-dir",
        action="append",
        type=Path,
        required=True,
        help="Arm OUT root containing <PDB>/result.csv (repeatable)",
    )
    ap.add_argument(
        "--panel",
        nargs="*",
        default=list(PILOT8),
        help="Target list (default pilot8)",
    )
    ap.add_argument("--json-out", type=Path, default=None)
    ap.add_argument(
        "--require-complete",
        action="store_true",
        help="Fail (exit 2) if any panel target missing on any arm",
    )
    args = ap.parse_args(argv)

    arms = [scan_arm(d, list(args.panel)) for d in args.arm_dir]
    # Pooled: a target counts once if any arm has it; success if any arm succeeds
    by_pdb: Dict[str, Dict[str, Any]] = {p: {"arms": []} for p in args.panel}
    for a in arms:
        for t in a["targets"]:
            pid = t["pdb_id"]
            by_pdb.setdefault(pid, {"arms": []})["arms"].append(t)

    n_targets = 0
    n_schema = 0
    n_bcr = 0
    n_s10 = 0
    for pid, blob in by_pdb.items():
        done = [t for t in blob["arms"] if t.get("status") in ("ok", "schema_incomplete")]
        if not done:
            continue
        n_targets += 1
        if any(t.get("schema_ok") for t in done):
            n_schema += 1
        if any(t.get("bcr_success") for t in done):
            n_bcr += 1
        if any(t.get("s_top10_success") for t in done):
            n_s10 += 1

    schema_ok = n_schema == n_targets and n_targets >= 1 and n_schema >= 1
    # Fail-closed: if any done row lacks schema, schema_ok false when incomplete schemas exist
    any_incomplete = any(
        t.get("status") == "schema_incomplete"
        for a in arms
        for t in a["targets"]
    )
    if any_incomplete:
        schema_ok = False

    missing = [
        pid
        for pid, blob in by_pdb.items()
        if not any(t.get("status") in ("ok", "schema_incomplete") for t in blob["arms"])
    ]

    # Incomplete pilot8 panel must not open P4 via a premature P3 pass:
    # evaluate_p3_pilot returns hold when science_hold is set.
    incomplete_panel = len(missing) > 0 or n_targets < len(args.panel)
    summary: Dict[str, Any] = {
        "schema": "pilot8_summary/v1",
        "panel": list(args.panel),
        "n_targets": n_targets,
        "n_panel": len(args.panel),
        "n_missing": len(missing),
        "missing_pdbs": missing,
        "schema_ok": bool(schema_ok),
        "mode_rmsd_present": bool(n_schema > 0),
        "bcr_success": n_bcr,
        "s_top10_success": n_s10,
        "S_top10_success": n_s10,
        "BCR_success": n_bcr,
        "arms": arms,
        "science_hold": bool(incomplete_panel),
        "SCIENCE_HOLD": bool(incomplete_panel),
        "full85_authorized": False,
        "note": (
            "Real OUT scan only; reconstruction A/B labels are outside this summary. "
            "n_targets counts panel members with ≥1 arm result.csv. "
            "science_hold=true while panel incomplete (blocks premature P3 pass / P4)."
        ),
    }

    text = json.dumps(summary, indent=2) + "\n"
    if args.json_out:
        args.json_out.expanduser().parent.mkdir(parents=True, exist_ok=True)
        args.json_out.expanduser().write_text(text, encoding="utf-8")
        print(f"wrote {args.json_out}", file=sys.stderr)
    print(text, end="")

    if args.require_complete and missing:
        return 2
    if n_targets < 1:
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main())
