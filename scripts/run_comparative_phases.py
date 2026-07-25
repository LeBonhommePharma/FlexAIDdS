#!/usr/bin/env python3
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""CLI: comparative campaign phases P0–P5 (goal methodology).

Usage:
  python3 scripts/run_comparative_phases.py --phase P0
  python3 scripts/run_comparative_phases.py --pipeline-dry
  python3 scripts/run_comparative_phases.py --pipeline-dry --force-p2-pass --force-p3-pass
  python3 scripts/run_comparative_phases.py --phase P4 --check-only

See docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

# Allow running without install
_REPO = Path(__file__).resolve().parents[1]
_PY = _REPO / "python"
if str(_PY) not in sys.path:
    sys.path.insert(0, str(_PY))

from flexaidds.comparative_phases.gates import (  # noqa: E402
    can_run_p4,
    load_phase_state,
)
from flexaidds.comparative_phases.pipeline import (  # noqa: E402
    default_state_path,
    run_phase,
    run_pipeline_dry,
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Comparative phases P0–P5")
    ap.add_argument("--phase", choices=["P0", "P1", "P2", "P3", "P4", "P5"])
    ap.add_argument("--pipeline-dry", action="store_true", help="Run dry P0–P5 gate path")
    ap.add_argument("--local-root", default=None)
    ap.add_argument("--campaign", default="comparative_pilot8")
    ap.add_argument("--allow-reconstruction", action="store_true")
    ap.add_argument("--force-p2-pass", action="store_true")
    ap.add_argument("--force-p3-pass", action="store_true")
    ap.add_argument("--check-only", action="store_true", help="For P4: print ALLOW/BLOCK only")
    ap.add_argument("--oracle-json", default=None, help="Path to oracle_status.json for P2")
    ap.add_argument("--pilot-json", default=None, help="Path to pilot summary for P3")
    ap.add_argument("--json-out", default=None)
    args = ap.parse_args(argv)

    if args.pipeline_dry:
        out = run_pipeline_dry(
            args.local_root,
            allow_reconstruction=True,
            campaign=args.campaign,
            force_p2_pass=args.force_p2_pass,
            force_p3_pass=args.force_p3_pass,
        )
        print(json.dumps(out, indent=2, default=str))
        if args.json_out:
            Path(args.json_out).write_text(json.dumps(out, indent=2, default=str) + "\n")
        # Always exit 0 for dry pipeline that completed; inspect steps for holds
        print("PHASE=PIPELINE_DRY status=done")
        return 0

    if not args.phase:
        ap.error("need --phase or --pipeline-dry")

    if args.check_only and args.phase == "P4":
        sp = default_state_path(args.local_root)
        st = load_phase_state(sp)
        allow = can_run_p4(st)
        print("ALLOW" if allow else "BLOCK")
        print(f"PHASE=P4 check can_run_p4={allow} state={ {k: st.get(k) for k in ['P2','P3']} }")
        return 0 if allow else 3

    oracle = None
    if args.oracle_json:
        oracle = json.loads(Path(args.oracle_json).read_text())
    pilot = None
    if args.pilot_json:
        pilot = json.loads(Path(args.pilot_json).read_text())

    result = run_phase(
        args.phase,
        local_root_path=args.local_root,
        allow_reconstruction=args.allow_reconstruction,
        campaign=args.campaign,
        oracle_status=oracle,
        pilot_summary=pilot,
        dry_run=True,
    )
    print(json.dumps(result, indent=2, default=str))
    print(f"PHASE={result.get('phase')} status={result.get('status')}")
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(result, indent=2, default=str) + "\n")
    if result.get("status") == "fail" and result.get("blocked"):
        return 3
    if result.get("status") == "fail":
        return 2
    if result.get("status") == "hold":
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
