#!/usr/bin/env python3
# Copyright 2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
# SPDX-License-Identifier: Apache-2.0
"""CLI for FlexAIDdS comparative pipeline serial phase gates (P2/P3/P4).

Fail-closed: P4 requires P2==pass AND P3==pass. Prints ALLOW or BLOCK.

Usage
-----
  # Inspect / gate check
  python3 scripts/comparative_phase_gate.py --state-file state.json --check-p4
  python3 scripts/comparative_phase_gate.py --state-file state.json --check-p3

  # Record oracle verdict
  python3 scripts/comparative_phase_gate.py --state-file state.json --set-p2 pass
  python3 scripts/comparative_phase_gate.py --state-file state.json --set-p3 hold

  # Evaluate from oracle / pilot JSON payloads
  python3 scripts/comparative_phase_gate.py --state-file state.json \\
      --eval-p2-json oracle_status.json
  python3 scripts/comparative_phase_gate.py --state-file state.json \\
      --eval-p3-json pilot_summary.json

  # Built-in synthetic self-test
  python3 scripts/comparative_phase_gate.py --dry-run

Exit codes
----------
  0  ALLOW / dry-run OK / successful write / eval pass
  1  BLOCK on --check-p3/--check-p4; dry-run failure; eval non-pass
  2  usage / I/O error
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

# Resolve package without requiring an editable install.
_REPO = Path(__file__).resolve().parents[1]
_PY = _REPO / "python"
if str(_PY) not in sys.path:
    sys.path.insert(0, str(_PY))

from flexaidds.comparative_phases.gates import (  # noqa: E402
    MATRIX_MD5_PIN,
    PHASES,
    can_run_p3,
    can_run_p4,
    empty_phase_state,
    evaluate_p2_oracle,
    evaluate_p3_pilot,
    load_phase_state,
    next_allowed_phase,
    save_phase_state,
    self_test,
    set_phase,
)


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"{path}: expected JSON object")
    # Panel oracle may wrap summary under "summary".
    if "summary" in data and isinstance(data["summary"], dict):
        # Prefer top-level ranking_allowed if present; else merge summary.
        merged = dict(data["summary"])
        for k, v in data.items():
            if k != "summary" and k not in merged:
                merged[k] = v
        return merged
    return data


def _print_state(state: Dict[str, Any]) -> None:
    print(f"matrix_md5_pin={state.get('matrix_md5_pin', MATRIX_MD5_PIN)}")
    for ph in PHASES:
        print(f"  {ph}={state.get(ph, 'pending')}")
    nxt = next_allowed_phase(state)
    print(f"next_allowed_phase={nxt}")
    print(f"can_run_p3={'yes' if can_run_p3(state) else 'no'}")
    print(f"can_run_p4={'yes' if can_run_p4(state) else 'no'}")


def run_dry_run() -> int:
    """Synthetic self-test with temporary state file round-trip."""
    errs = self_test()
    with tempfile.TemporaryDirectory(prefix="cmp_phase_gate_") as td:
        path = Path(td) / "phase_state.json"
        st = empty_phase_state()
        set_phase(st, "P0", "pass")
        set_phase(st, "P1", "pass")
        set_phase(st, "P2", "pass")
        set_phase(st, "P3", "pass")
        save_phase_state(path, st)
        loaded = load_phase_state(path)
        if not can_run_p4(loaded):
            errs.append("dry-run: can_run_p4 False after synthetic pass state")
        if next_allowed_phase(loaded) != "P4":
            errs.append(
                f"dry-run: next_allowed_phase expected P4, got {next_allowed_phase(loaded)}"
            )
        set_phase(loaded, "P4", "pass")
        if next_allowed_phase(loaded) != "P5":
            errs.append(
                f"dry-run: next after P4 pass expected P5, got {next_allowed_phase(loaded)}"
            )
        # Fail-closed: wipe P3 → BLOCK
        set_phase(loaded, "P3", "hold")
        save_phase_state(path, loaded)
        blocked = load_phase_state(path)
        if can_run_p4(blocked):
            errs.append("dry-run: can_run_p4 True with P3=hold (should BLOCK)")

    if errs:
        print("DRY-RUN FAIL", file=sys.stderr)
        for e in errs:
            print(f"  - {e}", file=sys.stderr)
        return 1
    print("DRY-RUN OK")
    print(f"MATRIX_MD5_PIN={MATRIX_MD5_PIN}")
    print(f"PHASES={','.join(PHASES)}")
    print("serial gates: P3 requires P2=pass; P4 requires P2=pass AND P3=pass")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Fail-closed serial gates for comparative pipeline phases P2/P3/P4. "
            "P4 prints ALLOW only when P2 and P3 are both pass."
        )
    )
    p.add_argument(
        "--state-file",
        type=Path,
        default=None,
        help="JSON phase state file (created if missing)",
    )
    p.add_argument(
        "--set-p2",
        choices=("pass", "hold", "fail", "pending"),
        default=None,
        help="Set P2 verdict and save state",
    )
    p.add_argument(
        "--set-p3",
        choices=("pass", "hold", "fail", "pending"),
        default=None,
        help="Set P3 verdict and save state",
    )
    p.add_argument(
        "--set-phase",
        nargs=2,
        metavar=("PHASE", "VERDICT"),
        default=None,
        help="Set arbitrary phase (e.g. P0 pass)",
    )
    p.add_argument(
        "--check-p3",
        action="store_true",
        help="Print ALLOW or BLOCK for P3 (requires P2=pass)",
    )
    p.add_argument(
        "--check-p4",
        action="store_true",
        help="Print ALLOW or BLOCK for P4 (requires P2=pass and P3=pass)",
    )
    p.add_argument(
        "--eval-p2-json",
        type=Path,
        default=None,
        help="Evaluate oracle status JSON; optionally write P2 into --state-file",
    )
    p.add_argument(
        "--eval-p3-json",
        type=Path,
        default=None,
        help="Evaluate pilot results_summary JSON; optionally write P3 into --state-file",
    )
    p.add_argument(
        "--apply",
        action="store_true",
        help="When used with --eval-p2-json / --eval-p3-json, write verdict into state",
    )
    p.add_argument(
        "--show",
        action="store_true",
        help="Print current phase state",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Run synthetic self-test (no external files required)",
    )
    return p


def main(argv: Optional[list] = None) -> int:
    args = build_parser().parse_args(argv)

    if args.dry_run:
        return run_dry_run()

    if args.state_file is None and not (args.eval_p2_json or args.eval_p3_json):
        build_parser().print_help()
        print(
            "\nerror: --state-file required unless --dry-run or pure --eval-*-json",
            file=sys.stderr,
        )
        return 2

    state: Dict[str, Any]
    if args.state_file is not None:
        try:
            state = load_phase_state(args.state_file)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"error loading state: {exc}", file=sys.stderr)
            return 2
    else:
        state = empty_phase_state()

    mutated = False

    if args.set_phase is not None:
        phase, verdict = args.set_phase
        try:
            set_phase(state, phase, verdict)
        except ValueError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        mutated = True

    if args.set_p2 is not None:
        set_phase(state, "P2", args.set_p2)
        mutated = True

    if args.set_p3 is not None:
        set_phase(state, "P3", args.set_p3)
        mutated = True

    eval_exit: Optional[int] = None

    if args.eval_p2_json is not None:
        try:
            payload = _load_json(args.eval_p2_json)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"error loading P2 json: {exc}", file=sys.stderr)
            return 2
        verdict, reason = evaluate_p2_oracle(payload)
        print(f"P2={verdict}  # {reason}")
        if args.apply or args.state_file is not None and args.set_p2 is None:
            # Apply only when --apply; avoid surprise writes on pure eval.
            if args.apply:
                set_phase(state, "P2", verdict)
                mutated = True
        eval_exit = 0 if verdict == "pass" else 1

    if args.eval_p3_json is not None:
        try:
            payload = _load_json(args.eval_p3_json)
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            print(f"error loading P3 json: {exc}", file=sys.stderr)
            return 2
        verdict, reason = evaluate_p3_pilot(payload)
        print(f"P3={verdict}  # {reason}")
        if args.apply:
            set_phase(state, "P3", verdict)
            mutated = True
        eval_exit = 0 if verdict == "pass" else 1

    if mutated:
        if args.state_file is None:
            print("error: cannot save without --state-file", file=sys.stderr)
            return 2
        try:
            save_phase_state(args.state_file, state)
        except OSError as exc:
            print(f"error saving state: {exc}", file=sys.stderr)
            return 2
        print(f"wrote {args.state_file}", file=sys.stderr)

    if args.show or (
        mutated
        and not args.check_p3
        and not args.check_p4
        and args.eval_p2_json is None
        and args.eval_p3_json is None
    ):
        _print_state(state)

    # Gate checks: print ALLOW or BLOCK; exit 0=ALLOW, 1=BLOCK (fail-closed).
    if args.check_p3:
        ok = can_run_p3(state)
        print("ALLOW" if ok else "BLOCK")
        return 0 if ok else 1

    if args.check_p4:
        ok = can_run_p4(state)
        print("ALLOW" if ok else "BLOCK")
        return 0 if ok else 1

    if eval_exit is not None and not mutated:
        return eval_exit

    return 0


if __name__ == "__main__":
    sys.exit(main())
