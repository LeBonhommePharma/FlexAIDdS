# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
"""Integrated P0–P5 comparative campaign pipeline (fail-closed serial docks)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

from .gates import (
    can_run_p3,
    can_run_p4,
    evaluate_p2_oracle,
    evaluate_p3_pilot,
    load_phase_state,
    next_allowed_phase,
    save_phase_state,
)
from .p0_layout import local_root, run_p0
from .p1_binaries import run_p1
from .p5_analyze import run_p5


def default_state_path(local: Optional[str] = None) -> Path:
    return local_root(local) / "campaigns" / "three_engine" / "phase_state.json"


def run_phase(
    phase: str,
    *,
    local_root_path: Optional[str] = None,
    allow_reconstruction: bool = False,
    campaign: str = "comparative_pilot8",
    oracle_status: Optional[Dict[str, Any]] = None,
    pilot_summary: Optional[Dict[str, Any]] = None,
    dry_run: bool = True,
    state_path: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run one phase and update state file. Docking phases respect serial gates."""
    sp = Path(state_path) if state_path else default_state_path(local_root_path)
    state = load_phase_state(sp)
    phase = phase.upper()
    result: Dict[str, Any]

    if phase == "P0":
        result = run_p0(local_root_path)
    elif phase == "P1":
        if str(state.get("P0")) != "pass" and not dry_run:
            return {
                "phase": "P1",
                "status": "fail",
                "reason": "P0 not pass (run P0 first)",
                "blocked": True,
            }
        result = run_p1(
            local_root_path, allow_reconstruction=allow_reconstruction
        )
    elif phase == "P2":
        if str(state.get("P1")) != "pass":
            return {
                "phase": "P2",
                "status": "fail",
                "reason": "P1 not pass",
                "blocked": True,
            }
        osrc = oracle_status or {}
        # Load from disk if empty
        if not osrc:
            cand = (
                local_root(local_root_path)
                / "campaigns"
                / "three_engine"
                / f"{campaign}_oracle_status.json"
            )
            if cand.is_file():
                osrc = json.loads(cand.read_text(encoding="utf-8"))
            elif dry_run:
                # Dry-run self-test default: deferred hold (not silent pass)
                osrc = {"deferred": True, "status": "DEFERRED"}
        st, reason = evaluate_p2_oracle(osrc)
        result = {
            "phase": "P2",
            "status": st,
            "reason": reason,
            "oracle_status": osrc,
        }
    elif phase == "P3":
        if not can_run_p3(state):
            return {
                "phase": "P3",
                "status": "fail",
                "reason": "P4-serial rule: P3 blocked until P2=pass",
                "blocked": True,
                "can_run_p3": False,
            }
        summary = pilot_summary or {}
        if not summary and dry_run:
            # Dry-run: schema-ok pilot with n=8 and nonzero BCR is a pass example;
            # real runs must supply summary or result scan.
            summary = {
                "n_targets": 8,
                "schema_ok": True,
                "bcr_success": 0,
                "s_top10_success": 0,
            }
        st, reason = evaluate_p3_pilot(summary)
        result = {
            "phase": "P3",
            "status": st,
            "reason": reason,
            "pilot_summary": summary,
            "launcher": "scripts/run_3dsig_red_pair_serial.sh --only A|B",
            "dry_run": dry_run,
        }
    elif phase == "P4":
        if not can_run_p4(state):
            return {
                "phase": "P4",
                "status": "fail",
                "reason": "P4 blocked: requires P2=pass and P3=pass",
                "blocked": True,
                "can_run_p4": False,
            }
        result = {
            "phase": "P4",
            "status": "pass" if dry_run else "pending",
            "reason": (
                "dry-run: P4 gate open (would call run_3dsig_red_pair_full85.sh)"
                if dry_run
                else "invoke scripts/run_3dsig_red_pair_full85.sh serially"
            ),
            "launcher": "scripts/run_3dsig_red_pair_full85.sh",
            "dry_run": dry_run,
        }
    elif phase == "P5":
        result = run_p5(
            campaign,
            local_root_path=local_root_path,
            dry_run=dry_run,
        )
    else:
        return {"phase": phase, "status": "fail", "reason": f"unknown phase {phase}"}

    # Update state
    state[phase] = result.get("status", "fail")
    state[f"{phase}_reason"] = result.get("reason", "")
    state["last_phase"] = phase
    save_phase_state(sp, state)
    result["state_path"] = str(sp)
    result["phase_state"] = {k: state.get(k) for k in ("P0", "P1", "P2", "P3", "P4", "P5")}
    result["next_allowed"] = next_allowed_phase(state)
    return result


def run_pipeline_dry(
    local_root_path: Optional[str] = None,
    *,
    allow_reconstruction: bool = True,
    campaign: str = "comparative_pilot8_dry",
    force_p2_pass: bool = False,
    force_p3_pass: bool = False,
) -> Dict[str, Any]:
    """Run P0→P5 dry path enforcing gates; optional synthetic P2/P3 pass for P4 check."""
    steps: List[Dict[str, Any]] = []
    sp = default_state_path(local_root_path)
    if sp.is_file():
        sp.unlink()

    steps.append(
        run_phase("P0", local_root_path=local_root_path, state_path=sp)
    )
    steps.append(
        run_phase(
            "P1",
            local_root_path=local_root_path,
            allow_reconstruction=allow_reconstruction,
            state_path=sp,
        )
    )
    # P2: default deferred → hold; optionally inject pass
    oracle = (
        {"ranking_allowed": True, "status": "PASS"}
        if force_p2_pass
        else {"deferred": True, "status": "DEFERRED"}
    )
    steps.append(
        run_phase(
            "P2",
            local_root_path=local_root_path,
            oracle_status=oracle,
            campaign=campaign,
            state_path=sp,
        )
    )
    # Attempt P3 without P2 pass should block when P2 hold
    steps.append(
        run_phase(
            "P3",
            local_root_path=local_root_path,
            dry_run=True,
            pilot_summary=(
                {
                    "n_targets": 8,
                    "schema_ok": True,
                    "bcr_success": 1,
                    "s_top10_success": 1,
                }
                if force_p3_pass
                else None
            ),
            state_path=sp,
        )
    )
    steps.append(
        run_phase(
            "P4",
            local_root_path=local_root_path,
            dry_run=True,
            state_path=sp,
        )
    )
    # If we forced passes, re-set and open P4
    if force_p2_pass and force_p3_pass:
        st = load_phase_state(sp)
        st["P2"] = "pass"
        st["P3"] = "pass"
        save_phase_state(sp, st)
        steps.append(
            run_phase(
                "P4",
                local_root_path=local_root_path,
                dry_run=True,
                state_path=sp,
            )
        )
    steps.append(
        run_phase(
            "P5",
            local_root_path=local_root_path,
            campaign=campaign,
            dry_run=True,
            state_path=sp,
        )
    )
    return {
        "steps": steps,
        "state": load_phase_state(sp),
        "state_path": str(sp),
    }
