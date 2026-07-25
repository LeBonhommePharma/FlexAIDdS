# Copyright 2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
# SPDX-License-Identifier: Apache-2.0
"""Fail-closed serial gate logic for comparative pipeline phases P2/P3/P4.

Pure helpers (no I/O except JSON state load/save). Unit-testable without a
docking binary or live results tree.

Science order (COMPARATIVE_GOAL_METHODOLOGY.md §5)::

    P0 pins → P1 build → P2 oracle → P3 pilot8 → P4 full85 → P5 bootstrap

Fail-closed rules
-----------------
* P3 may run only when P2 == ``pass``.
* P4 may run only when P2 == ``pass`` **and** P3 == ``pass``.
* Any missing / unknown / hold state blocks later phases.
* SCIENCE_HOLD and dual-arm zero-success (BCR=0 and S_top10=0, N≥8) are
  **hold** (not pass) — do not advance to full85 claim tables.

Tests outline (if not yet in ``python/tests/``)
----------------------------------------------
1. ``evaluate_p2_oracle``: ranking_allowed True → pass; False → hold;
   status PASS / HOLD / SCIENCE_HOLD; native_cf_oracle PASS; empty {} +
   deferred → not hard-fail; missing keys without deferred → fail.
2. ``evaluate_p3_pilot``: schema_ok False → fail; science_hold → hold;
   n_targets≥8 and bcr_success=0 and s_top10_success=0 → hold;
   n_targets≥1 schema_ok True → pass.
3. ``can_run_p3`` / ``can_run_p4``: only exact ``"pass"`` opens the gate.
4. ``next_allowed_phase``: walks PHASES; first non-pass is next.
5. ``load_phase_state`` / ``save_phase_state``: round-trip JSON; missing file
   yields empty state (fail-closed defaults).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, List, Mapping, MutableMapping, Optional, Tuple, Union

# Production energy-matrix pin (Phase 0). Must match MC_st0r5.2_6.dat under
# local-first layout. See COMPARATIVE_GOAL_METHODOLOGY.md and
# benchmarks/protocols/three_engine_entropy_comparison.md.
MATRIX_MD5_PIN = "9dc93717dfed0698006d88dd6a9627bc"

# Serial phase order. Later phases require earlier ones to be ``pass``.
PHASES: Tuple[str, ...] = ("P0", "P1", "P2", "P3", "P4", "P5")

# Verdict strings (stable contract for state files + CLI).
GateVerdict = str  # "pass" | "hold" | "fail" | "pending"

_VALID_VERDICTS = frozenset({"pass", "hold", "fail", "pending"})

PathLike = Union[str, os.PathLike]


def empty_phase_state() -> Dict[str, Any]:
    """Return a fail-closed default state: every phase ``pending``."""
    return {
        "matrix_md5_pin": MATRIX_MD5_PIN,
        "phases": {p: "pending" for p in PHASES},
        # Convenience top-level mirrors (CLI / scripts often set these).
        **{p: "pending" for p in PHASES},
    }


def _normalize_verdict(value: Any) -> Optional[str]:
    if value is None:
        return None
    s = str(value).strip().lower()
    if s in _VALID_VERDICTS:
        return s
    # Accept common oracle spellings.
    if s in ("ok", "true", "yes", "allowed"):
        return "pass"
    if s in ("blocked", "block", "science_hold"):
        return "hold"
    if s in ("error", "false", "no", "forbidden"):
        return "fail"
    return None


def _phase_value(phase_state: Mapping[str, Any], phase: str) -> str:
    """Read phase status from top-level key or nested ``phases`` map."""
    if phase in phase_state:
        v = _normalize_verdict(phase_state[phase])
        if v is not None:
            return v
    nested = phase_state.get("phases")
    if isinstance(nested, Mapping) and phase in nested:
        v = _normalize_verdict(nested[phase])
        if v is not None:
            return v
    return "pending"


def load_phase_state(path: PathLike) -> Dict[str, Any]:
    """Load JSON phase state from *path*.

    Missing file → empty fail-closed state (all ``pending``).
    Invalid JSON → raises ``json.JSONDecodeError`` / ``OSError``.
    """
    p = Path(path)
    if not p.is_file():
        return empty_phase_state()
    with p.open("r", encoding="utf-8") as fh:
        raw = json.load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"phase state must be a JSON object, got {type(raw).__name__}")
    # Ensure pin + phase keys exist without wiping caller extras.
    state = empty_phase_state()
    state.update(raw)
    # Re-sync nested phases map from top-level if present.
    phases_map: Dict[str, str] = dict(state.get("phases") or {})
    for ph in PHASES:
        if ph in raw:
            nv = _normalize_verdict(raw[ph])
            if nv is not None:
                phases_map[ph] = nv
                state[ph] = nv
        elif ph not in phases_map:
            phases_map[ph] = "pending"
            state[ph] = "pending"
        else:
            nv = _normalize_verdict(phases_map[ph])
            state[ph] = nv if nv is not None else "pending"
            phases_map[ph] = state[ph]
    state["phases"] = phases_map
    if "matrix_md5_pin" not in state:
        state["matrix_md5_pin"] = MATRIX_MD5_PIN
    return state


def save_phase_state(path: PathLike, state: Mapping[str, Any]) -> None:
    """Atomically write *state* as pretty JSON to *path*."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    # Normalize known phase keys before write.
    out: Dict[str, Any] = dict(state)
    phases_map: Dict[str, str] = {}
    nested = out.get("phases")
    if isinstance(nested, Mapping):
        for k, v in nested.items():
            nv = _normalize_verdict(v)
            if nv is not None:
                phases_map[str(k)] = nv
    for ph in PHASES:
        if ph in out:
            nv = _normalize_verdict(out[ph])
            if nv is not None:
                phases_map[ph] = nv
                out[ph] = nv
        elif ph in phases_map:
            out[ph] = phases_map[ph]
        else:
            phases_map[ph] = "pending"
            out[ph] = "pending"
    out["phases"] = phases_map
    if "matrix_md5_pin" not in out:
        out["matrix_md5_pin"] = MATRIX_MD5_PIN

    payload = json.dumps(out, indent=2, sort_keys=True) + "\n"
    tmp = p.with_suffix(p.suffix + ".tmp")
    tmp.write_text(payload, encoding="utf-8")
    os.replace(tmp, p)


def _truthy_hold_flag(oracle_status: Mapping[str, Any]) -> bool:
    """True when SCIENCE_HOLD / science_hold is asserted."""
    if oracle_status.get("SCIENCE_HOLD") is True:
        return True
    if oracle_status.get("science_hold") is True:
        return True
    status = oracle_status.get("status")
    if isinstance(status, str) and status.strip().upper() in {
        "SCIENCE_HOLD",
        "HOLD",
    }:
        # HOLD handled by caller for status==HOLD; SCIENCE_HOLD always hold.
        if status.strip().upper() == "SCIENCE_HOLD":
            return True
    msg = oracle_status.get("message")
    if isinstance(msg, str) and "SCIENCE_HOLD" in msg.upper():
        return True
    return False


def _is_deferred_empty(oracle_status: Mapping[str, Any]) -> bool:
    """Empty or explicitly deferred status: do not hard-fail for missing keys."""
    if not oracle_status:
        return True
    if oracle_status.get("deferred") is True:
        return True
    if str(oracle_status.get("status", "")).strip().upper() == "DEFERRED":
        return True
    # Only bookkeeping keys, no evaluation signal.
    signal_keys = {
        "ranking_allowed",
        "status",
        "native_cf_oracle",
        "SCIENCE_HOLD",
        "science_hold",
        "n_pass",
        "n_fail_pathology",
        "pass_rate",
    }
    present_signals = [k for k in signal_keys if k in oracle_status]
    if not present_signals and oracle_status.get("deferred") is not False:
        # Treat as deferred/empty when no evaluation signal present.
        return True
    return False


def _has_p2_required_keys(oracle_status: Mapping[str, Any]) -> bool:
    """At least one primary oracle signal key must be present."""
    return any(
        k in oracle_status
        for k in ("ranking_allowed", "status", "native_cf_oracle")
    )


def evaluate_p2_oracle(oracle_status: dict) -> Tuple[str, str]:
    """Evaluate Phase-2 native CF oracle / mechanism gate.

    Returns
    -------
    (verdict, reason) where verdict is ``pass`` | ``hold`` | ``fail``.

    Rules (fail-closed)
    -------------------
    * **pass** if ``ranking_allowed is True`` OR ``status == "PASS"`` OR
      ``native_cf_oracle == "PASS"`` (and not SCIENCE_HOLD).
    * **hold** if SCIENCE_HOLD, or ``ranking_allowed is False``, or
      ``status == "HOLD"`` (or deferred empty).
    * **fail** if required signal keys are missing and the payload is not
      an empty/deferred deferral.
    """
    if not isinstance(oracle_status, Mapping):
        return "fail", "oracle_status must be a mapping"

    # SCIENCE_HOLD always blocks ranking experiments.
    if _truthy_hold_flag(oracle_status):
        return "hold", "SCIENCE_HOLD asserted — ranking / pilot advancement blocked"

    ranking = oracle_status.get("ranking_allowed", None)
    if ranking is False:
        return "hold", "ranking_allowed=false — Softβ / claim ranking forbidden"

    status_raw = oracle_status.get("status")
    status_u = str(status_raw).strip().upper() if status_raw is not None else ""

    if status_u == "HOLD":
        return "hold", "status=HOLD"

    native = oracle_status.get("native_cf_oracle")
    native_u = str(native).strip().upper() if native is not None else ""

    # Pass conditions (explicit allow signals).
    if ranking is True:
        return "pass", "ranking_allowed=true"
    if status_u == "PASS":
        return "pass", "status=PASS"
    if native_u == "PASS":
        return "pass", "native_cf_oracle=PASS"

    # Deferred / empty: hold (not hard fail) so operators can still stage P2.
    if _is_deferred_empty(oracle_status):
        return "hold", "oracle status empty or deferred — P2 not yet evaluated"

    if not _has_p2_required_keys(oracle_status):
        return (
            "fail",
            "missing required keys (need ranking_allowed and/or status "
            "and/or native_cf_oracle); not deferred",
        )

    # Keys present but no positive allow signal → fail-closed.
    return (
        "fail",
        "no pass signal (ranking_allowed/status/native_cf_oracle) and not HOLD",
    )


def _schema_ok_flag(results_summary: Mapping[str, Any]) -> Optional[bool]:
    """Resolve schema_ok from explicit flag or mode_rmsd presence marker."""
    if "schema_ok" in results_summary:
        return bool(results_summary["schema_ok"])
    if "mode_rmsd_present" in results_summary:
        return bool(results_summary["mode_rmsd_present"])
    # mode_rmsd as list/dict implies present.
    if "mode_rmsd" in results_summary:
        mr = results_summary["mode_rmsd"]
        if mr is False or mr is None:
            return False
        return True
    # mode_rmsd_0..9 style keys
    mode_keys = [k for k in results_summary if str(k).startswith("mode_rmsd")]
    if mode_keys:
        return True
    return None


def evaluate_p3_pilot(results_summary: dict) -> Tuple[str, str]:
    """Evaluate Phase-3 pilot8 exit gate.

    Returns
    -------
    (verdict, reason) where verdict is ``pass`` | ``hold`` | ``fail``.

    Rules
    -----
    * **fail** if ``schema_ok`` is False (mode_rmsd schema missing/broken).
    * **hold** if ``science_hold`` / SCIENCE_HOLD, OR both arms zero success
      with ``n_targets >= 8`` (``bcr_success==0`` and ``s_top10_success==0``).
    * **pass** if ``n_targets >= 1`` and schema_ok and not science_hold.
    """
    if not isinstance(results_summary, Mapping):
        return "fail", "results_summary must be a mapping"

    schema = _schema_ok_flag(results_summary)
    if schema is False:
        return "fail", "schema_ok=false — mode_rmsd columns missing or invalid"

    science_hold = bool(
        results_summary.get("science_hold")
        or results_summary.get("SCIENCE_HOLD")
        or str(results_summary.get("status", "")).strip().upper() == "SCIENCE_HOLD"
    )
    if science_hold:
        return "hold", "science_hold — pilot8 not claim-eligible; block full85"

    try:
        n_targets = int(results_summary.get("n_targets", 0) or 0)
    except (TypeError, ValueError):
        n_targets = 0

    def _as_int(key: str) -> Optional[int]:
        if key not in results_summary:
            return None
        try:
            return int(results_summary[key])
        except (TypeError, ValueError):
            return None

    bcr = _as_int("bcr_success")
    s10 = _as_int("s_top10_success")
    # Also accept alternate key spellings used in campaign summaries.
    if s10 is None:
        s10 = _as_int("S_top10_success")
    if bcr is None:
        bcr = _as_int("BCR_success")

    if (
        n_targets >= 8
        and bcr is not None
        and s10 is not None
        and bcr == 0
        and s10 == 0
    ):
        return (
            "hold",
            "both arms zero success (bcr_success=0 and s_top10_success=0) "
            f"with n_targets={n_targets} — sampling/prep failure pattern",
        )

    if schema is None:
        return (
            "fail",
            "schema_ok / mode_rmsd presence not reported — fail-closed",
        )

    if n_targets >= 1 and schema is True and not science_hold:
        return (
            "pass",
            f"n_targets={n_targets} schema_ok=true — pilot interpretable",
        )

    if n_targets < 1:
        return "fail", "n_targets < 1 — pilot incomplete"

    return "fail", "P3 criteria not met (fail-closed)"


def can_run_p3(phase_state: dict) -> bool:
    """True only when P2 has passed (fail-closed)."""
    if not isinstance(phase_state, Mapping):
        return False
    return _phase_value(phase_state, "P2") == "pass"


def can_run_p4(phase_state: dict) -> bool:
    """True only when P2 **and** P3 have passed (fail-closed)."""
    if not isinstance(phase_state, Mapping):
        return False
    return (
        _phase_value(phase_state, "P2") == "pass"
        and _phase_value(phase_state, "P3") == "pass"
    )


def next_allowed_phase(phase_state: dict) -> Optional[str]:
    """Return the first phase in ``PHASES`` that is not yet ``pass``.

    If all phases are ``pass``, returns ``None``.
    Fail-closed: ``hold`` / ``fail`` / ``pending`` all block later phases,
    so the first non-pass phase is the only one allowed to advance next.
    """
    if not isinstance(phase_state, Mapping):
        return PHASES[0]
    for ph in PHASES:
        if _phase_value(phase_state, ph) != "pass":
            return ph
    return None


def set_phase(
    phase_state: MutableMapping[str, Any],
    phase: str,
    verdict: str,
) -> Dict[str, Any]:
    """Mutate *phase_state* setting *phase* to *verdict*; return the mapping."""
    phase = str(phase).strip().upper()
    if phase not in PHASES:
        raise ValueError(f"unknown phase {phase!r}; expected one of {PHASES}")
    nv = _normalize_verdict(verdict)
    if nv is None or nv not in ("pass", "hold", "fail", "pending"):
        raise ValueError(f"invalid verdict {verdict!r}")
    phase_state[phase] = nv
    nested = phase_state.setdefault("phases", {})
    if isinstance(nested, dict):
        nested[phase] = nv
    return dict(phase_state)


# ---------------------------------------------------------------------------
# Lightweight self-test (also used by CLI --dry-run)
# ---------------------------------------------------------------------------

def self_test() -> List[str]:
    """Run synthetic checks; return list of failure messages (empty = OK)."""
    errors: List[str] = []

    def _check(cond: bool, msg: str) -> None:
        if not cond:
            errors.append(msg)

    # P2 pass signals
    v, _ = evaluate_p2_oracle({"ranking_allowed": True})
    _check(v == "pass", f"P2 ranking_allowed True → pass, got {v}")
    v, _ = evaluate_p2_oracle({"status": "PASS"})
    _check(v == "pass", f"P2 status PASS → pass, got {v}")
    v, _ = evaluate_p2_oracle({"native_cf_oracle": "PASS"})
    _check(v == "pass", f"P2 native_cf_oracle PASS → pass, got {v}")

    # P2 hold
    v, _ = evaluate_p2_oracle({"ranking_allowed": False})
    _check(v == "hold", f"P2 ranking_allowed False → hold, got {v}")
    v, _ = evaluate_p2_oracle({"status": "HOLD"})
    _check(v == "hold", f"P2 status HOLD → hold, got {v}")
    v, _ = evaluate_p2_oracle({"SCIENCE_HOLD": True, "ranking_allowed": True})
    _check(v == "hold", f"P2 SCIENCE_HOLD overrides pass → hold, got {v}")
    v, _ = evaluate_p2_oracle({})
    _check(v == "hold", f"P2 empty deferred → hold, got {v}")

    # P2 fail (signal keys missing, not deferred)
    v, _ = evaluate_p2_oracle({"n_targets": 2, "deferred": False})
    _check(v == "fail", f"P2 missing required keys → fail, got {v}")

    # P3
    v, _ = evaluate_p3_pilot(
        {"n_targets": 8, "schema_ok": True, "bcr_success": 1, "s_top10_success": 2}
    )
    _check(v == "pass", f"P3 healthy pilot → pass, got {v}")
    v, _ = evaluate_p3_pilot(
        {"n_targets": 8, "schema_ok": True, "bcr_success": 0, "s_top10_success": 0}
    )
    _check(v == "hold", f"P3 both-zero pattern → hold, got {v}")
    v, _ = evaluate_p3_pilot({"n_targets": 8, "schema_ok": False})
    _check(v == "fail", f"P3 schema_ok False → fail, got {v}")
    v, _ = evaluate_p3_pilot(
        {"n_targets": 8, "schema_ok": True, "science_hold": True}
    )
    _check(v == "hold", f"P3 science_hold → hold, got {v}")
    v, _ = evaluate_p3_pilot(
        {"n_targets": 1, "mode_rmsd_present": True}
    )
    _check(v == "pass", f"P3 mode_rmsd_present → pass, got {v}")

    # Serial gates
    st = empty_phase_state()
    _check(not can_run_p3(st), "can_run_p3 pending → False")
    _check(not can_run_p4(st), "can_run_p4 pending → False")
    set_phase(st, "P2", "pass")
    _check(can_run_p3(st), "can_run_p3 after P2 pass → True")
    _check(not can_run_p4(st), "can_run_p4 without P3 → False")
    set_phase(st, "P3", "hold")
    _check(not can_run_p4(st), "can_run_p4 with P3 hold → False")
    set_phase(st, "P3", "pass")
    _check(can_run_p4(st), "can_run_p4 P2+P3 pass → True")

    # next_allowed_phase
    st2 = empty_phase_state()
    _check(next_allowed_phase(st2) == "P0", "next on empty → P0")
    for ph in ("P0", "P1", "P2"):
        set_phase(st2, ph, "pass")
    _check(next_allowed_phase(st2) == "P3", "next after P0-P2 → P3")

    _check(MATRIX_MD5_PIN == "9dc93717dfed0698006d88dd6a9627bc", "matrix pin")
    _check(list(PHASES) == ["P0", "P1", "P2", "P3", "P4", "P5"], "phase order")

    return errors
