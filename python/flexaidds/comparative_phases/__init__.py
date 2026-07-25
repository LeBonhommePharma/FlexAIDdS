# Copyright 2026 Louis-Philippe Morency / NRGlab, Universite de Montreal
# SPDX-License-Identifier: Apache-2.0
"""Fail-closed serial phase gates for the FlexAIDdS comparative pipeline.

Phases (serial, fail-closed): P0 → P1 → P2 → P3 → P4 → P5

- P2  mechanism / native CF oracle
- P3  pilot8 dual-arm schema + zero-success hold
- P4  full85 (requires P2 pass AND P3 pass)

Authoritative matrix pin matches Phase 0 / COMPARATIVE_GOAL_METHODOLOGY.md.
"""

from .gates import (
    MATRIX_MD5_PIN,
    PHASES,
    GateVerdict,
    can_run_p3,
    can_run_p4,
    empty_phase_state,
    evaluate_p2_oracle,
    evaluate_p3_pilot,
    load_phase_state,
    next_allowed_phase,
    save_phase_state,
)

__all__ = [
    "MATRIX_MD5_PIN",
    "PHASES",
    "GateVerdict",
    "can_run_p3",
    "can_run_p4",
    "empty_phase_state",
    "evaluate_p2_oracle",
    "evaluate_p3_pilot",
    "load_phase_state",
    "next_allowed_phase",
    "save_phase_state",
]
