#!/usr/bin/env python3
"""CLI: Phase P5 — bootstrap S_top10 + COMPARATIVE_TABLE.md + thin iCloud sync.

Implements docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md Phase 5.

Uses the **shipped** ``scripts/bootstrap_3dsig_s_top10.py`` via subprocess when
arm dirs contain ``result.csv``. Does not reimplement the bootstrap median.

Usage:
  python3 scripts/comparative_p5_analyze.py --campaign comparative_pilot8 --dry-run
  python3 scripts/comparative_p5_analyze.py --campaign comparative_full85 --arms A,B,C
  python3 scripts/comparative_p5_analyze.py --campaign pilot8 --local-root ~/flexaidds_results

Exit codes:
  0  PHASE=P5 status=pass (table written; bootstrap invokable)
  1  PHASE=P5 status=fail
  2  usage / argument error

Copyright 2026 Le Bonhomme Pharma
SPDX-License-Identifier: Apache-2.0
"""
from __future__ import annotations

import sys
from pathlib import Path

# Prefer in-repo package without requiring install.
_REPO = Path(__file__).resolve().parents[1]
_PY = _REPO / "python"
if _PY.is_dir() and str(_PY) not in sys.path:
    sys.path.insert(0, str(_PY))

from flexaidds.comparative_phases.p5_analyze import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
