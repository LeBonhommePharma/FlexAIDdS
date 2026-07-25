#!/usr/bin/env python3
"""CLI: Phase P1 — pin comparative arm binaries (A/B/C) and write receipts.

Loads ``docs/implementation/arm_pins.json``, inspects
``$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/{A,B,C}/``, SHA256s
present Mach-Os, writes
``$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/receipts/arm_{X}_binary.json``.

Does **not** compile FlexAID. Reconstruction stubs require ``--allow-reconstruction``
and never invent binary digests for missing files.

Usage:
  python3 scripts/comparative_p1_pin_binaries.py
  python3 scripts/comparative_p1_pin_binaries.py --local-root ~/flexaidds_results
  python3 scripts/comparative_p1_pin_binaries.py --allow-reconstruction

Exit codes:
  0  A and B both present with different SHA256; or reconstruction receipts OK
  1  A and B both present but SHA256 identical (claim split fail)
  2  A/B missing without --allow-reconstruction, or fatal pin load error

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

from flexaidds.comparative_phases.p1_binaries import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
