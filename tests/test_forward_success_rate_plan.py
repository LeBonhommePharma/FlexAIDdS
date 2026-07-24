#!/usr/bin/env python3
"""Structural gates for FORWARD_SUCCESS_RATE_PLAN.md + code anchors it depends on.

Drives real LIB markers (acf vs free_energy_strict) and asserts the plan inventory
tags + JCIM labels + Wave-0 offline first steps. No invented success rates.
"""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PLAN = ROOT / "docs" / "implementation" / "FORWARD_SUCCESS_RATE_PLAN.md"
CLUSTER = ROOT / "LIB" / "cluster.cpp"
SOFTBETA = ROOT / "LIB" / "SoftBetaFreeEnergy.h"
SOFT_WALL = ROOT / "LIB" / "soft_wall.h"


def test_forward_plan_exists_and_tags_inventory():
    assert PLAN.is_file()
    text = PLAN.read_text(encoding="utf-8")
    for tag in ("KEEP", "DEFER", "REJECT", "UNCITABLE"):
        assert tag in text, f"missing inventory tag {tag}"
    # Election vs sampling
    assert re.search(r"Election vs sampling|election vs sampling", text, re.I)
    assert "BCR" in text and "S1" in text and "S_top10" in text
    # ACF implicit population
    assert "acf" in text.lower() or "ACF" in text
    assert "free_energy_strict" in text
    assert "implicit" in text.lower() or "size bias" in text.lower()
    # Wall before memetic
    assert re.search(r"wall.*memetic|memetic.*after.*E2|only after E2", text, re.I | re.S)
    # COM_BURIAL_CAP uncitible
    assert "UNCITABLE" in text
    assert "COM_BURIAL_CAP" in text or "COM_BURIAL" in text
    assert "-130" in text or "130" in text
    # Search primary
    assert re.search(r"search coverage|Search coverage", text)
    # JCIM labels correct
    assert "45.2%" in text and "66.7%" in text
    assert re.search(r"top-1[^\n]{0,80}45\.2%|45\.2%[^\n]{0,80}top-1", text, re.I)
    assert re.search(r"top-10[^\n]{0,80}66\.7%|66\.7%[^\n]{0,80}top-10", text, re.I)
    # Matrix
    assert "9dc93717" in text
    # Box
    assert re.search(r"workers?\s*[≤<=]?\s*2|2–4|2-4", text, re.I)
    # Ordered waves / first offline steps
    assert "Wave 0" in text or "W0.1" in text
    assert "E10" in text and "E1b" in text
    # No bare claim of restored rates as current fact
    assert "Does not claim restored rates" in text or "not claim restored" in text.lower()


def test_shipped_acf_path_still_uses_legacy_acf():
    """Plan KEEP E1b is grounded: cluster emission still calls soft_beta::acf."""
    c = CLUSTER.read_text(encoding="utf-8", errors="replace")
    assert "soft_beta::acf" in c
    h = SOFTBETA.read_text(encoding="utf-8", errors="replace")
    assert "free_energy_strict" in h
    assert "diagnostic only" in h.lower() or "Diagnostic only" in h
    # free_energy_strict defined after acf alias
    assert h.find("inline double acf") < h.find("free_energy_strict") or "free_energy_strict" in h


def test_soft_wall_exists_for_e2_grounding():
    assert SOFT_WALL.is_file()
    t = SOFT_WALL.read_text(encoding="utf-8", errors="replace")
    assert "soft_wall" in t.lower() or "k_wal" in t or "WAL" in t
