"""Structural + evidence tests for campaign methodology gates (B1/B3, BOOM, pb_clash).

These drive real repo artifacts and shipped script entry points — not reimplemented
oracle math.
"""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]


def test_b1_boom_inject_requires_interval_and_fraction():
    """gaboom.cpp must require both interval>0 and fraction>0 (B1)."""
    src = (REPO / "LIB" / "gaboom.cpp").read_text(encoding="utf-8", errors="replace")
    assert "boom_inject_interval > 0" in src
    assert "boom_inject_fraction > 0.0" in src
    # AND in the same condition block
    m = re.search(
        r"boom_inject_interval\s*>\s*0\s*&&\s*[^\n]*boom_inject_fraction\s*>\s*0",
        src,
    )
    assert m, "expected AND of interval and fraction on inject guard"


def test_claim_path_emits_boom_fraction_zero_deliberately():
    """DatasetRunner claim JSON hardcodes boom_inject_fraction 0.0 (anti-collapse)."""
    src = (REPO / "LIB" / "DatasetRunner.cpp").read_text(encoding="utf-8", errors="replace")
    assert "boom_inject_fraction" in src
    # literal emit of 0.0 after the catastrophic-regression comment
    assert re.search(r"boom_inject_fraction.*\n(?:.*\n){0,20}.*<<\s*0\.0", src)


def test_boom_frac_env_override_after_json():
    """config_parser reads FLEXAIDDS_BOOM_FRAC after JSON defaults."""
    src = (REPO / "LIB" / "config_parser.cpp").read_text(encoding="utf-8", errors="replace")
    i_json = src.find('boom_inject_fraction')
    i_env = src.find("FLEXAIDDS_BOOM_FRAC")
    assert i_json >= 0 and i_env > i_json


def test_pb_clash_env_override_present():
    src = (REPO / "LIB" / "config_parser.cpp").read_text(encoding="utf-8", errors="replace")
    assert "FLEXAIDDS_PB_CLASH_WEIGHT" in src
    # env after JSON default assignment
    i_json = src.find("pb_clash_weight")
    i_env = src.find("FLEXAIDDS_PB_CLASH_WEIGHT")
    assert i_env > i_json


def test_pb_clash_burial_oracle_script_exists_and_one_variable():
    script = REPO / "scripts" / "pb_clash_burial_oracle.py"
    assert script.is_file()
    text = script.read_text(encoding="utf-8")
    assert "FLEXAIDDS_PB_CLASH_WEIGHT" in text
    assert "diagnostic/probe_config" in text  # refuse non-production
    assert "CLEAN_PANEL" in text
    # parseable python
    ast.parse(text)


def test_campaign_gate_summary_labels_b1_b3_instrumentation():
    md = (REPO / "workorders" / "CAMPAIGN_GATE_SUMMARY.md").read_text(encoding="utf-8")
    assert "structurally unpassable" in md.lower() or "STRUCTURAL" in md
    assert "B3" in md and "B1" in md
    assert "scientifically invalid" in md.lower() or "INVALID as BOOM" in md
    assert "not" in md.lower() and "docking-quality" in md.lower() or "instrumentation" in md.lower()
    assert "WALL_PILOT_PASS" in md
    assert "full-85" in md.lower() or "Full-85" in md


def test_methodology_points_to_pb_clash_replacement():
    md = (
        REPO / "docs" / "implementation" / "CAMPAIGN_METHODOLOGY_for_Grok.md"
    ).read_text(encoding="utf-8")
    assert "structurally unpassable" in md
    assert "pb_clash" in md.lower() or "PB_CLASH" in md
    assert "BOOM_FRAC" in md or "boom_inject_fraction" in md


def test_pb_clash_oracle_workorder_present():
    p = REPO / "workorders" / "PB_CLASH_ORACLE.md"
    assert p.is_file()
    t = p.read_text(encoding="utf-8")
    assert "FLEXAIDDS_PB_CLASH_WEIGHT" in t
    assert "PASS" in t or "FAIL" in t


def test_inversion_map_script_and_workorder():
    script = REPO / "scripts" / "native_elected_cf_inversion_map.py"
    assert script.is_file()
    text = script.read_text(encoding="utf-8")
    assert "SCORING-LOCKED" in text and "SEARCH-MISS" in text
    assert "probe_cf" in text and "dock_config" in text
    wo = REPO / "workorders" / "INVERSION_MAP.md"
    assert wo.is_file()
    wt = wo.read_text(encoding="utf-8")
    assert "SEARCH-MISS" in wt
    # must not claim full-85
    assert "full-85" in wt.lower() or "Full-85" in wt or "not claimed" in wt.lower()


def test_next_campaign_step_is_single_lever_and_blocks_exhausted_paths():
    """NEXT_CAMPAIGN_STEP must name one primary experiment and ban closed gates."""
    p = REPO / "workorders" / "NEXT_CAMPAIGN_STEP.md"
    assert p.is_file(), "missing workorders/NEXT_CAMPAIGN_STEP.md"
    t = p.read_text(encoding="utf-8")
    assert "inversion" in t.lower() or "Native–Elected" in t or "Native-Elected" in t
    assert "SCORING-LOCKED" in t and "SEARCH-MISS" in t
    # one primary, not full-85 / WAL re-run / interval-only BOOM
    assert "full-85" in t.lower() or "Full-85" in t
    assert "Rejected" in t or "rejected" in t
    assert "WAL_COERCIVE" in t or "WAL" in t
    assert "BOOM_INTERVAL" in t or "interval-only" in t.lower()
    # separates liveness / docking claims
    assert "Not claimed" in t or "not claimed" in t.lower() or "Instrumentation" in t
    # grounds in existing workorders
    assert "CAMPAIGN_GATE_SUMMARY" in t or "E10" in t
    assert "probe_cf" in t and "dock_config" in t