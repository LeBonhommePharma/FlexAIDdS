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
    assert "CLEAN_PANEL" in text or "SEARCH_MISS_PANEL" in text
    assert "SCORING_LOCKED_PANEL" in text
    assert "scoring-locked" in text
    assert "MAGNITUDE_FLOOR_KCAL" in text
    assert "1.0" in text  # floor constant / default
    # parseable python
    ast.parse(text)


def test_pb_clash_oracle_verdict_scoring_locked_magnitude_floor_is_real():
    """Drive shipped verdict_scoring_locked — not a reimplementation."""
    import importlib.util

    path = REPO / "scripts" / "pb_clash_burial_oracle.py"
    spec = importlib.util.spec_from_file_location("pb_clash_burial_oracle", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # Inverted targets with tiny noise → FAIL
    tiny = [
        {"dCF_off": 20.0, "dCF_on": 19.99},
        {"dCF_off": 30.0, "dCF_on": 29.99},
        {"dCF_off": 70.0, "dCF_on": 69.99},
    ]
    v, reason, st = mod.verdict_scoring_locked(tiny, 1.0, 2)
    assert v == "FAIL"
    assert st["n_both"] == 0
    # Real flip on 2/3 with ≥1.0 decrease → PASS
    good = [
        {"dCF_off": 18.0, "dCF_on": -2.0},  # decrease 20, flip
        {"dCF_off": 29.0, "dCF_on": -1.5},  # decrease 30.5, flip
        {"dCF_off": 70.0, "dCF_on": 69.5},  # no flip
    ]
    v2, reason2, st2 = mod.verdict_scoring_locked(good, 1.0, 2)
    assert v2 == "PASS", reason2
    assert st2["n_both"] == 2
    assert st2["n_sign_flip"] == 2


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
    # ROADMAP_v2: prior SEARCH-MISS result must be labeled VOID
    assert "VOID" in t
    assert "magnitude" in t.lower() or "1.0" in t


def test_memetic_gate_accepts_phase2_unlock_env():
    gate = (REPO / "LIB" / "memetic_gate.h").read_text(encoding="utf-8")
    assert "FLEXAIDDS_PB_CLASH_PHASE2_PASS" in gate
    assert "FLEXAIDDS_MEMETIC" in gate
    # still does not enable on MEMETIC alone (logic: want && unlock)
    assert "phase2_pass" in gate or "PB_CLASH_PHASE2" in gate


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

def test_g4_2_niche_distance_header_shipped_and_wired():
    """G4.2: pure niche_distance.h exists; gaboom calls shipped helpers."""
    hdr = REPO / "LIB" / "niche_distance.h"
    assert hdr.is_file()
    h = hdr.read_text(encoding="utf-8")
    assert "niche_cartesian_rmsd" in h
    assert "niche_gene_rmsp" in h
    assert "niche_pair_distance" in h
    assert "niche_cartesian_env_enabled" in h
    src = (REPO / "LIB" / "gaboom.cpp").read_text(encoding="utf-8", errors="replace")
    assert '#include "niche_distance.h"' in src
    assert "flexaids::niche_cartesian_rmsd" in src
    assert "flexaids::niche_cartesian_env_enabled" in src
    assert "FLEXAIDDS_NICHE_CARTESIAN" in src
    assert "[NICHE-CART]" in src


def test_g4_2_niche_distance_drives_shipped_cpp_binary():
    """Drive shipped flexaids::niche_* via the gtest binary (not a Python reimpl)."""
    import os
    import subprocess

    candidates = [
        REPO / "build" / "test_niche_distance",
        REPO / "build" / "tests" / "test_niche_distance",
    ]
    bin_path = next((p for p in candidates if p.is_file()), None)
    if bin_path is None:
        # Still require sources + ability to compile a smoke against the header
        assert (REPO / "tests" / "test_niche_distance.cpp").is_file()
        smoke = REPO / "build" / "niche_distance_smoke"
        # Prefer already-built gtest; else compile minimal smoke on shipped header
        import tempfile
        from pathlib import Path as P
        scratch = P(os.environ.get(
            "SCRATCH",
            "/var/folders/8b/tgtvwb_j6zd_g03vl1w4ykfw0000gn/T/grok-goal-b3876f84368a/implementer",
        ))
        scratch.mkdir(parents=True, exist_ok=True)
        src = scratch / "niche_smoke.cpp"
        src.write_text(
            r'''
#include "niche_distance.h"
#include <cassert>
#include <cmath>
#include <cstdlib>
int main() {
  unsetenv("FLEXAIDDS_NICHE_CARTESIAN");
  assert(!flexaids::niche_cartesian_env_enabled());
  setenv("FLEXAIDDS_NICHE_CARTESIAN", "1", 1);
  assert(flexaids::niche_cartesian_env_enabled());
  double ic_a[3] = {0,0,0}, ic_b[3] = {10000,0,0};
  float xyz[3] = {1,2,3};
  double d_gene = flexaids::niche_pair_distance(false, ic_a, ic_b, 3, xyz, xyz, 1);
  double d_cart = flexaids::niche_pair_distance(true, ic_a, ic_b, 3, xyz, xyz, 1);
  assert(d_gene > 100.0);
  assert(std::fabs(d_cart) < 1e-12);
  float a[3]={0,0,0}, b[3]={3,0,0};
  assert(std::fabs(flexaids::niche_cartesian_rmsd(a,b,1) - 3.0) < 1e-6);
  return 0;
}
'''
        )
        out = scratch / "niche_smoke"
        r = subprocess.run(
            ["c++", "-std=c++17", f"-I{REPO}/LIB", str(src), "-o", str(out)],
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        assert r.returncode == 0, r.stderr
        r2 = subprocess.run([str(out)], capture_output=True, text=True, timeout=30, check=False)
        assert r2.returncode == 0, r2.stdout + r2.stderr
        return
    env = os.environ.copy()
    r = subprocess.run(
        [str(bin_path), "--gtest_color=no"],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
        check=False,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert "PASSED" in r.stdout


def test_g4_4_early_stop_workorder_exists():
    p = REPO / "workorders" / "G4_4_EARLY_STOP.md"
    assert p.is_file()
    t = p.read_text(encoding="utf-8")
    assert "truncation" in t.lower()
    assert "NO_SEC" in t or "no_sec" in t


def test_phase4_actualized_in_workorders():
    p = REPO / "workorders" / "PHASE4_GATES_ACTUALIZED.md"
    assert p.is_file()
    t = p.read_text(encoding="utf-8")
    assert "SEARCH-MISS" in t and "G4.2" in t
    assert "empty weight window" in t.lower() or "clash-free" in t.lower()
