#!/usr/bin/env python3
"""Structural + shipped-code checks for the three-arm comparative methodology.

Drives real modules (bootstrap_3dsig_s_top10, generate_flexaid_inp ARM_SPEC) and
asserts fairness axes match docs/implementation/COMPARATIVE_BENCHMARK_METHODOLOGY.md
and the 3Dsig / METHODOLOGY contracts. No hard-coded success rates.
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
METHOD = ROOT / "docs" / "implementation" / "COMPARATIVE_BENCHMARK_METHODOLOGY.md"
MATRIX_PIN = "9dc93717dfed0698006d88dd6a9627bc"
BOOTSTRAP = ROOT / "scripts" / "bootstrap_3dsig_s_top10.py"
GEN_INP = ROOT / "scripts" / "generate_flexaid_inp.py"
RED_PAIR = ROOT / "scripts" / "run_3dsig_red_pair_serial.sh"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_methodology_document_exists_and_names_three_arms():
    assert METHOD.is_file(), f"missing {METHOD}"
    text = METHOD.read_text(encoding="utf-8")
    # Arms A/B/C named
    assert re.search(r"Arm A.*FlexAID 2015|JCIM-era CF", text, re.I | re.S)
    assert re.search(r"Arm B.*soft-β|Shannon", text, re.I | re.S)
    assert re.search(r"Arm C.*FlexAIDdS|current FlexAIDdS", text, re.I | re.S)
    # Fairness axes
    assert "1000" in text and "2000" in text
    assert "2 000 000" in text or "2000000" in text or "2,000,000" in text
    assert MATRIX_PIN[:8] in text
    assert "S_top10" in text
    assert "10 000" in text or "10000" in text
    # Metric non-mixing + correct JCIM Table 2 labels (skeptic-verified)
    # top-1 = 45.2%, top-10 = 66.7% — never inverted
    assert "45.2%" in text and "66.7%" in text
    # 45.2% must appear with top-1 context, not only as top-10
    assert re.search(
        r"top-1[^\n]{0,80}45\.2%|45\.2%[^\n]{0,80}[Tt]op-1",
        text,
    ), "45.2% must be labeled as JCIM top-1 (Table 2)"
    assert re.search(
        r"top-10[^\n]{0,80}66\.7%|66\.7%[^\n]{0,80}[Tt]op-10",
        text,
    ), "66.7% must be labeled as JCIM top-10 (Table 2)"
    # Forbidden inverted protocol (historical bug)
    assert not re.search(
        r"45\.2%[^\n]{0,60}[Tt]op-10 over 10",
        text,
    ), "do not label 45.2% as top-10-over-10-runs"
    assert "Forbidden" in text
    # Arm B entropy definition
    assert "SoftBetaFreeEnergy" in text or "soft_beta" in text
    assert "TEMPER" in text and "FO" in text


def test_bootstrap_shipped_s_top10_contract():
    """Drive real bootstrap module — inclusive ≤2.0, top-10, 10k default."""
    m = _load(BOOTSTRAP, "bootstrap_3dsig_s_top10")
    assert m.DEFAULT_THRESH == 2.0
    assert m.DEFAULT_BOOTSTRAPS == 10_000
    assert m.TOP_N == 10
    # Inclusive threshold (claim contract)
    assert m.s_top10([2.0]) is True
    assert m.s_top10([2.0001]) is False
    assert m.s_top10([3.0, 1.9, 4.0]) is True
    assert m.s_top10([3.0] * 10) is False
    # Only first 10 ranks count
    assert m.s_top10([3.0] * 10 + [0.5]) is False
    # Sentinels / missing ignored (callers must pass Optional[float], not raw CSV strings)
    assert m.s_top10([-1.0, 1.5]) is True
    assert m.s_top10([None, None, 2.5]) is False
    assert m._finite_rmsd("NA") is None
    assert m._finite_rmsd(-1.0) is None


def test_arm_spec_a_b_science_variables():
    """generate_flexaid_inp ARM_SPEC must implement A CF-only vs B soft-β engine path."""
    gen = _load(GEN_INP, "generate_flexaid_inp")
    a = gen.ARM_SPEC["A"]
    b = gen.ARM_SPEC["B"]
    b0 = gen.ARM_SPEC["B0"]
    assert a["temper"] == 0
    assert a["clusta"] == "CF"
    assert b["temper"] == 21
    assert b["clusta"] == "FO"
    assert b.get("fo_minpts_policy") == "single_literature"
    # B0 is CF control on master, not entropy
    assert b0["temper"] == 0
    assert b0["clusta"] == "CF"
    # Budget constants used by prep (shipped defaults)
    assert int(gen.DEFAULT_POP) == 1000
    assert int(gen.DEFAULT_GEN) == 2000
    assert int(gen.DEFAULT_POP) * int(gen.DEFAULT_GEN) == 2_000_000


def test_red_pair_script_freezes_budget_and_matrix_pin():
    text = RED_PAIR.read_text(encoding="utf-8")
    assert 'FLEXAID_POP="${FLEXAID_POP:-1000}"' in text
    assert 'FLEXAID_GEN="${FLEXAID_GEN:-2000}"' in text
    assert 'FLEXAID_RESTARTS="${FLEXAID_RESTARTS:-10}"' in text
    assert MATRIX_PIN in text
    assert "--dry-run" in text


def test_matrix_file_pin_at_repo_root():
    mat = ROOT / "MC_st0r5.2_6.dat"
    if not mat.is_file():
        mat = ROOT / "WRK" / "MC_st0r5.2_6.dat"
    if not mat.is_file():
        pytest.skip("matrix file not present in checkout")
    import hashlib

    md5 = hashlib.md5(mat.read_bytes()).hexdigest()
    assert md5 == MATRIX_PIN, f"matrix md5 {md5} != production pin {MATRIX_PIN}"


def test_softbeta_identity_header_exists():
    h = ROOT / "LIB" / "SoftBetaFreeEnergy.h"
    assert h.is_file()
    text = h.read_text(encoding="utf-8", errors="replace")
    assert "free_energy" in text
    assert "E_min" in text or "Emin" in text or "ln" in text
