#!/usr/bin/env python3
"""Production PSHARE defaults for classic FlexAID ga.inp generation.

Pilot8 incorrectly wrote SHARESCL 0.20 (ADAPTKCO k4 copy-paste) and SHAREPEK 3.
Engine/config_defaults and JCIM fair logs use SHARESCL 10 / SHAREPEK 5.
"""

from __future__ import annotations

import importlib.util
import sys
import tempfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
GEN = ROOT / "scripts" / "generate_flexaid_inp.py"


def _load_gen():
    spec = importlib.util.spec_from_file_location("generate_flexaid_inp", GEN)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def test_default_sharescl_is_ten():
    mod = _load_gen()
    assert mod.DEFAULT_SHARESCL == 10.0
    assert mod.DEFAULT_SHAREPEK == 5.0
    assert mod.DEFAULT_SHAREALF == 4.0


def test_write_ga_emits_production_pshare(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    mod = _load_gen()
    monkeypatch.delenv("FLEXAIDDS_GA_SHARESCL", raising=False)
    monkeypatch.delenv("FLEXAIDDS_GA_SHAREPEK", raising=False)
    monkeypatch.delenv("FLEXAIDDS_GA_SHAREALF", raising=False)
    path = tmp_path / "ga.inp"
    mod.write_ga(path, pop=100, gen=50, seed=1)
    text = path.read_text()
    assert "SHARESCL 10.00" in text
    assert "SHAREPEK 5.00" in text
    assert "SHAREALF 4.00" in text
    # Must not reintroduce the pilot typo
    assert "SHARESCL 0.20" not in text
    assert "SHAREPEK 3.00" not in text


def test_env_override_sharescl(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    mod = _load_gen()
    monkeypatch.setenv("FLEXAIDDS_GA_SHARESCL", "10")
    monkeypatch.setenv("FLEXAIDDS_GA_SHAREPEK", "5")
    alf, pek, scl = mod.resolve_pshare_knobs()
    assert scl == 10.0
    assert pek == 5.0
    path = tmp_path / "ga.inp"
    mod.write_ga(path, pop=10, gen=5, seed=2, sharealf=alf, sharepek=pek, sharescl=scl)
    assert "SHARESCL 10.00" in path.read_text()


def test_read_lig_source_has_inclusive_latm():
    """Static guard: read_lig must assign latm = atm_cnt (not atm_cnt - 1)."""
    src = (ROOT / "LIB" / "read_lig.cpp").read_text(encoding="utf-8", errors="ignore")
    assert "latm[0] = FA->atm_cnt;" in src or "latm[0]=FA->atm_cnt;" in src
    # Bug line must not return as the live assignment
    assert "latm[0] = FA->atm_cnt - 1;" not in src
