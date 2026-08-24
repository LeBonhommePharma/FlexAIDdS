#!/usr/bin/env python3
"""Wave 4: five-way RMSD cross-check named in METHODOLOGY.md §0."""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "rmsd_five_way_crosscheck.py"


def _load():
    spec = importlib.util.spec_from_file_location("rmsd_five_way_crosscheck", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_methodology_lists_exactly_the_five_implementations():
    text = (ROOT / "METHODOLOGY.md").read_text(encoding="utf-8")
    assert "python/flexaidds/benchmark.py::compute_rmsd" in text
    assert "benchmarks/astex_repro/score_reference.py" in text
    assert "benchmarks/astex_repro/score_offline.py" in text
    assert "LIB/calc_rmsd.cpp::calc_Hungarian_RMSD" in text
    assert "dataset::hungarian_rmsd" in text
    numbered = re.findall(r"^\s+(\d)\.\s+\*\*", text, flags=re.M)
    assert numbered[:5] == ["1", "2", "3", "4", "5"]


def test_registry_matches_methodology_five():
    mod = _load()
    assert len(mod.FIVE_METHODS) == 5
    assert [m["id"] for m in mod.FIVE_METHODS] == [1, 2, 3, 4, 5]
    assert mod.CLAIM_CUTOFF_A == 2.0
    assert mod.claim_success(2.0) is True
    assert mod.claim_success(2.0001) is False
    assert mod.claim_success(1.99) is True


def test_python_methods_agree_on_inplace_claim_cutoff():
    mod = _load()
    hit = mod.python_methods_agree_on_claim(1.5)
    assert hit["agreed_claim_success"] is True
    assert hit["method1_inplace"] == 1.5
    assert hit["method1_superposed"] < 1e-6
    assert hit["superposed_must_not_be_claim_metric"] is True

    edge = mod.python_methods_agree_on_claim(2.0)
    assert edge["agreed_claim_success"] is True

    miss = mod.python_methods_agree_on_claim(2.5)
    assert miss["agreed_claim_success"] is False


def test_cxx_claim_cutoff_tests_exist():
    text = (ROOT / "tests" / "test_dataset_runner.cpp").read_text(encoding="utf-8")
    assert "TEST(RmsdClaimCutoff, Rank0InPlaceHitAt1p5A)" in text
    assert "TEST(RmsdClaimCutoff, Rank0InPlaceInclusiveTwoAngstrom)" in text
    assert "TEST(RmsdClaimCutoff, Rank0InPlaceMissAt2p5A)" in text
    assert "TEST(RmsdCrossCheck, AlignedTyping_SolversAgreeOnASharedAssignment)" in text
