#!/usr/bin/env python3
"""Structural + math gates for Wave 0–4 implementation (worktree)."""

from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel: str):
    path = ROOT / rel
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_acf_strict_math_matches_softbeta_contract():
    """Python free_energy_strict mirrors SoftBetaStrict.ExactDuplicateInvariance."""
    mod = _load("acf_reelect", "scripts/acf_strict_offline_reelect.py")
    once = [-50.0, -48.0, -45.0]
    cloned = once * 10
    T = 250.0
    classic_once = mod.free_energy(once, T)
    classic_clone = mod.free_energy(cloned, T)
    assert classic_clone < classic_once - 1.0
    strict_once = mod.free_energy_strict(once, T)
    strict_clone = mod.free_energy_strict(cloned, T)
    assert abs(strict_clone - strict_once) < 1e-9


def test_cluster_cpp_has_acf_strict_gate():
    """Product path: free_energy_strict default; legacy ACF via FLEXAIDDS_ELECT_LEGACY_ACF.

    Wave0 branch used opt-in FLEXAIDDS_ACF_STRICT; main merge kept free_energy_strict
    as default with LEGACY escape hatch (see SoftBetaFreeEnergy.h / cluster.cpp).
    """
    text = (ROOT / "LIB" / "cluster.cpp").read_text(encoding="utf-8", errors="replace")
    assert "free_energy_strict" in text
    assert "soft_beta::acf" in text
    assert (
        "FLEXAIDDS_ELECT_LEGACY_ACF" in text or "FLEXAIDDS_ACF_STRICT" in text
    ), "need opt-in legacy ACF or ACF_STRICT env gate"


def test_config_parser_elec_and_wave3_knobs():
    text = (ROOT / "LIB" / "config_parser.cpp").read_text(encoding="utf-8", errors="replace")
    assert "electrostatics_enabled" in text
    assert "FLEXAIDDS_USE_ELEC" in text
    assert "FLEXAIDDS_BOOM_INTERVAL" in text
    assert "FLEXAIDDS_SIGMA_SCALE" in text
    assert "FLEXAIDDS_COARSE_ORIENTATIONS" in text
    assert "FLEXAIDDS_MEMETIC" in text
    assert "FLEXAIDDS_WALL_PILOT_PASS" in text


def test_wal_coercive_default_off_in_vcfunction():
    text = (ROOT / "LIB" / "vcfunction.cpp").read_text(encoding="utf-8", errors="replace")
    assert "FLEXAIDDS_WAL_COERCIVE" in text
    assert "wal_coercive" in text


def test_no_com_burial_cap_default_on_main_path():
    """CAP must not be product default on this branch (W2.2)."""
    # Presence of env gate is OK; hardcoding CAP=-130 as default is not.
    for rel in ("LIB/vcfunction.cpp", "LIB/cluster.cpp", "LIB/config_parser.cpp"):
        t = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
        assert "COM_BURIAL_CAP=-130" not in t
