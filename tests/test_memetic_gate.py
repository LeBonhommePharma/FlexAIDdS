#!/usr/bin/env python3
"""Structural test: use_memetic is a real FA field gated by WALL_PILOT_PASS."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_use_memetic_field_and_gate_in_parser():
    h = (ROOT / "LIB" / "flexaid.h").read_text(encoding="utf-8", errors="replace")
    assert "use_memetic" in h
    cpp = (ROOT / "LIB" / "config_parser.cpp").read_text(encoding="utf-8", errors="replace")
    assert "FA->use_memetic" in cpp
    assert "FLEXAIDDS_WALL_PILOT_PASS" in cpp
    # Must set enable path, not only fprintf
    assert "use_memetic = 1" in cpp or "use_memetic=1" in cpp
    gab = (ROOT / "LIB" / "gaboom.cpp").read_text(encoding="utf-8", errors="replace")
    assert "use_memetic" in gab  # consumed in GA post-path
