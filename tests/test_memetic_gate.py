#!/usr/bin/env python3
"""Structural test: use_memetic gated by Phase2 pb_clash PASS or legacy WALL."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_use_memetic_field_and_gate_in_parser():
    h = (ROOT / "LIB" / "flexaid.h").read_text(encoding="utf-8", errors="replace")
    assert "use_memetic" in h
    assert "memetic_armed_at_gen" in h
    cpp = (ROOT / "LIB" / "config_parser.cpp").read_text(encoding="utf-8", errors="replace")
    assert "FA->use_memetic" in cpp
    assert "FLEXAIDDS_WALL_PILOT_PASS" in cpp
    assert "FLEXAIDDS_PB_CLASH_PHASE2_PASS" in cpp
    assert "memetic_gate.h" in cpp
    assert "resolve_use_memetic_from_env" in cpp
    gab = (ROOT / "LIB" / "gaboom.cpp").read_text(encoding="utf-8", errors="replace")
    assert "use_memetic" in gab
    assert "memetic_armed_at_gen" in gab  # durable arm marker on GA path
    gate = (ROOT / "LIB" / "memetic_gate.h").read_text(encoding="utf-8", errors="replace")
    assert "resolve_use_memetic_from_env" in gate
    assert "FLEXAIDDS_PB_CLASH_PHASE2_PASS" in gate
    assert "FLEXAIDDS_WALL_PILOT_PASS" in gate


def test_memetic_gate_header_logic_via_subprocess_cpp():
    """Drive shipped resolve_use_memetic_from_env through the C++ unit binary if built."""
    import os
    import subprocess

    bin_path = ROOT / "build_wave0" / "test_memetic_gate"
    if not bin_path.is_file():
        # still assert source exists for structural gate
        assert (ROOT / "tests" / "test_memetic_gate.cpp").is_file()
        return
    env = os.environ.copy()
    r = subprocess.run(
        [str(bin_path)],
        capture_output=True,
        text=True,
        env=env,
        timeout=60,
        check=False,
    )
    assert r.returncode == 0, r.stdout + r.stderr
    assert "PASSED" in r.stdout
