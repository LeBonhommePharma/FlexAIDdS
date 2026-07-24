"""Regression: claim path always emits coarse_init.enabled=true + boom_frac=0.

Background (2026-07-24): mode UNSET left coarse_init false while seed_fraction=0.
Tight pockets (1M2Z) then started gen-0 with only random orientations; every
chromosome hit Vcontacts CLASH_THRESHOLD=1e4 (CF=10000 all gens). Native still
scored fine via SCORE_NATIVE coordinate override. boom_frac=1.0 for UNSET also
wiped any partial progress every 100 gens.

Source of truth: LIB/DatasetRunner.cpp dock_config emission.
"""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
SRC = (ROOT / "LIB" / "DatasetRunner.cpp").read_text()


def test_coarse_init_enabled_unconditionally():
    # Literal true in the dock_config emission (no mode ternary).
    assert '<< "    \\"enabled\\": true,\\n"' in SRC
    block_start = SRC.find("Coarse pocket scan")
    assert block_start > 0
    block = SRC[block_start : block_start + 900]
    assert "BenchmarkMode::AUTONOMOUS" not in block, (
        "coarse_init must not be mode-gated; UNSET claim canaries need it ON"
    )
    assert "CLASH_THRESHOLD" in block  # documents the 1M2Z pathology


def test_boom_frac_hard_zero_on_claim_path():
    block_start = SRC.find("No-seed modes must preserve accumulated GA progress")
    assert block_start > 0
    block = SRC[block_start : block_start + 1200]
    # Hard-coded 0.0 emission (not mode-ternary / effective_boom_frac).
    assert re.search(r"<<\s*0\.0\s*\n", block), block[-200:]
    # Must not call effective_boom_frac() for the emission (comment may name it).
    assert not re.search(r"<<\s*protocol_cfg_\.effective_boom_frac", block)


def test_unset_mode_forces_seed_elitism_off():
    # UNSET is in the force-off list for both receipt and selection override.
    assert re.search(
        r"DEFINED_CLEFT_REDOCK \|\|\s*\n\s*config\.mode == BenchmarkMode::AUTONOMOUS \|\|\s*\n\s*config\.mode == BenchmarkMode::UNSET",
        SRC,
    ) or (
        "BenchmarkMode::UNSET" in SRC
        and "receipt_seed_elitism = false" in SRC
    )


if __name__ == "__main__":
    test_coarse_init_enabled_unconditionally()
    test_boom_frac_hard_zero_on_claim_path()
    test_unset_mode_forces_seed_elitism_off()
    print("ALL PASS")
