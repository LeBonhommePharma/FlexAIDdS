"""Drive benchmark_self_eval contract validation + posteriori on synthetic arms."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "benchmark_self_eval.py"
CONTRACT = ROOT / "workorders" / "BENCHMARK_SELF_EVAL_CONTRACT.md"


def _mod():
    spec = importlib.util.spec_from_file_location("benchmark_self_eval", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(m)
    return m


def test_contract_doc_present():
    assert CONTRACT.is_file()
    text = CONTRACT.read_text()
    assert "a priori" in text.lower()
    assert "a posteriori" in text.lower()
    assert "PHASE4" in text
    assert "9dc9" in text or "9dc93717" in text


def test_apriori_rejects_wrong_near_miss_codes(tmp_path: Path):
    m = _mod()
    ap = {
        "one_variable": "x",
        "panel_class": "NEAR_MISS",
        "codes": ["1J3J", "1N1M"],
        "matrix_pin": "9dc9",
        "no_sec": True,
        "sol9": True,
        "matched_control": True,
        "magnitude_floor": "x",
        "report_tiers_separately": True,
    }
    errs = m.validate_apriori(ap)
    assert any("NEAR_MISS" in e for e in errs)


def test_apriori_ok_near_miss():
    m = _mod()
    ap = json.loads((ROOT / "workorders" / "G4_1_NEAR_MISS_APRIORI.json").read_text())
    assert m.validate_apriori(ap) == []


def test_posteriori_pass_magnitude(tmp_path: Path):
    m = _mod()

    def wr(out, code, bcr, elect):
        d = out / code
        d.mkdir(parents=True)
        (d / "result.csv").write_text(
            "pdb_id,best_cluster_rmsd,conditional_scanned_pool_ceiling,rmsd_to_crystal\n"
            f"{code},{bcr},{bcr},{elect}\n"
        )
        (d / "stdout.log").write_text("")
        (d / "stderr.log").write_text("")

    ctrl, tx = tmp_path / "c", tmp_path / "t"
    for code in ("1N1M", "1L7F"):
        wr(ctrl, code, 4.0, 5.0)
        wr(tx, code, 3.4, 5.0)
    # Production path: [BOOM] is on stderr, not stdout
    (tx / "1N1M" / "stderr.log").write_text(
        "[BOOM] injection #1 at gen 100: re-randomized worst 50/1000 chromosomes\n"
    )
    (tx / "1N1M" / "stdout.log").write_text("Generation:  1999\nDone.\n")
    rec = m.posteriori_g4_1_style(ctrl, {"frac010": tx}, ["1N1M", "1L7F"])
    assert rec["accept_magnitude"] is True
    assert rec["status"] == "PASS"
    assert rec["treatments"]["frac010"]["n_markers"] >= 1


def test_count_marker_reads_stderr_not_only_stdout(tmp_path: Path):
    """Regression: live engine puts [BOOM] on stderr.log only."""
    m = _mod()
    arm = tmp_path / "arm"
    (arm / "1N1M" / "r1").mkdir(parents=True)
    (arm / "1N1M" / "stdout.log").write_text("Generation: 100\n")
    (arm / "1N1M" / "stderr.log").write_text("[BOOM] injection #1\n")
    (arm / "1N1M" / "r1" / "stderr.log").write_text("[BOOM] injection #2\n")
    assert m.count_marker(arm, "[BOOM]") == 2
    assert m.count_marker(arm / "1N1M", "[BOOM]") == 2


def test_live_g4_1_l4_markers_on_stderr():
    """Pin against finished G4.1 OUT — stderr-only BOOM must be counted."""
    m = _mod()
    g4 = Path("/Users/lp.more/flexaidds_results/g4_1_boom_near_miss_20260726_200953")
    if not (g4 / "arm_frac010" / "1L7F" / "result.csv").is_file():
        pytest.skip("G4.1 OUT missing")
    n_ctrl = m.count_marker(g4 / "arm_control", "[BOOM]")
    n_tx = m.count_marker(g4 / "arm_frac010", "[BOOM]")
    assert n_ctrl == 0, f"control must have zero [BOOM], got {n_ctrl}"
    assert n_tx > 0, f"frac010 must have live [BOOM] on stderr, got {n_tx}"
    rec = m.posteriori_g4_1_style(
        g4 / "arm_control",
        {
            "frac005": g4 / "arm_frac005",
            "frac010": g4 / "arm_frac010",
            "frac020": g4 / "arm_frac020",
        },
        ["1N1M", "1L7F"],
    )
    # Magnitude null, but L4 liveness must NOT false-fail markers
    assert rec["control_markers"] == 0
    for name, tx in rec["treatments"].items():
        assert tx["n_markers"] > 0, f"{name} n_markers={tx['n_markers']}"
        # status is FAIL (null magnitude) or PASS_LIVENESS — never "no markers"
        assert tx["status"] in ("FAIL", "PASS_LIVENESS", "PASS")
        assert tx["status"] != "IN_PROGRESS"


def test_cli_validate_contract():
    m = _mod()
    assert m.main(["validate-contract-doc", "--path", str(CONTRACT)]) == 0
