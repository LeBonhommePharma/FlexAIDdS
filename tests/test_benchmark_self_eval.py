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


def test_validate_pins_fails_without_accept(tmp_path: Path):
    """S2: missing evidence/accept.txt must fail the shipped checker."""
    m = _mod()
    out = tmp_path / "out"
    arm = out / "arm_control"
    (arm / "bin").mkdir(parents=True)
    # Stamped binary present so SHA alone is not enough
    blob = b"fake-flexaidds-binary-for-pin-test-control"
    (arm / "bin" / "FlexAIDdS.stamped").write_bytes(blob)
    rep = m.validate_pins_report(out)
    assert rep["ok"] is False
    assert any("accept.txt" in e for e in rep["errors"])
    rc = m.main(["validate-pins", "--out", str(out)])
    assert rc == 2


def test_validate_pins_fails_without_arm_sha(tmp_path: Path):
    """S2: accept present but no per-arm SHA / stamped binary → fail."""
    m = _mod()
    out = tmp_path / "out"
    (out / "arm_control").mkdir(parents=True)
    (out / "evidence").mkdir(parents=True)
    (out / "evidence" / "accept.txt").write_text("ACCEPT_X=False\nstatus=FAIL\n")
    # No bin, no arm_pins.json
    rep = m.validate_pins_report(out, arms=["control"])
    assert rep["ok"] is False
    assert any("binary SHA256" in e for e in rep["errors"])
    assert m.main(["validate-pins", "--out", str(out), "--arms", "control"]) == 2


def test_validate_pins_ok_with_accept_and_arm_pins(tmp_path: Path):
    """S2: complete pin pack (accept + arm_pins SHA) passes shipped entry point."""
    m = _mod()
    out = tmp_path / "out"
    (out / "arm_control").mkdir(parents=True)
    (out / "arm_tx").mkdir(parents=True)
    (out / "evidence").mkdir(parents=True)
    (out / "evidence" / "accept.txt").write_text(
        "ACCEPT_DEMO=False\nstatus=PASS_LIVENESS\n"
    )
    sha_c = "a" * 64
    sha_t = "b" * 64
    pins = {
        "matrix_pin": "9dc93717dfed0698006d88dd6a9627bc",
        "shared_binary": False,
        "arms": {
            "control": {"binary_sha256": sha_c, "git_tip": "deadbeef"},
            "tx": {"binary_sha256": sha_t},
        },
    }
    (out / "evidence" / "arm_pins.json").write_text(json.dumps(pins, indent=2))
    rep = m.validate_pins_report(out)
    assert rep["ok"] is True, rep["errors"]
    assert rep["arm_sha256"]["control"] == sha_c
    assert rep["arm_sha256"]["tx"] == sha_t
    assert m.main(["validate-pins", "--out", str(out)]) == 0
    # CLI JSON path also drives shipped report
    assert m.main(["validate-pins", "--out", str(out), "--json"]) == 0


def test_validate_pins_ok_from_stamped_binary_hash(tmp_path: Path):
    """S2: hashing arm_*/bin/FlexAIDdS.stamped satisfies per-arm SHA."""
    import hashlib

    m = _mod()
    out = tmp_path / "out"
    (out / "evidence").mkdir(parents=True)
    (out / "evidence" / "accept.txt").write_text("ACCEPT_Y=True\nstatus=PASS\n")
    blob = b"stamped-binary-content-xyz"
    for arm in ("control", "mut_gran"):
        d = out / f"arm_{arm}" / "bin"
        d.mkdir(parents=True)
        (d / "FlexAIDdS.stamped").write_bytes(blob)
    expect = hashlib.sha256(blob).hexdigest()
    rep = m.validate_pins_report(out)
    assert rep["ok"] is True, rep["errors"]
    assert rep["arm_sha256"]["control"] == expect
    assert rep["arm_sha256"]["mut_gran"] == expect
    assert m.main(["validate-pins", "--out", str(out)]) == 0


def test_contract_doc_mentions_s2_pin_pack():
    text = CONTRACT.read_text()
    assert "accept.txt" in text
    assert "binary_sha256" in text
    assert "arm_pins" in text
    assert "validate-pins" in text
