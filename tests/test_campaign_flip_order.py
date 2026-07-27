"""Drive campaign_flip_order on synthetic result.csv trees + live DUMP_POP elect."""
from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "campaign_flip_order.py"
DUMP = Path("/Users/lp.more/flexaidds_results/dump_pop_search_miss_20260726_172356")


def _mod():
    spec = importlib.util.spec_from_file_location("campaign_flip_order", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(m)
    return m


def _write_result(path: Path, bcr: float, elect: float):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "pdb_id,best_cluster_rmsd,conditional_scanned_pool_ceiling,rmsd_to_crystal,success_rmsd\n"
        f"X,{bcr},{bcr},{elect},0\n"
    )


def test_g4_1_magnitude_promotes_boom(tmp_path: Path):
    m = _mod()
    ctrl = tmp_path / "ctrl"
    tx = tmp_path / "tx"
    for code, bcr in [("1N1M", 4.0), ("1L7F", 4.0)]:
        _write_result(ctrl / code / "result.csv", bcr, bcr + 1)
        _write_result(tx / code / "result.csv", bcr - 0.6, bcr + 1)  # mean dBCR -0.6
    # Production path: [BOOM] on stderr.log (stdout has gens only)
    (tx / "1N1M").mkdir(exist_ok=True)
    (tx / "1N1M" / "stdout.log").write_text("Generation:  100\n")
    (tx / "1N1M" / "stderr.log").write_text(
        "[BOOM] injection #1 at gen 100: re-randomized worst 50/1000 chromosomes\n"
    )
    (ctrl / "1N1M").mkdir(exist_ok=True)
    (ctrl / "1N1M" / "stdout.log").write_text("no boom here\n")
    (ctrl / "1N1M" / "stderr.log").write_text("no boom here\n")
    # L4 is applied in main() CLI path via boom_l4; unit-test boom_l4 + evaluate
    l4 = m.boom_l4(tx)
    assert l4["n_boom_markers"] >= 1
    assert m.boom_l4(ctrl)["n_boom_markers"] == 0
    rec = m.evaluate_g4_1(
        m.collect_arm(ctrl, ["1N1M", "1L7F"]),
        {"frac010": m.collect_arm(tx, ["1N1M", "1L7F"])},
        codes=["1N1M", "1L7F"],
    )
    assert rec["accept_g4_1"] is True
    assert rec["flip_order"]["rule"] == "G4.1_BOOM_hits_magnitude"


def test_boom_l4_reads_stderr_and_restart_logs(tmp_path: Path):
    m = _mod()
    arm = tmp_path / "arm_frac010"
    (arm / "1L7F" / "r1").mkdir(parents=True)
    (arm / "1L7F" / "stdout.log").write_text("Generation: 7999\n")
    (arm / "1L7F" / "stderr.log").write_text("[BOOM] injection #1\n[BOOM] injection #2\n")
    (arm / "1L7F" / "r1" / "stderr.log").write_text("[BOOM] injection #3\n")
    rec = m.boom_l4(arm)
    assert rec["n_boom_markers"] == 3
    assert rec["wipeout_flag"] is False


def test_live_g4_1_boom_l4_not_false_negative():
    """Finished G4.1: stdout-only would report 0; stderr scan must see markers."""
    m = _mod()
    g4 = Path("/Users/lp.more/flexaidds_results/g4_1_boom_near_miss_20260726_200953")
    if not (g4 / "arm_frac010" / "1N1M" / "result.csv").is_file():
        pytest.skip("G4.1 OUT missing")
    ctrl = m.boom_l4(g4 / "arm_control")
    tx = m.boom_l4(g4 / "arm_frac010")
    assert ctrl["n_boom_markers"] == 0
    assert tx["n_boom_markers"] > 0
    # CLI path must not l4_fail after fix
    import subprocess
    import json

    p = subprocess.run(
        [
            "python3",
            str(SCRIPT),
            "g4_1",
            "--control",
            str(g4 / "arm_control"),
            "--treatment",
            "frac010",
            str(g4 / "arm_frac010"),
            "--codes",
            "1N1M,1L7F",
            "--json",
        ],
        capture_output=True,
        text=True,
    )
    start = p.stdout.find("{")
    rec = json.loads(p.stdout[start:])
    assert rec["l4_control"]["n_boom_markers"] == 0
    assert rec["l4_treatments"]["frac010"]["n_boom_markers"] > 0
    assert "l4_fail" not in rec["treatments"]["frac010"]
    # magnitude still null on real OUT
    assert rec["accept_g4_1"] is False
    assert rec["flip_order"]["rule"] == "G4.1_null_or_l4_fail"


def test_g4_1_null_flips_to_election_priority(tmp_path: Path):
    m = _mod()
    ctrl = tmp_path / "ctrl"
    tx = tmp_path / "tx"
    for code in ("1N1M", "1L7F"):
        _write_result(ctrl / code / "result.csv", 4.0, 5.0)
        _write_result(tx / code / "result.csv", 3.9, 5.0)  # only -0.1
    rec = m.evaluate_g4_1(
        m.collect_arm(ctrl, ["1N1M", "1L7F"]),
        {"frac010": m.collect_arm(tx, ["1N1M", "1L7F"])},
        codes=["1N1M", "1L7F"],
    )
    assert rec["accept_g4_1"] is False
    assert "election" in rec["flip_order"]["action"].lower() or "G4.3" in rec["flip_order"]["action"]


def test_election_offline_1n1m_live():
    m = _mod()
    if not (DUMP / "1N1M" / "1N1M.pop.tsv").is_file():
        pytest.skip("DUMP_POP missing")
    rec = m.evaluate_election_offline(
        DUMP / "1N1M" / "1N1M.pop.tsv",
        DUMP / "1N1M" / "result.csv",
    )
    assert rec["pop_best_rmsd_sym"] < 3.0
    assert rec["actual_elect_rmsd"] > rec["pop_best_rmsd_sym"]
    assert rec["would_flip_to_election_P0"] is True
