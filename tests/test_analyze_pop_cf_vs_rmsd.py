"""Drive analyze_pop_cf_vs_rmsd on engine .pop.tsv schema."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_pop_cf_vs_rmsd.py"


def _mod():
    spec = importlib.util.spec_from_file_location("analyze_pop_cf_vs_rmsd", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(m)
    return m


def test_scoring_pull_when_best_cf_is_far(tmp_path: Path):
    m = _mod()
    # min rmsd at 12Å with mediocre CF; best CF at 25Å
    tsv = (
        "idx\tcluster\trmsd_to_head\trmsd_raw\trmsd_sym\tcf_total\tcf_com\tcf_wal\tpose_id\tis_elected\n"
        "0\t0\t0\t12.0\t12.0\t-10.0\t-10\t0\t-1\t0\n"
        "1\t0\t0\t25.0\t25.0\t-90.0\t-90\t0\t0\t1\n"
        "2\t0\t0\t13.0\t13.0\t-11.0\t-11\t0\t-1\t0\n"
        "3\t0\t0\t24.0\t24.0\t-80.0\t-80\t0\t-1\t0\n"
        "4\t0\t0\t14.0\t14.0\t-12.0\t-12\t0\t-1\t0\n"
        "5\t0\t0\t23.0\t23.0\t-70.0\t-70\t0\t-1\t0\n"
        "6\t0\t0\t15.0\t15.0\t-13.0\t-13\t0\t-1\t0\n"
        "7\t0\t0\t22.0\t22.0\t-60.0\t-60\t0\t-1\t0\n"
        "8\t0\t0\t16.0\t16.0\t-14.0\t-14\t0\t-1\t0\n"
        "9\t0\t0\t21.0\t21.0\t-50.0\t-50\t0\t-1\t0\n"
    )
    code = "1AAA"
    d = tmp_path / code
    d.mkdir()
    (d / f"{code}.pop.tsv").write_text(tsv)
    rec = m.diagnose(m.load_pop(d / f"{code}.pop.tsv"))
    assert rec["label"] == "SCORING_PULL"
    assert rec["rmsd_at_min_cf"] == pytest.approx(25.0)
    assert rec["min_rmsd_sym"] == pytest.approx(12.0)


def test_search_fail_when_best_cf_near_min_rmsd_but_far(tmp_path: Path):
    m = _mod()
    # all far; best CF at min rmsd
    lines = ["idx\tcluster\trmsd_to_head\trmsd_raw\trmsd_sym\tcf_total\tcf_com\tcf_wal\tpose_id\tis_elected"]
    for i in range(10):
        # rmsd 12+i*0.1, cf improves (more neg) as rmsd decreases slightly
        r = 12.0 + i * 0.1
        c = -100.0 + i  # best CF at i=0 rmsd=12
        lines.append(f"{i}\t0\t0\t{r}\t{r}\t{c}\t{c}\t0\t-1\t0")
    code = "1BBB"
    d = tmp_path / code
    d.mkdir()
    (d / f"{code}.pop.tsv").write_text("\n".join(lines) + "\n")
    rec = m.diagnose(m.load_pop(d / f"{code}.pop.tsv"))
    assert rec["label"] == "SEARCH_FAIL"
    assert rec["delta_rmsd_bestCF_minus_minRMSD"] <= 1.0


def test_cli_on_real_dump_if_present():
    m = _mod()
    live = Path("/Users/lp.more/flexaidds_results/dump_pop_search_miss_20260726_172356")
    if not (live / "1J3J" / "1J3J.pop.tsv").is_file():
        pytest.skip("live DUMP_POP OUT not present")
    rc = m.main([str(live), "--codes", "1J3J,1M2Z,1K3U"])
    assert rc == 0
