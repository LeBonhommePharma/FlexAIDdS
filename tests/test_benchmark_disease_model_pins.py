"""Pin campaign disease-model numbers to live DUMP_POP OUT (shipped .pop.tsv columns)."""
from __future__ import annotations

import csv
from pathlib import Path

import pytest

DUMP = Path("/Users/lp.more/flexaidds_results/dump_pop_search_miss_20260726_172356")
G42 = Path("/Users/lp.more/flexaidds_results/g4_2_r5_near_miss_20260726_175237")
NEAR = ("1N1M", "1L7F")
GROSS = ("1J3J", "1K3U", "1M2Z")


@pytest.fixture(scope="module")
def dump_ok():
    if not DUMP.is_dir():
        pytest.skip("DUMP_POP OUT missing")
    return DUMP


def test_sub2_count_zero_all_search_miss(dump_ok):
    for t in NEAR + GROSS:
        rows = list(csv.DictReader((dump_ok / t / f"{t}.pop.tsv").open(), delimiter="\t"))
        assert rows, t
        n2 = sum(1 for r in rows if float(r["rmsd_sym"]) <= 2.0)
        assert n2 == 0, f"{t} has sub-2 poses — retention hypothesis reopens"


def test_retention_gap_under_half_angstrom(dump_ok):
    for t in NEAR + GROSS:
        rows = list(csv.DictReader((dump_ok / t / f"{t}.pop.tsv").open(), delimiter="\t"))
        pop_best = min(float(r["rmsd_sym"]) for r in rows)
        elects = [r for r in rows if str(r.get("is_elected", "0")) in ("1", "true", "True")]
        emit_best = min(float(r["rmsd_sym"]) for r in elects)
        assert emit_best - pop_best <= 0.5, (t, emit_best, pop_best)


def test_two_tier_split_pop_bests(dump_ok):
    def pop_best(t):
        rows = list(csv.DictReader((dump_ok / t / f"{t}.pop.tsv").open(), delimiter="\t"))
        return min(float(r["rmsd_sym"]) for r in rows)

    assert pop_best("1N1M") < 3.0
    assert pop_best("1L7F") < 5.0
    assert pop_best("1J3J") > 20.0
    assert pop_best("1K3U") > 10.0
    assert pop_best("1M2Z") > 10.0


def test_1j3j_best_cf_farther_than_min_rmsd(dump_ok):
    rows = list(csv.DictReader((dump_ok / "1J3J" / "1J3J.pop.tsv").open(), delimiter="\t"))
    i_rm = min(range(len(rows)), key=lambda i: float(rows[i]["rmsd_sym"]))
    i_cf = min(range(len(rows)), key=lambda i: float(rows[i]["cf_total"]))
    assert float(rows[i_cf]["rmsd_sym"]) - float(rows[i_rm]["rmsd_sym"]) > 20.0


def test_g4_2_r5_near_miss_no_sub2_if_present():
    if not G42.is_dir():
        pytest.skip("G4.2 R5 OUT missing")
    for t in NEAR:
        r = next(csv.DictReader((G42 / t / "result.csv").open()))
        assert float(r["best_cluster_rmsd"]) > 2.0
        assert float(r["conditional_scanned_pool_ceiling"]) == pytest.approx(
            float(r["best_cluster_rmsd"]), abs=1e-3
        )
