"""Unit tests for ranking-bias diagnostics (pure stdlib; no package root import)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_MOD_PATH = Path(__file__).resolve().parents[1] / "flexaidds" / "diagnostics" / "ranking_bias.py"


def _load():
    """Load ranking_bias without importing flexaidds package (avoids numpy)."""
    spec = importlib.util.spec_from_file_location("ranking_bias_mod", _MOD_PATH)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    # register so dataclasses can resolve module
    import sys
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def rb():
    return _load()


def test_rank_by_cf_lower_is_better(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("a", cf=-10.0, rmsd=3.0),
        PoseCF("b", cf=-50.0, rmsd=5.0),
        PoseCF("c", cf=-20.0, rmsd=1.0),
    ]
    ranked = rb.rank_by_cf(poses)
    assert [p.name for p in ranked] == ["b", "c", "a"]


def test_near_native_missed_when_good_rmsd_not_top1(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("decoy", cf=-50.0, rmsd=8.0),
        PoseCF("nativeish", cf=-20.0, rmsd=1.2),
        PoseCF("mid", cf=-15.0, rmsd=4.0),
    ]
    assert rb.near_native_missed_by_top1(poses) is True
    assert rb.rank_of_best_rmsd(poses) == 2


def test_near_native_not_missed_when_top1_is_good(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("good", cf=-50.0, rmsd=1.0),
        PoseCF("decoy", cf=-20.0, rmsd=6.0),
    ]
    assert rb.near_native_missed_by_top1(poses) is False


def test_no_oracle_near_native_returns_false(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("a", cf=-50.0, rmsd=5.0),
        PoseCF("b", cf=-40.0, rmsd=4.0),
    ]
    assert rb.near_native_missed_by_top1(poses) is False


def test_scoring_pathology_when_decoy_beats_seed(rb):
    assert rb.scoring_pathology_gap(cf_top1=-189.0, cf_seed=-54.0, pathology_cut=5.0) is True
    assert rb.scoring_pathology_gap(cf_top1=-29.0, cf_seed=-83.0, pathology_cut=5.0) is False


def test_search_never_beats_seed(rb):
    PoseCF = rb.PoseCF
    poses = [PoseCF("a", cf=-29.0), PoseCF("b", cf=-10.0)]
    assert rb.search_never_beats_seed(poses, cf_seed=-83.0) is True
    assert rb.search_never_beats_seed(poses, cf_seed=-5.0) is False


def test_spearman_positive_when_cf_tracks_rmsd(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("n", cf=-50.0, rmsd=1.0),
        PoseCF("m", cf=-30.0, rmsd=3.0),
        PoseCF("f", cf=-10.0, rmsd=8.0),
    ]
    rho = rb.spearman_cf_rmsd(poses)
    assert rho is not None and rho > 0.9


def test_spearman_negative_when_anti_correlated(rb):
    PoseCF = rb.PoseCF
    poses = [
        PoseCF("decoy", cf=-100.0, rmsd=12.0),
        PoseCF("mid", cf=-50.0, rmsd=6.0),
        PoseCF("near", cf=-10.0, rmsd=1.5),
    ]
    rho = rb.spearman_cf_rmsd(poses)
    assert rho is not None and rho < -0.9


def test_wal_over_abs_com(rb):
    assert rb.wal_over_abs_com(444.0, -218.0) == pytest.approx(444.0 / 218.0)
    assert rb.wal_over_abs_com(None, -10.0) is None
