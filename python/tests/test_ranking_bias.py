"""Unit tests for ranking-bias diagnostics (pure functions on CF/RMSD tables)."""

from __future__ import annotations

import pytest

from flexaidds.diagnostics.ranking_bias import (
    PoseCF,
    near_native_missed_by_top1,
    rank_by_cf,
    rank_of_best_rmsd,
    scoring_pathology_gap,
    search_never_beats_seed,
    spearman_cf_rmsd,
    wal_over_abs_com,
)


def test_rank_by_cf_lower_is_better():
    poses = [
        PoseCF("a", cf=-10.0, rmsd=3.0),
        PoseCF("b", cf=-50.0, rmsd=5.0),
        PoseCF("c", cf=-20.0, rmsd=1.0),
    ]
    ranked = rank_by_cf(poses)
    assert [p.name for p in ranked] == ["b", "c", "a"]


def test_near_native_missed_when_good_rmsd_not_top1():
    # Near-native at CF=-20, decoy top1 at CF=-50 with high RMSD
    poses = [
        PoseCF("decoy", cf=-50.0, rmsd=8.0),
        PoseCF("nativeish", cf=-20.0, rmsd=1.2),
        PoseCF("mid", cf=-15.0, rmsd=4.0),
    ]
    assert near_native_missed_by_top1(poses) is True
    assert rank_of_best_rmsd(poses) == 2


def test_near_native_not_missed_when_top1_is_good():
    poses = [
        PoseCF("good", cf=-50.0, rmsd=1.0),
        PoseCF("decoy", cf=-20.0, rmsd=6.0),
    ]
    assert near_native_missed_by_top1(poses) is False


def test_no_oracle_near_native_returns_false():
    poses = [
        PoseCF("a", cf=-50.0, rmsd=5.0),
        PoseCF("b", cf=-40.0, rmsd=4.0),
    ]
    assert near_native_missed_by_top1(poses) is False


def test_scoring_pathology_when_decoy_beats_seed():
    # top1 much better than seed (gap = -100)
    assert scoring_pathology_gap(cf_top1=-189.0, cf_seed=-54.0, pathology_cut=5.0) is True
    # search loses to seed (gap positive)
    assert scoring_pathology_gap(cf_top1=-29.0, cf_seed=-83.0, pathology_cut=5.0) is False


def test_search_never_beats_seed():
    poses = [
        PoseCF("a", cf=-29.0),
        PoseCF("b", cf=-10.0),
    ]
    assert search_never_beats_seed(poses, cf_seed=-83.0) is True
    assert search_never_beats_seed(poses, cf_seed=-5.0) is False


def test_spearman_positive_when_cf_tracks_rmsd():
    # better CF (lower) with lower RMSD → positive Spearman
    poses = [
        PoseCF("n", cf=-50.0, rmsd=1.0),
        PoseCF("m", cf=-30.0, rmsd=3.0),
        PoseCF("f", cf=-10.0, rmsd=8.0),
    ]
    rho = spearman_cf_rmsd(poses)
    assert rho is not None
    assert rho > 0.9


def test_spearman_negative_when_anti_correlated():
    # better CF with higher RMSD → ranking anti-correlated with geometry
    poses = [
        PoseCF("decoy", cf=-100.0, rmsd=12.0),
        PoseCF("mid", cf=-50.0, rmsd=6.0),
        PoseCF("near", cf=-10.0, rmsd=1.5),
    ]
    rho = spearman_cf_rmsd(poses)
    assert rho is not None
    assert rho < -0.9


def test_wal_over_abs_com():
    assert wal_over_abs_com(444.0, -218.0) == pytest.approx(444.0 / 218.0)
    assert wal_over_abs_com(None, -10.0) is None
