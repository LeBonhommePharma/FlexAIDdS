"""sampling_power: library / S_top10 success, independent of CF election.

This is the sampling half of the Tier-1 split gate. docking_power asks
whether the *elected* top-N pose is near-native. sampling_power asks
whether *any* pose in the emitted library is near-native.
"""

from flexaidds.dataset_runner.metrics import (
    PoseScore,
    compute_all_metrics,
    docking_power,
    sampling_power,
)


def _pose(target, rmsd, score, rank=1):
    return PoseScore(
        target_id=target,
        ligand_id=target + "_lig",
        pose_rank=rank,
        rmsd=rmsd,
        enthalpy_score=score,
        entropy_correction=0.0,
        total_score=score,
        is_active=True,
    )


def test_sampling_succeeds_when_near_native_is_not_rank_one():
    """The 35016255184 phenotype: rank-1 is a decoy, a later pose is ≤2 Å."""
    poses = [
        _pose("1gpk", 4.71, -22108.0, rank=1),
        _pose("1gpk", 1.94, -17977.0, rank=6),
    ]
    assert docking_power(poses, top_n=1) == 0.0
    assert sampling_power(poses) == 1.0


def test_sampling_failed_targets_stay_in_the_denominator():
    poses = [
        _pose("HIT", 1.2, -10.0),
        _pose("MISS", 5.0, -9.0),
    ]
    assert sampling_power(poses, n_targets=4) == 1.0 / 4.0


def test_sampling_sentinels_never_succeed():
    for sentinel in (-1.0, 999.0):
        assert sampling_power([_pose("T1", sentinel, -10.0)]) == 0.0


def test_sampling_empty_roster_with_pinned_n_is_zero_not_absent():
    """Silent top-1=0 from empty results must measure sampling_power=0.0."""
    assert sampling_power([], n_targets=4) == 0.0


def test_compute_all_metrics_emits_sampling_power():
    poses = [
        _pose("A", 4.0, -20.0, rank=1),
        _pose("A", 1.5, -5.0, rank=2),
        _pose("B", 3.5, -15.0, rank=1),
    ]
    r = compute_all_metrics(poses, requested=["sampling_power", "docking_power_top1"])
    assert r["sampling_power"] == 0.5
    assert r["docking_power_top1"] == 0.0
