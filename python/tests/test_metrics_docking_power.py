"""Docking-power denominator and sentinel-RMSD integrity.

These pin the two ways a reported success rate can be silently inflated:
dropping targets that produced no usable pose from the denominator, and
letting a sentinel RMSD (-1 = not computed, 999 = no pose) satisfy the
``rmsd < 2.0`` comparison.
"""

from flexaidds.dataset_runner.metrics import PoseScore, docking_power


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


def test_sentinel_rmsd_is_never_a_success():
    for sentinel in (-1.0, 999.0):
        assert docking_power([_pose("T1", sentinel, -10.0)]) == 0.0


def test_failed_targets_stay_in_the_denominator():
    poses = [
        _pose("HIT", 1.2, -10.0),
        _pose("MISS", -1.0, -9.0),   # timed out / no RMSD computed
        _pose("NOPOSE", 999.0, -8.0),
    ]
    # 1 of 3 attempted, not 1 of 1 with a usable pose.
    assert docking_power(poses) == 1.0 / 3.0


def test_n_targets_pins_the_full_dataset_size():
    # Only two targets reported anything at all; the rest never ran.
    poses = [_pose("HIT", 1.2, -10.0), _pose("MISS", 5.0, -9.0)]
    assert docking_power(poses, n_targets=85) == 1.0 / 85.0


def test_sentinel_pose_still_occupies_its_rank():
    # A sentinel that outranks the near-native pose must block the top-1 hit
    # rather than being filtered out and promoting the good pose.
    poses = [_pose("T1", 999.0, -20.0, rank=1), _pose("T1", 1.0, -5.0, rank=2)]
    assert docking_power(poses, top_n=1) == 0.0
    assert docking_power(poses, top_n=2) == 1.0
