"""Pose-accuracy aggregates: mean_rmsd / median_rmsd.

These metrics are declared as baselines by pose-reproduction datasets
(astex_diverse, hap2, ...) but had no implementation in ``compute_all_metrics``
until #341 — so a run that produced real per-pose RMSDs still reported them as
"never measured" and the #326 completeness gate failed on a metric the code
simply never computed.

Definition (mirrors ``benchmark.BenchmarkSummary``): per target take the
best (minimum) valid pose RMSD, then mean / median across targets.  A target
with no valid RMSD contributes nothing, and when *no* target yields a valid
best-pose RMSD the metric is left ABSENT so the gate reports it unmeasured
rather than emitting a misleading 0.0.
"""

from flexaidds.dataset_runner.metrics import PoseScore, compute_all_metrics

_RMSD_KEYS = ["mean_rmsd", "median_rmsd"]


def _pose(target, rmsd, rank=1):
    return PoseScore(
        target_id=target,
        ligand_id=f"{target}_lig",
        pose_rank=rank,
        rmsd=rmsd,
        enthalpy_score=-10.0,
        entropy_correction=0.0,
        total_score=-10.0,
        is_active=True,
    )


def test_best_pose_per_target_then_mean_and_median():
    poses = [
        _pose("A", 1.0, rank=1), _pose("A", 2.5, rank=2),  # best = 1.0
        _pose("B", 3.0, rank=1), _pose("B", 2.0, rank=2),  # best = 2.0
        _pose("C", 3.0, rank=1),                            # best = 3.0
    ]
    r = compute_all_metrics(poses, requested=_RMSD_KEYS)
    assert r["mean_rmsd"] == (1.0 + 2.0 + 3.0) / 3.0
    assert r["median_rmsd"] == 2.0


def test_even_target_count_median_is_mean_of_two_middle():
    poses = [_pose("A", 1.0), _pose("B", 2.0), _pose("C", 3.0), _pose("D", 6.0)]
    r = compute_all_metrics(poses, requested=_RMSD_KEYS)
    assert r["median_rmsd"] == (2.0 + 3.0) / 2.0
    assert r["mean_rmsd"] == (1.0 + 2.0 + 3.0 + 6.0) / 4.0


def test_sentinel_rmsds_are_excluded_from_the_aggregate():
    # -1 (not computed) and 999 (no pose) must not enter the best-pose set.
    poses = [
        _pose("A", -1.0, rank=1), _pose("A", 1.5, rank=2),  # best = 1.5, not -1
        _pose("B", 999.0),                                   # no valid RMSD -> dropped
    ]
    r = compute_all_metrics(poses, requested=_RMSD_KEYS)
    assert r["mean_rmsd"] == 1.5
    assert r["median_rmsd"] == 1.5


def test_no_valid_rmsd_leaves_metrics_absent_for_the_gate():
    # Every pose has a sentinel RMSD -> the completeness gate must see the
    # metric as unmeasured (absent), NOT as 0.0.
    poses = [_pose("A", -1.0), _pose("B", 999.0)]
    r = compute_all_metrics(poses, requested=_RMSD_KEYS)
    assert "mean_rmsd" not in r
    assert "median_rmsd" not in r


def test_single_target_reproduces_the_1gpk_run():
    # The first real Tier-1 run: one target (1gpk), 11 near-identical poses at
    # ~3.156 A.  best = min; mean == median == that best for a single target.
    rmsds = [
        3.155735, 3.155741, 3.155754, 3.155757, 3.155762, 3.155792,
        3.155814, 3.155819, 3.155820, 3.155822, 3.155878,
    ]
    poses = [_pose("1gpk", v, rank=i + 1) for i, v in enumerate(rmsds)]
    r = compute_all_metrics(poses, requested=_RMSD_KEYS)
    best = min(rmsds)
    assert r["mean_rmsd"] == best
    assert r["median_rmsd"] == best
    # And it is now PRESENT — the exact metric the gate reported "never
    # measured" on run 30680006724.
    assert set(_RMSD_KEYS) <= r.keys()
