"""Baseline tolerance must be additive slack, not a scaling of the baseline.

`baseline * (1 ± tol)` is correct only while every baseline is positive.
erds_specificity declares ``target_specificity_zscore: -2.50`` (more negative
= better binding), and scaling a negative baseline moves the threshold AWAY
from zero -- demanding -2.625, which is stricter than the target itself.  The
metric could then only pass by beating its baseline by 5%, and sitting exactly
on target flagged as a regression.
"""
import pytest

from flexaidds.dataset_runner.runner import (
    DatasetResult,
    _HIGHER_IS_BETTER,
    _LOWER_IS_BETTER,
)


class _Cfg:
    def __init__(self, baselines, tol=0.05):
        self.expected_baselines = baselines
        self.baseline_tolerance = tol


def _flags(metrics, baselines):
    dr = DatasetResult.__new__(DatasetResult)
    dr.metrics = metrics
    dr.config = _Cfg(baselines)
    dr.regression_flags = {}
    dr.inconclusive_metrics = []
    return dr.check_regressions()


# `mean_rmsd` is lower-is-better; `docking_power_top1` is higher-is-better.
# Both are produced and registered, so they exercise the real lookup.
_LOWER = "mean_rmsd"
_HIGHER = "docking_power_top1"


def test_registry_membership_assumed_by_this_module_holds():
    """Guard the fixtures themselves: a re-file would silently void the tests."""
    assert _LOWER in _LOWER_IS_BETTER
    assert _HIGHER in _HIGHER_IS_BETTER


@pytest.mark.parametrize(
    "measured,flagged",
    [
        (-3.00, False),  # better than target
        (-2.50, False),  # EXACTLY on target -- scaling form flagged this
        (-2.40, False),  # 4% worse, inside the 5% slack
        (-2.30, True),   # 8% worse, genuinely outside
        (0.00, True),    # no signal at all
    ],
)
def test_negative_baseline_lower_is_better(measured, flagged):
    """Slack runs toward zero for a negative lower-is-better baseline."""
    assert _flags({_LOWER: measured}, {_LOWER: -2.50})[_LOWER] is flagged


@pytest.mark.parametrize(
    "measured,flagged",
    [
        (-2.00, False),  # above a negative baseline
        (-2.50, False),  # exactly on target
        (-3.00, True),   # below the slack floor
    ],
)
def test_negative_baseline_higher_is_better(measured, flagged):
    """The mirror branch: -5.0 must not read as a regression from -2.50."""
    assert _flags({_HIGHER: measured}, {_HIGHER: -2.50})[_HIGHER] is flagged


@pytest.mark.parametrize("baseline", [0.70, 2.10, 2.30, 0.85, 0.30])
def test_positive_baselines_are_unchanged(baseline):
    """The fix is a no-op wherever the baseline is positive.

    Every baseline in the repo today is positive except one, so this pins the
    property that made the change safe to land ahead of anything else.
    """
    tol = 0.05
    for metric, scaled in (
        (_LOWER, baseline * (1 + tol)),
        (_HIGHER, baseline * (1 - tol)),
    ):
        # Just inside the old threshold must stay unflagged; just outside must
        # flag -- identical behaviour to the scaling form.
        inside = scaled - 0.001 if metric is _LOWER else scaled + 0.001
        outside = scaled + 0.001 if metric is _LOWER else scaled - 0.001
        assert _flags({metric: inside}, {metric: baseline})[metric] is False
        assert _flags({metric: outside}, {metric: baseline})[metric] is True


def test_zero_baseline_has_no_slack():
    """abs(0) * tol == 0: a zero baseline permits nothing, which is honest."""
    assert _flags({_LOWER: 0.001}, {_LOWER: 0.0})[_LOWER] is True
    assert _flags({_LOWER: 0.0}, {_LOWER: 0.0})[_LOWER] is False
