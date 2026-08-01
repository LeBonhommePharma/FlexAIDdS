"""Regression-direction registry: distance metrics must flag when they WORSEN.

Guards the bug where check_regressions inferred direction from the substrings
"rmse"/"mae", so mean_rmsd (spelled rmsD) took the higher-is-better branch and
a 37% worse RMSD scored as no-regression.
"""
import yaml
from pathlib import Path

from flexaidds.dataset_runner.runner import (
    DatasetResult,
    _HIGHER_IS_BETTER,
    _LOWER_IS_BETTER,
)


class _Cfg:
    def __init__(self, baselines, tol=0.05):
        self.expected_baselines = baselines
        self.baseline_tolerance = tol


def _result(metrics, baselines):
    dr = DatasetResult.__new__(DatasetResult)
    dr.metrics = metrics
    dr.config = _Cfg(baselines)
    dr.regression_flags = {}
    dr.inconclusive_metrics = []
    return dr


def test_mean_rmsd_worse_than_baseline_is_flagged():
    """The exact 1gpk case: 3.1557 measured against a 2.30 baseline."""
    dr = _result({"mean_rmsd": 3.1557}, {"mean_rmsd": 2.30})
    assert dr.check_regressions()["mean_rmsd"] is True


def test_mean_rmsd_better_than_baseline_is_not_flagged():
    dr = _result({"mean_rmsd": 1.80}, {"mean_rmsd": 2.30})
    assert dr.check_regressions()["mean_rmsd"] is False


def test_median_rmsd_and_crossdock_variants_flag_on_increase():
    for name, base in (
        ("median_rmsd", 1.95),
        ("crossdock_mean_rmsd", 2.0),
        ("crossdock_median_rmsd", 2.0),
        ("selectivity_log_error", 0.5),
    ):
        dr = _result({name: base * 2}, {name: base})
        assert dr.check_regressions()[name] is True, name


def test_higher_is_better_metric_flags_on_decrease():
    dr = _result({"docking_power_top1": 0.0}, {"docking_power_top1": 0.70})
    assert dr.check_regressions()["docking_power_top1"] is True


def test_registries_are_disjoint():
    assert not (_LOWER_IS_BETTER & _HIGHER_IS_BETTER)


def test_every_declared_baseline_has_a_direction():
    """Any metric a dataset declares must be in exactly one registry.

    shannon_energy_collapse is the known exception: direction undetermined,
    warns at runtime rather than silently assuming.
    """
    known_gap = {"shannon_energy_collapse"}
    repo = Path(__file__).resolve().parents[2]
    # BOTH dataset trees. benchmarks/datasets is canonical (CANONICAL.md), but
    # the package copy is a declared mirror that has drifted in 10 of 12 shared
    # pairs (#337) -- reading only one tree would pass by luck of the drift not
    # currently being in metric names.
    roots = [
        repo / "benchmarks" / "datasets",
        repo / "python" / "flexaidds" / "dataset_runner" / "datasets",
    ]
    roots = [r for r in roots if r.is_dir()]
    assert roots, "no dataset directory found in either tree"
    undeclared = set()
    for root in roots:
        for f in sorted(root.glob("*.yaml")):
            raw = yaml.safe_load(f.read_text()) or {}
            for metric in (raw.get("expected_baselines") or {}):
                if metric not in _LOWER_IS_BETTER and metric not in _HIGHER_IS_BETTER:
                    undeclared.add(f"{metric} ({root.name} tree)")
    assert not (
        {u.split(" (")[0] for u in undeclared} - known_gap
    ), f"metrics with no declared direction: {sorted(undeclared)}"
