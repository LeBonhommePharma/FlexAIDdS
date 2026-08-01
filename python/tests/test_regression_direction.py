"""Regression-direction registry: distance metrics must flag when they WORSEN.

Guards the bug where check_regressions inferred direction from the substrings
"rmse"/"mae", so mean_rmsd (spelled rmsD) took the higher-is-better branch and
a 37% worse RMSD scored as no-regression.
"""
import re

import yaml
from pathlib import Path

from flexaidds.dataset_runner import metrics
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


def test_median_rmsd_flags_on_increase():
    """Only produced distance metrics are asserted here.

    This test used to cover crossdock_mean_rmsd, crossdock_median_rmsd and
    selectivity_log_error as well. Nothing computes any of those -- the
    assertions were exercising a registry entry that had never met a value, so
    what they pinned was a guess, not a behaviour. They are removed with their
    entries and come back with their producers.
    """
    dr = _result({"median_rmsd": 3.90}, {"median_rmsd": 1.95})
    assert dr.check_regressions()["median_rmsd"] is True


def test_higher_is_better_metric_flags_on_decrease():
    dr = _result({"docking_power_top1": 0.0}, {"docking_power_top1": 0.70})
    assert dr.check_regressions()["docking_power_top1"] is True


def test_registries_are_disjoint():
    assert not (_LOWER_IS_BETTER & _HIGHER_IS_BETTER)


def _produced_metric_names() -> set[str]:
    """Metric names `compute_all_metrics` can actually write into its results.

    Derived by reading the source rather than kept as a hand-maintained list.
    A list would be another table with no upstream -- exactly the object this
    module exists to shrink -- and it would go stale silently the moment someone
    dispatched a metric without editing it. This cannot go stale: the commit
    that adds a producer is the same commit that creates the obligation below.

    Crude on purpose: a regex over source text, not an import-and-introspect.
    `results["x"] = ...` assignments are not visible without executing the
    function, and executing it needs poses. Crude-and-derived beats
    clean-and-hand-maintained here; the failure mode of the regex (missing a
    producer written some other way) is a metric exempted that should not be,
    which is the same state as today rather than a new one.
    """
    src = (
        Path(metrics.__file__).read_text()
        if hasattr(metrics, "__file__")
        else ""
    )
    return set(re.findall(r'results\["([a-z0-9_]+)"\]', src))


def test_every_declared_baseline_has_a_direction():
    """A metric a dataset declares AND something produces must have a direction.

    Declaration alone is not enough. Eleven names are declared in dataset YAMLs
    with nothing computing them -- requiring a direction for those would red
    `main` for whoever next touches a Python file, and the thing they would be
    asked to resolve is the direction of a metric they are not implementing:
    a red with nobody present who can answer it.

    So the obligation is keyed on the producer, not the declaration. The moment
    a producer appears the exemption evaporates and this test demands an entry
    -- in the same edit, from the person writing the assignment, who is the only
    one who knows what the stored value means.

    shannon_energy_collapse stays a named exception: it IS produced, and its
    direction is genuinely undetermined pending someone who knows the physics.
    """
    known_gap = {"shannon_energy_collapse"}
    produced = _produced_metric_names()
    assert produced, "producer scan found nothing -- the regex or the file moved"
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
                if metric not in produced:
                    continue  # nothing can contradict a direction for it yet
                if metric not in _LOWER_IS_BETTER and metric not in _HIGHER_IS_BETTER:
                    undeclared.add(f"{metric} ({root.name} tree)")
    assert not (
        {u.split(" (")[0] for u in undeclared} - known_gap
    ), f"metrics with no declared direction: {sorted(undeclared)}"


def test_registry_holds_only_metrics_something_can_produce():
    """The registry must not describe a metric nothing computes.

    An entry with no producer cannot be contradicted by any behaviour, so it is
    not a claim that can be wrong -- it is one that instructs. This is the guard
    that keeps the registry to what is falsifiable; it is the same criterion as
    the exemption above, applied from the other side.
    """
    produced = _produced_metric_names()
    registered = set(_LOWER_IS_BETTER) | set(_HIGHER_IS_BETTER)
    assert not (registered - produced), (
        "registry entries with no producer (a direction nothing can contradict): "
        f"{sorted(registered - produced)}"
    )
