"""#326 benchmark verdict gate: a crashed / empty / unmeasured run must be
INCONCLUSIVE, never a silent pass.

These test the liveness / productivity / completeness gates (1-3) directly,
without docking anything — the check Honey noted "would have caught this in
2019": feed the gate a run where the binary exited non-zero and assert it does
not read as "no regression".
"""

import os
import types

import pytest

from flexaidds.dataset_runner.cli import _benchmark_inconclusive_reasons


def _dr(slug, *, completed=(), failed=(), crashes=0, total_poses=0,
        exit_codes=None, missing_metrics=()):
    """Build a minimal DatasetResult-shaped object for the gate helper."""
    return types.SimpleNamespace(
        config=types.SimpleNamespace(slug=slug),
        targets_completed=list(completed),
        targets_failed=list(failed),
        flexaid_crashes=crashes,
        total_poses=total_poses,
        entry_exit_codes=exit_codes or {},
        inconclusive_metrics=list(missing_metrics),
    )


def test_crashed_binary_is_inconclusive():
    # The exact production bug: every entry crashed, zero poses, no metrics.
    dr = _dr(
        "astex_diverse",
        failed=["1gpk", "1mq6", "1n2j", "1t46", "1t9b"],
        crashes=5,
        total_poses=0,
        exit_codes={"1gpk/lig": -6},
        missing_metrics=["docking_power_top1", "mean_rmsd"],
    )
    reasons = _benchmark_inconclusive_reasons([dr])
    assert reasons, "a 5/5 crashed run must be INCONCLUSIVE, not a pass"
    assert any("liveness" in r for r in reasons)
    assert any("productivity" in r for r in reasons)
    assert any("completeness" in r for r in reasons)


def test_zero_poses_without_crash_is_inconclusive():
    # Binary exited 0 but produced nothing across all attempted targets.
    dr = _dr("astex_diverse", completed=["a", "b"], crashes=0, total_poses=0)
    reasons = _benchmark_inconclusive_reasons([dr])
    assert any("productivity" in r for r in reasons)


def test_unmeasured_metric_is_inconclusive():
    dr = _dr("astex_diverse", completed=["a"], total_poses=3,
             missing_metrics=["docking_power_top1"])
    reasons = _benchmark_inconclusive_reasons([dr])
    assert reasons and all("completeness" in r for r in reasons)


def test_healthy_run_passes_gates_1_3():
    dr = _dr("astex_diverse", completed=["a", "b"], total_poses=42)
    assert _benchmark_inconclusive_reasons([dr]) == []


def test_no_attempts_is_not_flagged_productivity():
    # An empty dataset (nothing attempted) is not a productivity failure.
    dr = _dr("empty", completed=[], failed=[], total_poses=0)
    assert _benchmark_inconclusive_reasons([dr]) == []


def test_allowlist_bypasses_productivity(monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_BENCH_ALLOW_EMPTY", "hard_dataset,other")
    dr = _dr("hard_dataset", completed=["a"], total_poses=0)
    assert _benchmark_inconclusive_reasons([dr]) == []
    # ...but a crash is never allowlistable.
    dr2 = _dr("hard_dataset", failed=["a"], crashes=1, total_poses=0)
    assert any("liveness" in r for r in _benchmark_inconclusive_reasons([dr2]))
