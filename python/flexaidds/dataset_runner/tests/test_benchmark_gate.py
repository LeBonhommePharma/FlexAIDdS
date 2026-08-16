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
        exit_codes=None, missing_metrics=(), newly_executed=None, resumed=0,
        timeouts=None):
    """Build a minimal DatasetResult-shaped object for the gate helper.

    ``newly_executed`` defaults to attempted (a fresh run); pass 0 to model a
    run that ran nothing this session. ``resumed`` is the count of targets
    loaded from --resume checkpoints — the field that separates a legitimate
    resume from a run that scheduled nothing at all. ``timeouts`` mirrors
    DatasetResult.entry_timeouts (entries killed at the wall-time cap).
    """
    completed = list(completed)
    failed = list(failed)
    if newly_executed is None:
        newly_executed = len(completed) + len(failed)
    return types.SimpleNamespace(
        config=types.SimpleNamespace(slug=slug),
        targets_completed=completed,
        targets_failed=failed,
        flexaid_crashes=crashes,
        total_poses=total_poses,
        entry_exit_codes=exit_codes or {},
        entry_timeouts=timeouts or {},
        inconclusive_metrics=list(missing_metrics),
        newly_executed=newly_executed,
        resumed=resumed,
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


def test_fully_resumed_run_is_not_inconclusive():
    # Regression for the review finding (Bumble + Honey): a pure --resume run
    # executes no fresh targets — its poses come from checkpoints, not all_poses
    # — so 0 poses + unmeasured metrics must NOT read as a dead binary.
    dr = _dr(
        "astex_diverse",
        completed=["a", "b", "c"],   # all seeded from already_completed
        newly_executed=0,            # nothing dispatched this run
        resumed=3,                   # ...but three came from checkpoints
        total_poses=0,               # checkpoint poses not reloaded
        missing_metrics=["docking_power_top1", "mean_rmsd"],
    )
    assert _benchmark_inconclusive_reasons([dr]) == []


def test_run_that_scheduled_nothing_is_inconclusive():
    # Bumble's re-review finding: newly_executed==0 AND resumed==0 means the run
    # scheduled no work at all (empty tier / bad filter). Zero output there IS
    # meaningful — it must not pass silently, the same inversion one size down.
    dr = _dr(
        "astex_diverse",
        completed=[],
        newly_executed=0,
        resumed=0,
        total_poses=0,
        missing_metrics=["docking_power_top1", "mean_rmsd"],
    )
    reasons = _benchmark_inconclusive_reasons([dr])
    assert reasons and any("completeness" in r for r in reasons)


def test_partial_resume_still_judges_fresh_work():
    # Some targets resumed, some freshly run and crashed → still INCONCLUSIVE.
    dr = _dr(
        "astex_diverse",
        completed=["resumed_a", "resumed_b"],
        failed=["fresh_c"],
        newly_executed=1,            # one fresh target ran
        crashes=1,
        total_poses=0,
        exit_codes={"fresh_c/lig": -6},
    )
    reasons = _benchmark_inconclusive_reasons([dr])
    assert any("liveness" in r for r in reasons)


def test_allowlist_bypasses_productivity(monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_BENCH_ALLOW_EMPTY", "hard_dataset,other")
    dr = _dr("hard_dataset", completed=["a"], total_poses=0)
    assert _benchmark_inconclusive_reasons([dr]) == []
    # ...but a crash is never allowlistable.
    dr2 = _dr("hard_dataset", failed=["a"], crashes=1, total_poses=0)
    assert any("liveness" in r for r in _benchmark_inconclusive_reasons([dr2]))


# --- WO-4: the liveness reason names the failure class ---------------------
# A timeout and a crash are the SAME science verdict (both count as crashes —
# that must not change), but different triage: a crash implicates the PR, a
# timeout at the cap usually implicates the job budget. The reason string must
# let a reader tell them apart without pulling the log blob.


def test_timeout_reason_names_cap():
    dr = _dr(
        "astex_diverse",
        failed=["2bys"],
        crashes=1,
        total_poses=0,
        exit_codes={"2bys/lig": None},
        timeouts={"2bys/lig": 3600},
    )
    (reason,) = [r for r in _benchmark_inconclusive_reasons([dr]) if "liveness" in r]
    assert "TIMED OUT at 3600 s" in reason
    assert "2bys/lig" in reason


def test_crash_reason_names_exit_code():
    dr = _dr(
        "astex_diverse",
        failed=["1gpk"],
        crashes=1,
        total_poses=0,
        exit_codes={"1gpk/lig": 134},
    )
    (reason,) = [r for r in _benchmark_inconclusive_reasons([dr]) if "liveness" in r]
    assert "CRASHED with code 134" in reason


def test_exec_failure_reason_when_no_exit_code_and_no_timeout():
    # None exit code WITHOUT a timeout record = the binary never executed.
    dr = _dr(
        "astex_diverse",
        failed=["1mq6"],
        crashes=1,
        total_poses=0,
        exit_codes={"1mq6/lig": None},
        timeouts={},
    )
    (reason,) = [r for r in _benchmark_inconclusive_reasons([dr]) if "liveness" in r]
    assert "EXEC FAILURE" in reason


def test_timeout_still_counts_as_crash_for_the_verdict():
    # The class distinction is triage-only: a timed-out entry keeps the run
    # INCONCLUSIVE exactly like a crash. Exit 3 is the gate working.
    dr = _dr(
        "astex_diverse",
        failed=["2bys"],
        crashes=1,
        exit_codes={"2bys/lig": None},
        timeouts={"2bys/lig": 5400},
    )
    reasons = _benchmark_inconclusive_reasons([dr])
    assert any("liveness" in r for r in reasons)


def test_reason_survives_results_without_timeout_field():
    # Back-compat: gate objects predating entry_timeouts (SimpleNamespace
    # mocks, old pickles) must still classify via getattr default.
    dr = _dr("astex_diverse", failed=["a"], crashes=1, exit_codes={"a/lig": -6})
    del dr.entry_timeouts
    (reason,) = [r for r in _benchmark_inconclusive_reasons([dr]) if "liveness" in r]
    assert "CRASHED with code -6" in reason
