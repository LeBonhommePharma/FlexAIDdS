"""Unit tests for new DatasetRunner / EntryTaskManager features (per-entry, cost-aware, manager)."""

import json
import tempfile
from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import (
    DRY_RUN_METRICS_NOTE,
    DatasetConfig,
    DatasetRunner,
    EntryTaskManager,
    _filter_requested_metrics_for_dry_run,
    _is_docking_power_metric,
    _strip_docking_power_metrics,
)


def test_entry_manager_basic_local():
    work = [("t1", "holo"), ("t2", "holo")]
    mgr = EntryTaskManager(work, n_workers=2)

    def fake(item):
        return (*item, [], 0.1, "")

    res = mgr.run(fake)
    assert len(res) == 2


def test_entry_manager_cost_hints_sorting():
    work = [("expensive", "holo"), ("cheap", "holo")]
    hints = {"cheap_holo": 1.0, "expensive_holo": 100.0}
    mgr = EntryTaskManager(work, cost_hints=hints)
    # Cheapest should be first after __init__ sorting
    assert mgr.work_items[0] == ("cheap", "holo")


def test_entry_manager_load_cost_hints_from_manifest():
    with tempfile.TemporaryDirectory() as tmp:
        man = Path(tmp) / "_entry_manifest.json"
        man.write_text(json.dumps({
            "timings": {
                "per_entry_cost_cpu_seconds": {"1a30_holo": 4.2}
            }
        }))
        hints = EntryTaskManager.load_cost_hints_from_manifest(man)
        assert hints["1a30_holo"] == 4.2


def test_entry_manager_hybrid_pool():
    mgr = EntryTaskManager([("t", "holo")], n_workers=2)

    def fake(item):
        return (*item, [], 0.05, "")

    res = mgr.run(fake)
    assert len(res) == 1


# ---------------------------------------------------------------------------
# Dry-run must not sell synthetic docking_power as real docking rates
# ---------------------------------------------------------------------------


def test_strip_docking_power_helpers():
    assert _is_docking_power_metric("docking_power_top1")
    assert _is_docking_power_metric("docking_power_top3")
    assert not _is_docking_power_metric("scoring_power_pearson_r")
    assert not _is_docking_power_metric("entropy_rescue_rate")

    raw = {
        "docking_power_top1": 0.99,
        "docking_power_top3": 1.0,
        "entropy_rescue_rate": 0.4,
        "scoring_power_pearson_r": 0.5,
    }
    cleaned = _strip_docking_power_metrics(raw)
    assert "docking_power_top1" not in cleaned
    assert "docking_power_top3" not in cleaned
    assert cleaned["entropy_rescue_rate"] == 0.4
    assert cleaned["scoring_power_pearson_r"] == 0.5

    assert _filter_requested_metrics_for_dry_run(None) is None
    assert _filter_requested_metrics_for_dry_run(
        ["docking_power_top1", "entropy_rescue_rate"]
    ) == ["entropy_rescue_rate"]


def test_dry_run_omits_docking_power_metrics(tmp_path):
    """Dry-run synthetic poses must never report docking_power_* as real rates.

    Does not hardcode fake success rates as expected production values —
    only asserts that docking_power keys are absent and dry_run provenance
    is explicit.
    """
    cfg = DatasetConfig(
        slug="tiny_fake_dry_run",
        name="Tiny Fake Dry-Run",
        description="Minimal synthetic config for dry-run metric gating",
        targets=["tA", "tB", "tC"],
        tier1_subset_size=3,
        metrics=[
            "docking_power_top1",
            "docking_power_top3",
            "scoring_power_pearson_r",
            "scoring_power_rmse",
            "entropy_rescue_rate",
        ],
        # Baselines exist so a buggy dry-run path could incorrectly treat
        # synthetic docking_power as a production regression signal.
        expected_baselines={
            "docking_power_top1": 0.70,
            "docking_power_top3": 0.85,
        },
    )

    runner = DatasetRunner(results_dir=tmp_path, dry_run=True, n_workers=1)
    dr = runner.run_dataset(cfg, tier=1)

    assert dr.dry_run is True
    assert dr.metrics_note
    assert "synthetic" in dr.metrics_note.lower()
    assert "docking_power" in dr.metrics_note.lower()
    assert dr.metrics_note == DRY_RUN_METRICS_NOTE

    docking_keys = [k for k in dr.metrics if k.startswith("docking_power")]
    assert docking_keys == [], (
        f"dry-run must not report docking_power_* as real rates; got {docking_keys} "
        f"with values {[dr.metrics[k] for k in docking_keys]}"
    )

    # Regression checks are skipped in dry-run — synthetic rates must not
    # be compared to production baselines.
    assert dr.regression_flags == {}

    payload = dr.to_dict()
    assert payload["dry_run"] is True
    assert "metrics_note" in payload
    assert "docking_power_top1" not in payload["metrics"]
    assert "docking_power_top3" not in payload["metrics"]
    for key in payload["metrics"]:
        assert not key.startswith("docking_power_"), key


def test_dry_run_default_all_metrics_still_strips_docking_power(tmp_path):
    """When metrics list is empty (compute-all), dry-run still strips docking_power_*."""
    cfg = DatasetConfig(
        slug="tiny_all_metrics_dry_run",
        name="Tiny All Metrics Dry-Run",
        description="Empty metrics list → compute_all_metrics default suite",
        targets=["x1", "x2"],
        tier1_subset_size=2,
        metrics=[],  # None path via empty → falls through to compute-all
    )
    runner = DatasetRunner(results_dir=tmp_path, dry_run=True, n_workers=1)
    # Explicit metric_subset=None so runner uses compute-all defaults
    dr = runner.run_dataset(cfg, tier=1, metric_subset=None)

    assert dr.dry_run is True
    assert not any(k.startswith("docking_power_") for k in dr.metrics)
    assert not any(k.startswith("docking_power_") for k in dr.ci_95)


# ── grand_xi serialization regression (PR #311, Codex review) ────────────────
#
# tr.grand_xi = g.log_Xi (bare, no parens) assigned a BOUND METHOD, which then
# failed to JSON-serialize in _save_target_result. That path is wrapped in
# `except Exception: logger.debug`, so it failed *quietly* — no red CI, no
# raised error, just a silently missing number. The full suite did not catch
# it; Grok found it by reading. These tests lock that path.

def test_grand_xi_serializes_as_number_through_save_path(tmp_path):
    """The exact quiet-failure path: grand_xi must reach JSON as a number."""
    import json
    from flexaidds.grand_canonical import compute_grand_partition

    g = compute_grand_partition([("lig1", -12.5, 1e-6)], temperature_K=298.0)

    # Computed exactly as dataset_runner does at the production call site.
    grand_xi = g.log_Xi()

    payload = {"target_id": "lig1", "grand_xi": grand_xi}
    text = json.dumps(payload)  # would raise TypeError on a bound method
    assert isinstance(json.loads(text)["grand_xi"], (int, float))


def test_log_Xi_is_a_method_not_a_property():
    """Pins the upstream contract the bug depended on.

    If GrandPartitionFunction.log_Xi ever becomes a @property, bare-attribute
    access starts being correct and this test flips — which is the signal to
    revisit every `log_Xi()` call site rather than discovering it in a receipt.
    """
    from flexaidds.grand_canonical import compute_grand_partition

    g = compute_grand_partition([("lig1", -12.5, 1e-6)], temperature_K=298.0)
    attr = type(g).__dict__.get("log_Xi")
    assert not isinstance(attr, property), "log_Xi became a property — audit all call sites"
    assert callable(g.log_Xi)
    assert isinstance(g.log_Xi(), float)
