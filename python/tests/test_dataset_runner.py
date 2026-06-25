"""Unit tests for new DatasetRunner / EntryTaskManager features (per-entry, cost-aware, manager)."""

import json
import tempfile
from pathlib import Path

import pytest

from flexaidds.dataset_runner.runner import (
    EntryTaskManager,
    DatasetRunner,
    DatasetConfig,
    load_entry_manifest,
    completed_targets_from_manifest,
    plan_runtime,
    KNOWN_LARGE_DATASETS,
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


def test_load_entry_manifest_missing():
    assert load_entry_manifest("/nonexistent/_entry_manifest.json") is None


def test_completed_targets_from_manifest_per_entry_status():
    targets = [f"t{i:04d}" for i in range(1200)]
    status = {f"{tid}_holo": {"success": True, "duration_seconds": 10.0} for tid in targets[:800]}
    manifest = {"per_entry_status": status, "completed": targets[:800]}
    done = completed_targets_from_manifest(manifest, targets, ["holo"])
    assert done is not None
    assert len(done) == 800


def test_completed_targets_from_manifest_legacy_single_state():
    manifest = {
        "completed": ["1a30", "1a31"],
        "timings": {"per_entry_wall_seconds": {"1a30_holo": 12.0}},
    }
    done = completed_targets_from_manifest(manifest, ["1a30", "1a31", "1a32"], ["holo"])
    assert done == {"1a30", "1a31"}


def test_discover_completed_manifest_fast_path(tmp_path):
    config = DatasetConfig(
        slug="astex_nonnative",
        name="test",
        description="",
        targets=[f"t{i}" for i in range(50)],
        structural_states=["holo"],
    )
    tier_dir = tmp_path / "astex_nonnative" / "tier1"
    tier_dir.mkdir(parents=True)
    status = {f"t{i}_holo": {"success": True, "duration_seconds": 5.0} for i in range(30)}
    manifest = {
        "completed": [f"t{i}" for i in range(30)],
        "per_entry_status": status,
        "timings": {"summary": {"mean_entry_seconds": 5.0}},
    }
    (tier_dir / "_entry_manifest.json").write_text(json.dumps(manifest))

    runner = DatasetRunner(results_dir=tmp_path, dry_run=True, resume=True)
    done = runner._discover_completed_targets(config, tier=1)
    assert len(done) == 30


def test_plan_runtime_writes_estimates(tmp_path):
    tier_dir = tmp_path / "astex_nonnative" / "tier2"
    tier_dir.mkdir(parents=True)
    (tier_dir / "_entry_manifest.json").write_text(json.dumps({
        "timings": {
            "summary": {
                "mean_entry_seconds": 120.0,
                "median_entry_seconds": 100.0,
            }
        }
    }))
    out = tmp_path / "runtime_plan.txt"
    plan = plan_runtime(
        results_dir=tmp_path,
        workers=4,
        omp_threads=2,
        output_path=out,
    )
    assert out.is_file()
    text = out.read_text()
    assert "astex_nonnative" in text
    assert "posex" in text
    assert plan["datasets"]["astex_nonnative"]["estimated_wall_seconds_mean"] > 0
    assert KNOWN_LARGE_DATASETS["astex_nonnative"] == 1113
