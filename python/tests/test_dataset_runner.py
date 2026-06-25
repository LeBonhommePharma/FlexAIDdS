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
    load_timing_priors,
    load_large_dataset_catalog,
    sanitize_entry_manifest,
    collect_per_entry_fields,
    _iter_entry_result_jsons,
    KNOWN_LARGE_DATASETS,
)
from flexaidds.dataset_runner.launch_queue import (
    build_run_status,
    build_launch_plan,
    LARGE_N_ENTRIES,
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
    targets = [f"t{i}" for i in range(10)]
    config = DatasetConfig(
        slug="astex_nonnative",
        name="test",
        description="",
        targets=targets,
        tier1_subset_size=5,
        structural_states=["holo"],
    )
    tier_dir = tmp_path / "astex_nonnative" / "tier1"
    tier_dir.mkdir(parents=True)
    scheduled = config.scheduled_targets(1)
    status = {f"{tid}_holo": {"success": True, "duration_seconds": 5.0} for tid in scheduled}
    manifest = {
        "completed": list(scheduled),
        "per_entry_status": status,
        "timings": {"summary": {"mean_entry_seconds": 5.0}},
    }
    (tier_dir / "_entry_manifest.json").write_text(json.dumps(manifest))

    runner = DatasetRunner(results_dir=tmp_path, dry_run=True, resume=True)
    done = runner._discover_completed_targets(config, tier=1, scheduled_targets=scheduled)
    assert len(done) == len(scheduled)


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
    assert plan["datasets"]["astex_nonnative"]["timing_source"] != "default_prior"


def test_iter_entry_result_jsons_excludes_dotfiles(tmp_path):
    d = tmp_path / "tier1"
    d.mkdir()
    (d / "ACE_holo.json").write_text('{"target_id":"ACE","structural_state":"holo"}')
    (d / ".cost_history.json").write_text('{"ACE_holo": 1.0}')
    (d / "_entry_manifest.json").write_text("{}")
    names = {p.name for p in _iter_entry_result_jsons(d)}
    assert names == {"ACE_holo.json"}


def test_collect_per_entry_fields_ignores_cost_history(tmp_path):
    d = tmp_path / "tier1"
    d.mkdir()
    (d / "ACE_holo.json").write_text(json.dumps({
        "target_id": "ACE",
        "structural_state": "holo",
        "duration_seconds": 4.2,
        "success": True,
        "poses": [],
    }))
    (d / ".cost_history.json").write_text('{"ACE_holo": 1.0}')
    (d / "_entry_manifest.json").write_text("{}")
    status, wall, cost, durations = collect_per_entry_fields(d, omp_threads=2)
    assert set(status.keys()) == {"ACE_holo"}
    assert "None_None" not in status
    assert wall["ACE_holo"] == 4.2
    assert durations == [4.2]


def test_write_entry_manifest_no_none_none_on_disk(tmp_path):
    tier_dir = tmp_path / "astex_nonnative" / "tier1"
    tier_dir.mkdir(parents=True)
    (tier_dir / "ACE_holo.json").write_text(json.dumps({
        "target_id": "ACE",
        "structural_state": "holo",
        "duration_seconds": 5.0,
        "success": True,
        "poses": [],
    }))
    (tier_dir / ".cost_history.json").write_text('{"ACE_holo": 1.0}')
    config = DatasetConfig(slug="astex_nonnative", name="t", description="", targets=["ACE"])
    runner = DatasetRunner(results_dir=tmp_path, dry_run=True)
    runner._write_entry_manifest(config, tier=1, completed=["ACE"], failed=[])
    raw = json.loads((tier_dir / "_entry_manifest.json").read_text())
    assert "None_None" not in raw.get("per_entry_status", {})
    assert "ACE_holo" in raw["per_entry_status"]
    assert len(raw["per_entry_status"]) == 1


def test_effective_entry_count_tier2_large_dataset():
    cfg = DatasetConfig(
        slug="astex_nonnative",
        name="x",
        description="",
        targets=["ACE"],
        structural_states=["holo", "apo", "alternative"],
    )
    assert cfg.effective_entry_count(2) == 1113
    assert cfg.effective_entry_count(1) == 3


def test_scheduled_work_items_tier2_astex_nonnative():
    cfg = DatasetConfig(slug="astex_nonnative", name="x", description="", targets=["ACE"])
    items = cfg.scheduled_work_items(tier=2)
    assert len(items) == 1113
    assert items[0][1] == "crossdock"


def test_scheduled_work_items_tier2_posex_cd():
    cfg = DatasetConfig(slug="posex_cd", name="x", description="", targets=["7FVX_K7C"])
    items = cfg.scheduled_work_items(tier=2)
    assert len(items) == 1312


def test_sanitize_manifest_removes_none_none():
    dirty = {
        "per_entry_status": {"None_None": {"success": False}, "ACE_holo": {"success": True}},
        "timings": {"per_entry_wall_seconds": {"None_None": 0.0, "ACE_holo": 1.0}},
    }
    clean = sanitize_entry_manifest(dirty)
    assert "None_None" not in clean["per_entry_status"]
    assert "ACE_holo" in clean["per_entry_status"]


def test_benchmark_report_no_utcnow_deprecation():
    import warnings
    from flexaidds.dataset_runner.runner import BenchmarkReport, _runner_info, _git_sha
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        BenchmarkReport(
            datasets=[],
            generated_at=__import__("flexaidds.dataset_runner.runner", fromlist=["_utc_now_iso"])._utc_now_iso(),
            git_sha=_git_sha(),
            host="test",
            runner_info=_runner_info(),
        )


def test_plan_runtime_uses_timing_priors_not_default():
    priors = load_timing_priors()
    assert "astex_nonnative" in priors
    plan = plan_runtime(results_dir="/nonexistent", workers=4)
    src = plan["datasets"]["astex_nonnative"]["timing_source"]
    assert src != "default_prior"
    assert plan["datasets"]["astex_nonnative"]["mean_entry_seconds"] > 100


def test_build_run_status_waiting_has_full_n_launch_plan():
    status = build_run_status(214)
    assert status["status"] == "waiting_for_astex_diverse_siblings"
    assert status["launch_plan"]["astex_nonnative"]["n_entries"] == LARGE_N_ENTRIES["astex_nonnative"]
    assert status["launch_plan"]["posex_cd"]["n_pairs"] == LARGE_N_ENTRIES["posex_cd"]


def test_build_run_status_launched_full_has_commands():
    status = build_run_status(0)
    assert status["status"] == "launched_full"
    plan = build_launch_plan()
    assert plan["astex_nonnative"]["command"][2] == "astex_nonnative"
    assert plan["posex_cd"]["command"][0] == "benchmark_datasets"
    assert status["posex_cd"]["n_pairs"] == 1312
