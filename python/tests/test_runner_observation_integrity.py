"""Runner-level regressions for fixed N and scientifically equivalent resume.

All docking is replaced with scalar synthetic observations. These tests never
run the engine, PoseBusters or tENCoM and make no molecular accuracy claims.
"""

import json

import pytest

from flexaidds.dataset_runner.cli import _benchmark_inconclusive_reasons
from flexaidds.dataset_runner.metrics import PoseScore
from flexaidds.dataset_runner.runner import BenchmarkReport, DatasetConfig, DatasetRunner


def _config():
    return DatasetConfig(
        slug="integrity", name="Synthetic integrity fixture", description="",
        targets=["a", "b"], metrics=["docking_power_top1"],
        expected_baselines={"docking_power_top1": 0.75}, baseline_tolerance=0.0,
    )


def _runner(tmp_path, **kwargs):
    return DatasetRunner(
        results_dir=tmp_path, binary="not-executed-test-engine",
        n_workers=1, omp_threads=1, **kwargs,
    )


def _pose(target, rmsd, state="holo"):
    return PoseScore(target, "lig", 1, rmsd, -10.0, 0.0, -10.0, False,
                     structural_state=state)


def _stub(runner, rmsds, calls):
    def dock(target, *args, structural_state="holo", **kwargs):
        calls.append(target)
        value = rmsds[target]
        return [] if value is None else [_pose(target, value, structural_state)]
    runner._dock_target = dock


def _verdict(result):
    return result.metrics, result.regression_flags, _benchmark_inconclusive_reasons([result])


def test_one_success_one_empty_keeps_fixed_denominator(tmp_path):
    runner = _runner(tmp_path)
    _stub(runner, {"a": 1.0, "b": None}, [])
    result = runner.run_dataset(_config())
    assert result.metrics["docking_power_top1"] == 0.5
    assert result.regression_flags["docking_power_top1"] is True
    assert result.targets_completed == ["a"]
    assert result.targets_failed == ["b"]
    assert result.flexaid_crashes == 0  # absent output is not an engine crash
    assert result.to_dict()["docking_power_denominator"] == 2
    assert result.observation_roster == [("a", "holo"), ("b", "holo")]
    assert result.observations_failed == [("b", "holo")]


def test_missing_input_is_an_observation_failure(tmp_path, monkeypatch):
    cfg = _config()
    cfg.data_dir = tmp_path / "inputs"
    runner = _runner(tmp_path / "results")
    monkeypatch.setattr(runner, "_resolve_entry_paths", lambda c, t, s: (
        (tmp_path / "receptor.pdb", [tmp_path / "lig.mol2"]) if t == "a" else (None, [])))
    calls = []
    _stub(runner, {"a": 1.0}, calls)
    result = runner.run_dataset(cfg)
    assert calls == ["a"]
    assert result.metrics["docking_power_top1"] == 0.5
    assert result.targets_failed == ["b"]
    assert result.flexaid_crashes == 0


@pytest.mark.parametrize("rmsds", [{"a": 5.0, "b": 5.0}, {"a": 1.0, "b": None}, {"a": None, "b": None}, {"a": 1.0, "b": 1.0}])
def test_full_resume_reconstructs_exact_quality_verdict(tmp_path, rmsds):
    runner = _runner(tmp_path)
    calls = []
    _stub(runner, rmsds, calls)
    fresh = runner.run_dataset(_config())
    runner.resume = True
    calls.clear()
    resumed = runner.run_dataset(_config())
    assert calls == []
    assert resumed.newly_executed == 0
    assert resumed.resumed == 2
    assert resumed.total_poses == fresh.total_poses
    assert resumed.targets_failed == fresh.targets_failed
    assert _verdict(resumed) == _verdict(fresh)


def test_partial_resume_includes_cached_miss(tmp_path):
    runner = _runner(tmp_path)
    _stub(runner, {"a": 5.0, "b": 1.0}, [])
    fresh = runner.run_dataset(_config())
    runner._target_result_path(_config(), 2, "b", "holo").unlink()
    runner.resume = True
    calls = []
    _stub(runner, {"b": 1.0}, calls)
    resumed = runner.run_dataset(_config())
    assert calls == ["b"]
    assert resumed.resumed == 1
    assert resumed.newly_executed == 1
    assert resumed.metrics["docking_power_top1"] == 0.5
    assert _verdict(resumed) == _verdict(fresh)


def test_resume_preserves_execution_and_termination_evidence(tmp_path):
    runner = _runner(tmp_path)
    def dock(target, *args, **kwargs):
        if target == "a":
            runner._entry_exit_codes["a/lig"] = -6
            runner._flexaid_crashes += 1
            runner._entry_early_termination["a/holo"] = {"terminated_early": True, "reason": "test", "generation": 10}
            return []
        return [_pose(target, 1.0)]
    runner._dock_target = dock
    fresh = runner.run_dataset(_config())
    runner.resume = True
    resumed = runner.run_dataset(_config())
    assert resumed.entry_exit_codes == fresh.entry_exit_codes == {"a/lig": -6}
    assert resumed.early_terminations == fresh.early_terminations
    assert resumed.flexaid_crashes == fresh.flexaid_crashes == 1
    assert _verdict(resumed) == _verdict(fresh)


@pytest.mark.parametrize("mutation", ["old_schema", "wrong_target", "missing_pose_field", "nonfinite", "wrong_score", "wrong_protocol"])
def test_incompatible_cache_fails_before_any_work(tmp_path, mutation):
    runner = _runner(tmp_path)
    _stub(runner, {"a": 1.0, "b": 1.0}, [])
    runner.run_dataset(_config())
    path = runner._target_result_path(_config(), 2, "a", "holo")
    data = json.loads(path.read_text())
    if mutation == "old_schema":
        data.pop("schema_version")
    elif mutation == "wrong_target":
        data["poses"][0]["target_id"] = "other"
    elif mutation == "missing_pose_field":
        data["poses"][0].pop("entropy")
    elif mutation == "nonfinite":
        data["poses"][0]["rmsd"] = float("nan")
    elif mutation == "wrong_score":
        data["poses"][0]["total_score"] = -100.0
    else:
        data["provenance"]["protocol"]["temperature"] = 500.0
    path.write_text(json.dumps(data))
    calls = []
    runner.resume = True
    _stub(runner, {}, calls)
    with pytest.raises(ValueError, match="Refusing unverified resume"):
        runner.run_dataset(_config())
    assert calls == []


def test_changed_ranking_objective_rejects_cache(tmp_path):
    runner = _runner(tmp_path)
    _stub(runner, {"a": 1.0, "b": 1.0}, [])
    runner.run_dataset(_config())
    legacy = _runner(tmp_path, resume=True, ranking_objective="legacy_cf_minus_s")
    with pytest.raises(ValueError, match="protocol differs"):
        legacy.run_dataset(_config())


def test_mixed_state_pool_refused_and_explicit_single_state_honored(tmp_path):
    cfg = _config()
    cfg.structural_states = ["holo", "apo"]
    runner = _runner(tmp_path)
    calls = []
    _stub(runner, {"a": 1.0, "b": None}, calls)
    with pytest.raises(ValueError, match="Mixed structural states"):
        runner.run_dataset(cfg)
    assert calls == []
    result = runner.run_dataset(cfg, structural_states=["apo"])
    assert result.observation_roster == [("a", "apo"), ("b", "apo")]
    assert result.metrics["docking_power_top1"] == 0.5


def test_run_all_partitions_states_with_fixed_denominators_and_resume(tmp_path, monkeypatch):
    cfg = _config()
    cfg.structural_states = ["holo", "apo"]
    runner = _runner(tmp_path)
    monkeypatch.setattr(runner, "discover_datasets", lambda: [cfg])
    calls = []

    def dock(target, *args, structural_state="holo", **kwargs):
        calls.append((target, structural_state))
        if target == "b":
            return []
        return [_pose(target, 1.0 if structural_state == "holo" else 5.0, structural_state)]

    runner._dock_target = dock
    fresh = runner.run_all()
    assert runner.results_dir == tmp_path
    assert cfg.structural_states == ["holo", "apo"]
    assert [dr.config.structural_states for dr in fresh.datasets] == [["holo"], ["apo"]]
    assert [dr.metrics["docking_power_top1"] for dr in fresh.datasets] == [0.5, 0.0]
    for dr in fresh.datasets:
        state = dr.config.structural_states[0]
        assert dr.targets_attempted == ["a", "b"]
        assert dr.observation_roster == [("a", state), ("b", state)]
        assert dr.to_dict()["structural_states"] == [state]
        output = tmp_path / "states" / state
        assert (output / "integrity_tier2.json").is_file()
        assert (output / "integrity/tier2" / f"a_{state}.json").is_file()
        manifest = BenchmarkReport._load_entry_manifest_summary(dr)
        assert manifest["observation_roster"] == [["a", state], ["b", state]]

    calls.clear()
    runner.resume = True
    resumed = runner.run_all()
    assert calls == []
    assert [dr.resumed for dr in resumed.datasets] == [2, 2]
    assert [_verdict(dr) for dr in resumed.datasets] == [_verdict(dr) for dr in fresh.datasets]
    saved = tmp_path / "report.json"
    saved.write_text(fresh.to_json())
    loaded = BenchmarkReport.load(saved)
    assert [dr.config.structural_states for dr in loaded.datasets] == [["holo"], ["apo"]]
    assert [dr.output_dir for dr in loaded.datasets] == [dr.output_dir for dr in fresh.datasets]


def test_run_all_single_state_preserves_existing_namespace(tmp_path, monkeypatch):
    cfg = _config()
    runner = _runner(tmp_path)
    monkeypatch.setattr(runner, "discover_datasets", lambda: [cfg])
    _stub(runner, {"a": 1.0, "b": None}, [])
    report = runner.run_all()
    assert len(report.datasets) == 1
    assert report.datasets[0].config is cfg
    assert (tmp_path / "integrity_tier2.json").is_file()
    assert (tmp_path / "integrity/tier2/a_holo.json").is_file()
    assert not (tmp_path / "states").exists()


def test_run_all_restores_namespace_on_state_failure(tmp_path, monkeypatch):
    cfg = _config()
    cfg.structural_states = ["holo", "apo"]
    runner = _runner(tmp_path)
    monkeypatch.setattr(runner, "discover_datasets", lambda: [cfg])

    def fail(*args, **kwargs):
        raise ValueError("invalid checkpoint")

    monkeypatch.setattr(runner, "run_dataset", fail)
    with pytest.raises(ValueError, match="invalid checkpoint"):
        runner.run_all()
    assert runner.results_dir == tmp_path


def test_crossdock_catalog_respects_declared_and_explicit_states(tmp_path, monkeypatch):
    catalog = [
        {"entry_id": "pair", "family": "ACE", "state": "crossdock"},
        {"entry_id": "pair", "family": "ACE", "state": "apo"},
    ]
    monkeypatch.setattr("flexaidds.dataset_runner.runner.load_large_dataset_catalog", lambda slug: catalog)
    cfg = DatasetConfig(slug="astex_nonnative", name="catalog fixture", description="",
                        targets=["ACE"], structural_states=["crossdock"],
                        metrics=["docking_power_top1"])
    assert cfg.scheduled_work_items(2) == [("pair", "crossdock")]
    assert cfg.scheduled_targets(2) == ["pair"]
    cfg.structural_states = ["crossdock", "apo"]
    runner = _runner(tmp_path)
    _stub(runner, {"pair": 1.0}, [])
    result = runner.run_dataset(cfg, structural_states=["apo"])
    assert result.observation_roster == [("pair", "apo")]
    assert result.metrics["docking_power_top1"] == 1.0


def test_bootstrap_samples_full_roster_including_empty_target(tmp_path, monkeypatch):
    runner = _runner(tmp_path, bootstrap_ci=True, n_bootstrap=20)
    samples = []
    def bootstrap(fn, data, **kwargs):
        samples.append(data)
        assert len(data) == 2
        assert sorted(len(bucket) for bucket in data) == [0, 1]
        value = fn(data)
        assert value == 0.5
        # Duplicate draws must count twice, not collapse into one target.
        assert fn([data[0], data[0]]) == 1.0
        return value, value
    monkeypatch.setattr("flexaidds.dataset_runner.runner.bootstrap_ci", bootstrap)
    _stub(runner, {"a": 1.0, "b": None}, [])
    result = runner.run_dataset(_config())
    assert samples
    assert result.ci_95["docking_power_top1"] == (0.5, 0.5)


def test_frozen_85_manifest_keeps_omitted_target_in_denominator(tmp_path):
    from pathlib import Path
    manifest = Path(__file__).resolve().parents[2] / "benchmarks/protocols/astex85_target_manifest.json"
    frozen = json.loads(manifest.read_text())
    assert frozen["N"] == 85
    scheduled = [target for target in frozen["targets"] if target != "2HR7"]
    assert len(scheduled) == 84
    cfg = DatasetConfig(slug="astex_diverse", name="Scalar frozen-roster fixture", description="",
                        targets=scheduled, metrics=["docking_power_top1"])
    runner = _runner(tmp_path, expected_target_manifest=manifest)
    calls = []
    _stub(runner, {target: 1.0 for target in scheduled}, calls)
    result = runner.run_dataset(cfg)
    assert len(calls) == 84
    assert result.to_dict()["docking_power_denominator"] == 85
    assert result.metrics["docking_power_top1"] == 84 / 85
    assert result.targets_failed == ["2HR7"]
    assert result.observations_failed == [("2HR7", "holo")]
    assert result.expected_target_manifest["sha256"] == "d207666a414422b2eea9f92ebe67a72aad9f3ac87ad029ed56609870aa414326"
    assert result.expected_target_manifest["sha256_of_sorted_codes"] == "da89650afd791e27924afbe403efdae5c63400ffce3b4c52a4c42351f4102afc"
    runner.resume = True
    calls.clear()
    resumed = runner.run_dataset(cfg)
    assert calls == []
    assert _verdict(resumed) == _verdict(result)


def test_changed_baseline_reuses_evidence_but_recomputes_verdict(tmp_path):
    cfg = _config()
    runner = _runner(tmp_path)
    _stub(runner, {"a": 1.0, "b": None}, [])
    fresh = runner.run_dataset(cfg)
    assert fresh.regression_flags["docking_power_top1"] is True
    cfg.expected_baselines["docking_power_top1"] = 0.25
    runner.resume = True
    resumed = runner.run_dataset(cfg)
    assert resumed.metrics == fresh.metrics
    assert resumed.regression_flags["docking_power_top1"] is False


def test_canonical_manifest_keeps_scheduled_checkpoint_case(tmp_path, monkeypatch):
    from pathlib import Path
    manifest = Path(__file__).resolve().parents[2] / "benchmarks/protocols/astex85_target_manifest.json"
    cfg = DatasetConfig(slug="astex_diverse", name="Case-sensitive fixture", description="",
                        targets=["1g9v"], metrics=["docking_power_top1"])
    runner = _runner(tmp_path, expected_target_manifest=manifest)
    _stub(runner, {"1g9v": 1.0}, [])
    actual_path = runner._target_result_path
    def case_sensitive_path(config, tier, target_id, state):
        # This assertion catches the Linux failure even on case-insensitive APFS.
        assert target_id == "1g9v", "canonical claim label was used as checkpoint filename"
        return actual_path(config, tier, target_id, state)
    monkeypatch.setattr(runner, "_target_result_path", case_sensitive_path)
    result = runner.run_dataset(cfg)
    assert result.targets_completed == ["1G9V"]
    saved_manifest = json.loads((tmp_path / "astex_diverse/tier2/_entry_manifest.json").read_text())
    assert "1G9V_holo" in saved_manifest["timings"]["per_entry_wall_seconds"]
    assert saved_manifest["checkpoint_files"] == {"1G9V_holo": "1g9v_holo.json"}
    runner.resume = True
    resumed = runner.run_dataset(cfg)
    assert resumed.metrics == result.metrics
