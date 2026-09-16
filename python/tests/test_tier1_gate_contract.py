"""Tier-1 CI gate contract: sampling vs ranking, receipt replay, matrix pin.

These tests reconstruct the 2026-09-15 Astex CI miss (Actions run
35016255184, seed 20260816, roster 1gpk/1mq6/1xm6/2cet) without docking.
They pin:

* ranking top-1 = 0.0 is advisory and must not fail CI
* sampling / min-RMSD still hard-fail a real library collapse
* empty results still hard-fail (productivity + completeness + sampling 0)
* seed 20260816 still draws the locked roster
* the repo CF matrix is the 9dc9 pin (this PR must not retune it)
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

from flexaidds.dataset_runner.cli import _benchmark_inconclusive_reasons
from flexaidds.dataset_runner.metrics import PoseScore, compute_all_metrics
from flexaidds.dataset_runner.runner import DatasetConfig, DatasetResult
from flexaidds.dataset_runner.tier1_contract import (
    MATRIX_PIN_MD5,
    PINNED_ROSTER,
    PINNED_TIER1_SEED,
    RECEIPT_SCHEMA,
    assert_election_objective,
    assert_matrix_pin,
    assert_tier1_roster_lock,
    assert_unregistered_ranking_is_advisory,
    build_receipt,
    md5_file,
    quality_exit_code,
    split_regression_flags,
    validate_receipt,
)

REPO = Path(__file__).resolve().parents[2]


def _pose(target, rmsd, score, rank=1):
    return PoseScore(
        target_id=target,
        ligand_id=f"{target}_lig",
        pose_rank=rank,
        rmsd=rmsd,
        enthalpy_score=score,
        entropy_correction=0.0,
        total_score=score,
        is_active=False,
    )


def _astex_cfg(**kw) -> DatasetConfig:
    path = REPO / "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml"
    cfg = DatasetConfig.from_yaml(path)
    for key, val in kw.items():
        setattr(cfg, key, val)
    return cfg


def _result(cfg, metrics, poses_n=40, completed=None):
    dr = DatasetResult(config=cfg, tier=1)
    dr.metrics = dict(metrics)
    dr.total_poses = poses_n
    dr.targets_attempted = list(PINNED_ROSTER)
    dr.targets_completed = list(completed if completed is not None else PINNED_ROSTER)
    dr.targets_failed = []
    dr.flexaid_crashes = 0
    dr.newly_executed = len(dr.targets_completed)
    dr.check_regressions()
    return dr


# Rank-1 / min-RMSD from the attached 35016255184 entry JSONs, plus rank-2/3
# decoys so docking_power_top3 stays 0 (the real libraries elected three
# non-natives before any near-native).
_FAILED_RUN_POSES = [
    _pose("1gpk", 4.7138, -22108.57, 1),
    _pose("1gpk", 5.0408, -20643.61, 2),
    _pose("1gpk", 5.1751, -20224.34, 3),
    _pose("1gpk", 1.9388, -17977.91, 6),
    _pose("1mq6", 5.4031, -26092.09, 1),
    _pose("1mq6", 5.1668, -24254.26, 2),
    _pose("1mq6", 4.7592, -23201.49, 3),
    _pose("1mq6", 0.7846, -7725.55, 10),
    _pose("1xm6", 6.8063, -24690.31, 1),
    _pose("1xm6", 7.5739, -24832.66, 2),
    _pose("1xm6", 5.8168, -22982.44, 3),
    _pose("1xm6", 3.2020, -17786.61, 6),  # no ≤2 Å pose
    _pose("2cet", 3.6110, -30171.05, 1),
    _pose("2cet", 3.7267, -28185.96, 2),
    _pose("2cet", 4.9496, -25881.04, 3),
    _pose("2cet", 1.5240, -20999.57, 10),
]


def test_failed_run_sampling_passes_ranking_is_advisory():
    """The exact CI-miss phenotype must be a green hard gate + advisory ranking."""
    cfg = _astex_cfg()
    metrics = compute_all_metrics(
        _FAILED_RUN_POSES,
        requested=list(cfg.metrics),
        n_targets=4,
    )
    assert metrics["docking_power_top1"] == 0.0
    assert metrics["docking_power_top3"] == 0.0
    assert metrics["sampling_power"] == 0.75
    assert metrics["mean_rmsd"] == pytest.approx(
        (1.9388 + 0.7846 + 3.2020 + 1.5240) / 4.0
    )

    dr = _result(cfg, metrics)
    blocking, advisory = split_regression_flags(cfg, dr.regression_flags)
    assert "docking_power_top1" in advisory
    assert "docking_power_top3" in advisory
    assert "entropy_rescue_rate" in advisory
    assert "sampling_power" not in blocking
    assert "mean_rmsd" not in blocking
    assert "median_rmsd" not in blocking
    assert quality_exit_code(
        inconclusive=_benchmark_inconclusive_reasons([dr]),
        blocking=blocking,
    ) == 0


def test_sampling_collapse_still_fails_ci():
    """All far poses: ranking 0 is no longer the only signal — sampling hard-fails."""
    cfg = _astex_cfg()
    poses = [_pose(t, 8.0, -10.0, 1) for t in PINNED_ROSTER]
    metrics = compute_all_metrics(poses, requested=list(cfg.metrics), n_targets=4)
    assert metrics["sampling_power"] == 0.0
    dr = _result(cfg, metrics)
    blocking, _advisory = split_regression_flags(cfg, dr.regression_flags)
    assert blocking.get("sampling_power") is True
    assert quality_exit_code(inconclusive=[], blocking=blocking) == 1


def test_empty_results_cannot_skip_as_ranking_zero():
    """Zero poses: #326 productivity + completeness, not an advisory ranking miss."""
    cfg = _astex_cfg()
    dr = _result(cfg, {}, poses_n=0, completed=["1gpk", "1mq6", "1xm6", "2cet"])
    dr.metrics = {}
    dr.check_regressions()
    reasons = _benchmark_inconclusive_reasons([dr])
    assert any("productivity" in r for r in reasons)
    assert any("completeness" in r for r in reasons)
    blocking, _advisory = split_regression_flags(cfg, dr.regression_flags)
    assert quality_exit_code(inconclusive=reasons, blocking=blocking) == 3


def test_other_datasets_default_to_hard_ranking_gates():
    """ci_gate_class omitted ⇒ docking_power_top1 remains a blocking regression."""
    cfg = DatasetConfig(
        slug="other",
        name="n",
        description="d",
        expected_baselines={"docking_power_top1": 0.70},
        baseline_tolerance=0.05,
    )
    blocking, advisory = split_regression_flags(
        cfg, {"docking_power_top1": True}
    )
    assert blocking == {"docking_power_top1": True}
    assert advisory == {}


def test_yaml_copies_agree_on_the_contract():
    copies = [
        REPO / "benchmarks/datasets/astex_diverse.yaml",
        REPO / "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml",
    ]
    seen = []
    for path in copies:
        raw = yaml.safe_load(path.read_text())
        seen.append(
            (
                raw.get("tier1_seed"),
                tuple(raw.get("tier1_expected_roster") or ()),
                raw.get("matrix_pin_md5"),
                raw.get("ci_ranking_status"),
                raw.get("ci_election_objective"),
                dict(raw.get("ci_gate_class") or {}),
                raw.get("expected_baselines"),
            )
        )
    assert seen[0] == seen[1]
    assert seen[0][0] == PINNED_TIER1_SEED
    assert seen[0][1] == PINNED_ROSTER
    assert seen[0][2] == MATRIX_PIN_MD5
    assert seen[0][3] == "unregistered"
    assert seen[0][4] == "cf_minus_ts"
    assert seen[0][5]["docking_power_top1"] == "advisory"
    assert seen[0][5]["sampling_power"] == "hard"
    assert seen[0][6]["sampling_power"] == 0.50
    assert seen[0][6]["docking_power_top1"] == 0.70


def test_every_astex_baseline_has_an_explicit_gate_class():
    cfg = _astex_cfg()
    for metric in cfg.expected_baselines:
        assert metric in cfg.ci_gate_class, metric
    assert cfg.ci_ranking_status == "unregistered"
    assert cfg.ci_gate_class["docking_power_top1"] == "advisory"


def test_unregistered_ranking_cannot_block(monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", str(PINNED_TIER1_SEED))
    cfg = _astex_cfg()
    assert cfg.tier1_targets() == list(PINNED_ROSTER)
    assert cfg.ci_ranking_status == "unregistered"
    blocking, advisory = split_regression_flags(
        cfg,
        {
            "docking_power_top1": True,
            "sampling_power": False,
            "mean_rmsd": False,
        },
    )
    assert "docking_power_top1" in advisory
    assert blocking == {}


def test_roster_lock_refuses_a_different_draw(monkeypatch):
    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", str(PINNED_TIER1_SEED))
    cfg = _astex_cfg()
    assert_tier1_roster_lock(cfg, list(PINNED_ROSTER), tier=1)
    with pytest.raises(RuntimeError, match="roster replay"):
        assert_tier1_roster_lock(cfg, ["1gpk", "1mq6", "1n2j", "1t46"], tier=1)


def test_matrix_pin_matches_repo_file():
    """This PR must not retune MC_st0r5.2_6.dat. The 9dc9 pin is the gate."""
    matrix = REPO / "MC_st0r5.2_6.dat"
    assert matrix.is_file()
    assert md5_file(matrix) == MATRIX_PIN_MD5
    assert assert_matrix_pin(repo_root=REPO) == MATRIX_PIN_MD5


def test_receipt_roundtrip_replays_the_draw(tmp_path):
    payload = build_receipt(
        seed=PINNED_TIER1_SEED,
        roster=PINNED_ROSTER,
        matrix_md5=MATRIX_PIN_MD5,
        election_objective="cf_minus_ts",
        ci_ranking_status="unregistered",
        hard_gates={"sampling_power": {"measured": 0.75, "baseline": 0.50, "regressed": False}},
        advisory_ranking={
            "docking_power_top1": {"measured": 0.0, "baseline": 0.70, "regressed": True}
        },
        exit_code=0,
    )
    assert payload["schema"] == RECEIPT_SCHEMA
    assert payload["verdict"] == "pass"
    assert payload["drawn"] == ["1xm6", "2cet"]
    dest = tmp_path / "TIER1_RECEIPT.json"
    dest.write_text(json.dumps(payload))
    assert validate_receipt(json.loads(dest.read_text())) == []


def test_receipt_rejects_wrong_matrix_or_roster():
    payload = build_receipt(
        seed=PINNED_TIER1_SEED,
        roster=PINNED_ROSTER,
        matrix_md5=MATRIX_PIN_MD5,
    )
    bad_matrix = dict(payload, matrix_md5="72d7c7396702331d96ff12d18f831796", matrix_pin_ok=False)
    assert any("matrix_md5" in e for e in validate_receipt(bad_matrix))
    bad_roster = dict(payload, roster=["1gpk", "1mq6"])
    assert any("roster" in e for e in validate_receipt(bad_roster))


def test_dry_run_skips_matrix_pin():
    # Synthetic CI must not require a staged binary-adjacent matrix.
    assert assert_matrix_pin(binary="/no/such/flexaid", dry_run=True) == MATRIX_PIN_MD5


def test_quality_exit_advisory_is_zero():
    assert quality_exit_code(inconclusive=[], blocking={}) == 0
    assert quality_exit_code(inconclusive=[], blocking={"sampling_power": True}) == 1
    assert quality_exit_code(inconclusive=["completeness"], blocking={}) == 3
    assert quality_exit_code(inconclusive=[], blocking={}, dry_run=True) == 0


def test_invalid_gate_class_is_rejected():
    from flexaidds.dataset_runner.tier1_contract import parse_ci_gate_class

    with pytest.raises(ValueError, match="ci_gate_class"):
        parse_ci_gate_class({"docking_power_top1": "ignore"})


def test_election_objective_lock_refuses_a_different_proxy():
    cfg = _astex_cfg()
    assert_election_objective(cfg, "cf_minus_ts")
    with pytest.raises(RuntimeError, match="election objective"):
        assert_election_objective(cfg, "legacy_cf_minus_s")


def test_unregistered_ranking_cannot_be_rehardened():
    cfg = _astex_cfg()
    assert_unregistered_ranking_is_advisory(cfg)
    cfg.ci_gate_class["docking_power_top1"] = "hard"
    with pytest.raises(ValueError, match="unregistered"):
        assert_unregistered_ranking_is_advisory(cfg)


def test_from_yaml_rejects_unregistered_hard_ranking(tmp_path):
    src = REPO / "python/flexaidds/dataset_runner/datasets/astex_diverse.yaml"
    raw = src.read_text()
    poisoned = raw.replace(
        "  docking_power_top1: advisory",
        "  docking_power_top1: hard",
    )
    dest = tmp_path / "astex_diverse.yaml"
    dest.write_text(poisoned)
    with pytest.raises(ValueError, match="unregistered"):
        DatasetConfig.from_yaml(dest)


def test_missing_sampling_power_is_completeness_not_advisory_pass():
    """A measured top-1 of 0.0 must not hide an unmeasured sampling floor."""
    cfg = _astex_cfg()
    metrics = {
        "docking_power_top1": 0.0,
        "docking_power_top3": 0.0,
        "mean_rmsd": 1.80,
        "median_rmsd": 1.70,
        "entropy_rescue_rate": 0.0,
    }
    dr = _result(cfg, metrics, poses_n=40)
    assert "sampling_power" in dr.inconclusive_metrics
    reasons = _benchmark_inconclusive_reasons([dr])
    assert any("completeness" in r and "sampling_power" in r for r in reasons)
    blocking, _advisory = split_regression_flags(cfg, dr.regression_flags)
    assert quality_exit_code(inconclusive=reasons, blocking=blocking) == 3


def test_receipt_is_written_for_astex_tier1(tmp_path, monkeypatch):
    from types import SimpleNamespace

    from flexaidds.dataset_runner.cli import _maybe_write_tier1_receipt

    monkeypatch.setenv("FLEXAIDDS_TIER1_SEED", str(PINNED_TIER1_SEED))
    cfg = _astex_cfg()
    metrics = compute_all_metrics(
        _FAILED_RUN_POSES,
        requested=list(cfg.metrics),
        n_targets=4,
    )
    dr = _result(cfg, metrics)
    dr.tier = 1
    dr.ranking_objective = "cf_minus_ts"
    report = SimpleNamespace(datasets=[dr])
    runner = SimpleNamespace(binary="", dry_run=False)
    _maybe_write_tier1_receipt(
        report,
        runner,
        results_dir=str(tmp_path),
        inconclusive=[],
        blocking={},
        advisory={"astex_diverse:docking_power_top1": True},
        exit_code=0,
    )
    path = tmp_path / "TIER1_RECEIPT.json"
    payload = json.loads(path.read_text())
    assert validate_receipt(payload) == []
    assert payload["verdict"] == "pass"
    assert payload["advisory_ranking"]["docking_power_top1"]["regressed"] is True


def test_receipt_is_not_written_for_other_datasets(tmp_path):
    from types import SimpleNamespace

    from flexaidds.dataset_runner.cli import _maybe_write_tier1_receipt

    cfg = DatasetConfig(slug="casf2016", name="n", description="d")
    dr = DatasetResult(config=cfg, tier=1)
    _maybe_write_tier1_receipt(
        SimpleNamespace(datasets=[dr]),
        SimpleNamespace(binary="", dry_run=False),
        results_dir=str(tmp_path),
        inconclusive=[],
        blocking={},
        advisory={},
        exit_code=0,
    )
    assert not (tmp_path / "TIER1_RECEIPT.json").exists()
