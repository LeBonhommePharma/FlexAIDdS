"""Tests for the flexaidds.dataset_runner subpackage."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")

from flexaidds.dataset_runner import (  # noqa: E402
    DatasetRunner,
    DatasetConfig,
    DatasetResult,
    BenchmarkReport,
    PoseScore,
)
from flexaidds.dataset_runner.metrics import (  # noqa: E402
    entropy_rescue_rate,
    enrichment_factor,
    log_auc,
    scoring_power,
    docking_power,
    target_specificity_zscore,
    hit_rate_top_n,
    bootstrap_ci,
    compute_all_metrics,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DATASETS_DIR = Path(__file__).resolve().parent.parent / "flexaidds" / "dataset_runner" / "datasets"


def _make_pose(
    target_id: str = "1abc",
    ligand_id: str = "lig1",
    pose_rank: int = 1,
    rmsd: float = 1.0,
    enthalpy_score: float = -8.0,
    entropy_correction: float = 1.5,
    total_score: float = -9.5,
    is_active: bool = True,
    exp_affinity: float | None = -7.0,
    structural_state: str = "holo",
) -> PoseScore:
    return PoseScore(
        target_id=target_id,
        ligand_id=ligand_id,
        pose_rank=pose_rank,
        rmsd=rmsd,
        enthalpy_score=enthalpy_score,
        entropy_correction=entropy_correction,
        total_score=total_score,
        is_active=is_active,
        exp_affinity=exp_affinity,
        structural_state=structural_state,
    )


@pytest.fixture
def sample_poses() -> list[PoseScore]:
    """A small set of poses spanning two targets for testing metrics."""
    poses = []
    # Target A: crystal pose (low RMSD) ranked poorly by enthalpy but well by total
    poses.append(_make_pose("A", "lig1", 1, rmsd=0.8, enthalpy_score=-5.0,
                            entropy_correction=3.0, total_score=-8.0,
                            is_active=True, exp_affinity=-8.0))
    poses.append(_make_pose("A", "lig1", 2, rmsd=3.0, enthalpy_score=-7.0,
                            entropy_correction=0.5, total_score=-7.5,
                            is_active=False, exp_affinity=-6.0))
    poses.append(_make_pose("A", "lig1", 3, rmsd=4.0, enthalpy_score=-8.0,
                            entropy_correction=0.2, total_score=-8.2,
                            is_active=False, exp_affinity=-5.0))
    poses.append(_make_pose("A", "lig1", 4, rmsd=5.0, enthalpy_score=-6.0,
                            entropy_correction=0.1, total_score=-6.1,
                            is_active=True, exp_affinity=-7.5))

    # Target B: crystal pose ranked well by enthalpy (no rescue needed)
    poses.append(_make_pose("B", "lig2", 1, rmsd=0.5, enthalpy_score=-10.0,
                            entropy_correction=1.0, total_score=-11.0,
                            is_active=True, exp_affinity=-10.0))
    poses.append(_make_pose("B", "lig2", 2, rmsd=3.5, enthalpy_score=-6.0,
                            entropy_correction=0.3, total_score=-6.3,
                            is_active=False, exp_affinity=-5.5))
    return poses


# ---------------------------------------------------------------------------
# DatasetConfig tests
# ---------------------------------------------------------------------------


class TestDatasetConfig:
    def test_construction(self):
        cfg = DatasetConfig(slug="test", name="Test Dataset", description="A test")
        assert cfg.slug == "test"
        assert cfg.tier == 2
        assert cfg.tier1_subset_size == 5

    def test_tier1_targets(self):
        cfg = DatasetConfig(
            slug="test", name="Test", description="",
            targets=["a", "b", "c", "d", "e", "f"],
            tier1_subset_size=3,
        )
        assert cfg.tier1_targets() == ["a", "b", "c"]

    def test_tier1_targets_when_fewer(self):
        cfg = DatasetConfig(
            slug="test", name="Test", description="",
            targets=["a", "b"],
            tier1_subset_size=5,
        )
        assert cfg.tier1_targets() == ["a", "b"]

    @pytest.mark.skipif(not DATASETS_DIR.is_dir(), reason="datasets dir missing")
    def test_from_yaml(self):
        yaml_files = list(DATASETS_DIR.glob("*.yaml"))
        assert len(yaml_files) > 0, "No YAML configs found in shipped datasets"
        try:
            import yaml  # noqa: F401
        except ImportError:
            pytest.skip("pyyaml not installed")
        cfg = DatasetConfig.from_yaml(yaml_files[0])
        assert cfg.slug
        assert cfg.name
        assert isinstance(cfg.targets, list)


# ---------------------------------------------------------------------------
# Metrics tests
# ---------------------------------------------------------------------------


class TestMetrics:
    def test_enrichment_factor_perfect(self):
        # All actives at the top
        scores = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        labels = [True, True, False, False, False, False, False, False, False, False]
        ef = enrichment_factor(scores, labels, fraction=0.2)
        assert ef == pytest.approx(5.0)

    def test_enrichment_factor_empty(self):
        assert enrichment_factor([], [], fraction=0.01) == 0.0

    def test_enrichment_factor_no_actives(self):
        assert enrichment_factor([1, 2, 3], [False, False, False]) == 0.0

    def test_log_auc_empty(self):
        assert log_auc([], []) == 0.0

    def test_log_auc_returns_float(self):
        scores = list(range(100))
        labels = [i < 10 for i in range(100)]
        result = log_auc(scores, labels)
        assert isinstance(result, float)
        assert 0.0 <= result <= 1.0

    def test_scoring_power_perfect_correlation(self):
        pred = [1.0, 2.0, 3.0, 4.0, 5.0]
        exp = [1.0, 2.0, 3.0, 4.0, 5.0]
        sp = scoring_power(pred, exp)
        assert sp["pearson_r"] == pytest.approx(1.0)
        assert sp["rmse"] == pytest.approx(0.0, abs=1e-10)

    def test_scoring_power_too_few(self):
        sp = scoring_power([1.0], [2.0])
        assert math.isnan(sp["pearson_r"])

    def test_docking_power(self, sample_poses):
        dp = docking_power(sample_poses, rmsd_threshold=2.0, top_n=1)
        assert 0.0 <= dp <= 1.0

    def test_entropy_rescue_rate_type(self, sample_poses):
        rate = entropy_rescue_rate(sample_poses)
        assert isinstance(rate, float)
        assert 0.0 <= rate <= 1.0

    def test_target_specificity_zscore(self):
        target = [1.0, 2.0, 3.0]
        background = [5.0, 6.0, 7.0]
        z = target_specificity_zscore(target, background)
        assert z < 0  # targets score lower (better) than background

    def test_target_specificity_zscore_zero_variance(self):
        z = target_specificity_zscore([1.0], [3.0])
        assert math.isnan(z)

    def test_hit_rate_top_n(self):
        scores = [1.0, 2.0, 3.0, 4.0, 5.0]
        labels = [True, False, True, False, False]
        rate = hit_rate_top_n(scores, labels, n=2)
        assert rate == pytest.approx(0.5)

    def test_bootstrap_ci_shape(self):
        data = list(range(50))
        lo, hi = bootstrap_ci(lambda d: sum(d) / len(d), data, n_resamples=100)
        assert lo <= hi
        assert not math.isnan(lo)

    def test_compute_all_metrics(self, sample_poses):
        metrics = compute_all_metrics(sample_poses)
        assert isinstance(metrics, dict)
        assert "entropy_rescue_rate" in metrics
        assert "docking_power_top1" in metrics


# ---------------------------------------------------------------------------
# DatasetResult tests
# ---------------------------------------------------------------------------


class TestDatasetResult:
    def test_to_dict_keys(self):
        cfg = DatasetConfig(slug="test", name="Test", description="")
        dr = DatasetResult(config=cfg, tier=1, metrics={"docking_power_top1": 0.5})
        d = dr.to_dict()
        assert d["dataset"] == "test"
        assert d["tier"] == 1
        assert "metrics" in d

    def test_check_regressions_passes(self):
        cfg = DatasetConfig(
            slug="test", name="Test", description="",
            expected_baselines={"docking_power_top1": 0.5},
            baseline_tolerance=0.05,
        )
        dr = DatasetResult(config=cfg, tier=1, metrics={"docking_power_top1": 0.6})
        flags = dr.check_regressions()
        assert flags.get("docking_power_top1") is False

    def test_check_regressions_flags(self):
        cfg = DatasetConfig(
            slug="test", name="Test", description="",
            expected_baselines={"docking_power_top1": 0.8},
            baseline_tolerance=0.05,
        )
        dr = DatasetResult(config=cfg, tier=1, metrics={"docking_power_top1": 0.5})
        flags = dr.check_regressions()
        assert flags.get("docking_power_top1") is True


# ---------------------------------------------------------------------------
# BenchmarkReport tests
# ---------------------------------------------------------------------------


class TestBenchmarkReport:
    def test_to_json(self):
        report = BenchmarkReport(
            generated_at="2026-01-01T00:00:00Z",
            git_sha="abc1234",
            host="testhost",
        )
        j = report.to_json()
        data = json.loads(j)
        assert data["git_sha"] == "abc1234"

    def test_to_markdown(self):
        cfg = DatasetConfig(slug="test", name="Test Dataset", description="")
        dr = DatasetResult(config=cfg, tier=1, metrics={"docking_power_top1": 0.75})
        report = BenchmarkReport(
            datasets=[dr],
            generated_at="2026-01-01T00:00:00Z",
        )
        md = report.to_markdown()
        assert "# FlexAIDdS Benchmark Report" in md
        assert "Test Dataset" in md

    def test_save_roundtrip(self, tmp_path):
        cfg = DatasetConfig(slug="rt", name="Roundtrip", description="")
        dr = DatasetResult(config=cfg, tier=1, metrics={"ef_1pct": 3.5})
        report = BenchmarkReport(
            datasets=[dr],
            generated_at="2026-01-01T00:00:00Z",
            git_sha="deadbeef",
            host="testhost",
        )
        prefix = tmp_path / "report"
        json_path, md_path = report.save(prefix)
        assert json_path.is_file()
        assert md_path.is_file()

        loaded = BenchmarkReport.load(json_path)
        assert len(loaded.datasets) == 1
        assert loaded.datasets[0].config.slug == "rt"
        assert loaded.git_sha == "deadbeef"


# ---------------------------------------------------------------------------
# DatasetRunner tests
# ---------------------------------------------------------------------------


class TestDatasetRunner:
    def test_discover_datasets(self):
        """Shipped YAML configs should be discoverable."""
        try:
            import yaml  # noqa: F401
        except ImportError:
            pytest.skip("pyyaml not installed")
        runner = DatasetRunner(dry_run=True, results_dir="/tmp/flexaidds_test_results")
        configs = runner.discover_datasets()
        assert len(configs) >= 9, f"Expected >=9 dataset configs, found {len(configs)}"
        slugs = {c.slug for c in configs}
        assert "casf2016" in slugs
        assert {"astex_diverse", "astex_nonnative", "hap2"} <= slugs

    def test_flexaid_2015_datasets_are_first(self):
        """FlexAID 2015 validation datasets should lead all-dataset runs."""
        try:
            import yaml  # noqa: F401
        except ImportError:
            pytest.skip("pyyaml not installed")
        runner = DatasetRunner(dry_run=True, results_dir="/tmp/flexaidds_test_results")
        configs = runner.discover_datasets()
        assert [cfg.slug for cfg in configs[:3]] == [
            "astex_diverse",
            "astex_nonnative",
            "hap2",
        ]

    def test_run_single_dry_run(self, tmp_path):
        """Dry-run a single dataset at tier 1."""
        try:
            import yaml  # noqa: F401
        except ImportError:
            pytest.skip("pyyaml not installed")
        runner = DatasetRunner(dry_run=True, results_dir=str(tmp_path))
        dr = runner.run_single("casf2016", tier=1)
        assert dr.config.slug == "casf2016"
        assert dr.tier == 1
        assert len(dr.targets_attempted) > 0
        assert isinstance(dr.metrics, dict)

    def test_run_single_nonexistent(self, tmp_path):
        runner = DatasetRunner(dry_run=True, results_dir=str(tmp_path))
        with pytest.raises(FileNotFoundError):
            runner.run_single("nonexistent_dataset_xyz")


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------


class TestCLI:
    def test_help(self):
        from flexaidds.dataset_runner.cli import main
        with pytest.raises(SystemExit) as exc_info:
            main(["--help"])
        assert exc_info.value.code == 0

    def test_dry_run_tier1(self, tmp_path):
        try:
            import yaml  # noqa: F401
        except ImportError:
            pytest.skip("pyyaml not installed")
        from flexaidds.dataset_runner.cli import main
        ret = main([
            "--dataset", "casf2016",
            "--tier", "1",
            "--dry-run",
            "--results-dir", str(tmp_path),
            "--report-prefix", str(tmp_path / "report"),
        ])
        assert ret == 0 or ret == 1  # 1 if regression flags triggered
