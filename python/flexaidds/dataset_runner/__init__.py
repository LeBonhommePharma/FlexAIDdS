"""FlexAIDdS distributed benchmarking suite.

Import the main orchestrator::

    from flexaidds.dataset_runner import DatasetRunner
    runner = DatasetRunner(results_dir="results/")
    report = runner.run_all(tier=1)
    print(report.to_markdown())

Or use the CLI::

    python -m flexaidds.dataset_runner --dataset casf2016 --tier 1
    python -m flexaidds.dataset_runner --all --tier 2 --distributed --nodes 4
"""

from .runner import (
    DatasetConfig,
    DatasetResult,
    TargetResult,
    BenchmarkReport,
    DatasetRunner,
    load_entry_manifest,
    completed_targets_from_manifest,
    plan_runtime,
    KNOWN_LARGE_DATASETS,
)
from .metrics import (
    PoseScore,
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

__all__ = [
    # Orchestrator
    "DatasetRunner",
    "DatasetConfig",
    "DatasetResult",
    "TargetResult",
    "BenchmarkReport",
    "load_entry_manifest",
    "completed_targets_from_manifest",
    "plan_runtime",
    "KNOWN_LARGE_DATASETS",
    # Metrics
    "PoseScore",
    "entropy_rescue_rate",
    "enrichment_factor",
    "log_auc",
    "scoring_power",
    "docking_power",
    "target_specificity_zscore",
    "hit_rate_top_n",
    "bootstrap_ci",
    "compute_all_metrics",
]
