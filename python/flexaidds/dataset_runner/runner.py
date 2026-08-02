"""DatasetRunner — distributed benchmarking orchestrator for FlexAIDdS.

Discovers dataset configs, distributes work across MPI nodes or local
processes, runs FlexAIDdS docking, computes all metrics, and produces
structured JSON + Markdown reports.

THERMO/ITC MODE (bulletproof, always active for itc187/bindingdb_itc/scorpio):
- Temperature is auto-forced to 298.0 K.
- Metrics are restricted to scoring_power_* + entropy_rescue_rate only
  (docking_power_* deliberately excluded).
- Full provenance is ALWAYS emitted: git_sha, binary, host, temperature,
  full_command (even in dry-run, for any dataset/method/sequence).
- Works for --dry-run, real runs, --all, comma sequences, custom YAMLs.

Typical usage
-------------
Library::

    from flexaidds.dataset_runner import DatasetRunner

    runner = DatasetRunner(results_dir="results/benchmark_run")
    report = runner.run_all(tier=1)
    json_path, md_path = report.save("results/benchmark_run/report")

CLI::

    python -m flexaidds.dataset_runner --dataset itc187 --tier 1 --temperature 298 --dry-run
    python -m flexaidds.dataset_runner --dataset itc187,bindingdb_itc --tier 1 --dry-run

MPI distributed run::

    mpirun -n 8 python -m flexaidds.dataset_runner --all --distributed
"""

from __future__ import annotations

import datetime
import json
import csv
import logging
import os
import re
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

try:
    import yaml  # PyYAML
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False

import numpy as np

from .data_paths import (
    _CF_APP_REMARK_RE,
    _CF_REMARK_RE,
    _ENTROPY_REMARK_RE,
    _RMSD_REMARK_RE,
    resolve_benchmark_paths,
)
from .metrics import (
    PoseScore,
    bootstrap_ci,
    compute_all_metrics,
    entropy_rescue_rate,
    docking_power,
)

logger = logging.getLogger(__name__)

# Tier-2 full-N datasets use committed entry catalogs (C++ parity).
KNOWN_LARGE_DATASETS: Dict[str, int] = {
    "astex_nonnative": 1113,
    "posex_cd": 1312,
}
_CATALOG_PATH = (
    Path(__file__).resolve().parents[3]
    / "benchmarks"
    / "m3pro"
    / "large_dataset_entry_catalogs.json"
)

# ---------------------------------------------------------------------------
# Thermo/ITC focus (bulletproof per handoff + validation prompt)
# ITC/thermo datasets ALWAYS use 298 K, only scoring_power_* + entropy_rescue_rate
# metrics, and full provenance (git_sha + binary + host + temperature + full_command).
# This applies for dry-run and real runs, sequences, any kernel/method.
# ---------------------------------------------------------------------------

THERMO_DATASETS: set[str] = {"itc187", "bindingdb_itc", "scorpio"}
THERMO_DEFAULT_TEMP: float = 298.0
THERMO_METRICS: set[str] = {
    "scoring_power_pearson_r",
    "scoring_power_spearman_r",
    "scoring_power_rmse",
    "scoring_power_mae",
    "entropy_rescue_rate",
}

# Dry-run uses synthetic poses for pipeline smoke tests only. Never report
# docking_power_* as if it were a real docking success rate.
DRY_RUN_METRICS_NOTE: str = (
    "dry_run=True: poses are synthetic (pipeline smoke test only). "
    "docking_power_* metrics are omitted because they would not reflect real "
    "docking success rates. Any remaining metrics (scoring_power_*, "
    "entropy_rescue_rate, etc.) are also computed on synthetic data and must "
    "not be reported as production results."
)


def _is_docking_power_metric(name: str) -> bool:
    """True for docking_power_* keys (including bootstrap CI suffixes)."""
    return name.startswith("docking_power_")


def _strip_docking_power_metrics(metrics: Dict[str, Any]) -> Dict[str, Any]:
    """Drop docking_power_* entries so synthetic dry-run rates are never sold as real."""
    return {k: v for k, v in metrics.items() if not _is_docking_power_metric(k)}


def _filter_requested_metrics_for_dry_run(
    requested: Optional[List[str]],
) -> Optional[List[str]]:
    """Remove docking_power_* from a requested metric list for dry-run runs."""
    if requested is None:
        return None
    filtered = [m for m in requested if not _is_docking_power_metric(m)]
    return filtered


# ---------------------------------------------------------------------------
# Configuration dataclasses
# ---------------------------------------------------------------------------


@dataclass
class DatasetConfig:
    """Declarative specification for one benchmark dataset.

    Loaded from a YAML file in ``benchmarks/datasets/``.

    Attributes:
        slug:                  Filesystem-safe short identifier (e.g. ``casf2016``).
        name:                  Human-readable dataset name.
        description:           One-paragraph description of the dataset.
        zenodo_doi:            Zenodo DOI for the canonical download (may be empty).
        download_url:          Alternative download URL.
        tier:                  Minimum tier required to run the full dataset (1 or 2).
        tier1_subset_size:     Number of targets to use for tier-1 (PR sanity) runs.
        benchmark_order:       Stable run/display order; lower values run first.
        targets:               Full list of target identifiers.
        structural_states:     Receptor states available (``holo``, ``apo``, ``af2``).
        docking_mode:          Normative semantics: self_docking | cross_docking |
                               affinity_scoring | virtual_screening | specialized.
        metrics:               Names of metrics to compute (must exist in metrics.py).
        expected_baselines:    ``{metric: value}`` reference values for regression checks.
        published_baselines:   ``{metric: value}`` published reference values for comparisons.
        published_source:      Human-readable citation for the published baselines.
        baseline_tolerance:    Fractional tolerance for regression detection (default 0.05).
        data_dir:              Local path to dataset files (None = not yet downloaded).
        data_format:           File format: ``"pdb"`` or ``"mol2"``.
        active_label_field:    Field in target metadata that encodes active/decoy status.
    """

    slug: str
    name: str
    description: str
    zenodo_doi: str = ""
    download_url: str = ""
    tier: int = 2
    tier1_subset_size: int = 5
    benchmark_order: int = 1000
    targets: List[str] = field(default_factory=list)
    structural_states: List[str] = field(default_factory=lambda: ["holo"])
    docking_mode: str = "self_docking"
    metrics: List[str] = field(default_factory=list)
    expected_baselines: Dict[str, float] = field(default_factory=dict)
    published_baselines: Dict[str, float] = field(default_factory=dict)
    published_source: str = ""
    baseline_tolerance: float = 0.05
    data_dir: Optional[Path] = None
    data_format: str = "pdb"
    active_label_field: str = "is_active"
    default_conc_M: float = 1.0  # P3: for grand canonical per-ligand or default
    ligand_concs: Dict[str, float] = field(default_factory=dict)  # P3: per-ligand_id -> conc_M from yaml

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "DatasetConfig":
        """Load a DatasetConfig from a YAML file."""
        if not _HAS_YAML:
            raise RuntimeError(
                "PyYAML is required to load dataset configs: pip install pyyaml"
            )
        with open(yaml_path) as fh:
            raw: dict = yaml.safe_load(fh)

        # Fail closed on self vs cross-docking contradictions before scheduling work.
        # parents: runner.py → dataset_runner → flexaidds → python → repo root
        skill_validator = (
            Path(__file__).resolve().parents[3]
            / ".grok"
            / "skills"
            / "flexaidds"
            / "scripts"
            / "validate_dataset_semantics.py"
        )
        if skill_validator.is_file():
            import importlib.util

            spec = importlib.util.spec_from_file_location(
                "validate_dataset_semantics", skill_validator
            )
            if spec and spec.loader:
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
                sem_errs = mod.validate_config(dict(raw), path=yaml_path)
                if sem_errs:
                    raise ValueError(
                        "Dataset docking_mode semantics failed:\n  - "
                        + "\n  - ".join(sem_errs)
                    )

        data_dir_raw = raw.pop("data_dir", None)
        # Extra YAML metadata keys are allowed; ignore unknown fields not on DatasetConfig.
        for extra in (
            "canonical_data_root",
            "checksum_manifest",
            "canonical_doc",
            "structure_recovery_policy",
            "claude_reevaluated_baselines",
            "claude_reevaluated_rationale",
        ):
            raw.pop(extra, None)
        docking_mode = str(raw.pop("docking_mode", "self_docking")).strip().lower()
        config = cls(
            slug=raw.pop("slug", yaml_path.stem),
            name=raw.pop("name", yaml_path.stem),
            description=raw.pop("description", ""),
            zenodo_doi=raw.pop("zenodo_doi", ""),
            download_url=raw.pop("download_url", ""),
            tier=int(raw.pop("tier", 2)),
            tier1_subset_size=int(raw.pop("tier1_subset_size", 5)),
            benchmark_order=int(raw.pop("benchmark_order", 1000)),
            targets=list(raw.pop("targets", [])),
            structural_states=list(raw.pop("structural_states", ["holo"])),
            docking_mode=docking_mode,
            metrics=list(raw.pop("metrics", [])),
            expected_baselines=dict(raw.pop("expected_baselines", {})),
            published_baselines=dict(raw.pop("published_baselines", {})),
            published_source=str(raw.pop("published_source", "")),
            baseline_tolerance=float(raw.pop("baseline_tolerance", 0.05)),
            data_format=str(raw.pop("data_format", "pdb")),
            active_label_field=str(raw.pop("active_label_field", "is_active")),
            default_conc_M=float(raw.pop("default_conc_M", 1.0)),
            ligand_concs=DatasetConfig._parse_ligand_concs(raw),
        )
        if data_dir_raw:
            # Absolute at construction: the engine runs with cwd=<per-entry tmp
            # dir> (DatasetRunner._run_flexaid), so a relative data_dir resolves
            # against that tmp dir for the child and finds nothing. It is never
            # correct to hand a subprocess a relative directory, regardless of
            # where the runner itself was invoked from.
            config.data_dir = Path(data_dir_raw).expanduser().resolve()
        return config

    def tier1_targets(self) -> List[str]:
        """Return the subset of targets used for tier-1 (fast) runs."""
        return self.targets[: self.tier1_subset_size]

    def _uses_crossdock_catalog(self, tier: int) -> bool:
        """Full-N crossdock catalog applies only when YAML requests crossdock state."""
        return (
            tier == 2
            and self.slug in KNOWN_LARGE_DATASETS
            and "crossdock" in self.structural_states
        )

    def scheduled_targets(self, tier: int) -> List[str]:
        """Targets actually scheduled for this tier."""
        if self._uses_crossdock_catalog(tier):
            return [e["entry_id"] for e in load_large_dataset_catalog(self.slug)]
        return self.tier1_targets() if tier == 1 else self.targets

    def scheduled_work_items(self, tier: int) -> List[Tuple[str, str]]:
        """(entry_id, state) pairs scheduled for this tier."""
        if self._uses_crossdock_catalog(tier):
            catalog = load_large_dataset_catalog(self.slug)
            families = {t.upper() for t in self.targets} if self.targets else None
            items: List[Tuple[str, str]] = []
            for entry in catalog:
                if families and entry.get("family", "").upper() not in families:
                    continue
                items.append((entry["entry_id"], entry.get("state", "crossdock")))
            if items:
                return items
            return [(e["entry_id"], e.get("state", "crossdock")) for e in catalog]
        targets = self.scheduled_targets(tier)
        return [(tid, st) for tid in targets for st in self.structural_states]

    def effective_entry_count(self, tier: int) -> int:
        """Entry count for runtime planning."""
        if self._uses_crossdock_catalog(tier):
            items = self.scheduled_work_items(tier)
            if items:
                return len(items)
            return KNOWN_LARGE_DATASETS[self.slug]
        return len(self.scheduled_work_items(tier))

    @staticmethod
    def _parse_ligand_concs(raw: dict) -> Dict[str, float]:
        """P3: parse per-ligand conc_M from 'ligands' list or competition_sets in yaml."""
        concs: Dict[str, float] = {}
        # direct 'ligands': [{ligand_id: , conc_M: }, ...]
        for lig in raw.get("ligands", []) or []:
            if isinstance(lig, dict):
                lid = lig.get("ligand_id") or lig.get("name")
                if lid and "conc_M" in lig:
                    concs[str(lid)] = float(lig["conc_M"])
        # competition_sets style from example yamls
        for cs in raw.get("competition_sets", []) or []:
            for lig in cs.get("ligands", []) or []:
                if isinstance(lig, dict):
                    lid = lig.get("ligand_id") or lig.get("name")
                    if lid and "conc_M" in lig:
                        concs[str(lid)] = float(lig["conc_M"])
        return concs


def load_large_dataset_catalog(
    slug: str,
    catalog_path: Optional[Union[str, Path]] = None,
) -> List[Dict[str, Any]]:
    """Load tier-2 work-item catalog for large datasets (1113+ entries)."""
    p = Path(catalog_path) if catalog_path else _CATALOG_PATH
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text())
        return list(data.get("datasets", {}).get(slug, {}).get("entries", []))
    except Exception:
        return []


# ---------------------------------------------------------------------------
# Result dataclasses
# ---------------------------------------------------------------------------


@dataclass
class TargetResult:
    """Docking results for a single target.

    Attributes:
        target_id:          Target identifier.
        structural_state:   Receptor state used (``holo`` / ``apo`` / ``af2``).
        poses:              All scored poses from this docking run.
        duration_seconds:   Wall-clock time for this target.
        error:              Non-empty if docking failed.
    """

    target_id: str
    structural_state: str
    poses: List[PoseScore]
    duration_seconds: float = 0.0
    error: str = ""
    grand_log_Z: Optional[float] = None  # P3: ensemble log_Z for this ligand (for grand Xi)
    conc_M: float = 1.0  # P3: the conc used for this ligand
    grand_xi: Optional[float] = None  # P3: computed log_Xi if grand_log_Z present

    @property
    def success(self) -> bool:
        return not self.error and bool(self.poses)


# ---------------------------------------------------------------------------
# Regression direction registry
# ---------------------------------------------------------------------------
# check_regressions used to infer "lower is better" from the substring "rmse"
# or "mae". That silently inverted every distance metric whose name is spelled
# differently -- mean_rmsd, median_rmsd, crossdock_*_rmsd, selectivity_log_error
# -- so a 37% WORSE RMSD scored as "no regression". Direction is a property of
# the metric, not of its spelling: declare it here, beside nothing else.
#
# Every metric named in any benchmarks/datasets/*.yaml expected_baselines block
# AND produced by compute_all_metrics must appear in exactly one of these sets.
# An unlisted metric logs a warning and falls back to higher-is-better rather
# than failing the run.
#
# A direction is only registered for a metric something can PRODUCE. A metric
# that nothing computes has no behaviour to contradict its entry, so an entry
# for it is not a record -- it is an instruction: the next person to implement
# the metric reads the registry and writes code to match, and the guess becomes
# true by recruiting the person who could have falsified it. Eleven such names
# were removed here (they remain declared in dataset YAMLs and will be
# re-registered by whoever dispatches them, in the same diff as the producer,
# where the sign convention and the assignment are visible together).
_LOWER_IS_BETTER = frozenset({
    "mean_rmsd",
    "median_rmsd",
    "scoring_power_rmse",
})

_HIGHER_IS_BETTER = frozenset({
    "docking_power_top1",
    "docking_power_top3",
    "ef_1pct",
    "ef_5pct",
    "entropy_rescue_rate",
    "hit_rate_top10",
    "log_auc",
    "scoring_power_pearson_r",
    "scoring_power_spearman_r",
})

# Deliberately in NEITHER set: shannon_energy_collapse. Its preferred direction
# is not determinable from the name or from any docstring in this repo, so it
# warns rather than silently taking a side. Declare it once someone who knows
# the physics says which way is better.


@dataclass
class DatasetResult:
    """Aggregated results for one complete dataset run.

    Attributes:
        config:             The dataset config that was run.
        tier:               Tier at which this run was executed.
        metrics:            ``{metric_name: scalar_value}`` computed across all targets.
        ci_95:              ``{metric_name: (lower, upper)}`` bootstrap CIs.
        regression_flags:   ``{metric_name: True}`` when metric regressed vs baseline.
        targets_attempted:  All target IDs that were scheduled.
        targets_completed:  Target IDs that produced at least one pose.
        targets_failed:     Target IDs where docking failed or produced no poses.
        duration_seconds:   Total wall-clock time for this dataset.
        timestamp:          ISO-8601 timestamp of run completion.
        git_sha:            Git commit SHA of the FlexAIDdS repo at run time.
        host:               Hostname of the machine that ran this dataset.
        grand_summary:      P3: per-ligand grand info (log_Z, conc, Xi, p_bind etc) for grand canonical emission.
    """

    config: DatasetConfig
    tier: int
    metrics: Dict[str, float] = field(default_factory=dict)
    ci_95: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    regression_flags: Dict[str, bool] = field(default_factory=dict)
    targets_attempted: List[str] = field(default_factory=list)
    targets_completed: List[str] = field(default_factory=list)
    targets_failed: List[str] = field(default_factory=list)
    duration_seconds: float = 0.0
    timestamp: str = ""
    git_sha: str = ""
    host: str = ""
    grand_summary: Dict[str, Dict[str, float]] = field(default_factory=dict)  # P3 grand
    binary: str = ""
    temperature: float = 300.0
    full_command: str = ""
    dry_run: bool = False
    metrics_note: str = ""
    # Benchmark liveness/productivity gate (#326): a run where the binary
    # crashed or produced nothing must not read as "no regression". These
    # fields let cli.main distinguish PASS / FAIL / INCONCLUSIVE.
    flexaid_crashes: int = 0
    total_poses: int = 0
    # value is a real exit/return code, or None when the engine never executed
    # (exec failure) — mirrors DatasetRunner._entry_exit_codes.
    entry_exit_codes: Dict[str, Optional[int]] = field(default_factory=dict)
    inconclusive_metrics: List[str] = field(default_factory=list)
    # Targets actually executed this run (excludes --resume checkpoints, whose
    # poses are not reloaded into all_poses). Gates 2-3 only judge a run that
    # did fresh work — a pure resume/package-regen run must not false-fail.
    newly_executed: int = 0
    # Targets loaded from --resume checkpoints (not re-executed). Distinguishes
    # a legitimate resume (skip gates 2-3) from a run that scheduled nothing at
    # all (newly==0 AND resumed==0 → still INCONCLUSIVE, not a silent pass).
    resumed: int = 0

    def check_regressions(self) -> Dict[str, bool]:
        """Flag metrics that regressed below baseline − tolerance.

        Returns dict of ``{metric: True}`` for regressed metrics.
        """
        flags: Dict[str, bool] = {}
        missing: List[str] = []
        tol = self.config.baseline_tolerance
        for metric, baseline in self.config.expected_baselines.items():
            measured = self.metrics.get(metric)
            if measured is None or np.isnan(measured):
                # Completeness gate (#326): a metric that was never measured is
                # NOT "no regression" — it is an inconclusive run. Record it so
                # cli.main can fail the run rather than silently pass it.
                missing.append(metric)
                continue
            # Direction is looked up, never inferred from the name. See
            # _LOWER_IS_BETTER: substring-sniffing silently inverted
            # mean_rmsd/median_rmsd because "rmsd" is not "rmse".
            # Tolerance is ADDITIVE slack around the baseline, not a scaling
            # of it.  `baseline * (1 + tol)` happens to equal
            # `baseline + abs(baseline) * tol` for every POSITIVE baseline, so
            # the scaling form was correct until a negative one appeared:
            # erds_specificity declares target_specificity_zscore: -2.50, and
            # scaling moves a negative baseline AWAY from zero -- demanding
            # -2.625, stricter than the target itself, so exactly-on-baseline
            # flagged as a regression.  The additive form is sign-correct in
            # both directions and byte-identical on all 143 positive baselines
            # currently in the repo.
            #
            # Degenerate case, stated so the next person finds it rather than
            # discovers it: at baseline == 0 the slack is abs(0) * tol == 0, so
            # ANY non-zero measurement flags.  Arguably correct -- 5% of zero is
            # zero -- but "additive slack" stops being slack there.  The scaling
            # form has the identical hole; no dataset declares a zero baseline
            # today (checked across both trees).
            slack = abs(baseline) * tol
            if metric in _LOWER_IS_BETTER:
                threshold = baseline + slack
                flags[metric] = bool(measured > threshold)
            else:
                if metric not in _HIGHER_IS_BETTER:
                    logger.warning(
                        "Metric %r has no declared direction in "
                        "_LOWER_IS_BETTER/_HIGHER_IS_BETTER; assuming "
                        "higher-is-better. Declare it before trusting this "
                        "regression verdict.",
                        metric,
                    )
                threshold = baseline - slack
                flags[metric] = bool(measured < threshold)
        self.regression_flags = flags
        self.inconclusive_metrics = missing
        return flags

    def to_dict(self) -> dict:
        """Serialise to a JSON-compatible dict."""
        payload = {
            "dataset": self.config.slug,
            "tier": self.tier,
            "timestamp": self.timestamp,
            "git_sha": self.git_sha,
            "host": self.host,
            "binary": self.binary,
            "temperature": self.temperature,
            "full_command": self.full_command,
            "duration_seconds": self.duration_seconds,
            "targets_attempted": len(self.targets_attempted),
            "targets_completed": len(self.targets_completed),
            "targets_failed": len(self.targets_failed),
            "metrics": self.metrics,
            "ci_95": {k: list(v) for k, v in self.ci_95.items()},
            "regression_flags": self.regression_flags,
            "expected_baselines": self.config.expected_baselines,
            "published_baselines": self.config.published_baselines,
            "published_source": self.config.published_source,
            "dry_run": self.dry_run,
            "grand_summary": self.grand_summary,  # P3
            # #326 gate provenance
            "flexaid_crashes": self.flexaid_crashes,
            "total_poses": self.total_poses,
            "entry_exit_codes": self.entry_exit_codes,
            "inconclusive_metrics": self.inconclusive_metrics,
            "newly_executed": self.newly_executed,
            "resumed": self.resumed,
        }
        if self.metrics_note:
            payload["metrics_note"] = self.metrics_note
        return payload


# ---------------------------------------------------------------------------
# Benchmark report
# ---------------------------------------------------------------------------


@dataclass
class BenchmarkReport:
    """Aggregated report across all datasets in a run.

    Attributes:
        datasets:      One :class:`DatasetResult` per dataset run.
        generated_at:  ISO-8601 timestamp.
        git_sha:       Repo commit SHA.
        host:          Hostname.
        runner_info:   Dict of environment/runtime metadata.
    """

    datasets: List[DatasetResult] = field(default_factory=list)
    generated_at: str = ""
    git_sha: str = ""
    host: str = ""
    runner_info: Dict[str, str] = field(default_factory=dict)
    temperature: float = 300.0
    binary: str = ""
    full_command: str = ""

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        return {
            "generated_at": self.generated_at,
            "git_sha": self.git_sha,
            "host": self.host,
            "runner_info": self.runner_info,
            "temperature": self.temperature,
            "binary": self.binary,
            "full_command": self.full_command,
            "datasets": [d.to_dict() for d in self.datasets],
        }

    def to_json(self, indent: int = 2) -> str:
        return json.dumps(self.to_dict(), indent=indent)

    def to_markdown(self) -> str:
        """Generate a human-readable Markdown summary."""
        lines = [
            "# FlexAIDdS Benchmark Report",
            "",
            f"**Generated**: {self.generated_at}  ",
            f"**Commit**: `{self.git_sha or 'unknown'}`  ",
            f"**Host**: {self.host}  ",
            f"**Temperature**: {getattr(self, 'temperature', 300.0)} K  ",
            f"**Binary**: {getattr(self, 'binary', 'unknown')}  ",
            f"**Command**: {getattr(self, 'full_command', '')[:200]}  ",
            "",
            "---",
            "",
        ]

        for dr in self.datasets:
            lines += self._dataset_section(dr)

        return "\n".join(lines)

    @staticmethod
    def _load_entry_manifest_summary(dr: DatasetResult) -> Optional[dict]:
        """Try to locate and load the rich per-entry manifest generated by EntryTaskManager."""
        candidates = [
            Path(f"results/{dr.config.slug}/tier{dr.tier}/_entry_manifest.json"),
            Path("results/benchmarks") / dr.config.slug / f"tier{dr.tier}" / "_entry_manifest.json",
            Path.cwd() / "results" / dr.config.slug / f"tier{dr.tier}" / "_entry_manifest.json",
        ]
        for p in candidates:
            if p.is_file():
                try:
                    return json.loads(p.read_text())
                except Exception:
                    continue
        return None

    @staticmethod
    def _dataset_section(dr: DatasetResult) -> List[str]:
        lines = [
            f"## {dr.config.name}",
            "",
            f"*Tier {dr.tier} · {len(dr.targets_completed)}/{len(dr.targets_attempted)} targets · "
            f"{dr.duration_seconds:.1f}s*",
            "",
        ]

        if dr.dry_run:
            note = dr.metrics_note or DRY_RUN_METRICS_NOTE
            lines += [
                f"> **DRY-RUN** — {note}",
                "",
            ]

        if dr.metrics:
            lines += ["| Metric | Value | 95% CI | Baseline | Published | Regressed? |",
                      "|--------|-------|--------|----------|-----------|------------|"]
            for metric, value in sorted(dr.metrics.items()):
                ci = dr.ci_95.get(metric)
                ci_str = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "—"
                baseline = dr.config.expected_baselines.get(metric)
                baseline_str = f"{baseline:.3f}" if baseline is not None else "—"
                published = dr.config.published_baselines.get(metric)
                published_str = f"{published:.3f}" if published is not None else "—"
                regressed = dr.regression_flags.get(metric, False)
                flag = "⚠ YES" if regressed else "OK"
                lines.append(
                    f"| {metric} | {value:.4f} | {ci_str} | {baseline_str} | {published_str} | {flag} |"
                )
            lines.append("")

        if dr.config.published_source:
            lines += [
                f"**Published comparator**: {dr.config.published_source}",
                "",
            ]

        # P3: grand canonical summary emission (Ξ, p_bind etc) if present
        if dr.grand_summary:
            lines += ["### Grand Canonical Summary (P3)", ""]
            lines += ["| Ligand | log_Z | conc_M | log_Xi | p_bind |", "|--------|-------|--------|--------|--------|"]
            for lid, info in sorted(dr.grand_summary.items()):
                lines.append(f"| {lid} | {info.get('log_Z', 0):.4f} | {info.get('conc_M', 1):.2e} | {info.get('log_Xi', 0):.4f} | {info.get('p_bind', 0):.4f} |")
            lines.append("")

        if dr.targets_failed:
            lines += [
                f"**Failed targets** ({len(dr.targets_failed)}): "
                + ", ".join(dr.targets_failed[:10])
                + ("…" if len(dr.targets_failed) > 10 else ""),
                "",
            ]

        # === Richer per-entry cost exposure ===
        manifest = BenchmarkReport._load_entry_manifest_summary(dr)
        manifest_rel = f"results/{dr.config.slug}/tier{dr.tier}/_entry_manifest.json"

        if manifest and "timings" in manifest:
            t = manifest["timings"]
            summary = t.get("summary", {})
            per_cost = t.get("per_entry_cost_cpu_seconds", {})

            lines.append("### Per-Entry Cost & Timing Summary (from EntryTaskManager)")
            lines.append(f"**Manifest**: `{manifest_rel}`")
            lines.append("")

            lines.append(f"- Entries timed: **{summary.get('num_timed_entries', 0)}**")
            lines.append(f"- Total estimated cost: **{summary.get('estimated_total_cost_cpu_seconds', 0):.1f}** CPU-seconds")
            lines.append(f"- Mean / Median / P90: {summary.get('mean_entry_seconds', 0):.2f}s / "
                         f"{summary.get('median_entry_seconds', 0):.2f}s / "
                         f"{summary.get('p90_entry_seconds', 0):.2f}s")
            lines.append("")

            if per_cost:
                sorted_cost = sorted(per_cost.items(), key=lambda x: -x[1])[:5]
                lines.append("**Top 5 most expensive entries:**")
                for name, cost in sorted_cost:
                    lines.append(f"- `{name}`: {cost:.2f} CPU-s")
                lines.append("")

            lines.append(f"Full per-entry breakdown (wall time + costs) is in the manifest.")
            lines.append("")
        else:
            lines.append(f"**Per-entry artifacts**: `{manifest_rel}`")
            lines.append("Run with the DatasetRunner (especially `--resume`) to generate detailed per-entry timing and cost data.")
            lines.append("")

        return lines

    def save(self, prefix: Union[str, Path]) -> Tuple[Path, Path]:
        """Save JSON and Markdown reports.

        Args:
            prefix: File path prefix (without extension).

        Returns:
            ``(json_path, md_path)``
        """
        prefix = Path(prefix)
        prefix.parent.mkdir(parents=True, exist_ok=True)
        json_path = prefix.with_suffix(".json")
        md_path = prefix.with_suffix(".md")
        json_path.write_text(self.to_json())
        md_path.write_text(self.to_markdown())
        logger.info("Report saved: %s, %s", json_path, md_path)
        return json_path, md_path

    @classmethod
    def load(cls, json_path: Union[str, Path]) -> "BenchmarkReport":
        """Load a previously saved JSON report (metadata only, no raw poses)."""
        data = json.loads(Path(json_path).read_text())
        # Reconstruct lightweight DatasetResult objects (no config, no poses)
        ds_results = []
        for d in data.get("datasets", []):
            stub_config = DatasetConfig(
                slug=d["dataset"],
                name=d["dataset"],
                description="",
                expected_baselines=d.get("expected_baselines", {}),
            )
            dr = DatasetResult(
                config=stub_config,
                tier=d["tier"],
                metrics=d.get("metrics", {}),
                ci_95={k: tuple(v) for k, v in d.get("ci_95", {}).items()},  # type: ignore[misc]
                regression_flags=d.get("regression_flags", {}),
                timestamp=d.get("timestamp", ""),
                git_sha=d.get("git_sha", ""),
                host=d.get("host", ""),
                duration_seconds=d.get("duration_seconds", 0.0),
                dry_run=bool(d.get("dry_run", False)),
                metrics_note=str(d.get("metrics_note", "") or ""),
                grand_summary=d.get("grand_summary", {}) or {},
            )
            ds_results.append(dr)
        return cls(
            datasets=ds_results,
            generated_at=data.get("generated_at", ""),
            git_sha=data.get("git_sha", ""),
            host=data.get("host", ""),
            runner_info=data.get("runner_info", {}),
        )


# ---------------------------------------------------------------------------
# Helper: Git SHA + env info
# ---------------------------------------------------------------------------


def _git_sha(repo_root: Optional[Path] = None) -> str:
    try:
        r = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True,
            cwd=repo_root or Path.cwd(),
            timeout=5,
        )
        return r.stdout.strip() if r.returncode == 0 else "unknown"
    except Exception:
        return "unknown"


def _runner_info() -> Dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "platform": sys.platform,
        "cpu_count": str(os.cpu_count() or 1),
        "flexaidds_version": _flexaidds_version(),
    }


def _flexaidds_version() -> str:
    try:
        import importlib.metadata
        return importlib.metadata.version("flexaidds")
    except Exception:
        try:
            from flexaidds import __version__
            return __version__
        except Exception:
            return "unknown"


# ---------------------------------------------------------------------------
# MPI / multiprocessing helpers
# ---------------------------------------------------------------------------


def _mpi_context():
    """Return (rank, size, is_root) for the current MPI context.

    Falls back to (0, 1, True) when mpi4py is unavailable.
    """
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        return comm.Get_rank(), comm.Get_size(), comm.Get_rank() == 0, comm
    except ImportError:
        return 0, 1, True, None


def _split_targets(targets: List[Any], rank: int, size: int) -> List[Any]:
    """Distribute any list of work items across MPI ranks (round-robin).
    Works for str target_ids or (target_id, state) tuples.
    """
    return [t for i, t in enumerate(targets) if i % size == rank]


# ---------------------------------------------------------------------------
# EntryTaskManager — the Master coordinator for individual dataset entries
# ---------------------------------------------------------------------------

class EntryTaskManager:
    """Master manager that automates allocation of individual dataset entries
    (target + structural state) to workers.

    Responsibilities:
    - Owns the queue of fine-grained work items (per-entry).
    - Dispatches to local workers (ThreadPoolExecutor today; easy to extend to ProcessPool or MPI master-worker).
    - Tracks completion for resource accounting and resume.
    - Provides a clean hook for future dynamic resource management (CPU pinning, memory-aware scheduling, heterogeneous target costs, etc.).

    This is the central automation the user requested for "how a master datasetrunner manager manages the allocation of resources for individual entries/workers to work together in parallel."
    """

    def __init__(
        self,
        work_items: List[Tuple[str, str]],
        n_workers: int = 1,
        omp_threads: int = 4,
        mpi_comm: Any = None,
        mpi_rank: int = 0,
        mpi_size: int = 1,
        mpi_root: bool = True,
        cost_hints: Optional[Dict[str, float]] = None,
    ):
        self.work_items = list(work_items)
        self.n_workers = max(1, int(n_workers))
        self.omp_threads = max(1, int(omp_threads))
        self.completed: List[Tuple[str, str, List[PoseScore], float, str]] = []

        # MPI context for stronger master-worker mode
        self._mpi_comm = mpi_comm
        self._mpi_rank = mpi_rank
        self._mpi_size = mpi_size
        self._mpi_root = mpi_root

        # Cost-aware scheduling hints (target_state -> estimated cost, e.g. from previous manifest)
        self.cost_hints = cost_hints or {}

        # Apply cost-aware ordering immediately if hints provided (cheaper first)
        if self.cost_hints:
            def _cost_key(item):
                key = f"{item[0]}_{item[1]}"
                return self.cost_hints.get(key, 999999.0)
            self.work_items.sort(key=_cost_key)

    @classmethod
    def load_cost_hints_from_manifest(cls, manifest_path: Union[str, Path]) -> Dict[str, float]:
        """Load per-entry cost hints from a previous _entry_manifest.json (timings section)."""
        p = Path(manifest_path)
        if not p.is_file():
            return {}
        try:
            data = json.loads(p.read_text())
            timings = data.get("timings", {}).get("per_entry_cost_cpu_seconds", {})
            # Also fall back to wall time if cost not present
            if not timings:
                timings = data.get("timings", {}).get("per_entry_wall_seconds", {})
            return {k: float(v) for k, v in timings.items()}
        except Exception:
            return {}

    def run(
        self,
        processor: Callable[[Tuple[str, str]], Tuple[str, str, List[PoseScore], float, str]],
    ) -> List[Tuple[str, str, List[PoseScore], float, str]]:
        """Execute all work items using the processor function.

        Dispatch strategy (in priority order):
        1. MPI multi-rank: dynamic master-worker (root farms tasks on demand).
        2. Local parallel: ThreadPoolExecutor with per-future exception isolation.
        3. Serial: simple loop (n_workers=1 or single item).

        Cost-aware ordering is applied in __init__ when hints are provided.
        """
        if not self.work_items:
            return []

        if self.cost_hints and (self._mpi_root or self._mpi_size <= 1):
            logger.info(
                "EntryTaskManager: cost-aware scheduling enabled (%d hints loaded)",
                len(self.cost_hints),
            )

        if self._mpi_size > 1 and self._mpi_comm is not None:
            return self._run_mpi_master_worker(processor)

        if self.n_workers <= 1 or len(self.work_items) == 1:
            for item in self.work_items:
                res = processor(item)
                self.completed.append(res)
            return list(self.completed)

        logger.info(
            "EntryTaskManager: dispatching %d entries across %d workers "
            "(OMP_NUM_THREADS=%d)",
            len(self.work_items), self.n_workers, self.omp_threads,
        )
        with ThreadPoolExecutor(max_workers=self.n_workers) as pool:
            future_to_item = {
                pool.submit(processor, item): item for item in self.work_items
            }
            for future in as_completed(future_to_item):
                try:
                    res = future.result()
                    self.completed.append(res)
                except Exception as exc:
                    item = future_to_item[future]
                    tid, st = item
                    logger.error("Worker exception on %s/%s: %s", tid, st, exc)
                    self.completed.append((tid, st, [], 0.0, str(exc)))

        return list(self.completed)

    def _run_mpi_master_worker(
        self,
        processor: Callable[[Tuple[str, str]], Tuple[str, str, List[PoseScore], float, str]],
    ) -> List[Tuple[str, str, List[PoseScore], float, str]]:
        """Dynamic MPI master-worker: root queues tasks, workers pull on demand."""
        comm = self._mpi_comm
        rank = self._mpi_rank
        root = 0

        try:
            from mpi4py import MPI as _MPI
            ANY_SRC = _MPI.ANY_SOURCE
            ANY_TAG = _MPI.ANY_TAG
        except Exception:
            ANY_SRC = 0
            ANY_TAG = 0

        TAG_REQ = 11
        TAG_TASK = 12
        TAG_RES = 13

        if rank == root:
            queue = list(self.work_items)
            results: List[Tuple[str, str, List[PoseScore], float, str]] = []
            active_tasks = 0
            total_expected = len(self.work_items)

            while queue or active_tasks > 0 or len(results) < total_expected:
                while comm.Iprobe(source=ANY_SRC, tag=TAG_RES):
                    res = comm.recv(source=ANY_SRC, tag=TAG_RES)
                    results.append(res)
                    active_tasks = max(0, active_tasks - 1)
                if comm.Iprobe(source=ANY_SRC, tag=TAG_REQ):
                    req_rank = comm.recv(source=ANY_SRC, tag=TAG_REQ)
                    if queue:
                        task = queue.pop(0)
                        comm.send(task, dest=req_rank, tag=TAG_TASK)
                        active_tasks += 1
                    else:
                        comm.send(None, dest=req_rank, tag=TAG_TASK)

            while len(results) < total_expected:
                results.append(comm.recv(source=ANY_SRC, tag=TAG_RES))

            self.completed = results
            return results
        else:
            if self.n_workers <= 1:
                while True:
                    comm.send(rank, dest=root, tag=TAG_REQ)
                    task = comm.recv(source=root, tag=ANY_TAG)
                    if task is None:
                        break
                    res = processor(task)
                    comm.send(res, dest=root, tag=TAG_RES)
                return []

            local_pool = ThreadPoolExecutor(max_workers=self.n_workers)
            futures: dict = {}

            def _request_more(n: int) -> None:
                for _ in range(n):
                    comm.send(rank, dest=root, tag=TAG_REQ)

            try:
                _request_more(self.n_workers)
                while True:
                    while comm.Iprobe(source=root, tag=ANY_TAG):
                        task = comm.recv(source=root, tag=ANY_TAG)
                        if task is None:
                            for fut in as_completed(list(futures.keys())):
                                try:
                                    comm.send(fut.result(), dest=root, tag=TAG_RES)
                                except Exception as exc:
                                    tid, st = futures[fut]
                                    comm.send((tid, st, [], 0.0, str(exc)),
                                              dest=root, tag=TAG_RES)
                            return []
                        fut = local_pool.submit(processor, task)
                        futures[fut] = task
                    done = []
                    for fut in list(futures.keys()):
                        if fut.done():
                            try:
                                comm.send(fut.result(), dest=root, tag=TAG_RES)
                            except Exception as exc:
                                tid, st = futures[fut]
                                comm.send((tid, st, [], 0.0, str(exc)),
                                          dest=root, tag=TAG_RES)
                            done.append(fut)
                            comm.send(rank, dest=root, tag=TAG_REQ)
                    for fut in done:
                        futures.pop(fut, None)
                    import time as _t
                    _t.sleep(0.01)
            finally:
                local_pool.shutdown(wait=True)


class CostHistory:
    """Persistent per-entry cost model with exponential moving average (EMA).

    Tracks historical docking costs so cost-aware scheduling improves
    automatically across repeated benchmark campaigns.
    """

    def __init__(self, history_path: Union[str, Path], alpha: float = 0.3):
        self.path = Path(history_path)
        self.alpha = max(0.05, min(0.8, alpha))
        self.data: Dict[str, float] = {}
        self._load()

    def _load(self) -> None:
        if self.path.is_file():
            try:
                self.data = json.loads(self.path.read_text())
            except Exception:
                self.data = {}

    def update(self, target_costs: Dict[str, float]) -> None:
        """Update historical costs using exponential moving average."""
        for key, cost in target_costs.items():
            if key in self.data and self.data[key] > 0:
                self.data[key] = self.alpha * cost + (1 - self.alpha) * self.data[key]
            else:
                self.data[key] = cost
        self._save()

    def get_hints(self) -> Dict[str, float]:
        return dict(self.data)

    def _save(self) -> None:
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            self.path.write_text(json.dumps(self.data, indent=2))
        except Exception:
            pass  # non-fatal

    @classmethod
    def for_dataset(cls, base_results_dir: Path, slug: str, tier: int) -> "CostHistory":
        history_path = base_results_dir / slug / f"tier{tier}" / ".cost_history.json"
        return cls(history_path)


# ---------------------------------------------------------------------------
# DatasetRunner
# ---------------------------------------------------------------------------


class DatasetRunner:
    """Orchestrates FlexAIDdS benchmarks across datasets and compute nodes.

    Discovers dataset YAML configs, dispatches docking jobs (locally or via
    MPI), aggregates metrics, and produces structured reports.

    Args:
        datasets_dir:   Directory containing ``*.yaml`` dataset configs.
        results_dir:    Directory for output files (created if missing).
        binary:         Path to the ``FlexAID`` executable.  Auto-detected
                        from ``FLEXAIDDS_BINARY`` env var or ``$PATH`` if None.
        temperature:    Simulation temperature in Kelvin (default 300).
        n_workers:      Local worker processes for parallel target evaluation
                        (ignored when ``use_mpi=True``).
        use_mpi:        Whether to use MPI for distributed execution.
        cache_dir:      Dataset cache directory; uses ``$FLEXAIDDS_BENCHMARK_DATA``
                        env var if None.
        bootstrap_ci:   Whether to compute 95% bootstrap CIs (slower).
        n_bootstrap:    Number of bootstrap resamples.
        dry_run:        Skip actual docking; use synthetic poses for pipeline
                        smoke tests. docking_power_* is never reported as a real
                        success rate in this mode.
        repo_root:      Root of the FlexAIDdS repository (for git SHA detection).
    """

    def __init__(
        self,
        datasets_dir: Union[str, Path, None] = None,
        results_dir: Union[str, Path] = "results/benchmarks",
        binary: Optional[str] = None,
        temperature: float = 300.0,
        n_workers: int = 1,
        omp_threads: Optional[int] = None,
        use_mpi: bool = False,
        cache_dir: Optional[Union[str, Path]] = None,
        bootstrap_ci: bool = False,
        n_bootstrap: int = 5_000,
        dry_run: bool = False,
        repo_root: Optional[Union[str, Path]] = None,
        resume: bool = False,
        command_line: Optional[str] = None,
        default_conc_M: float = 1.0,  # P3
    ) -> None:
        _default_datasets = Path(__file__).resolve().parent / "datasets"
        self.datasets_dir = Path(datasets_dir) if datasets_dir is not None else _default_datasets
        self.results_dir = Path(results_dir)
        self.binary = binary or os.environ.get("FLEXAIDDS_BINARY") or "FlexAID"
        # _resolve_binary already existed and was correct -- but it was wired
        # only into the *report* (provenance), while the actual invocation used
        # the raw string. Resolve once here so the engine and the report agree.
        # It deliberately does NOT blanket-absolutise: a bare name is left for
        # PATH (subprocess ignores cwd for those), a name containing a separator
        # is made absolute, and an unfindable name is left alone so subprocess
        # raises and the liveness gate counts it.
        self.binary = self._resolve_binary()
        self.temperature = temperature
        self.n_workers = max(1, int(n_workers))
        # OMP threads per worker: explicit > env > auto (aim for 2 on M3 Pro's 5 P-cores).
        if omp_threads is None:
            env_val = os.environ.get("FLEXAIDDS_OMP_THREADS")
            omp_threads = int(env_val) if env_val else max(1, 2 if self.n_workers > 1 else 4)
        self.omp_threads = max(1, int(omp_threads))
        self.use_mpi = use_mpi
        # Same reason as DatasetConfig.data_dir: this root is joined into paths
        # handed to a subprocess running in a tmp dir. Note the fallback literal
        # is relative, so the *default* is the trap, not just a bad override.
        self.cache_dir = Path(
            cache_dir or os.environ.get("FLEXAIDDS_BENCHMARK_DATA", "benchmark_data")
        ).expanduser().resolve()
        self.do_bootstrap = bootstrap_ci
        self.n_bootstrap = n_bootstrap
        self.dry_run = dry_run
        self.repo_root = Path(repo_root) if repo_root else Path.cwd()
        self.resume = bool(resume)
        self.command_line = command_line
        self.default_conc_M = float(default_conc_M)
        # Pose PDBs are written by the engine into a per-ligand TemporaryDirectory
        # and destroyed when that ligand's iteration ends -- before the target
        # finishes, let alone the job.  Every question that needs coordinates
        # rather than scalars (was the search collapsed, do these poses share a
        # conformer, what ΔS did the engine actually compute) is unanswerable
        # from a run's output for that reason alone, and no workflow change can
        # fix it: the CI upload step names RESULTS_DIR, which never held a PDB.
        # Opt-in, because keeping them costs disk on every entry and the default
        # behaviour of every existing run must not change.
        self.keep_poses = os.environ.get("FLEXAIDDS_KEEP_POSES", "").strip() not in ("", "0")

        self.results_dir.mkdir(parents=True, exist_ok=True)

        if use_mpi:
            self._mpi_rank, self._mpi_size, self._mpi_root, self._mpi_comm = (
                _mpi_context()
            )
        else:
            self._mpi_rank, self._mpi_size, self._mpi_root, self._mpi_comm = (
                0, 1, True, None
            )

        # Store original requested temp for reference; effective is per-dataset
        self._requested_temperature = temperature
        self.command_line: Optional[str] = None  # set by CLI for full provenance

        # #326 liveness gate: FlexAID non-zero exits are recorded here (under a
        # lock, since entries dispatch across a local thread pool) so a crashed
        # binary can never read as "no regression". Reset per dataset run.
        self._crash_lock = threading.Lock()
        self._flexaid_crashes = 0
        # int = a real process exit/return code; None = the engine never
        # executed (exec failure), which no completed subprocess.run can produce.
        self._entry_exit_codes: Dict[str, Optional[int]] = {}
        # Poses whose element list disagreed with the reference's.  Per
        # instance and reset per dataset, like the crash tally above it: a
        # process-global would report the previous dataset's refusals as
        # this one's, which is the misattribution this PR exists to close.
        self._element_mismatches: List[str] = []

    def _effective_temperature(self, slug: str) -> float:
        """Auto-force 298 K for ITC/thermo datasets (handoff requirement).

        Non-thermo datasets retain the configured temperature (default 300).
        Always returns the value that should be used for this run.
        """
        if slug in THERMO_DATASETS:
            return THERMO_DEFAULT_TEMP
        return self.temperature

    def _resolve_binary(self) -> str:
        """Return resolved absolute path to FlexAID binary when possible."""
        import shutil
        b = self.binary
        if os.path.isabs(b) and os.path.isfile(b):
            return str(Path(b).resolve())
        found = shutil.which(b)
        if found and os.path.isfile(found):
            return str(Path(found).resolve())
        return b

    # ------------------------------------------------------------------
    # Dataset discovery
    # ------------------------------------------------------------------

    def discover_datasets(self) -> List[DatasetConfig]:
        """Discover and load all ``*.yaml`` configs from ``datasets_dir``.

        Returns:
            List of :class:`DatasetConfig` objects sorted by benchmark order,
            then slug.
        """
        if not self.datasets_dir.is_dir():
            logger.warning("datasets_dir does not exist: %s", self.datasets_dir)
            return []

        configs: List[DatasetConfig] = []
        for yaml_path in sorted(self.datasets_dir.glob("*.yaml")):
            try:
                cfg = DatasetConfig.from_yaml(yaml_path)
                if cfg.data_dir is None:
                    cfg.data_dir = self.cache_dir / cfg.slug
                configs.append(cfg)
                logger.debug("Loaded dataset config: %s (%d targets)",
                             cfg.slug, len(cfg.targets))
            except Exception as exc:
                logger.error("Failed to load %s: %s", yaml_path, exc)

        return sorted(configs, key=lambda cfg: (cfg.benchmark_order, cfg.slug))

    def load_dataset_config(self, yaml_path: Union[str, Path]) -> DatasetConfig:
        """Load a single dataset config from an explicit path."""
        cfg = DatasetConfig.from_yaml(Path(yaml_path))
        if cfg.data_dir is None:
            cfg.data_dir = self.cache_dir / cfg.slug
        return cfg

    # ------------------------------------------------------------------
    # Core docking dispatch
    # ------------------------------------------------------------------

    def _resolve_entry_paths(
        self,
        config: DatasetConfig,
        target_id: str,
        state: str,
    ) -> Tuple[Optional[Path], List[Path]]:
        """Locate receptor PDB and ligand files for a benchmark entry."""
        if config.data_dir is None:
            return None, []
        receptor, ligands = resolve_benchmark_paths(
            config.slug, config.data_dir, target_id, state
        )
        if receptor is None:
            logger.warning(
                "No receptor found for %s (%s) in %s",
                target_id, state, config.data_dir,
            )
        return receptor, ligands

    def _dock_target(
        self,
        target_id: str,
        receptor_path: Path,
        ligand_paths: List[Path],
        structural_state: str = "holo",
        with_entropy: bool = True,
        conc_M: float = 1.0,  # P3
    ) -> List[PoseScore]:
        """Run FlexAIDdS on one target and return scored poses.

        In dry-run mode, returns synthetic poses for framework testing.

        Args:
            target_id:        Target identifier.
            receptor_path:    Receptor PDB file.
            ligand_paths:     Ligand files (Mol2 or SDF).
            structural_state: Receptor structural state.
            with_entropy:     Include TΔS correction.

        Returns:
            List of :class:`PoseScore` objects.
        """
        if self.dry_run:
            return self._synthetic_poses(target_id, ligand_paths, structural_state)

        try:
            return self._run_flexaid(
                target_id,
                receptor_path,
                ligand_paths,
                structural_state,
                with_entropy,
                conc_M=conc_M,
            )
        except Exception as exc:
            logger.error("Docking failed for %s: %s", target_id, exc)
            return []

    def _run_flexaid(
        self,
        target_id: str,
        receptor_path: Path,
        ligand_paths: List[Path],
        structural_state: str,
        with_entropy: bool,
        conc_M: float = 1.0,  # P3
    ) -> List[PoseScore]:
        """Invoke the FlexAID binary and parse output poses."""
        poses: List[PoseScore] = []

        for ligand_path in ligand_paths:
            ligand_id = ligand_path.stem

            with tempfile.TemporaryDirectory(prefix=f"flexaid_{target_id}_") as tmp:
                tmp_path = Path(tmp)

                sub_env = os.environ.copy()
                sub_env["OMP_NUM_THREADS"] = str(self.omp_threads)
                # Disable early termination (stagnation/entropy) during benchmarking
                # so the *full* configured generation budget is always consumed.
                # Spare generations after a would-be early exit are used for
                # additional conformational search (see GABOOM no_sec + exploration boost).
                sub_env["FLEXAIDDS_NO_SEC"] = "1"
                # Signal to core that this is a benchmark run needing equal search effort
                sub_env["FLEXAIDDS_BENCHMARK"] = "1"
                # Wrap ONLY the subprocess call: the exceptions here (timeout,
                # exec failure) are properties of *running the engine*. Pose
                # parsing lives below, deliberately outside this try, so an
                # OSError while reading result files is never miscounted as a
                # liveness crash — the engine ran; only the read failed.
                try:
                    # Direct CLI: <receptor> <ligand> -o prefix  (binary auto-detects files)
                    # Write outputs as flexaid_*.pdb so parser glob *.pdb catches them.
                    result = subprocess.run(
                        [self.binary, str(receptor_path), str(ligand_path), "-o", "flexaid", "--conc", str(conc_M)],
                        capture_output=True,
                        text=True,
                        timeout=3600,
                        cwd=tmp_path,
                        env=sub_env,
                    )
                except subprocess.TimeoutExpired:
                    # #326 liveness: a timeout is the same class as the OSError
                    # below — the engine produced no result — and the rationale
                    # written there applies verbatim: without recording it, "the
                    # run would look like 'executed, 0 poses' (productivity)
                    # rather than 'engine did not run' (liveness)."  That fix was
                    # applied to the OSError sibling and not to this one, so a
                    # timed-out dock left NOTHING in the artifact: no crash count,
                    # no exit code, only this log line.  Record it the same way.
                    #
                    # None, not a numeric sentinel, for the same reason the
                    # OSError branch gives: subprocess encodes signal death as a
                    # negative returncode, so -1 already means "ran, killed by
                    # signal 1".  None is a value no completed subprocess.run can
                    # yield, so in entry_exit_codes it means "did not complete"
                    # and only that.
                    with self._crash_lock:
                        self._flexaid_crashes += 1
                        self._entry_exit_codes[f"{target_id}/{ligand_id}"] = None
                    logger.error("Docking timed out: %s/%s", target_id, ligand_id)
                    continue
                except OSError as exc:
                    # #326 liveness: subprocess.run itself raising (the binary is
                    # missing or not executable -> FileNotFoundError; permission
                    # denied -> PermissionError; both OSError) means the engine
                    # NEVER RAN. That is exactly the case the liveness gate is
                    # named for, yet the returncode-based counter below never
                    # sees it — the exception would otherwise be swallowed
                    # per-entry by _process_one_item and the run would look like
                    # "executed, 0 poses" (productivity) rather than "engine did
                    # not run" (liveness). Record it as a crash so gate 1 fires
                    # even when productivity is relaxed (FLEXAIDDS_BENCH_ALLOW_EMPTY)
                    # or vacuous (a dataset declaring 0 baselines). Record None,
                    # not a numeric sentinel: the engine produced no exit code at
                    # all. A value like -1 would be ambiguous — subprocess encodes
                    # signal death as a negative returncode (SIGHUP -> -1), so -1
                    # already means "ran, killed by signal 1" in the returncode
                    # branch above. None is a value no completed subprocess.run
                    # can yield, so in entry_exit_codes it means "did not execute"
                    # and only that.
                    with self._crash_lock:
                        self._flexaid_crashes += 1
                        self._entry_exit_codes[f"{target_id}/{ligand_id}"] = None
                    logger.error(
                        "FlexAID failed to execute for %s/%s: %s",
                        target_id, ligand_id, exc,
                    )
                    continue

                # The engine ran (exec succeeded); everything below is parsing
                # and is NOT under the exec try — a read error here is not a
                # liveness crash.
                if result.returncode != 0:
                    # #326 liveness: a non-zero exit is a gate failure, not
                    # a debug aside. Record it so the run is scored
                    # INCONCLUSIVE instead of silently passing on 0 poses.
                    stderr_head = (result.stderr or "").strip().splitlines()
                    stderr_head = stderr_head[0][:200] if stderr_head else ""
                    with self._crash_lock:
                        self._flexaid_crashes += 1
                        self._entry_exit_codes[f"{target_id}/{ligand_id}"] = result.returncode
                    logger.warning(
                        "FlexAID non-zero exit %d for %s/%s: %s",
                        result.returncode, target_id, ligand_id, stderr_head,
                    )
                    if result.stderr:
                        logger.debug("stderr: %s", result.stderr[:500])
                parsed = self._parse_flexaid_output(
                    tmp_path,
                    target_id,
                    ligand_id,
                    structural_state,
                    reference_ligand=ligand_path,
                    mismatches=self._element_mismatches,
                )
                # Copy the coordinates out BEFORE the with-block closes and
                # takes them.  Deliberately after the parse, so a preservation
                # failure can never cost us the scalars we already have.
                if self.keep_poses:
                    self._preserve_pose_pdbs(
                        tmp_path, target_id, ligand_id, structural_state
                    )
                # P3: capture grand log_Z from [GRAND] stdout (emitted by C++ cluster hook when --conc used)
                grand_log_z = None
                if result.stdout:
                    for line in result.stdout.splitlines():
                        if '[GRAND]' in line and 'log_Z=' in line:
                            try:
                                grand_log_z = float(line.split('log_Z=')[1].split()[0])
                            except Exception:
                                pass
                if grand_log_z is not None:
                    for p in parsed:
                        p.ensemble_log_Z = grand_log_z
                poses.extend(parsed)

        return poses

    def pose_pdb_dir(self, target_id: str, ligand_id: str, structural_state: str) -> Path:
        """Destination for preserved pose PDBs, ALWAYS under ``results_dir``.

        The destination root is not a detail.  ``benchmark-tier1.yml`` uploads
        exactly ``${{ env.RESULTS_DIR }}/``, so files preserved anywhere else
        survive the run, satisfy any local check, and still reach nobody --
        the same empty-artifact outcome, relocated one step later and harder to
        see.  Tests assert against this method rather than a literal path so
        that "it persists" and "CI can reach it" cannot drift apart.
        """
        safe_target = re.sub(r"[^A-Za-z0-9_.-]", "_", target_id)
        safe_ligand = re.sub(r"[^A-Za-z0-9_.-]", "_", ligand_id)
        safe_state = re.sub(r"[^A-Za-z0-9_.-]", "_", structural_state)
        return self.results_dir / "poses" / safe_target / f"{safe_ligand}_{safe_state}"

    def _preserve_pose_pdbs(
        self,
        work_dir: Path,
        target_id: str,
        ligand_id: str,
        structural_state: str,
    ) -> List[Path]:
        """Copy pose PDBs out of the engine's temp dir into ``results_dir``.

        Returns the files written.  Never raises: preserving coordinates is a
        diagnostic convenience, and it must not be able to fail a docking run
        that has already produced its scores.
        """
        copied: List[Path] = []
        try:
            sources = sorted(work_dir.glob("*.pdb"))
            if not sources:
                return copied
            dest_dir = self.pose_pdb_dir(target_id, ligand_id, structural_state)
            dest_dir.mkdir(parents=True, exist_ok=True)
            for src in sources:
                dest = dest_dir / src.name
                shutil.copy2(src, dest)
                copied.append(dest)
        except OSError as exc:
            logger.warning(
                "Could not preserve pose PDBs for %s/%s: %s",
                target_id, ligand_id, exc,
            )
        return copied

    @staticmethod
    def _parse_flexaid_output(
        work_dir: Path,
        target_id: str,
        ligand_id: str,
        structural_state: str,
        reference_ligand: Optional[Path] = None,
        mismatches: Optional[List[str]] = None,
    ) -> List[PoseScore]:
        """Parse FlexAIDdS result PDB files from a completed docking run.

        Accepts both legacy ``REMARK CF_SCORE:`` headers and current
        ``REMARK CF=`` / ``REMARK <rmsd> RMSD to ref. structure`` output.
        Computes post-hoc RMSD vs the reference ligand when REMARK omits it.
        """
        poses: List[PoseScore] = []
        pdb_files = sorted(
            p for p in work_dir.glob("*.pdb")
            if (p.name.startswith("flexaid") or p.name.startswith("FlexAID"))
            and not _is_initial_pose_file(p)
        )
        if not pdb_files:
            pdb_files = sorted(
                p for p in work_dir.glob("*.pdb") if not _is_initial_pose_file(p)
            )

        # A contaminated reference is a refusal, not a docking failure, and it
        # must not reach the caller's broad `except Exception` looking like one.
        # Without this the raise escapes the per-file try below, degrades to
        # "this target produced nothing", and is indistinguishable from bad
        # docking -- the exact shape #366 exists to remove.
        ref_coords = None
        ref_elements = None
        if reference_ligand:
            try:
                ref_coords, ref_elements = _reference_ligand_atoms(
                    reference_ligand)
            except ValueError as exc:
                logger.warning(
                    "Reference ligand refused for %s/%s: %s  No RMSD will be "
                    "computed for this target; it is scored as a miss.  This "
                    "is a REFERENCE FILE problem, not a docking result.",
                    target_id, ligand_id, exc,
                )

        for rank, pdb_path in enumerate(pdb_files, start=1):
            rmsd = -1.0
            enthalpy_score = 0.0
            entropy_correction = 0.0
            total_score = 0.0
            is_active = False
            exp_affinity: Optional[float] = None

            try:
                lines = pdb_path.read_text().splitlines()
                for line in lines:
                    if not line.startswith("REMARK"):
                        continue
                    if "RMSD:" in line:
                        rmsd = _parse_remark_float(line, "RMSD:")
                    elif "CF_SCORE:" in line:
                        enthalpy_score = _parse_remark_float(line, "CF_SCORE:")
                    elif "ENTROPY:" in line:
                        entropy_correction = _parse_remark_float(line, "ENTROPY:")
                    elif "TOTAL_SCORE:" in line:
                        total_score = _parse_remark_float(line, "TOTAL_SCORE:")
                    elif "EXP_AFFINITY:" in line:
                        exp_affinity = _parse_remark_float(line, "EXP_AFFINITY:")
                    elif "ACTIVE:1" in line:
                        is_active = True
                    elif "ENSEMBLE_LOG_Z:" in line or "GRAND_LOG_Z:" in line:
                        try:
                            grand_log_z = float( line.split(":")[-1].strip() )
                        except:
                            pass
                    else:
                        m = _RMSD_REMARK_RE.search(line)
                        if m:
                            rmsd = float(m.group(1))
                        m = _CF_REMARK_RE.search(line)
                        if m and enthalpy_score == 0.0:
                            enthalpy_score = float(m.group(1))
                        m = _CF_APP_REMARK_RE.search(line)
                        if m and enthalpy_score == 0.0:
                            enthalpy_score = float(m.group(1))
                        m = _ENTROPY_REMARK_RE.search(line)
                        if m and entropy_correction == 0.0:
                            entropy_correction = float(m.group(1))

                if rmsd < 0.0 and ref_coords is not None:
                    rmsd = _pose_rmsd_vs_reference(
                        pdb_path, ref_coords, ref_elements, mismatches)

            except Exception as exc:
                logger.debug("Error parsing %s: %s", pdb_path, exc)
                continue

            if total_score == 0.0 and enthalpy_score != 0.0:
                total_score = enthalpy_score - entropy_correction

            poses.append(
                PoseScore(
                    target_id=target_id,
                    ligand_id=ligand_id,
                    pose_rank=rank,
                    rmsd=rmsd,
                    enthalpy_score=enthalpy_score,
                    entropy_correction=entropy_correction,
                    total_score=total_score,
                    is_active=is_active,
                    exp_affinity=exp_affinity,
                    structural_state=structural_state,
                )
            )

        return poses

    @staticmethod
    def _synthetic_poses(
        target_id: str,
        ligand_paths: List[Path],
        structural_state: str,
    ) -> List[PoseScore]:
        """Generate synthetic poses for dry-run / framework testing."""
        import random
        rng = random.Random(hash(target_id) & 0xFFFFFFFF)
        poses: List[PoseScore] = []
        for lig_path in ligand_paths:
            ligand_id = lig_path.stem
            for rank in range(1, 6):
                enthalpy = rng.uniform(-12, -4)
                entropy = rng.uniform(0.5, 3.0)
                total = enthalpy - entropy
                # Synthetic near-native pose at rank 2 sometimes
                rmsd = rng.uniform(0.5, 1.5) if rank == 2 else rng.uniform(1.0, 5.0)
                poses.append(
                    PoseScore(
                        target_id=target_id,
                        ligand_id=ligand_id,
                        pose_rank=rank,
                        rmsd=rmsd,
                        enthalpy_score=enthalpy,
                        entropy_correction=entropy,
                        total_score=total,
                        is_active=rng.random() < 0.3,
                        exp_affinity=rng.uniform(-12, -6),
                        ensemble_log_Z = rng.uniform(5, 15),  # synthetic for grand test
                        structural_state=structural_state,
                    )
                )
        return poses

    # ------------------------------------------------------------------
    # Per-dataset run
    # ------------------------------------------------------------------

    def run_dataset(
        self,
        config: DatasetConfig,
        tier: int = 2,
        metric_subset: Optional[List[str]] = None,
        structural_states: Optional[List[str]] = None,
    ) -> DatasetResult:
        """Run benchmarks for one dataset.

        Now uses automated per-entry (per-target) saving and processing.
        Individual TargetResult files are written atomically as soon as each
        (target, structural_state) finishes. This enables resume and fine-grained
        progress / resource monitoring.

        A lightweight EntryTaskManager (see below) coordinates the work items.
        """
        t0 = time.monotonic()
        targets = config.scheduled_targets(tier)
        states = structural_states or config.structural_states
        scheduled_items = config.scheduled_work_items(tier)
        slug = config.slug
        eff_temp = self._effective_temperature(slug)
        if slug in THERMO_DATASETS and eff_temp != self.temperature:
            logger.info(
                "ITC/thermo dataset %s: auto-forcing temperature=%.1f K "
                "(bulletproof per handoff; was %.1f)",
                slug, eff_temp, self.temperature,
            )

        requested_metrics = metric_subset or config.metrics or None
        if slug in THERMO_DATASETS:
            # Strictly limit to scoring + entropy_rescue (docking_power OFF for thermo)
            if requested_metrics:
                filtered = [
                    m for m in requested_metrics
                    if m in THERMO_METRICS or (m.startswith("scoring_power_") and "docking" not in m.lower())
                ]
                if set(filtered) != set(requested_metrics or []):
                    logger.info(
                        "Thermo dataset %s: restricting metrics to scoring_power_* + entropy_rescue_rate only",
                        slug,
                    )
                    requested_metrics = filtered or list(THERMO_METRICS)
            else:
                requested_metrics = list(THERMO_METRICS)

        if self.dry_run:
            # Synthetic poses must never produce docking_power success rates.
            before = list(requested_metrics) if requested_metrics is not None else None
            requested_metrics = _filter_requested_metrics_for_dry_run(requested_metrics)
            if before is not None and before != requested_metrics:
                logger.info(
                    "Dry-run: omitting docking_power_* metrics (synthetic poses are not real docking results)"
                )

        dr = DatasetResult(
            config=config,
            tier=tier,
            targets_attempted=list(targets),
            timestamp=datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(self.repo_root),
            host=socket.gethostname(),
            binary=self._resolve_binary(),
            temperature=eff_temp,
            full_command=getattr(self, "command_line", None) or " ".join(sys.argv),
            dry_run=bool(self.dry_run),
            metrics_note=DRY_RUN_METRICS_NOTE if self.dry_run else "",
        )

        if not targets:
            logger.warning("Dataset %s has no targets", config.slug)
            dr.duration_seconds = 0.0
            return dr

        # --- NEW: Per-entry resume + individual processing automation ---
        already_completed: set[str] = set()
        if self.resume:
            already_completed = self._discover_completed_targets(config, tier)
            if already_completed:
                logger.info(
                    "Resume mode: %d/%d targets already have complete per-entry results — skipping them",
                    len(already_completed), len(targets)
                )

        # Build work items from catalog (large-N tier-2) or targets × states
        all_work_items: List[Tuple[str, str]] = []
        for tid, st in scheduled_items:
            if tid in already_completed:
                continue
            all_work_items.append((tid, st))

        if config._uses_crossdock_catalog(tier):
            logger.info(
                "Large-N dataset %s tier-%d: %d scheduled crossdock work items",
                config.slug, tier, len(scheduled_items),
            )

        # For stronger dynamic MPI master-worker, root sees the full remaining list.
        # Non-roots will participate as workers when the manager detects MPI.
        if self.use_mpi and self._mpi_size > 1:
            work_for_manager = all_work_items if self._mpi_root else []
        else:
            work_for_manager = _split_targets(all_work_items, self._mpi_rank, self._mpi_size)

        logger.info(
            "[rank %d/%d] Dataset %s: preparing %d work items for EntryTaskManager (dynamic MPI master-worker when applicable)",
            self._mpi_rank, self._mpi_size, config.slug, len(work_for_manager),
        )

        # Master entry manager — with MPI context + persistent CostHistory (EMA)
        cost_hints = {}
        cost_history = None
        try:
            cost_history = CostHistory.for_dataset(self.results_dir, config.slug, tier)
            cost_hints = cost_history.get_hints()

            # Also layer in fresh manifest costs on resume (they take precedence)
            if self.resume:
                prev_manifest = self._entry_results_dir(config, tier) / "_entry_manifest.json"
                fresh = EntryTaskManager.load_cost_hints_from_manifest(prev_manifest)
                cost_hints.update(fresh)
        except Exception:
            pass

        entry_manager = EntryTaskManager(
            work_items=work_for_manager,
            n_workers=self.n_workers if self._mpi_size == 1 else 1,
            omp_threads=self.omp_threads,
            mpi_comm=self._mpi_comm,
            mpi_rank=self._mpi_rank,
            mpi_size=self._mpi_size,
            mpi_root=self._mpi_root,
            cost_hints=cost_hints,
        )

        all_poses: List[PoseScore] = []
        completed_targets: set[str] = set(already_completed)
        failed_targets: set[str] = set()

        # #326: reset per-dataset FlexAID crash tally before this run's entries.
        with self._crash_lock:
            self._flexaid_crashes = 0
            self._entry_exit_codes = {}
            self._element_mismatches = []

        def _process_one_item(item: Tuple[str, str]) -> Tuple[str, str, List[PoseScore], float, str]:
            """Process a single (target_id, structural_state) entry.

            All exceptions are caught and returned as an error string so that one
            bad target never kills the entire campaign (in both serial and parallel
            dispatch modes).
            """
            target_id, state = item
            t_start = time.monotonic()
            try:
                error = ""
                receptor = None
                ligands: List[Path] = []

                if config.data_dir and not self.dry_run:
                    receptor, ligands = self._resolve_entry_paths(config, target_id, state)
                    if receptor is None:
                        error = f"No receptor found for {target_id}/{state}"
                        return target_id, state, [], time.monotonic() - t_start, error
                    if not ligands:
                        error = f"No ligand found for {target_id}/{state}"
                        return target_id, state, [], time.monotonic() - t_start, error

                # P3: per-ligand conc from config.ligand_concs if present (from yaml)
                eff_conc = getattr(config, 'default_conc_M', self.default_conc_M)
                if ligands:
                    lid = ligands[0].stem
                    eff_conc = getattr(config, 'ligand_concs', {}).get(lid, eff_conc)
                poses = self._dock_target(
                    target_id,
                    receptor or Path("/dev/null"),
                    ligands or [Path(f"{target_id}.mol2")],
                    structural_state=state,
                    conc_M=eff_conc,
                )

                elapsed = time.monotonic() - t_start
                cost_cpu = elapsed * max(1, self.omp_threads)

                tr = TargetResult(
                    target_id=target_id,
                    structural_state=state,
                    poses=poses,
                    duration_seconds=elapsed,
                    error=error,
                    conc_M=eff_conc,
                )
                if poses and getattr(poses[0], 'ensemble_log_Z', None) is not None:
                    tr.grand_log_Z = poses[0].ensemble_log_Z
                    try:
                        from ..grand_canonical import compute_grand_partition
                        g = compute_grand_partition([(target_id, tr.grand_log_Z, tr.conc_M)], temperature_K=298.0)
                        tr.grand_xi = g.log_Xi()
                    except Exception as e:
                        logger.debug("P3 grand compute: %s", e)
                try:
                    saved_path = self._save_target_result(tr, config, tier, cost_cpu=cost_cpu)
                    logger.debug(
                        "Saved per-entry result: %s (cost~%.1fs CPU)", saved_path, cost_cpu
                    )
                except Exception as save_exc:
                    logger.warning(
                        "Failed to save per-entry result for %s/%s: %s",
                        target_id, state, save_exc,
                    )

                return target_id, state, poses, elapsed, error

            except Exception as exc:
                elapsed = time.monotonic() - t_start
                logger.error(
                    "Unhandled exception for %s/%s: %s", target_id, state, exc,
                    exc_info=True,
                )
                return target_id, state, [], elapsed, f"EXCEPTION: {exc}"

        # Dispatch via the EntryTaskManager (local ThreadPool or sequential)
        results = entry_manager.run(_process_one_item)

        for target_id, state, poses, elapsed, error in results:
            all_poses.extend(poses)
            if error:
                failed_targets.add(target_id)
                logger.warning("Failed %s/%s: %s", target_id, state, error)
            else:
                completed_targets.add(target_id)

            logger.debug(
                "Entry %s/%s: %d poses in %.1fs",
                target_id, state, len(poses), elapsed,
            )

        # P3: emit grand canonical summary for ligands that have grand_log_Z (e.g. from multi-ligand on receptor or grand datasets)
        # Group by target (ligand), use conc from config, compute Xi, p etc using the module.
        grand_summary: Dict[str, Dict[str, float]] = {}
        for target_id, state, poses, elapsed, error in results:
            if not error and poses:
                logz = getattr(poses[0], 'ensemble_log_Z', None)
                if logz is not None:
                    conc = getattr(config, 'ligand_concs', {}).get(target_id, getattr(config, 'default_conc_M', 1.0))
                    try:
                        from ..grand_canonical import compute_grand_partition
                        g = compute_grand_partition([(target_id, logz, conc)], temperature_K=298.0)
                        p_bind = g.binding_probability(target_id) if hasattr(g, 'binding_probability') else (math.exp(logz + math.log(conc)) / math.exp(g.log_Xi()) if conc > 0 else 0)
                        grand_summary[target_id] = {
                            'log_Z': logz,
                            'conc_M': conc,
                            'log_Xi': g.log_Xi(),
                            'p_bind': p_bind,
                        }
                    except Exception as e:
                        logger.debug("P3 grand summary compute skip for %s: %s", target_id, e)
        if grand_summary:
            dr.grand_summary = grand_summary
            logger.info("P3: emitted grand summary for %d ligands", len(grand_summary))

        completed = sorted(completed_targets)
        failed = sorted(failed_targets)

        # #326: snapshot this rank's FlexAID crash tally + count of targets
        # actually executed this run. `results` is the dispatched work list,
        # which excludes --resume checkpoints (skipped before dispatch), so it
        # is exactly "fresh work" — the denominator gates 2-3 should judge.
        with self._crash_lock:
            crashes = self._flexaid_crashes
            exit_codes = dict(self._entry_exit_codes)
        newly = len(results)
        resumed = len(already_completed)

        # Element-order refusals are scored as misses, so they must be visible
        # or a correspondence bug reads as bad docking.
        if self._element_mismatches:
            logger.warning(
                "%d pose(s) had an element list disagreeing with their "
                "reference and were scored as misses: %s%s",
                len(self._element_mismatches),
                ", ".join(
                    Path(p).name for p in self._element_mismatches[:5]),
                " ..." if len(self._element_mismatches) > 5 else "",
            )

        # MPI gather (still works at target granularity for final aggregation)
        if self._mpi_comm is not None:
            all_results_by_rank = self._mpi_comm.gather(
                (all_poses, completed, failed, crashes, exit_codes, newly, resumed), root=0
            )
            if self._mpi_root:
                all_poses = []
                completed = []
                failed = []
                crashes = 0
                exit_codes = {}
                newly = 0
                resumed = 0
                for poses_i, comp_i, fail_i, crash_i, codes_i, newly_i, resumed_i in (all_results_by_rank or []):
                    all_poses.extend(poses_i)
                    completed.extend(comp_i)
                    failed.extend(fail_i)
                    crashes += crash_i
                    exit_codes.update(codes_i)
                    newly += newly_i
                    # `resumed` is REPLICATED, not partitioned: every rank runs
                    # _discover_completed_targets over the same shared disk, so
                    # each reports the same count. Summing would give
                    # resumed × n_ranks. Take the common value (max == that
                    # value) so the report-JSON provenance stays honest. The
                    # gate only asks `resumed == 0`, which summing preserved,
                    # but a wrong number in the artifact is the exact thing this
                    # whole change exists to prevent.
                    resumed = max(resumed, resumed_i)

        # Root finalizes
        if self._mpi_root:
            dr.targets_completed = completed
            dr.targets_failed = failed
            dr.flexaid_crashes = crashes
            dr.entry_exit_codes = exit_codes
            dr.total_poses = len(all_poses)
            dr.newly_executed = newly
            dr.resumed = resumed

            if all_poses:
                metrics = compute_all_metrics(all_poses, requested=requested_metrics)
                if self.dry_run:
                    # Safety net: never surface docking_power_* from synthetic poses.
                    metrics = _strip_docking_power_metrics(metrics)
                dr.metrics = metrics

                if self.do_bootstrap:
                    cis = self._compute_bootstrap_cis(all_poses, requested_metrics)
                    if self.dry_run:
                        cis = _strip_docking_power_metrics(cis)
                    dr.ci_95 = cis

            if not self.dry_run:
                dr.check_regressions()
            else:
                logger.info(
                    "Dry-run mode — skipping regression checks; docking_power_* omitted "
                    "(synthetic poses are not production docking rates)"
                )

            # Also write a small per-dataset manifest of individual entry status (for audit + reproducibility)
            self._write_entry_manifest(config, tier, completed, failed)

            # Update persistent cost history (CostHistory + EMA)
            if cost_history is not None:
                try:
                    observed_costs = {}
                    for target_id, state, poses, elapsed, error in results:
                        if not error and elapsed > 0:
                            key = f"{target_id}_{state}"
                            observed_costs[key] = elapsed * max(1, self.omp_threads)
                    if observed_costs:
                        cost_history.update(observed_costs)
                        logger.info("CostHistory updated with %d new observations", len(observed_costs))
                except Exception:
                    pass

        dr.duration_seconds = time.monotonic() - t0
        return dr

    def _compute_bootstrap_cis(
        self,
        poses: List[PoseScore],
        requested: Optional[List[str]],
    ) -> Dict[str, Tuple[float, float]]:
        """Compute bootstrap CIs for scalar pose-level metrics."""
        cis: Dict[str, Tuple[float, float]] = {}

        def _rescue_fn(sample):
            return entropy_rescue_rate(sample)

        def _dock_fn(sample):
            return docking_power(sample, top_n=1)

        fns = {
            "entropy_rescue_rate": _rescue_fn,
            "docking_power_top1": _dock_fn,
        }
        for name, fn in fns.items():
            # Dry-run never bootstraps docking_power (synthetic poses ≠ real rates).
            if self.dry_run and _is_docking_power_metric(name):
                continue
            if requested is None or name in requested:
                lo, hi = bootstrap_ci(fn, poses, n_resamples=self.n_bootstrap)
                cis[name] = (lo, hi)

        return cis

    # ------------------------------------------------------------------
    # Run all datasets
    # ------------------------------------------------------------------

    def run_all(
        self,
        datasets: Optional[List[str]] = None,
        tier: int = 2,
        distributed: bool = False,
        n_nodes: int = 1,
        metric_subset: Optional[List[str]] = None,
    ) -> BenchmarkReport:
        """Run benchmarks across all (or selected) datasets.

        Args:
            datasets:      Dataset slugs to run; None = all discovered.
            tier:          Benchmark tier (1 = fast, 2 = full).
            distributed:   Enable MPI-distributed execution.
            n_nodes:       Number of MPI nodes (used for logging only; actual
                           distribution is controlled by ``mpirun``).
            metric_subset: Restrict metrics to this list.

        Returns:
            :class:`BenchmarkReport` containing all dataset results.
        """
        all_configs = self.discover_datasets()
        if datasets:
            slugs = set(datasets)
            all_configs = [c for c in all_configs if c.slug in slugs]
            missing = slugs - {c.slug for c in all_configs}
            if missing:
                logger.warning("Datasets not found: %s", ", ".join(sorted(missing)))

        if not all_configs:
            logger.error("No datasets to run.")
            return BenchmarkReport(
                generated_at=datetime.datetime.utcnow().isoformat() + "Z",
                git_sha=_git_sha(self.repo_root),
                host=socket.gethostname(),
                runner_info=_runner_info(),
            )

        if self._mpi_root:
            logger.info(
                "DatasetRunner: %d dataset(s) · tier %d · %d MPI rank(s) · temperature=%s K (exact for best BindingMode thermo ledger)",
                len(all_configs), tier, self._mpi_size, self.temperature,
            )

        results: List[DatasetResult] = []
        for config in all_configs:
            logger.info("Running dataset: %s", config.slug)
            dr = self.run_dataset(config, tier=tier, metric_subset=metric_subset)
            results.append(dr)

            # Save incremental results after each dataset
            if self._mpi_root:
                self._save_dataset_result(dr)

        report = BenchmarkReport(
            datasets=results,
            generated_at=datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(self.repo_root),
            host=socket.gethostname(),
            runner_info=_runner_info(),
            temperature=results[0].temperature if results else self.temperature,
            binary=results[0].binary if results else self._resolve_binary(),
            full_command=getattr(self, "command_line", None) or "",
        )
        return report

    def _save_dataset_result(self, dr: DatasetResult) -> None:
        """Write a per-dataset JSON result file as soon as it's ready.
        P3: also emit grand_summary.csv when grand canonical data present (Ξ, p_bind, concs etc).
        """
        out_path = self.results_dir / f"{dr.config.slug}_tier{dr.tier}.json"
        out_path.write_text(json.dumps(dr.to_dict(), indent=2))
        logger.info("Dataset result saved: %s", out_path)
        if getattr(dr, "grand_summary", None):
            csv_path = self.results_dir / f"{dr.config.slug}_tier{dr.tier}_grand_summary.csv"
            try:
                rows = []
                for lid, info in sorted(dr.grand_summary.items()):
                    row = {"ligand": lid}
                    row.update({k: info.get(k) for k in ("log_Z", "conc_M", "log_Xi", "p_bind") if k in info})
                    rows.append(row)
                if rows:
                    with open(csv_path, "w", newline="") as fh:
                        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
                        w.writeheader()
                        w.writerows(rows)
                    logger.info("P3 grand CSV saved: %s", csv_path)
            except Exception as e:
                logger.debug("P3 grand CSV emission skipped: %s", e)

    # ------------------------------------------------------------------
    # Per-entry (per-target) persistence — the new automated model
    # ------------------------------------------------------------------

    def _entry_results_dir(self, config: DatasetConfig, tier: int) -> Path:
        """Directory layout for individual entry results:
        results/<slug>/tier<tier>/
        """
        d = self.results_dir / config.slug / f"tier{tier}"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def _target_result_path(
        self, config: DatasetConfig, tier: int, target_id: str, structural_state: str
    ) -> Path:
        safe_id = target_id.replace("/", "_")
        return self._entry_results_dir(config, tier) / f"{safe_id}_{structural_state}.json"

    def _save_target_result(self, tr: TargetResult, config: DatasetConfig, tier: int, cost_cpu: float = 0.0) -> Path:
        """Atomic write of one TargetResult (crash-safe). Includes simple cost tracking."""
        path = self._target_result_path(config, tier, tr.target_id, tr.structural_state)
        tmp = path.with_suffix(".json.tmp")
        payload = {
            "target_id": tr.target_id,
            "structural_state": tr.structural_state,
            "duration_seconds": tr.duration_seconds,
            "cost_cpu_seconds": round(cost_cpu, 2),
            "error": tr.error,
            "poses": [p.to_dict() if hasattr(p, "to_dict") else p.__dict__ for p in tr.poses],
            "success": tr.success,
            "grand_log_Z": getattr(tr, 'grand_log_Z', None),
            "grand_xi": getattr(tr, 'grand_xi', None),
            "conc_M": getattr(tr, 'conc_M', 1.0),
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
        }
        tmp.write_text(json.dumps(payload, indent=2))
        os.replace(tmp, path)  # atomic on POSIX + Windows
        return path

    def _load_target_result(
        self, config: DatasetConfig, tier: int, target_id: str, structural_state: str
    ) -> Optional["TargetResult"]:
        path = self._target_result_path(config, tier, target_id, structural_state)
        if not path.is_file():
            return None
        try:
            data = json.loads(path.read_text())
            poses = []
            for p in data.get("poses", []):
                if isinstance(p, dict):
                    # Best-effort reconstruction; fall back to storing raw dict if fields mismatch
                    try:
                        poses.append(PoseScore(**p))  # type: ignore[call-arg]
                    except Exception:
                        poses.append(p)  # type: ignore[arg-type]
                else:
                    poses.append(p)
            return TargetResult(
                target_id=data["target_id"],
                structural_state=data["structural_state"],
                poses=poses,
                duration_seconds=data.get("duration_seconds", 0.0),
                error=data.get("error", ""),
                grand_log_Z=data.get("grand_log_Z"),
                conc_M=data.get("conc_M", 1.0),
                grand_xi=data.get("grand_xi"),
            )
        except Exception as e:
            logger.warning("Corrupt target result %s — will re-run: %s", path, e)
            return None

    def _discover_completed_targets(self, config: DatasetConfig, tier: int) -> set[str]:
        """Return set of target_ids that have at least one successful result for every requested state."""
        states = config.structural_states
        completed = set()
        for target_id in config.targets:
            all_states_ok = True
            for st in states:
                tr = self._load_target_result(config, tier, target_id, st)
                if not (tr and tr.success):
                    all_states_ok = False
                    break
            if all_states_ok:
                completed.add(target_id)
        return completed

    def _write_entry_manifest(
        self,
        config: DatasetConfig,
        tier: int,
        completed: List[str],
        failed: List[str],
        raw_results: Optional[List[Tuple[str, str, List[PoseScore], float, str]]] = None,
    ) -> None:
        """Write rich per-entry manifest with timing and cost tracking.

        This is the enhanced version delivering automated per-target timing/cost
        data for manifests (audit, planning, reproducibility).
        """
        manifest_path = self._entry_results_dir(config, tier) / "_entry_manifest.json"

        # Build timing/cost data by scanning the individual result files we just wrote.
        # This is robust across MPI (each rank wrote its own) and resume runs.
        entry_dir = self._entry_results_dir(config, tier)
        per_entry: Dict[str, float] = {}
        per_entry_cost: Dict[str, float] = {}
        durations: List[float] = []

        for jf in sorted(entry_dir.glob("*_*.json")):
            if jf.name.startswith("_"):
                continue
            try:
                data = json.loads(jf.read_text())
                key = f"{data.get('target_id')}_{data.get('structural_state')}"
                dur = float(data.get("duration_seconds", 0.0))
                cost = float(data.get("cost_cpu_seconds", dur * max(1, self.omp_threads)))
                per_entry[key] = round(dur, 2)
                per_entry_cost[key] = round(cost, 2)
                if dur > 0:
                    durations.append(dur)
            except Exception:
                continue

        # Compute summary statistics (stdlib only, no extra imports)
        n = len(durations)
        total_wall = sum(durations)
        mean = total_wall / n if n > 0 else 0.0
        median = sorted(durations)[n // 2] if n > 0 else 0.0
        p90 = sorted(durations)[int(n * 0.9)] if n > 0 else 0.0
        total_cost = sum(per_entry_cost.values())

        data = {
            "dataset": config.slug,
            "tier": tier,
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
            "omp_threads": self.omp_threads,
            "n_workers_used": self.n_workers,
            "completed": completed,
            "failed": failed,
            "total_attempted": len(completed) + len(failed),
            # --- NEW: per-target timing + cost tracking (user priority #2 first) ---
            "timings": {
                "per_entry_wall_seconds": per_entry,
                "per_entry_cost_cpu_seconds": per_entry_cost,
                "summary": {
                    "num_timed_entries": n,
                    "total_entry_wall_seconds": round(total_wall, 2),
                    "mean_entry_seconds": round(mean, 2),
                    "median_entry_seconds": round(median, 2),
                    "p90_entry_seconds": round(p90, 2),
                    "estimated_total_cost_cpu_seconds": round(total_cost, 2),
                    "cost_model": "wall_time * omp_threads (simple upper-bound CPU-second estimate)",
                },
            },
        }
        manifest_path.write_text(json.dumps(data, indent=2))
        logger.info("Per-entry manifest with timing/cost written: %s", manifest_path)

    # ------------------------------------------------------------------
    # Public convenience helpers
    # ------------------------------------------------------------------

    def run_single(
        self,
        dataset_slug: str,
        tier: int = 2,
        metric: Optional[str] = None,
    ) -> DatasetResult:
        """Run a single named dataset.

        Args:
            dataset_slug: Dataset slug (matches YAML filename stem).
            tier:         Benchmark tier.
            metric:       Run only this metric (None = all).

        Returns:
            :class:`DatasetResult`.

        Raises:
            FileNotFoundError: When the YAML config is not found.
        """
        yaml_path = self.datasets_dir / f"{dataset_slug}.yaml"
        if not yaml_path.is_file():
            raise FileNotFoundError(
                f"No config found for dataset '{dataset_slug}' at {yaml_path}"
            )
        config = self.load_dataset_config(yaml_path)
        metrics = [metric] if metric else None
        return self.run_dataset(config, tier=tier, metric_subset=metrics)

    def generate_report(self, results: List[DatasetResult]) -> Tuple[dict, str]:
        """Generate JSON dict and Markdown string from a list of results.

        Args:
            results: List of completed :class:`DatasetResult`.

        Returns:
            ``(json_dict, markdown_str)``
        """
        report = BenchmarkReport(
            datasets=results,
            generated_at=datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(self.repo_root),
            host=socket.gethostname(),
            runner_info=_runner_info(),
            temperature=results[0].temperature if results else self.temperature,
            binary=results[0].binary if results else self._resolve_binary(),
            full_command=getattr(self, "command_line", None) or "",
        )
        return report.to_dict(), report.to_markdown()


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def _parse_remark_float(line: str, key: str) -> float:
    """Extract a float value following ``key`` in a REMARK line."""
    try:
        idx = line.index(key) + len(key)
        return float(line[idx:].split()[0])
    except (ValueError, IndexError):
        return 0.0


def _is_initial_pose_file(pdb_path: Path) -> bool:
    """True for FlexAID's ``flexaid_INI.pdb`` starting structure.

    FlexAID writes the *input* ligand placement alongside the docked results as
    ``flexaid_INI.pdb`` (``REMARK initial structure``).  It is not a docking
    result and must never enter the pose list.

    On a self-docking benchmark the input IS the crystal pose, so this file
    scores ~0 A RMSD against the reference — measured at 0.0320 A on 1gpk, and
    identical in every sweep cell because it never depends on the search.  Two
    consequences if it is globbed in:

    * it silently drags ``mean_rmsd``/``median_rmsd`` down (1gpk F=130:
      3.841 A over the ten real poses vs 3.495 A once the 0.032 A input is
      counted as an eleventh), and the pollution is constant across cells, so
      the aggregate looks frozen while ``docking_power`` moves;
    * more seriously it is a latent FALSE POSITIVE for ``docking_power`` — a
      near-0 A entry that counts as a success whenever its score ranks inside
      top-N.  It stayed out of top-N in the runs measured so far only because
      its CF was worse than the docked poses, which is luck, not a guarantee.
    """
    return pdb_path.stem.upper().endswith("_INI")


def _reference_ligand_atoms(ligand_path: Path):
    """Heavy-atom coordinates AND elements from a reference ligand file.

    The elements are needed to type the symmetry-corrected RMSD.  Typing the
    reference with the *pose's* element list would assume the two files list
    atoms in the same order -- and #354 exists precisely because FlexAID pose
    order and SDF file order are not guaranteed to agree.
    """
    from flexaidds.benchmark import extract_ligand_atoms_from_pdb

    suffix = ligand_path.suffix.lower()
    if suffix == ".pdb":
        # Reference PDBs can carry cofactors too.  No count is available here
        # (this IS the count source), so the extractor must refuse to union
        # rather than silently merge a cofactor into the reference (#366).
        return extract_ligand_atoms_from_pdb(ligand_path)
    if suffix in {".sdf", ".mol"}:
        return _extract_ligand_atoms_from_sdf(ligand_path)
    if suffix == ".mol2":
        from flexaidds.dataset_adapters import parse_mol2_atoms

        atoms = parse_mol2_atoms(str(ligand_path))
        if not atoms:
            raise ValueError(f"No atoms in {ligand_path}")
        import numpy as np

        return (np.array([[a.x, a.y, a.z] for a in atoms], dtype=np.float64),
                [str(getattr(a, "element", "") or "").upper() for a in atoms])
    raise ValueError(f"Unsupported reference ligand format: {ligand_path}")


def _reference_ligand_coords(ligand_path: Path):
    """Coordinates only.  See :func:`_reference_ligand_atoms`."""
    return _reference_ligand_atoms(ligand_path)[0]


def _extract_ligand_atoms_from_sdf(sdf_path: Path):
    """Parse heavy-atom XYZ rows and element symbols from a V2000 SD file."""
    import numpy as np

    lines = sdf_path.read_text().splitlines()
    counts_line = None
    for i, line in enumerate(lines):
        if "V2000" in line or "V3000" in line:
            counts_line = i
            break
    if counts_line is None:
        raise ValueError(f"Cannot parse SDF atom block in {sdf_path}")

    n_atoms = int(lines[counts_line][0:3].strip())
    coords = []
    elements = []
    for line in lines[counts_line + 1: counts_line + 1 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            continue
        element = parts[3] if len(parts) > 3 else ""
        if element.upper() == "H":
            continue
        coords.append((float(parts[0]), float(parts[1]), float(parts[2])))
        elements.append(element.upper())
    if not coords:
        raise ValueError(f"No heavy atoms in {sdf_path}")
    return np.array(coords, dtype=np.float64), elements


def _extract_ligand_coords_from_sdf(sdf_path: Path):
    """Coordinates only.  See :func:`_extract_ligand_atoms_from_sdf`."""
    return _extract_ligand_atoms_from_sdf(sdf_path)[0]


def _pose_rmsd_vs_reference(
    pose_pdb: Path, ref_coords, ref_elements=None, mismatches=None
) -> float:
    """RMSD between docked pose ligand heavy atoms and reference coordinates."""
    from flexaidds.benchmark import compute_rmsd, extract_ligand_atoms_from_pdb

    try:
        pred, pred_elements = extract_ligand_atoms_from_pdb(
            pose_pdb, expected_n_atoms=len(ref_coords))
    except ValueError:
        return -1.0
    if pred.shape != ref_coords.shape:
        # Do NOT truncate.  The previous behaviour trimmed both arrays to the
        # shorter length and computed an RMSD anyway, which silently compared
        # atom i of one molecule against atom i+1 of the other whenever the
        # extra atom sat at the front.  On 1mq6 that reported 6.49 A where the
        # true best-of-10 is 7.18 A -- a plausible number, in the optimistic
        # direction, indistinguishable from a measurement.  A count mismatch
        # here is a bug in selection, not a condition to work around.
        return -1.0
    try:
        # Symmetry-corrected, as FlexAID's calc_rmsd Hungarian branch.
        #
        # Both element lists are read from their own file.  If they disagree,
        # the two files do not list the same atoms in the same order -- the
        # #354 correspondence bug, resurfacing on a new pair -- and every
        # number downstream is measured against a shifted correspondence.
        # Report the sentinel rather than a plausible RMSD.
        if ref_elements is not None and list(ref_elements) != list(pred_elements):
            # Fail closed, but LOUDLY.  -1.0 is the same sentinel a genuine
            # miss produces, and metrics.docking_power keeps the target in the
            # denominator either way -- so a systematic element-order
            # divergence would silently depress reported docking power with no
            # diagnostic anywhere.  A refusal nobody can see has the same shape
            # as the bug it replaced.
            logger.warning(
                "RMSD refused for %s: reference and pose list different "
                "elements (ref %s..., pose %s...).  The two files do not "
                "describe the same atoms in the same order; this target is "
                "scored as a miss and the count is reported in the run "
                "summary.",
                pose_pdb.name, list(ref_elements)[:6], list(pred_elements)[:6],
            )
            if mismatches is not None:
                mismatches.append(str(pose_pdb))
            return -1.0
        return compute_rmsd(pred, ref_coords, pred_elements)
    except ValueError:
        return -1.0
