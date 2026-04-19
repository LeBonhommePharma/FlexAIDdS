"""DatasetRunner — distributed benchmarking orchestrator for FlexAIDdS.

Discovers dataset configs, distributes work across MPI nodes or local
processes, runs FlexAIDdS docking, computes all metrics, and produces
structured JSON + Markdown reports.

Typical usage
-------------
Library::

    from flexaidds.dataset_runner import DatasetRunner

    runner = DatasetRunner(results_dir="results/benchmark_run")
    report = runner.run_all(tier=1)
    json_path, md_path = report.save("results/benchmark_run/report")

CLI::

    python -m flexaidds.dataset_runner --dataset casf2016 --tier 1
    python -m flexaidds.dataset_runner --all --distributed --nodes 4

MPI distributed run::

    mpirun -n 8 python -m flexaidds.dataset_runner --all --distributed
"""

from __future__ import annotations

import datetime
import json
import logging
import os
import socket
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

try:
    import yaml  # PyYAML
    _HAS_YAML = True
except ImportError:
    _HAS_YAML = False

import numpy as np

from .metrics import (
    PoseScore,
    bootstrap_ci,
    compute_all_metrics,
    entropy_rescue_rate,
    docking_power,
)

logger = logging.getLogger(__name__)


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
        targets:               Full list of target identifiers.
        structural_states:     Receptor states available (``holo``, ``apo``, ``af2``).
        metrics:               Names of metrics to compute (must exist in metrics.py).
        expected_baselines:    ``{metric: value}`` reference values for regression checks.
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
    targets: List[str] = field(default_factory=list)
    structural_states: List[str] = field(default_factory=lambda: ["holo"])
    metrics: List[str] = field(default_factory=list)
    expected_baselines: Dict[str, float] = field(default_factory=dict)
    baseline_tolerance: float = 0.05
    data_dir: Optional[Path] = None
    data_format: str = "pdb"
    active_label_field: str = "is_active"

    @classmethod
    def from_yaml(cls, yaml_path: Path) -> "DatasetConfig":
        """Load a DatasetConfig from a YAML file."""
        if not _HAS_YAML:
            raise RuntimeError(
                "PyYAML is required to load dataset configs: pip install pyyaml"
            )
        with open(yaml_path) as fh:
            raw: dict = yaml.safe_load(fh)

        data_dir_raw = raw.pop("data_dir", None)
        config = cls(
            slug=raw.pop("slug", yaml_path.stem),
            name=raw.pop("name", yaml_path.stem),
            description=raw.pop("description", ""),
            zenodo_doi=raw.pop("zenodo_doi", ""),
            download_url=raw.pop("download_url", ""),
            tier=int(raw.pop("tier", 2)),
            tier1_subset_size=int(raw.pop("tier1_subset_size", 5)),
            targets=list(raw.pop("targets", [])),
            structural_states=list(raw.pop("structural_states", ["holo"])),
            metrics=list(raw.pop("metrics", [])),
            expected_baselines=dict(raw.pop("expected_baselines", {})),
            baseline_tolerance=float(raw.pop("baseline_tolerance", 0.05)),
            data_format=str(raw.pop("data_format", "pdb")),
            active_label_field=str(raw.pop("active_label_field", "is_active")),
        )
        if data_dir_raw:
            config.data_dir = Path(data_dir_raw)
        return config

    def tier1_targets(self) -> List[str]:
        """Return the subset of targets used for tier-1 (fast) runs."""
        return self.targets[: self.tier1_subset_size]


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

    @property
    def success(self) -> bool:
        return not self.error and bool(self.poses)


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

    def check_regressions(self) -> Dict[str, bool]:
        """Flag metrics that regressed below baseline − tolerance.

        Returns dict of ``{metric: True}`` for regressed metrics.
        """
        flags: Dict[str, bool] = {}
        tol = self.config.baseline_tolerance
        for metric, baseline in self.config.expected_baselines.items():
            measured = self.metrics.get(metric)
            if measured is None or np.isnan(measured):
                continue
            # For error/RMSE metrics lower is better — allow increase
            if "rmse" in metric or "mae" in metric:
                threshold = baseline * (1 + tol)
                flags[metric] = bool(measured > threshold)
            else:
                # For all other metrics higher is better
                threshold = baseline * (1 - tol)
                flags[metric] = bool(measured < threshold)
        self.regression_flags = flags
        return flags

    def to_dict(self) -> dict:
        """Serialise to a JSON-compatible dict."""
        return {
            "dataset": self.config.slug,
            "tier": self.tier,
            "timestamp": self.timestamp,
            "git_sha": self.git_sha,
            "host": self.host,
            "duration_seconds": self.duration_seconds,
            "targets_attempted": len(self.targets_attempted),
            "targets_completed": len(self.targets_completed),
            "targets_failed": len(self.targets_failed),
            "metrics": self.metrics,
            "ci_95": {k: list(v) for k, v in self.ci_95.items()},
            "regression_flags": self.regression_flags,
            "expected_baselines": self.config.expected_baselines,
        }


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

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        return {
            "generated_at": self.generated_at,
            "git_sha": self.git_sha,
            "host": self.host,
            "runner_info": self.runner_info,
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
            "",
            "---",
            "",
        ]

        for dr in self.datasets:
            lines += self._dataset_section(dr)

        return "\n".join(lines)

    @staticmethod
    def _dataset_section(dr: DatasetResult) -> List[str]:
        lines = [
            f"## {dr.config.name}",
            "",
            f"*Tier {dr.tier} · {len(dr.targets_completed)}/{len(dr.targets_attempted)} targets · "
            f"{dr.duration_seconds:.1f}s*",
            "",
        ]

        if dr.metrics:
            lines += ["| Metric | Value | 95% CI | Baseline | Regressed? |",
                      "|--------|-------|--------|----------|------------|"]
            for metric, value in sorted(dr.metrics.items()):
                ci = dr.ci_95.get(metric)
                ci_str = f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "—"
                baseline = dr.config.expected_baselines.get(metric)
                baseline_str = f"{baseline:.3f}" if baseline is not None else "—"
                regressed = dr.regression_flags.get(metric, False)
                flag = "⚠ YES" if regressed else "OK"
                lines.append(
                    f"| {metric} | {value:.4f} | {ci_str} | {baseline_str} | {flag} |"
                )
            lines.append("")

        if dr.targets_failed:
            lines += [
                f"**Failed targets** ({len(dr.targets_failed)}): "
                + ", ".join(dr.targets_failed[:10])
                + ("…" if len(dr.targets_failed) > 10 else ""),
                "",
            ]

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


def _split_targets(targets: List[str], rank: int, size: int) -> List[str]:
    """Distribute target list across MPI ranks (round-robin)."""
    return [t for i, t in enumerate(targets) if i % size == rank]


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
        dry_run:        Skip actual docking; useful for testing the framework.
        repo_root:      Root of the FlexAIDdS repository (for git SHA detection).
    """

    def __init__(
        self,
        datasets_dir: Union[str, Path, None] = None,
        results_dir: Union[str, Path] = "results/benchmarks",
        binary: Optional[str] = None,
        temperature: float = 300.0,
        n_workers: int = 1,
        use_mpi: bool = False,
        cache_dir: Optional[Union[str, Path]] = None,
        bootstrap_ci: bool = False,
        n_bootstrap: int = 5_000,
        dry_run: bool = False,
        repo_root: Optional[Union[str, Path]] = None,
    ) -> None:
        _default_datasets = Path(__file__).resolve().parent / "datasets"
        self.datasets_dir = Path(datasets_dir) if datasets_dir is not None else _default_datasets
        self.results_dir = Path(results_dir)
        self.binary = binary or os.environ.get("FLEXAIDDS_BINARY") or "FlexAID"
        self.temperature = temperature
        self.n_workers = n_workers
        self.use_mpi = use_mpi
        self.cache_dir = Path(
            cache_dir or os.environ.get("FLEXAIDDS_BENCHMARK_DATA", "benchmark_data")
        )
        self.do_bootstrap = bootstrap_ci
        self.n_bootstrap = n_bootstrap
        self.dry_run = dry_run
        self.repo_root = Path(repo_root) if repo_root else Path.cwd()

        self.results_dir.mkdir(parents=True, exist_ok=True)

        if use_mpi:
            self._mpi_rank, self._mpi_size, self._mpi_root, self._mpi_comm = (
                _mpi_context()
            )
        else:
            self._mpi_rank, self._mpi_size, self._mpi_root, self._mpi_comm = (
                0, 1, True, None
            )

    # ------------------------------------------------------------------
    # Dataset discovery
    # ------------------------------------------------------------------

    def discover_datasets(self) -> List[DatasetConfig]:
        """Discover and load all ``*.yaml`` configs from ``datasets_dir``.

        Returns:
            List of :class:`DatasetConfig` objects sorted by slug.
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

        return configs

    def load_dataset_config(self, yaml_path: Union[str, Path]) -> DatasetConfig:
        """Load a single dataset config from an explicit path."""
        cfg = DatasetConfig.from_yaml(Path(yaml_path))
        if cfg.data_dir is None:
            cfg.data_dir = self.cache_dir / cfg.slug
        return cfg

    # ------------------------------------------------------------------
    # Core docking dispatch
    # ------------------------------------------------------------------

    def _find_receptor(self, target_id: str, data_dir: Path, state: str) -> Optional[Path]:
        """Locate receptor PDB for a given target and structural state."""
        candidates = [
            data_dir / target_id / f"{target_id}_{state}.pdb",
            data_dir / target_id / f"{target_id}_protein.pdb",
            data_dir / target_id / f"receptor.pdb",
            data_dir / f"{target_id}.pdb",
        ]
        for p in candidates:
            if p.is_file():
                return p
        logger.warning("No receptor found for %s (%s) in %s", target_id, state, data_dir)
        return None

    def _find_ligands(self, target_id: str, data_dir: Path) -> List[Path]:
        """Locate all ligand files for a target (Mol2 or SDF)."""
        target_dir = data_dir / target_id
        if not target_dir.is_dir():
            target_dir = data_dir
        ligands = (
            list(target_dir.glob("*.mol2"))
            + list(target_dir.glob("*.sdf"))
            + list((target_dir / "ligands").glob("*.mol2") if (target_dir / "ligands").is_dir() else [])
        )
        return ligands

    def _dock_target(
        self,
        target_id: str,
        receptor_path: Path,
        ligand_paths: List[Path],
        structural_state: str = "holo",
        with_entropy: bool = True,
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
                target_id, receptor_path, ligand_paths,
                structural_state, with_entropy,
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
    ) -> List[PoseScore]:
        """Invoke the FlexAID binary and parse output poses."""
        poses: List[PoseScore] = []

        for ligand_path in ligand_paths:
            ligand_id = ligand_path.stem

            with tempfile.TemporaryDirectory(prefix=f"flexaid_{target_id}_") as tmp:
                tmp_path = Path(tmp)
                cfg_path = tmp_path / "dock.inp"

                cfg_lines = [
                    f"PDBNAM {receptor_path}",
                    f"INPLIG {ligand_path}",
                    f"TEMPER {int(self.temperature)}",
                    "METOPT GA",
                    "COMPLF VCT",
                ]
                if not with_entropy:
                    cfg_lines.append("NOENTROPY 1")

                cfg_path.write_text("\n".join(cfg_lines) + "\n")

                try:
                    result = subprocess.run(
                        [self.binary, str(cfg_path)],
                        capture_output=True,
                        text=True,
                        timeout=3600,
                        cwd=tmp_path,
                    )
                    if result.returncode != 0:
                        logger.warning(
                            "FlexAID returned %d for %s/%s",
                            result.returncode, target_id, ligand_id,
                        )
                    else:
                        parsed = self._parse_flexaid_output(
                            tmp_path, target_id, ligand_id, structural_state
                        )
                        poses.extend(parsed)
                except subprocess.TimeoutExpired:
                    logger.error("Docking timed out: %s/%s", target_id, ligand_id)

        return poses

    @staticmethod
    def _parse_flexaid_output(
        work_dir: Path,
        target_id: str,
        ligand_id: str,
        structural_state: str,
    ) -> List[PoseScore]:
        """Parse FlexAIDdS result PDB files from a completed docking run.

        Reads REMARK lines for energy, entropy, RMSD metadata.
        """
        poses: List[PoseScore] = []
        pdb_files = sorted(work_dir.glob("*.pdb"))

        for rank, pdb_path in enumerate(pdb_files, start=1):
            rmsd = -1.0
            enthalpy_score = 0.0
            entropy_correction = 0.0
            total_score = 0.0
            is_active = False
            exp_affinity: Optional[float] = None

            try:
                for line in pdb_path.read_text().splitlines():
                    if not line.startswith("REMARK"):
                        continue
                    # Parse structured REMARK fields emitted by FlexAID
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

            except Exception as exc:
                logger.debug("Error parsing %s: %s", pdb_path, exc)
                continue

            # Fallback: total = enthalpy - entropy if not set
            if total_score == 0.0:
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

        Args:
            config:           Dataset configuration.
            tier:             1 = fast subset, 2 = full target list.
            metric_subset:    Restrict metrics to this list (None = all).
            structural_states: Override structural states from config.

        Returns:
            :class:`DatasetResult` with computed metrics.
        """
        t0 = time.monotonic()
        targets = config.tier1_targets() if tier == 1 else config.targets
        states = structural_states or config.structural_states
        requested_metrics = metric_subset or config.metrics or None

        dr = DatasetResult(
            config=config,
            tier=tier,
            targets_attempted=list(targets),
            timestamp=datetime.datetime.utcnow().isoformat() + "Z",
            git_sha=_git_sha(self.repo_root),
            host=socket.gethostname(),
        )

        if not targets:
            logger.warning("Dataset %s has no targets", config.slug)
            dr.duration_seconds = 0.0
            return dr

        # MPI distribution
        my_targets = _split_targets(targets, self._mpi_rank, self._mpi_size)
        logger.info(
            "[rank %d/%d] Dataset %s: running %d/%d targets",
            self._mpi_rank, self._mpi_size,
            config.slug, len(my_targets), len(targets),
        )

        all_poses: List[PoseScore] = []
        completed: List[str] = []
        failed: List[str] = []

        for target_id in my_targets:
            t_start = time.monotonic()
            target_poses: List[PoseScore] = []

            for state in states:
                receptor = None
                if config.data_dir and not self.dry_run:
                    receptor = self._find_receptor(target_id, config.data_dir, state)
                    if receptor is None and not self.dry_run:
                        logger.warning("No receptor for %s/%s — skipping", target_id, state)
                        continue

                ligands: List[Path] = []
                if config.data_dir and not self.dry_run:
                    ligands = self._find_ligands(target_id, config.data_dir)

                poses = self._dock_target(
                    target_id,
                    receptor or Path("/dev/null"),
                    ligands or [Path(f"{target_id}.mol2")],
                    structural_state=state,
                )
                target_poses.extend(poses)

            if target_poses:
                all_poses.extend(target_poses)
                completed.append(target_id)
            else:
                failed.append(target_id)
                logger.warning("No poses for target %s", target_id)

            logger.debug(
                "Target %s: %d poses in %.1fs",
                target_id, len(target_poses), time.monotonic() - t_start,
            )

        # MPI gather poses to root
        if self._mpi_comm is not None:
            all_results_by_rank = self._mpi_comm.gather(
                (all_poses, completed, failed), root=0
            )
            if self._mpi_root:
                all_poses = []
                completed = []
                failed = []
                for poses_i, comp_i, fail_i in (all_results_by_rank or []):
                    all_poses.extend(poses_i)
                    completed.extend(comp_i)
                    failed.extend(fail_i)

        # Only root computes metrics and writes report
        if self._mpi_root:
            dr.targets_completed = completed
            dr.targets_failed = failed

            if all_poses:
                metrics = compute_all_metrics(all_poses, requested=requested_metrics)
                dr.metrics = metrics

                if self.do_bootstrap:
                    dr.ci_95 = self._compute_bootstrap_cis(
                        all_poses, requested_metrics
                    )

            if not self.dry_run:
                dr.check_regressions()
            else:
                logger.info("Dry-run mode — skipping regression checks against baselines")

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
                "DatasetRunner: %d dataset(s) · tier %d · %d MPI rank(s)",
                len(all_configs), tier, self._mpi_size,
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
        )
        return report

    def _save_dataset_result(self, dr: DatasetResult) -> None:
        """Write a per-dataset JSON result file as soon as it's ready."""
        out_path = self.results_dir / f"{dr.config.slug}_tier{dr.tier}.json"
        out_path.write_text(json.dumps(dr.to_dict(), indent=2))
        logger.info("Dataset result saved: %s", out_path)

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
