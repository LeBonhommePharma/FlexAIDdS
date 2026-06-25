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
        benchmark_order:       Stable run/display order; lower values run first.
        targets:               Full list of target identifiers.
        structural_states:     Receptor states available (``holo``, ``apo``, ``af2``).
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
    metrics: List[str] = field(default_factory=list)
    expected_baselines: Dict[str, float] = field(default_factory=dict)
    published_baselines: Dict[str, float] = field(default_factory=dict)
    published_source: str = ""
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
            benchmark_order=int(raw.pop("benchmark_order", 1000)),
            targets=list(raw.pop("targets", [])),
            structural_states=list(raw.pop("structural_states", ["holo"])),
            metrics=list(raw.pop("metrics", [])),
            expected_baselines=dict(raw.pop("expected_baselines", {})),
            published_baselines=dict(raw.pop("published_baselines", {})),
            published_source=str(raw.pop("published_source", "")),
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

    def scheduled_targets(self, tier: int) -> List[str]:
        """Targets actually scheduled for this tier."""
        if tier == 2 and self.slug in KNOWN_LARGE_DATASETS:
            return [e["entry_id"] for e in load_large_dataset_catalog(self.slug)]
        return self.tier1_targets() if tier == 1 else self.targets

    def scheduled_work_items(self, tier: int) -> List[Tuple[str, str]]:
        """(entry_id, state) pairs scheduled for this tier — full scale for large datasets."""
        if tier == 2 and self.slug in KNOWN_LARGE_DATASETS:
            catalog = load_large_dataset_catalog(self.slug)
            return [(e["entry_id"], e.get("state", "crossdock")) for e in catalog]
        targets = self.scheduled_targets(tier)
        return [(tid, st) for tid in targets for st in self.structural_states]

    def effective_entry_count(self, tier: int) -> int:
        """Entry count for runtime planning (uses C++-parity scales for large datasets)."""
        if tier == 2 and self.slug in KNOWN_LARGE_DATASETS:
            return len(load_large_dataset_catalog(self.slug)) or KNOWN_LARGE_DATASETS[self.slug]
        return len(self.scheduled_work_items(tier))

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
            "published_baselines": self.config.published_baselines,
            "published_source": self.config.published_source,
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


# ---------------------------------------------------------------------------
# Entry manifest loader — O(1) resume discovery for 1000+ entry campaigns
# ---------------------------------------------------------------------------

KNOWN_LARGE_DATASETS: Dict[str, int] = {
    "astex_nonnative": 1113,
    "posex_cd": 1312,
    "posex": 1319,
}

_TIMING_PRIORS_PATH = (
    Path(__file__).resolve().parents[3] / "benchmarks" / "m3pro" / "large_dataset_timing_priors.json"
)
_CATALOG_PATH = (
    Path(__file__).resolve().parents[3] / "benchmarks" / "m3pro" / "large_dataset_entry_catalogs.json"
)


def _utc_now_iso() -> str:
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _iter_entry_result_jsons(entry_dir: Path):
    """Yield per-entry result JSON files, excluding manifests and dotfiles."""
    for jf in sorted(entry_dir.glob("*.json")):
        if jf.name.startswith(("_", ".")):
            continue
        if "_" not in jf.stem:
            continue
        yield jf


def load_large_dataset_catalog(slug: str, catalog_path: Optional[Union[str, Path]] = None) -> List[Dict[str, Any]]:
    """Load tier-2 work-item catalog for large datasets (1113+ entries)."""
    p = Path(catalog_path) if catalog_path else _CATALOG_PATH
    if not p.is_file():
        return []
    try:
        data = json.loads(p.read_text())
        return list(data.get("datasets", {}).get(slug, {}).get("entries", []))
    except Exception:
        return []


def sanitize_entry_manifest(manifest: dict) -> dict:
    """Strip legacy polluted keys (e.g. None_None from .cost_history.json glob)."""
    for section in ("per_entry_status",):
        block = manifest.get(section)
        if isinstance(block, dict):
            manifest[section] = {
                k: v for k, v in block.items()
                if k and "None" not in k and "_" in k
            }
    timings = manifest.get("timings", {})
    for key in ("per_entry_wall_seconds", "per_entry_cost_cpu_seconds"):
        block = timings.get(key)
        if isinstance(block, dict):
            timings[key] = {
                k: v for k, v in block.items()
                if k and "None" not in k and "_" in k
            }
    return manifest


def load_timing_priors(
    priors_path: Optional[Union[str, Path]] = None,
) -> Dict[str, Dict[str, Any]]:
    """Load committed large-dataset timing priors (real benchmark wall_time_s aggregates)."""
    p = Path(priors_path) if priors_path else _TIMING_PRIORS_PATH
    if not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text())
        return data.get("datasets", {})
    except Exception:
        return {}

def load_entry_manifest(manifest_path: Union[str, Path]) -> Optional[dict]:
    """Load _entry_manifest.json once; return None on absence or corruption."""
    p = Path(manifest_path)
    if not p.is_file():
        return None
    try:
        return sanitize_entry_manifest(json.loads(p.read_text()))
    except Exception as exc:
        logger.warning("Corrupt entry manifest %s — falling back to per-file scan: %s", p, exc)
        return None


def completed_targets_from_manifest(
    manifest: dict,
    targets: Sequence[str],
    states: Sequence[str],
) -> Optional[set[str]]:
    """Derive fully-completed target_ids from a manifest without per-entry JSON reads.

    Returns None when the manifest lacks enough structure for a trustworthy fast path.
    """
    per_entry_status = manifest.get("per_entry_status")
    if isinstance(per_entry_status, dict) and per_entry_status:
        completed: set[str] = set()
        for tid in targets:
            if all(
                bool(per_entry_status.get(f"{tid}_{st}", {}).get("success", False))
                for st in states
            ):
                completed.add(tid)
        return completed

    # Legacy manifests: single-state campaigns can trust the completed list directly.
    if len(states) == 1:
        legacy = set(manifest.get("completed", []))
        if legacy:
            return {tid for tid in targets if tid in legacy}

        wall = manifest.get("timings", {}).get("per_entry_wall_seconds", {})
        if isinstance(wall, dict) and wall:
            state = states[0]
            return {
                tid for tid in targets
                if f"{tid}_{state}" in wall and float(wall[f"{tid}_{state}"]) > 0.0
            }

    return None


def plan_runtime(
    *,
    results_dir: Union[str, Path],
    workers: int = 1,
    omp_threads: int = 4,
    datasets: Optional[Dict[str, int]] = None,
    output_path: Optional[Union[str, Path]] = None,
    default_mean_entry_s: float = 180.0,
) -> Dict[str, Any]:
    """Emit scaled wall-clock / CPU-second estimates for large benchmark campaigns.

    Reads timing summaries from existing ``_entry_manifest.json`` files and/or
    ``.cost_history.json`` EMA hints when present; otherwise uses
    ``default_mean_entry_s``.
    """
    results_dir = Path(results_dir)
    workers = max(1, int(workers))
    omp_threads = max(1, int(omp_threads))
    scales = dict(KNOWN_LARGE_DATASETS)
    if datasets:
        scales.update(datasets)

    timing_priors = load_timing_priors()
    lines: List[str] = [
        "FlexAIDdS DatasetRunner — Computational Runtime Plan",
        f"Generated: {_utc_now_iso()}",
        f"Workers: {workers}  |  OMP threads/worker: {omp_threads}",
        "",
    ]
    plan: Dict[str, Any] = {
        "generated_at": _utc_now_iso(),
        "workers": workers,
        "omp_threads": omp_threads,
        "datasets": {},
    }

    for slug, n_entries in sorted(scales.items(), key=lambda kv: kv[0]):
        mean_s = default_mean_entry_s
        median_s = default_mean_entry_s
        source = "default_prior"
        prior = timing_priors.get(slug, {})
        if prior:
            n_entries = int(prior.get("n_entries", n_entries))
            mean_s = float(prior.get("mean_entry_seconds", mean_s))
            median_s = float(prior.get("median_entry_seconds", median_s))
            source = str(prior.get("timing_source", "timing_priors_json"))

        tier_dir = results_dir / slug / "tier2"
        if not tier_dir.is_dir():
            tier_dir = results_dir / slug / "tier1"
        manifest = load_entry_manifest(tier_dir / "_entry_manifest.json")
        if manifest and "timings" in manifest:
            summary = manifest["timings"].get("summary", {})
            if summary.get("mean_entry_seconds", 0) > 0:
                mean_s = float(summary["mean_entry_seconds"])
                source = f"manifest:{tier_dir / '_entry_manifest.json'}"
            if summary.get("median_entry_seconds", 0) > 0:
                median_s = float(summary["median_entry_seconds"])
        elif source == "default_prior":
            history = tier_dir / ".cost_history.json"
            if history.is_file():
                try:
                    hints = json.loads(history.read_text())
                    if hints:
                        mean_s = float(sum(hints.values()) / len(hints))
                        median_s = mean_s
                        source = f"cost_history:{history}"
                except Exception:
                    pass

        est_wall_mean = mean_s * n_entries / workers
        est_wall_median = median_s * n_entries / workers
        est_cpu_mean = est_wall_mean * omp_threads

        entry = {
            "n_entries": n_entries,
            "mean_entry_seconds": round(mean_s, 2),
            "median_entry_seconds": round(median_s, 2),
            "timing_source": source,
            "estimated_wall_seconds_mean": round(est_wall_mean, 1),
            "estimated_wall_seconds_median": round(est_wall_median, 1),
            "estimated_cpu_seconds_mean": round(est_cpu_mean, 1),
            "estimated_wall_hours_mean": round(est_wall_mean / 3600.0, 2),
        }
        plan["datasets"][slug] = entry

        lines += [
            f"## {slug}",
            f"  N entries: {n_entries}",
            f"  Timing source: {source}",
            f"  Mean entry: {mean_s:.1f}s  |  Median entry: {median_s:.1f}s",
            f"  Est. wall (mean): {entry['estimated_wall_hours_mean']:.2f} h "
            f"({entry['estimated_wall_seconds_mean']:.0f} s)",
            f"  Est. wall (median): {entry['estimated_wall_seconds_median']:.0f} s",
            f"  Est. CPU-seconds (mean × OMP): {entry['estimated_cpu_seconds_mean']:.0f}",
            "",
        ]

    text = "\n".join(lines)
    plan["text_summary"] = text
    if output_path is not None:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text)
        plan_json = out.with_suffix(".json") if out.suffix == ".txt" else out.parent / f"{out.stem}.json"
        plan_json.write_text(json.dumps(plan, indent=2))
    return plan


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
        omp_threads: Optional[int] = None,
        use_mpi: bool = False,
        cache_dir: Optional[Union[str, Path]] = None,
        bootstrap_ci: bool = False,
        n_bootstrap: int = 5_000,
        dry_run: bool = False,
        repo_root: Optional[Union[str, Path]] = None,
        resume: bool = False,
    ) -> None:
        _default_datasets = Path(__file__).resolve().parent / "datasets"
        self.datasets_dir = Path(datasets_dir) if datasets_dir is not None else _default_datasets
        self.results_dir = Path(results_dir)
        self.binary = binary or os.environ.get("FLEXAIDDS_BINARY") or "FlexAID"
        self.temperature = temperature
        self.n_workers = max(1, int(n_workers))
        # OMP threads per worker: explicit > env > auto (aim for 2 on M3 Pro's 5 P-cores).
        if omp_threads is None:
            env_val = os.environ.get("FLEXAIDDS_OMP_THREADS")
            omp_threads = int(env_val) if env_val else max(1, 2 if self.n_workers > 1 else 4)
        self.omp_threads = max(1, int(omp_threads))
        self.use_mpi = use_mpi
        self.cache_dir = Path(
            cache_dir or os.environ.get("FLEXAIDDS_BENCHMARK_DATA", "benchmark_data")
        )
        self.do_bootstrap = bootstrap_ci
        self.n_bootstrap = n_bootstrap
        self.dry_run = dry_run
        self.repo_root = Path(repo_root) if repo_root else Path.cwd()
        self.resume = bool(resume)

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

    def _find_receptor(self, target_id: str, data_dir: Path, state: str) -> Optional[Path]:
        """Locate receptor PDB for a given target and structural state.

        Probes multiple naming conventions:
        - PDBbind/CASF style:  <id>/<id>_holo.pdb  or  <id>/<id>_protein.pdb
        - Astex Diverse style: <id>/<id>.pdb        (uppercase or lowercase)
        - Flat layout:         <id>.pdb
        """
        def _probes(tid: str) -> List[Path]:
            return [
                data_dir / tid / f"{tid}_{state}.pdb",
                data_dir / tid / f"{tid}_protein.pdb",
                data_dir / tid / f"{tid}.pdb",          # Astex Diverse bare-name convention
                data_dir / tid / "receptor.pdb",
                data_dir / f"{tid}.pdb",
            ]

        for p in _probes(target_id):
            if p.is_file():
                return p
        # Case-insensitive fallback: Astex Diverse stores PDB IDs as uppercase on disk
        upper = target_id.upper()
        if upper != target_id:
            for p in _probes(upper):
                if p.is_file():
                    return p
        logger.warning("No receptor found for %s (%s) in %s", target_id, state, data_dir)
        return None

    def _find_ligands(self, target_id: str, data_dir: Path) -> List[Path]:
        """Locate all ligand files for a target (Mol2 or SDF)."""
        target_dir = data_dir / target_id
        if not target_dir.is_dir():
            # Try uppercase variant (Astex Diverse stores IDs as uppercase)
            target_dir = data_dir / target_id.upper()
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

                sub_env = os.environ.copy()
                sub_env["OMP_NUM_THREADS"] = str(self.omp_threads)
                try:
                    result = subprocess.run(
                        [self.binary, str(cfg_path)],
                        capture_output=True,
                        text=True,
                        timeout=3600,
                        cwd=tmp_path,
                        env=sub_env,
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

        Now uses automated per-entry (per-target) saving and processing.
        Individual TargetResult files are written atomically as soon as each
        (target, structural_state) finishes. This enables resume and fine-grained
        progress / resource monitoring.

        A lightweight EntryTaskManager (see below) coordinates the work items.
        """
        t0 = time.monotonic()
        targets = config.scheduled_targets(tier)
        states = structural_states or config.structural_states
        requested_metrics = metric_subset or config.metrics or None
        n_entries_scale = config.effective_entry_count(tier)
        scheduled_items = config.scheduled_work_items(tier)

        dr = DatasetResult(
            config=config,
            tier=tier,
            targets_attempted=list(targets),
            timestamp=_utc_now_iso(),
            git_sha=_git_sha(self.repo_root),
            host=socket.gethostname(),
        )

        if not targets:
            logger.warning("Dataset %s has no targets", config.slug)
            dr.duration_seconds = 0.0
            return dr

        # --- NEW: Per-entry resume + individual processing automation ---
        already_completed: set[str] = set()
        if self.resume:
            already_completed = self._discover_completed_targets(config, tier, targets)
            if already_completed:
                logger.info(
                    "Resume mode: %d/%d targets already have complete per-entry results — skipping them "
                    "(tier-%d scale: %d total entries for %s)",
                    len(already_completed), len(targets), tier, n_entries_scale, config.slug,
                )

        # Build work items from full catalog (tier-2 large-N) or targets × states
        all_work_items: List[Tuple[str, str]] = []
        for tid, st in scheduled_items:
            if tid in already_completed:
                continue
            all_work_items.append((tid, st))

        if tier == 2 and config.slug in KNOWN_LARGE_DATASETS:
            logger.info(
                "Large-N dataset %s tier-%d: %d scheduled work items (catalog scale)",
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
                    receptor = self._find_receptor(target_id, config.data_dir, state)
                    if receptor is None:
                        error = f"No receptor found for {target_id}/{state}"
                        return target_id, state, [], time.monotonic() - t_start, error
                    ligands = self._find_ligands(target_id, config.data_dir) or []

                poses = self._dock_target(
                    target_id,
                    receptor or Path("/dev/null"),
                    ligands or [Path(f"{target_id}.mol2")],
                    structural_state=state,
                )

                elapsed = time.monotonic() - t_start
                cost_cpu = elapsed * max(1, self.omp_threads)

                tr = TargetResult(
                    target_id=target_id,
                    structural_state=state,
                    poses=poses,
                    duration_seconds=elapsed,
                    error=error,
                )
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

        completed = sorted(completed_targets)
        failed = sorted(failed_targets)

        # MPI gather (still works at target granularity for final aggregation)
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

        # Root finalizes
        if self._mpi_root:
            dr.targets_completed = completed
            dr.targets_failed = failed

            if all_poses:
                metrics = compute_all_metrics(all_poses, requested=requested_metrics)
                dr.metrics = metrics

                if self.do_bootstrap:
                    dr.ci_95 = self._compute_bootstrap_cis(all_poses, requested_metrics)

            if not self.dry_run:
                dr.check_regressions()
            else:
                logger.info("Dry-run mode — skipping regression checks")

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
                generated_at=_utc_now_iso(),
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
            generated_at=_utc_now_iso(),
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
            "timestamp": _utc_now_iso(),
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
            )
        except Exception as e:
            logger.warning("Corrupt target result %s — will re-run: %s", path, e)
            return None

    def _probe_entry_success(
        self, config: DatasetConfig, tier: int, target_id: str, structural_state: str
    ) -> Optional[bool]:
        """Cheap success probe — reads JSON metadata only, skips pose parsing."""
        path = self._target_result_path(config, tier, target_id, structural_state)
        if not path.is_file():
            return None
        try:
            data = json.loads(path.read_text())
            if "success" in data:
                return bool(data["success"])
            return not bool(data.get("error")) and bool(data.get("poses"))
        except Exception as exc:
            logger.warning("Corrupt target result %s — will re-run: %s", path, exc)
            return False

    def _discover_completed_targets(
        self,
        config: DatasetConfig,
        tier: int,
        scheduled_targets: Optional[List[str]] = None,
    ) -> set[str]:
        """Return target_ids with successful results for every requested state.

        Uses manifest-first fast path (single JSON read) for 1000+ entry campaigns;
        falls back to per-entry success probes without full pose parsing.
        """
        targets = scheduled_targets or config.scheduled_targets(tier)
        if tier == 2 and config.slug in KNOWN_LARGE_DATASETS:
            states = ["crossdock"]
        else:
            states = config.structural_states
        manifest_path = self._entry_results_dir(config, tier) / "_entry_manifest.json"
        manifest = load_entry_manifest(manifest_path)
        if manifest is not None:
            fast = completed_targets_from_manifest(manifest, targets, states)
            if fast is not None:
                logger.info(
                    "Manifest-first resume: %d/%d targets complete "
                    "(1 manifest read, skipped %d per-entry JSON parses)",
                    len(fast),
                    len(targets),
                    max(0, len(targets) - len(fast)),
                )
                return fast

        completed: set[str] = set()
        for target_id in targets:
            all_states_ok = True
            for st in states:
                ok = self._probe_entry_success(config, tier, target_id, st)
                if ok is not True:
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
        per_entry_status: Dict[str, Dict[str, Any]] = {}
        durations: List[float] = []

        for jf in sorted(entry_dir.glob("*_*.json")):
            if jf.name.startswith("_"):
                continue
            try:
                data = json.loads(jf.read_text())
                key = f"{data.get('target_id')}_{data.get('structural_state')}"
                dur = float(data.get("duration_seconds", 0.0))
                cost = float(data.get("cost_cpu_seconds", dur * max(1, self.omp_threads)))
                success = bool(data.get("success", not data.get("error") and data.get("poses")))
                per_entry[key] = round(dur, 2)
                per_entry_cost[key] = round(cost, 2)
                per_entry_status[key] = {
                    "success": success,
                    "duration_seconds": round(dur, 2),
                }
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
            "timestamp": _utc_now_iso(),
            "omp_threads": self.omp_threads,
            "n_workers_used": self.n_workers,
            "completed": completed,
            "failed": failed,
            "total_attempted": len(completed) + len(failed),
            "per_entry_status": per_entry_status,
            # --- per-target timing + cost tracking ---
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
            generated_at=_utc_now_iso(),
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
