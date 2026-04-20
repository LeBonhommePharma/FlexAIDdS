"""Result loader for FlexAID∆S docking output directories.

The primary public entry point is :func:`load_results`, which scans a
directory of PDB output files produced by the FlexAID∆S C++ engine and
assembles them into a :class:`~flexaidds.models.DockingResult` hierarchy.

Typical usage::

    from flexaidds import load_results

    result = load_results("/path/to/docking/output")
    top = result.top_mode()
    print(f"Best ΔG: {top.free_energy:.2f} kcal/mol  (mode {top.mode_id})")
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .io import parse_pose_result
from .models import BindingModeResult, DockingResult, PoseResult

_PDB_SUFFIXES = {".pdb", ".ent"}


def _collect_pose_files(root: Path) -> List[Path]:
    """Recursively collect all PDB-like files under *root*.

    Args:
        root: Top-level directory to search.

    Returns:
        Sorted list of :class:`pathlib.Path` objects with ``.pdb`` or
        ``.ent`` extensions.
    """
    files: List[Path] = []
    for path in root.rglob("*"):
        if path.is_file() and path.suffix.lower() in _PDB_SUFFIXES:
            files.append(path)
    return sorted(files)


def _mode_temperature(poses: Sequence[PoseResult]) -> Optional[float]:
    """Return the first non-None temperature found across *poses*."""
    for pose in poses:
        if pose.temperature is not None:
            return pose.temperature
    return None


def _mode_metric(poses: Sequence[PoseResult], name: str) -> Optional[float]:
    """Return the first non-None value of attribute *name* across *poses*."""
    for pose in poses:
        value = getattr(pose, name)
        if value is not None:
            return value
    return None


def _mode_frequency(poses: Sequence[PoseResult]) -> Optional[int]:
    """Infer the cluster size (pose count) for a binding mode.

    Checks several common REMARK key names in priority order
    (``frequency``, ``nposes``, ``population``, ``cluster_size``, ``size``).
    Falls back to the length of *poses* if none are found.

    Args:
        poses: All poses belonging to one binding mode.

    Returns:
        Integer cluster size, or ``None`` if *poses* is empty.
    """
    for pose in poses:
        for key in ("frequency", "nposes", "population", "cluster_size", "size"):
            value = pose.remarks.get(key)
            if isinstance(value, int):
                return value
    return len(poses) if poses else None


def _mode_metadata(poses: Sequence[PoseResult]) -> Dict[str, object]:
    """Collect REMARK fields that are identical across all *poses*.

    Only key/value pairs present in every pose with the same value are
    included, so the resulting dictionary contains only globally shared
    metadata (e.g. receptor name).

    Args:
        poses: All poses belonging to one binding mode.

    Returns:
        Dictionary of shared metadata fields, possibly empty.
    """
    meta: Dict[str, object] = {}
    if not poses:
        return meta
    shared_keys = set(poses[0].remarks)
    for pose in poses[1:]:
        shared_keys &= set(pose.remarks)
    for key in shared_keys:
        values = {pose.remarks.get(key) for pose in poses}
        if len(values) == 1:
            meta[key] = poses[0].remarks[key]
    return meta


def _build_mode(mode_id: int, poses: Sequence[PoseResult], rank: int = 1) -> BindingModeResult:
    """Construct a :class:`BindingModeResult` from a flat list of poses.

    Poses are sorted by ``(pose_rank, filename)`` and thermodynamic
    aggregates are extracted from the first pose that carries each field.

    Args:
        mode_id: Numeric binding-mode identifier.
        poses: All :class:`~flexaidds.models.PoseResult` objects for this mode.
        rank: Rank of this mode among all modes (1-based, where 1 = best).

    Returns:
        Fully assembled :class:`~flexaidds.models.BindingModeResult`.
    """
    ordered = sorted(poses, key=lambda p: (p.pose_rank, p.path.name))
    best_pose = None
    scored = [p for p in ordered if p.cf is not None]
    if scored:
        best_pose = min(scored, key=lambda p: p.cf)
    elif ordered:
        best_pose = ordered[0]

    return BindingModeResult(
        mode_id=mode_id,
        rank=rank,
        poses=list(ordered),
        free_energy=_mode_metric(ordered, "free_energy"),
        enthalpy=_mode_metric(ordered, "enthalpy"),
        entropy=_mode_metric(ordered, "entropy"),
        heat_capacity=_mode_metric(ordered, "heat_capacity"),
        std_energy=_mode_metric(ordered, "std_energy"),
        best_cf=best_pose.cf if best_pose else None,
        frequency=_mode_frequency(ordered),
        temperature=_mode_temperature(ordered),
        metadata=_mode_metadata(ordered),
    )


def load_results(path: str | Path) -> DockingResult:
    """Load all docking results from a FlexAID∆S output directory.

    Recursively scans *path* for PDB files, parses their REMARK headers to
    extract thermodynamic quantities and cluster assignments, then assembles
    the results into a :class:`~flexaidds.models.DockingResult`.

    Args:
        path: Path to the directory containing docking result PDB files.
            May be a string or :class:`pathlib.Path`; ``~`` is expanded.

    Returns:
        :class:`~flexaidds.models.DockingResult` containing all discovered
        binding modes, sorted by ``mode_id``.

    Raises:
        FileNotFoundError: If *path* does not exist.
        NotADirectoryError: If *path* points to a file rather than a directory.

    Example::

        result = load_results("output/run1")
        for mode in result.binding_modes:
            print(mode.mode_id, mode.free_energy, mode.n_poses)
    """
    root = Path(path).expanduser().resolve()
    if not root.exists():
        raise FileNotFoundError(f"Docking results directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Expected a results directory, got file: {root}")

    pose_files = _collect_pose_files(root)
    if not pose_files:
        return DockingResult(
            source_dir=root,
            binding_modes=[],
            temperature=None,
            metadata={"n_pose_files": 0},
        )

    grouped: Dict[int, List[PoseResult]] = defaultdict(list)
    for pose_file in pose_files:
        pose = parse_pose_result(pose_file)
        grouped[pose.mode_id].append(pose)

    # Build modes (initially with rank=1, will be corrected after sorting)
    modes = [_build_mode(mode_id, poses, rank=1) for mode_id, poses in sorted(grouped.items())]

    # Sort modes by free energy (best = lowest ΔG), fallback to alphabetical mode_id
    modes_sorted = sorted(
        modes,
        key=lambda m: (m.free_energy is None, m.free_energy, m.mode_id)
    )

    # Assign correct ranks (1-based, 1 = best)
    modes_with_ranks = [
        BindingModeResult(
            mode_id=m.mode_id,
            rank=rank,
            poses=m.poses,
            free_energy=m.free_energy,
            enthalpy=m.enthalpy,
            entropy=m.entropy,
            heat_capacity=m.heat_capacity,
            std_energy=m.std_energy,
            best_cf=m.best_cf,
            frequency=m.frequency,
            temperature=m.temperature,
            metadata=m.metadata,
        )
        for rank, m in enumerate(modes_sorted, start=1)
    ]

    temperature = next((m.temperature for m in modes_with_ranks if m.temperature is not None), None)

    return DockingResult(
        source_dir=root,
        binding_modes=modes_with_ranks,
        temperature=temperature,
        metadata={"n_pose_files": len(pose_files)},
    )
