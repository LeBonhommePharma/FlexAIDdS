"""Data models for FlexAID∆S docking results.

This module defines frozen dataclasses that represent docking output at three
levels of granularity:

- :class:`PoseResult` – a single docked pose (one PDB file).
- :class:`BindingModeResult` – a cluster of poses sharing a binding geometry.
- :class:`DockingResult` – the top-level container returned by
  :func:`~flexaidds.results.load_results`.

All three classes are immutable (``frozen=True``) so they can be safely shared
across threads and used as dictionary keys.
"""

from __future__ import annotations

import ast
import csv
import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .thermodynamics import (
    ClaimValidity,
    ScientificProvenance,
    ThermodynamicBreakdown,
)

# Grand-canonical / competitive models (P2 additive, behind HAS_GRAND_BINDINGS in __init__)
# LigandSpec describes inputs for multi-ligand grand canonical calculations.
# Extended fields on DockingResult are optional (None/empty for legacy single-ligand paths).


def _as_breakdown(value):
    """Coerce a serialized thermodynamics dict into a ThermodynamicBreakdown.

    Robust round-trip helper: dict -> ThermodynamicBreakdown via from_dict();
    anything else (already a breakdown, or None) passes through unchanged.
    """
    if isinstance(value, dict):
        return ThermodynamicBreakdown.from_dict(value)
    return value


def _as_provenance(value: Any) -> ScientificProvenance:
    """Coerce serialized provenance without trusting claimed validity."""
    if isinstance(value, ScientificProvenance):
        return value
    if isinstance(value, dict):
        return ScientificProvenance.from_dict(value)
    return ScientificProvenance()


@dataclass(frozen=True)
class LigandSpec:
    """Specification for one ligand in a grand-canonical / competitive binding context.

    Used by compute_grand_partition and future multi-ligand loaders/campaigns.
    concentration_M is in molar units (standard state = 1.0 M).

    Additive only — existing single-ligand flows are unaffected.
    """

    name: str
    concentration_M: float = 1.0
    # Optional metadata for manifests, sidecars, reports
    smiles: Optional[str] = None
    ligand_id: Optional[str] = None
    results_dir: Optional[str] = None  # path for per-ligand result dir (multi-ligand layout)

    def __post_init__(self):
        if self.concentration_M <= 0.0:
            raise ValueError(f"concentration_M must be > 0 (got {self.concentration_M})")
        if self.concentration_M > 1000.0:
            raise ValueError(
                f"concentration_M > 1000 M — convert µM/nM to M first (got {self.concentration_M})"
            )


@dataclass(frozen=True)
class PoseResult:
    """A single docked pose read from one FlexAID∆S output PDB file.

    Attributes:
        path: Absolute path to the PDB file on disk.
        mode_id: Binding-mode (cluster) index this pose belongs to.
        pose_rank: Rank of this pose within its binding mode (1-based).
        cf: CF/contact-function scoring proxy (Voronoi CF). Lower is better.
            Not a free energy or experimental ΔG.
        cf_app: Apparent CF scoring proxy after grid-approximation correction.
        rmsd_raw: RMSD to reference structure without symmetry correction (Å).
        rmsd_sym: Symmetry-corrected RMSD to reference structure (Å).
        free_energy: Legacy ensemble transform. Its units and interpretation
            are determined exclusively by ``scientific_provenance``; current
            docking output is a CF-domain proxy, not experimental ΔG_bind.
        proxy_free_energy: Explicit schema-v2 name for the CF-domain transform.
        soft_beta_G: Mode-election objective emitted by the docking engine.
        enthalpy: Ensemble mean in the declared energy domain.
        entropy: Ensemble S-like diagnostic in the declared domain per kelvin.
        heat_capacity: Ensemble Cv-like diagnostic in the declared domain.
        std_energy: Standard deviation of ensemble energies σ_E in the
            declared energy domain (kcal/mol only under calibrated provenance).
        temperature: Simulation temperature (K) parsed from REMARK section.
        remarks: Raw key→value mapping of all ``REMARK`` fields parsed from the
            PDB header.
    """

    path: Path
    mode_id: int
    pose_rank: int
    cf: Optional[float] = None
    cf_app: Optional[float] = None
    rmsd_raw: Optional[float] = None
    rmsd_sym: Optional[float] = None
    free_energy: Optional[float] = None
    proxy_free_energy: Optional[float] = None
    soft_beta_G: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    heat_capacity: Optional[float] = None
    std_energy: Optional[float] = None
    temperature: Optional[float] = None
    scientific_provenance: ScientificProvenance = field(
        default_factory=ScientificProvenance
    )
    remarks: Dict[str, Any] = field(default_factory=dict)

    @property
    def claim_validity(self) -> ClaimValidity:
        """Strongest interpretation authorized by parsed evidence."""
        return self.scientific_provenance.claim_validity

    def __repr__(self) -> str:
        score = self.cf if self.cf is not None else self.cf_app
        parts = [f"mode={self.mode_id}", f"rank={self.pose_rank}"]
        if score is not None:
            parts.append(f"cf={score:.2f}")
        parts.append(f"path={self.path.name!r}")
        return f"<PoseResult {' '.join(parts)}>"

    def __lt__(self, other: "PoseResult") -> bool:
        """Sort by CF score (lower is better). Falls back to cf_app."""
        s = self.cf if self.cf is not None else (self.cf_app if self.cf_app is not None else float("inf"))
        o = other.cf if other.cf is not None else (other.cf_app if other.cf_app is not None else float("inf"))
        return s < o

    def __le__(self, other: "PoseResult") -> bool:
        return self == other or self < other

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PoseResult":
        """Reconstruct a PoseResult from a dictionary.

        Accepts both the internal field names (``cf``, ``rmsd_raw``) and the
        serialised names produced by :meth:`to_records`-style output
        (``best_pose_path``).  Unknown keys are silently ignored.

        Args:
            data: Dictionary with PoseResult field values.

        Returns:
            A new :class:`PoseResult` instance.
        """
        path = data.get("path", data.get("best_pose_path", ""))
        return cls(
            path=Path(path) if not isinstance(path, Path) else path,
            mode_id=data.get("mode_id", 0),
            pose_rank=data.get("pose_rank", 0),
            cf=data.get("cf"),
            cf_app=data.get("cf_app"),
            rmsd_raw=data.get("rmsd_raw"),
            rmsd_sym=data.get("rmsd_sym"),
            free_energy=data.get("free_energy"),
            proxy_free_energy=data.get("proxy_free_energy"),
            soft_beta_G=data.get("soft_beta_G"),
            enthalpy=data.get("enthalpy"),
            entropy=data.get("entropy"),
            heat_capacity=data.get("heat_capacity"),
            std_energy=data.get("std_energy"),
            temperature=data.get("temperature"),
            scientific_provenance=_as_provenance(
                data.get("scientific_provenance")
            ),
            remarks=data.get("remarks", {}),
        )


@dataclass(frozen=True)
class BindingModeResult:
    """A cluster of docked poses that share a common binding geometry.

    A binding mode aggregates :class:`PoseResult` objects that were grouped
    together by the OPTICS/DBSCAN clustering step inside the FlexAID∆S C++
    engine.  Thermodynamic quantities stored here are mode-level aggregates
    derived from the statistical mechanics engine (Helmholtz free energy,
    configurational entropy, heat capacity, etc.).

    Attributes:
        mode_id: Unique integer identifier for this binding mode.
        rank: Rank of this mode among all modes (1 = best ensemble F estimate).
        poses: Ordered list of individual poses belonging to this mode.
        free_energy: Legacy ensemble transform in the domain declared by
            ``scientific_provenance``. Current docking output is proxy-only.
        proxy_free_energy: Explicit schema-v2 CF-domain ensemble transform.
        soft_beta_G: Engine-emitted mode-election objective; lower ranks first.
        enthalpy: Ensemble mean in the declared energy domain.
        entropy: S-like diagnostic in the declared domain per kelvin.
        heat_capacity: Cv-like diagnostic in the declared domain.
        std_energy: Standard deviation of ensemble energies σ_E in the
            declared energy domain (kcal/mol only under calibrated provenance).
        best_cf: Lowest (most favourable) individual CF/contact-function scoring
            proxy within the mode (not free energy).
        frequency: Number of GA chromosomes assigned to this mode; proportional
            to Boltzmann population weight.
        temperature: Simulation temperature (K) associated with this mode.
        metadata: Arbitrary extra fields shared across all poses in the mode
            (e.g. receptor name, ligand SMILES).
    """

    mode_id: int
    rank: int
    poses: List[PoseResult]
    free_energy: Optional[float] = None
    proxy_free_energy: Optional[float] = None
    soft_beta_G: Optional[float] = None
    enthalpy: Optional[float] = None
    entropy: Optional[float] = None
    heat_capacity: Optional[float] = None
    std_energy: Optional[float] = None
    best_cf: Optional[float] = None
    frequency: Optional[int] = None
    temperature: Optional[float] = None
    scientific_provenance: ScientificProvenance = field(
        default_factory=ScientificProvenance
    )
    # Full audited ledger (rich dataclass when available from engine or pure fallback).
    # Legacy scalar fields above preserved for backward compat.
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Receptor-bound ions/cofactors in the complex that influenced this mode.
    # Format: "RESNAME:CHAIN:RESNUM" (e.g. "MG:A:101", "ZN:B:202").
    cofactors: List[str] = field(default_factory=list)
    thermodynamics: Optional[ThermodynamicBreakdown] = None

    @property
    def claim_validity(self) -> ClaimValidity:
        """Strongest interpretation authorized by mode-level evidence."""
        return self.scientific_provenance.claim_validity

    @property
    def n_poses(self) -> int:
        """Number of poses in this binding mode."""
        return len(self.poses)

    def best_pose(self) -> Optional[PoseResult]:
        """Return the pose with the lowest CF (or cf_app) score.

        Selection priority:

        1. Pose with the lowest ``cf`` value.
        2. Pose with the lowest ``cf_app`` value (if no ``cf`` is available).
        3. First pose in :attr:`poses` (fallback when no scores are present).

        Returns:
            The best-scored :class:`PoseResult`, or ``None`` if the mode is
            empty.
        """
        scored = [p for p in self.poses if p.cf is not None]
        if scored:
            return min(scored, key=lambda p: p.cf)
        scored = [p for p in self.poses if p.cf_app is not None]
        if scored:
            return min(scored, key=lambda p: p.cf_app)
        return self.poses[0] if self.poses else None

    def __repr__(self) -> str:
        parts = [f"mode_id={self.mode_id}", f"n_poses={self.n_poses}"]
        if self.free_energy is not None:
            parts.append(f"F={self.free_energy:.2f}")
        if self.best_cf is not None:
            parts.append(f"best_cf={self.best_cf:.2f}")
        return f"<BindingModeResult {' '.join(parts)}>"

    def __lt__(self, other: "BindingModeResult") -> bool:
        """Sort by emitted election objective, then legacy transform/rank."""
        if self.soft_beta_G is not None and other.soft_beta_G is not None:
            return self.soft_beta_G < other.soft_beta_G
        if self.free_energy is not None and other.free_energy is not None:
            return self.free_energy < other.free_energy
        return self.rank < other.rank

    def __le__(self, other: "BindingModeResult") -> bool:
        return self == other or self < other

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "BindingModeResult":
        """Reconstruct a BindingModeResult from a dictionary.

        Pose entries under the ``"poses"`` key are deserialised via
        :meth:`PoseResult.from_dict`.  If ``"poses"`` is absent an empty list
        is used.

        Args:
            data: Dictionary with BindingModeResult field values.

        Returns:
            A new :class:`BindingModeResult` instance.
        """
        poses = [PoseResult.from_dict(p) for p in data.get("poses", [])]
        thermo_data = data.get("thermodynamics")
        thermodynamics = (
            ThermodynamicBreakdown.from_dict(thermo_data)
            if isinstance(thermo_data, dict)
            else None
        )
        return cls(
            mode_id=data.get("mode_id", 0),
            rank=data.get("rank", 0),
            poses=poses,
            free_energy=data.get("free_energy"),
            proxy_free_energy=data.get("proxy_free_energy"),
            soft_beta_G=data.get("soft_beta_G"),
            enthalpy=data.get("enthalpy"),
            entropy=data.get("entropy"),
            heat_capacity=data.get("heat_capacity"),
            std_energy=data.get("std_energy"),
            best_cf=data.get("best_cf"),
            frequency=data.get("frequency"),
            temperature=data.get("temperature"),
            scientific_provenance=_as_provenance(
                data.get("scientific_provenance")
            ),
            metadata=data.get("metadata", {}),
            cofactors=data.get("cofactors", []),
            thermodynamics=thermodynamics,
        )


@dataclass(frozen=True)
class DockingResult:
    """Top-level container for a complete FlexAID∆S docking run.

    Returned by :func:`~flexaidds.results.load_results` after scanning a
    docking output directory.  Provides convenience methods for ranking,
    serialisation, and optional pandas integration.

    Attributes:
        source_dir: Absolute path to the directory that was scanned.
        binding_modes: List of :class:`BindingModeResult` objects, sorted by
            ascending ``mode_id``.
        temperature: Simulation temperature (K) inferred from the output files,
            or ``None`` if not available.
        metadata: Arbitrary extra information collected during loading (e.g.
            ``n_pose_files``).
    """

    source_dir: Path
    binding_modes: List[BindingModeResult]
    temperature: Optional[float] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    # ── Grand canonical / competitive augmentation (P2, additive, default=None/empty) ──
    # Present when results were produced/loaded under a competitive TargetServer/GPF context.
    # Legacy single-ligand load_results paths leave these at defaults.
    grand_log_xi: Optional[float] = None          # ln(Ξ)
    ligand_occupancies: Dict[str, float] = field(default_factory=dict)  # name -> p(bound)
    selectivities: Dict[str, float] = field(default_factory=dict)       # e.g. "A/B": ratio or log
    per_ligand_results: Dict[str, Any] = field(default_factory=dict)    # name -> summary dict or sub-result
    empty_probability: Optional[float] = None
    mean_occupancy: Optional[float] = None

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DockingResult":
        """Reconstruct a DockingResult from a dictionary.

        Accepts the format produced by :meth:`to_json` (flat
        ``binding_modes`` records from :meth:`to_records`) as well as
        nested structures where each mode contains a ``"poses"`` list.

        Args:
            data: Dictionary with DockingResult field values.

        Returns:
            A new :class:`DockingResult` instance.
        """
        raw_modes = data.get("binding_modes", [])
        modes: List[BindingModeResult] = []
        for i, m in enumerate(raw_modes):
            if "poses" in m:
                modes.append(BindingModeResult.from_dict(m))
            else:
                # Flat record from to_records(): wrap into a BindingModeResult
                modes.append(BindingModeResult(
                    mode_id=m.get("mode_id", i),
                    rank=m.get("rank", i + 1),
                    poses=[],
                    free_energy=m.get("free_energy"),
                    proxy_free_energy=m.get("proxy_free_energy"),
                    soft_beta_G=m.get("soft_beta_G"),
                    enthalpy=m.get("enthalpy"),
                    entropy=m.get("entropy"),
                    heat_capacity=m.get("heat_capacity"),
                    std_energy=m.get("std_energy"),
                    best_cf=m.get("best_cf"),
                    temperature=m.get("temperature"),
                    scientific_provenance=_as_provenance(
                        m.get("scientific_provenance")
                    ),
                    thermodynamics=_as_breakdown(m.get("thermodynamics")),
                ))
        return cls(
            source_dir=Path(data.get("source_dir", ".")),
            binding_modes=modes,
            temperature=data.get("temperature"),
            metadata=data.get("metadata", {}),
            # grand fields (additive, tolerant)
            grand_log_xi=data.get("grand_log_xi"),
            ligand_occupancies=data.get("ligand_occupancies", {}) or {},
            selectivities=data.get("selectivities", {}) or {},
            per_ligand_results=data.get("per_ligand_results", {}) or {},
            empty_probability=data.get("empty_probability"),
            mean_occupancy=data.get("mean_occupancy"),
        )

    @property
    def n_modes(self) -> int:
        """Number of binding modes in this result."""
        return len(self.binding_modes)

    def __repr__(self) -> str:
        parts = [f"n_modes={self.n_modes}"]
        if self.temperature is not None:
            parts.append(f"T={self.temperature:.0f}K")
        parts.append(f"source={self.source_dir.name!r}")
        return f"<DockingResult {' '.join(parts)}>"

    def top_mode(self) -> Optional[BindingModeResult]:
        """Return the binding mode with the lowest emitted election objective.

        Falls back to the legacy ensemble transform and then stored rank when
        schema-v2 ``soft_beta_G`` is absent.

        Returns:
            Best :class:`BindingModeResult`, or ``None`` if there are no modes.
        """
        if not self.binding_modes:
            return None
        # Prefer the engine-emitted objective; old results retain their former
        # free-energy/rank behavior.
        def _score(m):
            election = m.soft_beta_G
            fe = m.free_energy if m.free_energy is not None else float('inf')
            npos = m.n_poses if m.n_poses > 0 else 0
            tmatch = 0
            if getattr(self, 'temperature', None) and m.temperature:
                tmatch = 1 if abs(m.temperature - self.temperature) < 0.1 else -1
            return (
                0 if election is not None else 1,
                election if election is not None else fe,
                -npos,
                -tmatch,
                m.rank,
            )
        sane = [
            m for m in self.binding_modes
            if m.n_poses > 0 or m.free_energy is not None or m.soft_beta_G is not None
        ]
        if sane:
            return min(sane, key=_score)
        return min(self.binding_modes, key=lambda m: m.rank)

    def to_records(self) -> List[Dict[str, Any]]:
        """Serialise all binding modes to a list of flat dictionaries.

        Each dictionary contains mode-level scalar fields plus the path to the
        best pose.  Suitable for direct conversion to a
        :class:`pandas.DataFrame` via :meth:`to_dataframe`.

        Returns:
            List of dictionaries, one per binding mode, with keys:
            ``mode_id``, ``rank``, ``n_poses``, ``free_energy``,
            ``proxy_free_energy``, ``soft_beta_G``, ``enthalpy``, ``entropy``,
            ``heat_capacity``, ``std_energy``, ``best_cf``, ``temperature``,
            ``scientific_provenance``, ``thermodynamics``, ``best_pose_path``.

            Every value is JSON/CSV-serialisable: the provenance and ledger
            objects are emitted as plain mappings, never as dataclass objects.
        """
        records: List[Dict[str, Any]] = []
        for mode in self.binding_modes:
            best_pose = mode.best_pose()
            records.append(
                {
                    "mode_id": mode.mode_id,
                    "rank": mode.rank,
                    "n_poses": mode.n_poses,
                    "free_energy": mode.free_energy,
                    "proxy_free_energy": mode.proxy_free_energy,
                    "soft_beta_G": mode.soft_beta_G,
                    "enthalpy": mode.enthalpy,
                    "entropy": mode.entropy,
                    "heat_capacity": mode.heat_capacity,
                    "std_energy": mode.std_energy,
                    "best_cf": mode.best_cf,
                    "temperature": mode.temperature,
                    "scientific_provenance": mode.scientific_provenance.to_dict(),
                    # Serialisable mapping, not the dataclass object: records
                    # feed JSON/CSV/pandas and must not leak live objects.
                    "thermodynamics": (
                        mode.thermodynamics.to_dict()
                        if mode.thermodynamics is not None
                        else None
                    ),
                    "best_pose_path": str(best_pose.path) if best_pose else None,
                }
            )
        return records

    @staticmethod
    def _binding_mode_json_record(mode: BindingModeResult) -> Dict[str, Any]:
        """Return a JSON record with legacy flat fields plus new nested data."""
        best_pose = mode.best_pose()
        record: Dict[str, Any] = {
            "mode_id": mode.mode_id,
            "rank": mode.rank,
            "n_poses": mode.n_poses,
            "free_energy": mode.free_energy,
            "proxy_free_energy": mode.proxy_free_energy,
            "soft_beta_G": mode.soft_beta_G,
            "enthalpy": mode.enthalpy,
            "entropy": mode.entropy,
            "heat_capacity": mode.heat_capacity,
            "std_energy": mode.std_energy,
            "best_cf": mode.best_cf,
            "temperature": mode.temperature,
            "scientific_provenance": mode.scientific_provenance.to_dict(),
            "best_pose_path": str(best_pose.path) if best_pose else None,
        }
        if mode.thermodynamics is not None:
            record["thermodynamics"] = mode.thermodynamics.to_dict()
        return record

    def to_dataframe(self):
        """Convert binding-mode results to a :class:`pandas.DataFrame`.

        Each row corresponds to one binding mode.  Columns match the fields
        returned by :meth:`to_records`.

        Raises:
            ImportError: If ``pandas`` is not installed.  Use
                :meth:`to_records` for a dependency-free alternative.

        Returns:
            :class:`pandas.DataFrame` with one row per binding mode.
        """
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "pandas is required for DockingResult.to_dataframe(); use to_records() instead."
            ) from exc
        return pd.DataFrame(self.to_records())

    def to_json(self, path: Union[str, Path, None] = None, **kwargs) -> Optional[str]:
        """Serialise docking results to JSON.

        The output includes the source directory, temperature, metadata, and
        a ``binding_modes`` array produced by :meth:`to_records`.

        Args:
            path: Destination file path.  When *None* the JSON text is returned
                  as a string instead of being written to disk.
            **kwargs: Extra keyword arguments forwarded to :func:`json.dumps`
                (e.g. ``indent``, ``sort_keys``).

        Returns:
            JSON text when *path* is ``None``, otherwise ``None``.
        """
        payload = {
            "source_dir": str(self.source_dir),
            "temperature": self.temperature,
            "n_modes": self.n_modes,
            "metadata": self.metadata,
            "binding_modes": [
                self._binding_mode_json_record(mode)
                for mode in self.binding_modes
            ],
            # Grand canonical fields (only emitted if populated — additive, non-breaking)
            "grand_log_xi": self.grand_log_xi,
            "ligand_occupancies": self.ligand_occupancies,
            "selectivities": self.selectivities,
            "per_ligand_results": self.per_ligand_results,
            "empty_probability": self.empty_probability,
            "mean_occupancy": self.mean_occupancy,
        }
        kwargs.setdefault("indent", 2)
        text = json.dumps(payload, **kwargs)

        if path is None:
            return text

        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.write("\n")
        return None

    @classmethod
    def from_json(
        cls, source: Union[str, Path], *, source_dir: Union[str, Path, None] = None
    ) -> "DockingResult":
        """Load a :class:`DockingResult` from JSON produced by :meth:`to_json`.

        Accepts either a file path or a raw JSON string.  The binding-mode
        records are reconstructed into :class:`BindingModeResult` objects
        (each with a single :class:`PoseResult` placeholder pointing to the
        best pose path, when available).

        Args:
            source: Path to a JSON file, or a JSON string.
            source_dir: Override the ``source_dir`` stored in the JSON payload.
                Useful when the original output directory has moved.

        Returns:
            Reconstructed :class:`DockingResult`.

        Raises:
            json.JSONDecodeError: If the input is not valid JSON.
            KeyError: If required fields are missing from the JSON payload.
        """
        # Robust path-vs-content detection: avoid Path(<long-json-str>) on macOS
        # (and other platforms) which raises ENAMETOOLONG for >255 char components.
        text = None
        s = str(source).strip()
        looks_like_json = (s.startswith("{") or s.startswith("[")) and len(s) > 2
        if looks_like_json:
            text = str(source)
        else:
            try:
                source_path = Path(source)
                if source_path.is_file():
                    text = source_path.read_text(encoding="utf-8")
                else:
                    text = str(source)
            except (OSError, FileNotFoundError, RuntimeError):
                text = str(source)
        if text is None:
            text = str(source)

        payload = json.loads(text)
        resolved_dir = Path(source_dir) if source_dir else Path(payload["source_dir"])

        modes: List[BindingModeResult] = []
        for rec in payload.get("binding_modes", []):
            best_path = rec.get("best_pose_path")
            thermo_data = rec.get("thermodynamics")
            thermodynamics = (
                ThermodynamicBreakdown.from_dict(thermo_data)
                if isinstance(thermo_data, dict)
                else None
            )
            # Rebuilt from the record's own evidence fields; a serialized
            # ``claim_validity`` is discarded by _as_provenance.
            provenance = _as_provenance(rec.get("scientific_provenance"))
            poses: List[PoseResult] = []
            if best_path is not None:
                poses.append(
                    PoseResult(
                        path=Path(best_path),
                        mode_id=rec["mode_id"],
                        pose_rank=1,
                        cf=rec.get("best_cf"),
                        free_energy=rec.get("free_energy"),
                        proxy_free_energy=rec.get("proxy_free_energy"),
                        soft_beta_G=rec.get("soft_beta_G"),
                        enthalpy=rec.get("enthalpy"),
                        entropy=rec.get("entropy"),
                        heat_capacity=rec.get("heat_capacity"),
                        std_energy=rec.get("std_energy"),
                        temperature=rec.get("temperature"),
                        scientific_provenance=provenance,
                    )
                )
            modes.append(
                BindingModeResult(
                    mode_id=rec["mode_id"],
                    rank=rec["rank"],
                    poses=poses,
                    free_energy=rec.get("free_energy"),
                    proxy_free_energy=rec.get("proxy_free_energy"),
                    soft_beta_G=rec.get("soft_beta_G"),
                    enthalpy=rec.get("enthalpy"),
                    entropy=rec.get("entropy"),
                    heat_capacity=rec.get("heat_capacity"),
                    std_energy=rec.get("std_energy"),
                    best_cf=rec.get("best_cf"),
                    temperature=rec.get("temperature"),
                    scientific_provenance=provenance,
                    thermodynamics=thermodynamics,
                )
            )

        return cls(
            source_dir=resolved_dir,
            binding_modes=modes,
            temperature=payload.get("temperature"),
            metadata=payload.get("metadata", {}),
            # grand (tolerant defaults for old JSONs)
            grand_log_xi=payload.get("grand_log_xi"),
            ligand_occupancies=payload.get("ligand_occupancies") or {},
            selectivities=payload.get("selectivities") or {},
            per_ligand_results=payload.get("per_ligand_results") or {},
            empty_probability=payload.get("empty_probability"),
            mean_occupancy=payload.get("mean_occupancy"),
        )

    def to_csv(self, path: Union[str, Path, None] = None) -> Optional[str]:
        """Write binding mode summary to CSV.

        Args:
            path: Destination file path.  When *None* the CSV text is returned
                  as a string instead of being written to disk.

        Returns:
            CSV text when *path* is ``None``, otherwise ``None``.
        """
        records = self.to_records()
        if not records:
            fieldnames: List[str] = []
        else:
            fieldnames = list(records[0].keys())

        if path is None:
            buf = io.StringIO()
            writer = csv.DictWriter(buf, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
            return buf.getvalue()

        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        return None

    @classmethod
    def from_csv(cls, source: Union[str, Path]) -> "DockingResult":
        """Load a DockingResult from a CSV file or string.

        Accepts either a file path or raw CSV text.  The CSV format is the
        one produced by :meth:`to_csv` (flat binding-mode records).

        Numeric fields are coerced from their string representation; empty
        strings and the literal ``"None"`` are treated as ``None``.

        Args:
            source: Path to a ``.csv`` file, or a CSV-encoded string.

        Returns:
            A new :class:`DockingResult` instance.
        """
        # Path-vs-content detection must not call Path().exists() on raw CSV
        # text: a serialized provenance/ledger column easily exceeds the OS
        # path-component limit and raises ENAMETOOLONG (see from_json).
        if isinstance(source, Path):
            text = source.read_text(encoding="utf-8") if source.is_file() else str(source)
        else:
            raw = str(source)
            looks_like_csv_text = "\n" in raw or "\r" in raw or len(raw) > 200
            text = raw
            if not looks_like_csv_text:
                try:
                    candidate = Path(raw)
                    if candidate.is_file():
                        text = candidate.read_text(encoding="utf-8")
                except (OSError, ValueError):
                    text = raw

        reader = csv.DictReader(io.StringIO(text))
        records = []
        for row in reader:
            coerced: Dict[str, Any] = {}
            for key, value in row.items():
                coerced[key] = cls._coerce_csv_value(key, value)
            records.append(coerced)

        modes: List[BindingModeResult] = []
        for i, rec in enumerate(records):
            modes.append(BindingModeResult(
                mode_id=rec.get("mode_id", i),
                rank=rec.get("rank", i + 1),
                poses=[],
                free_energy=rec.get("free_energy"),
                proxy_free_energy=rec.get("proxy_free_energy"),
                soft_beta_G=rec.get("soft_beta_G"),
                enthalpy=rec.get("enthalpy"),
                entropy=rec.get("entropy"),
                heat_capacity=rec.get("heat_capacity"),
                std_energy=rec.get("std_energy"),
                best_cf=rec.get("best_cf"),
                temperature=rec.get("temperature"),
                # _as_provenance re-derives validity from the evidence fields;
                # a ``claim_validity`` column can never authorize a claim.
                scientific_provenance=_as_provenance(
                    rec.get("scientific_provenance")
                ),
                thermodynamics=_as_breakdown(rec.get("thermodynamics")),
            ))

        return cls(
            source_dir=Path("."),
            binding_modes=modes,
        )

    @staticmethod
    def _coerce_csv_value(key: str, value: str) -> Any:
        """Coerce a CSV string value to the appropriate Python type."""
        if value is None or value == "" or value == "None":
            return None
        _int_keys = {"mode_id", "rank", "n_poses"}
        if key in _int_keys:
            try:
                return int(float(value))
            except (ValueError, TypeError):
                return value
        _float_keys = {
            "free_energy", "proxy_free_energy", "soft_beta_G",
            "enthalpy", "entropy", "heat_capacity",
            "std_energy", "best_cf", "temperature",
        }
        if key in _float_keys:
            try:
                return float(value)
            except (ValueError, TypeError):
                return value
        _mapping_keys = {"scientific_provenance", "thermodynamics"}
        if key in _mapping_keys:
            # csv stores these columns as a repr/JSON blob. Parse with a
            # literal-only evaluator (never eval) and accept a mapping only;
            # anything else falls through and is rejected downstream.
            for parse in (ast.literal_eval, json.loads):
                try:
                    parsed = parse(value)
                except (ValueError, SyntaxError, TypeError, MemoryError,
                        RecursionError):
                    continue
                if isinstance(parsed, dict):
                    return parsed
            return None
        return value
