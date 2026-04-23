"""Comparative benchmarking of FlexAIDdS vs Boltz-2 for protein-ligand prediction.

Provides a systematic framework to evaluate pose accuracy, affinity correlation,
and enrichment metrics across a dataset of protein-ligand systems with
experimental ground truth.
"""

from __future__ import annotations

import io as _io
import json
import math
import os
import re
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np


# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

_R_kcal = 0.001987206  # kcal mol⁻¹ K⁻¹


# ---------------------------------------------------------------------------
# Hardware-aware worker detection
# ---------------------------------------------------------------------------


def auto_workers() -> int:
    """Detect optimal worker count for parallel benchmark execution.

    Heuristic: ``min(cpu_count // 2, mem_gb // 2.5)``, clamped to ``[1, 16]``.

    On an M3 Pro (12 cores, 18 GB) this returns 6, reserving efficiency cores
    for the OS and nested C++/OpenMP threads inside each docking subprocess.
    """
    cpu = os.cpu_count() or 4
    try:
        mem_gb = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES") / (1024 ** 3)
    except (AttributeError, ValueError, OSError):
        mem_gb = 16.0
    by_cpu = cpu // 2
    by_mem = int(mem_gb / 2.5)
    return max(1, min(by_cpu, by_mem, 16))


# ---------------------------------------------------------------------------
# Affinity conversions
# ---------------------------------------------------------------------------


def ki_to_dg(ki_nM: float, temperature_K: float = 298.15) -> float:
    """Convert Ki (nM) to binding free energy ΔG (kcal/mol).

    ΔG = RT ln(Ki), where Ki is in molar.
    """
    ki_M = ki_nM * 1e-9
    return _R_kcal * temperature_K * math.log(ki_M)


def ic50_to_dg(ic50_nM: float, temperature_K: float = 298.15) -> float:
    """Convert IC50 (nM) to approximate ΔG (kcal/mol).

    Uses the Cheng-Prusoff approximation: Ki ≈ IC50 / 2.
    """
    return ki_to_dg(ic50_nM / 2.0, temperature_K)


def pic50_to_dg(pic50: float, temperature_K: float = 298.15) -> float:
    """Convert pIC50 to ΔG (kcal/mol).

    pIC50 = -log10(IC50_M), so IC50_M = 10^(-pIC50).
    Uses Cheng-Prusoff: Ki ≈ IC50 / 2.
    """
    ic50_M = 10.0 ** (-pic50)
    ic50_nM = ic50_M * 1e9
    return ic50_to_dg(ic50_nM, temperature_K)


# ---------------------------------------------------------------------------
# Coordinate extraction
# ---------------------------------------------------------------------------


def extract_ligand_coords_from_pdb(pdb_path: Path) -> np.ndarray:
    """Extract ligand heavy-atom coordinates from a PDB file.

    Selects HETATM records excluding HOH, returns (N, 3) array sorted
    by atom name for deterministic ordering.
    """
    atoms = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("HETATM"):
                continue
            resname = line[17:20].strip()
            if resname == "HOH":
                continue
            element = line[76:78].strip() if len(line) > 77 else ""
            if element == "H":
                continue
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append((atom_name, x, y, z))

    if not atoms:
        raise ValueError(f"No ligand heavy atoms found in {pdb_path}")

    atoms.sort(key=lambda a: a[0])
    coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=np.float64)
    return coords


def extract_ligand_coords_from_mmcif(mmcif_string: str, ligand_id: str = "") -> np.ndarray:
    """Extract ligand heavy-atom coordinates from a Boltz-2 mmCIF string.

    Parses the ``_atom_site`` loop to find HETATM records for the ligand
    component.  Returns (N, 3) array sorted by atom name.
    """
    lines = mmcif_string.splitlines()
    # Find _atom_site loop
    in_atom_site = False
    columns: List[str] = []
    data_rows: List[List[str]] = []

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "loop_":
            # Check if next lines are _atom_site.* headers
            j = i + 1
            headers: List[str] = []
            while j < len(lines) and lines[j].strip().startswith("_atom_site."):
                headers.append(lines[j].strip())
                j += 1
            if headers:
                columns = [h.split(".")[1] for h in headers]
                in_atom_site = True
                i = j
                continue
        if in_atom_site:
            if line.startswith("_") or line.startswith("loop_") or line == "#":
                break
            if line:
                tokens = _tokenize_mmcif_line(line)
                if len(tokens) == len(columns):
                    data_rows.append(tokens)
        i += 1

    if not columns:
        raise ValueError("No _atom_site loop found in mmCIF data.")

    col_map = {name: idx for idx, name in enumerate(columns)}

    # Required columns
    for req in ("Cartn_x", "Cartn_y", "Cartn_z"):
        if req not in col_map:
            raise ValueError(f"Missing {req} column in _atom_site.")

    group_idx = col_map.get("group_PDB")
    comp_idx = col_map.get("label_comp_id")
    asym_idx = col_map.get("label_asym_id")
    type_idx = col_map.get("type_symbol")
    atom_idx = col_map.get("label_atom_id", col_map.get("auth_atom_id"))
    x_idx = col_map["Cartn_x"]
    y_idx = col_map["Cartn_y"]
    z_idx = col_map["Cartn_z"]

    atoms = []
    for row in data_rows:
        # Filter to HETATM or ligand chain
        if group_idx is not None and row[group_idx] != "HETATM":
            continue
        if comp_idx is not None and row[comp_idx] == "HOH":
            continue
        if ligand_id and asym_idx is not None and row[asym_idx] != ligand_id:
            continue
        # Skip hydrogens
        if type_idx is not None and row[type_idx] == "H":
            continue

        atom_name = row[atom_idx] if atom_idx is not None else ""
        x = float(row[x_idx])
        y = float(row[y_idx])
        z = float(row[z_idx])
        atoms.append((atom_name, x, y, z))

    if not atoms:
        raise ValueError("No ligand heavy atoms found in mmCIF data.")

    atoms.sort(key=lambda a: a[0])
    return np.array([[a[1], a[2], a[3]] for a in atoms], dtype=np.float64)


def _tokenize_mmcif_line(line: str) -> List[str]:
    """Tokenize an mmCIF data line, handling quoted strings."""
    tokens: List[str] = []
    i = 0
    while i < len(line):
        if line[i] in (" ", "\t"):
            i += 1
            continue
        if line[i] in ("'", '"'):
            quote = line[i]
            j = i + 1
            while j < len(line) and line[j] != quote:
                j += 1
            tokens.append(line[i + 1 : j])
            i = j + 1
        else:
            j = i
            while j < len(line) and line[j] not in (" ", "\t"):
                j += 1
            tokens.append(line[i:j])
            i = j
    return tokens


# ---------------------------------------------------------------------------
# RMSD computation
# ---------------------------------------------------------------------------


def compute_rmsd(coords_pred: np.ndarray, coords_ref: np.ndarray) -> float:
    """Compute RMSD after optimal Kabsch alignment.

    Both arrays must be (N, 3).  Returns RMSD in Angstrom.
    """
    if coords_pred.shape != coords_ref.shape:
        raise ValueError(
            f"Shape mismatch: predicted {coords_pred.shape} vs "
            f"reference {coords_ref.shape}."
        )
    n = coords_pred.shape[0]
    if n == 0:
        return 0.0

    # Center both
    center_pred = coords_pred.mean(axis=0)
    center_ref = coords_ref.mean(axis=0)
    p = coords_pred - center_pred
    q = coords_ref - center_ref

    # Kabsch: find optimal rotation via SVD
    H = p.T @ q
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1.0, 1.0, d])
    R = Vt.T @ sign_matrix @ U.T

    p_rotated = p @ R.T
    diff = p_rotated - q
    rmsd = float(np.sqrt((diff * diff).sum() / n))
    return rmsd


# ---------------------------------------------------------------------------
# Statistical metrics (stdlib-only)
# ---------------------------------------------------------------------------


def _ranks(values: List[float]) -> List[float]:
    """Assign ranks to values (1-based, average ties)."""
    indexed = sorted(enumerate(values), key=lambda x: x[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        avg_rank = (i + j + 1) / 2.0  # 1-based average
        for k in range(i, j):
            ranks[indexed[k][0]] = avg_rank
        i = j
    return ranks


def _pearson_r(x: List[float], y: List[float]) -> float:
    """Pearson correlation coefficient."""
    n = len(x)
    if n < 2:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    cov = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    sx = math.sqrt(sum((xi - mx) ** 2 for xi in x))
    sy = math.sqrt(sum((yi - my) ** 2 for yi in y))
    if sx == 0.0 or sy == 0.0:
        return 0.0
    return cov / (sx * sy)


def spearman_rho(x: List[float], y: List[float]) -> float:
    """Spearman rank correlation coefficient."""
    if len(x) != len(y) or len(x) < 2:
        return 0.0
    return _pearson_r(_ranks(x), _ranks(y))


def kendall_tau(x: List[float], y: List[float]) -> float:
    """Kendall's tau-b rank correlation."""
    n = len(x)
    if n < 2:
        return 0.0
    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            prod = dx * dy
            if prod > 0:
                concordant += 1
            elif prod < 0:
                discordant += 1
    denom = concordant + discordant
    if denom == 0:
        return 0.0
    return (concordant - discordant) / denom


def r_squared(x: List[float], y: List[float]) -> float:
    """Coefficient of determination R²."""
    r = _pearson_r(x, y)
    return r * r


def roc_auc(
    scores: List[float], labels: List[bool], higher_is_better: bool = True
) -> float:
    """ROC-AUC via trapezoidal integration."""
    n = len(scores)
    if n == 0:
        return 0.0

    paired = list(zip(scores, labels))
    paired.sort(key=lambda p: p[0], reverse=higher_is_better)

    n_pos = sum(1 for _, lbl in paired if lbl)
    n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return 0.0

    tp = 0
    fp = 0
    prev_tp = 0
    prev_fp = 0
    auc = 0.0
    prev_score = None

    for score, label in paired:
        if score != prev_score and prev_score is not None:
            auc += (fp - prev_fp) * (tp + prev_tp) / 2.0
            prev_tp = tp
            prev_fp = fp
        if label:
            tp += 1
        else:
            fp += 1
        prev_score = score

    auc += (fp - prev_fp) * (tp + prev_tp) / 2.0
    return auc / (n_pos * n_neg)


def enrichment_factor(
    scores: List[float],
    labels: List[bool],
    fraction: float = 0.01,
    higher_is_better: bool = True,
) -> float:
    """Enrichment factor at a given fraction of the ranked list."""
    n = len(scores)
    if n == 0:
        return 0.0

    n_pos = sum(1 for lbl in labels if lbl)
    if n_pos == 0:
        return 0.0

    paired = list(zip(scores, labels))
    paired.sort(key=lambda p: p[0], reverse=higher_is_better)

    cutoff = max(1, int(math.ceil(n * fraction)))
    hits = sum(1 for _, lbl in paired[:cutoff] if lbl)
    expected = n_pos * fraction
    if expected == 0.0:
        return 0.0
    return hits / expected


# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BenchmarkSystem:
    """One protein-ligand system in the benchmark dataset.

    Attributes:
        system_id: Unique identifier (e.g. ``"1HVR_indinavir"``).
        protein_pdb_path: Path to receptor PDB file (FlexAIDdS input).
        protein_sequence: Amino acid sequence (Boltz-2 input).
        ligand_mol2_path: Path to ligand MOL2 file (FlexAIDdS input).
        ligand_smiles: SMILES string (Boltz-2 input).
        reference_pose_pdb_path: Path to experimental reference pose PDB.
        experimental_dg_kcal_mol: Experimental ΔG (kcal/mol).
        experimental_ki_nM: Experimental Ki (nM).
        experimental_ic50_nM: Experimental IC50 (nM).
        is_active: Whether the ligand is an active binder (for enrichment).
        pocket_residues: 1-based residue indices for pocket constraints.
        metadata: Arbitrary extra fields.
    """

    system_id: str
    protein_pdb_path: Path
    protein_sequence: str
    ligand_mol2_path: Path
    ligand_smiles: str
    reference_pose_pdb_path: Path
    experimental_dg_kcal_mol: Optional[float] = None
    experimental_ki_nM: Optional[float] = None
    experimental_ic50_nM: Optional[float] = None
    is_active: Optional[bool] = None
    pocket_residues: Optional[Tuple[int, ...]] = None
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "system_id": self.system_id,
            "protein_pdb_path": str(self.protein_pdb_path),
            "protein_sequence": self.protein_sequence,
            "ligand_mol2_path": str(self.ligand_mol2_path),
            "ligand_smiles": self.ligand_smiles,
            "reference_pose_pdb_path": str(self.reference_pose_pdb_path),
        }
        if self.experimental_dg_kcal_mol is not None:
            d["experimental_dg_kcal_mol"] = self.experimental_dg_kcal_mol
        if self.experimental_ki_nM is not None:
            d["experimental_ki_nM"] = self.experimental_ki_nM
        if self.experimental_ic50_nM is not None:
            d["experimental_ic50_nM"] = self.experimental_ic50_nM
        if self.is_active is not None:
            d["is_active"] = self.is_active
        if self.pocket_residues is not None:
            d["pocket_residues"] = list(self.pocket_residues)
        if self.metadata:
            d["metadata"] = self.metadata
        return d

    @classmethod
    def from_dict(cls, data: Dict[str, Any], base_dir: Optional[Path] = None) -> "BenchmarkSystem":
        """Reconstruct from a dictionary.  Paths resolved relative to *base_dir*."""
        bd = base_dir or Path(".")
        pocket = data.get("pocket_residues")
        return cls(
            system_id=data["system_id"],
            protein_pdb_path=bd / data["protein_pdb_path"],
            protein_sequence=data["protein_sequence"],
            ligand_mol2_path=bd / data["ligand_mol2_path"],
            ligand_smiles=data["ligand_smiles"],
            reference_pose_pdb_path=bd / data["reference_pose_pdb_path"],
            experimental_dg_kcal_mol=data.get("experimental_dg_kcal_mol"),
            experimental_ki_nM=data.get("experimental_ki_nM"),
            experimental_ic50_nM=data.get("experimental_ic50_nM"),
            is_active=data.get("is_active"),
            pocket_residues=tuple(pocket) if pocket else None,
            metadata=data.get("metadata", {}),
        )


@dataclass(frozen=True)
class MethodResult:
    """Result from one method on one benchmark system.

    Attributes:
        method: ``"flexaidds"`` or ``"boltz2"``.
        system_id: Links to :class:`BenchmarkSystem`.
        best_pose_rmsd_angstrom: RMSD of best pose vs experimental reference.
        predicted_dg_kcal_mol: Predicted ΔG on a common scale.
        predicted_score: Raw score in method-native units.
        n_poses: Number of poses or diffusion samples generated.
        wall_time_seconds: Wall-clock time for this system.
    """

    method: str
    system_id: str
    best_pose_rmsd_angstrom: Optional[float] = None
    predicted_dg_kcal_mol: Optional[float] = None
    predicted_score: Optional[float] = None
    n_poses: int = 0
    wall_time_seconds: float = 0.0


@dataclass(frozen=True)
class SystemBenchmarkResult:
    """Comparison of both methods on a single system."""

    system: BenchmarkSystem
    flexaidds_result: Optional[MethodResult] = None
    boltz2_result: Optional[MethodResult] = None


@dataclass(frozen=True)
class BenchmarkSummary:
    """Aggregate statistics from a benchmark run.

    Attributes:
        n_systems: Total number of systems.
        flexaidds_success_rate: Fraction with RMSD < threshold.
        boltz2_success_rate: Fraction with RMSD < threshold.
        flexaidds_median_rmsd_angstrom: Median best-pose RMSD.
        boltz2_median_rmsd_angstrom: Median best-pose RMSD.
        rank_correlation_spearman: Spearman rho between methods' score rankings.
        rank_correlation_kendall: Kendall tau between methods' score rankings.
        flexaidds_dg_r_squared: R² of FlexAIDdS ΔG vs experimental ΔG.
        boltz2_dg_r_squared: R² of Boltz-2 ΔG vs experimental ΔG.
        flexaidds_roc_auc: ROC-AUC for active/decoy discrimination.
        boltz2_roc_auc: ROC-AUC for active/decoy discrimination.
        flexaidds_enrichment_factor_1pct: Enrichment factor at 1%.
        boltz2_enrichment_factor_1pct: Enrichment factor at 1%.
        flexaidds_mean_time_seconds: Mean wall-clock per system.
        boltz2_mean_time_seconds: Mean wall-clock per system.
    """

    n_systems: int = 0
    flexaidds_success_rate: Optional[float] = None
    boltz2_success_rate: Optional[float] = None
    flexaidds_median_rmsd_angstrom: Optional[float] = None
    boltz2_median_rmsd_angstrom: Optional[float] = None
    rank_correlation_spearman: Optional[float] = None
    rank_correlation_kendall: Optional[float] = None
    flexaidds_dg_r_squared: Optional[float] = None
    boltz2_dg_r_squared: Optional[float] = None
    flexaidds_roc_auc: Optional[float] = None
    boltz2_roc_auc: Optional[float] = None
    flexaidds_enrichment_factor_1pct: Optional[float] = None
    boltz2_enrichment_factor_1pct: Optional[float] = None
    flexaidds_mean_time_seconds: Optional[float] = None
    boltz2_mean_time_seconds: Optional[float] = None


@dataclass(frozen=True)
class BenchmarkResult:
    """Aggregate benchmark results across all systems.

    Attributes:
        systems: Per-system paired results.
        rmsd_threshold_angstrom: Success threshold for pose accuracy.
        metadata: Run metadata (date, versions, hardware, etc.).
    """

    systems: Tuple[SystemBenchmarkResult, ...]
    rmsd_threshold_angstrom: float = 2.0
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def n_systems(self) -> int:
        return len(self.systems)

    def summary(self) -> BenchmarkSummary:
        """Compute aggregate metrics across all systems."""
        fa_rmsds: List[float] = []
        b2_rmsds: List[float] = []
        fa_dgs: List[float] = []
        b2_dgs: List[float] = []
        exp_dgs_fa: List[float] = []
        exp_dgs_b2: List[float] = []
        cross_fa: List[float] = []
        cross_b2: List[float] = []
        fa_scores_enrich: List[float] = []
        b2_scores_enrich: List[float] = []
        labels_fa: List[bool] = []
        labels_b2: List[bool] = []
        fa_times: List[float] = []
        b2_times: List[float] = []

        for sr in self.systems:
            fa = sr.flexaidds_result
            b2 = sr.boltz2_result
            exp_dg = sr.system.experimental_dg_kcal_mol

            if fa is not None:
                fa_times.append(fa.wall_time_seconds)
                if fa.best_pose_rmsd_angstrom is not None:
                    fa_rmsds.append(fa.best_pose_rmsd_angstrom)
                if fa.predicted_dg_kcal_mol is not None and exp_dg is not None:
                    fa_dgs.append(fa.predicted_dg_kcal_mol)
                    exp_dgs_fa.append(exp_dg)
                if sr.system.is_active is not None and fa.predicted_score is not None:
                    # FlexAID: lower CF = better → negate for ROC
                    fa_scores_enrich.append(-fa.predicted_score)
                    labels_fa.append(sr.system.is_active)

            if b2 is not None:
                b2_times.append(b2.wall_time_seconds)
                if b2.best_pose_rmsd_angstrom is not None:
                    b2_rmsds.append(b2.best_pose_rmsd_angstrom)
                if b2.predicted_dg_kcal_mol is not None and exp_dg is not None:
                    b2_dgs.append(b2.predicted_dg_kcal_mol)
                    exp_dgs_b2.append(exp_dg)
                if sr.system.is_active is not None and b2.predicted_score is not None:
                    b2_scores_enrich.append(b2.predicted_score)
                    labels_b2.append(sr.system.is_active)

            # Cross-method rank correlation (systems where both have scores)
            if (fa is not None and b2 is not None
                    and fa.predicted_score is not None
                    and b2.predicted_score is not None):
                cross_fa.append(fa.predicted_score)
                cross_b2.append(b2.predicted_score)

        threshold = self.rmsd_threshold_angstrom

        def _success(rmsds: List[float]) -> Optional[float]:
            return sum(1 for r in rmsds if r < threshold) / len(rmsds) if rmsds else None

        def _median(vals: List[float]) -> Optional[float]:
            if not vals:
                return None
            s = sorted(vals)
            mid = len(s) // 2
            return (s[mid] + s[mid - 1]) / 2.0 if len(s) % 2 == 0 else s[mid]

        def _mean(vals: List[float]) -> Optional[float]:
            return sum(vals) / len(vals) if vals else None

        # For cross-method Spearman: negate FlexAID so both higher=better
        cross_fa_neg = [-v for v in cross_fa]

        return BenchmarkSummary(
            n_systems=len(self.systems),
            flexaidds_success_rate=_success(fa_rmsds),
            boltz2_success_rate=_success(b2_rmsds),
            flexaidds_median_rmsd_angstrom=_median(fa_rmsds),
            boltz2_median_rmsd_angstrom=_median(b2_rmsds),
            rank_correlation_spearman=(
                spearman_rho(cross_fa_neg, cross_b2) if len(cross_fa) >= 2 else None
            ),
            rank_correlation_kendall=(
                kendall_tau(cross_fa_neg, cross_b2) if len(cross_fa) >= 2 else None
            ),
            flexaidds_dg_r_squared=(
                r_squared(fa_dgs, exp_dgs_fa) if len(fa_dgs) >= 2 else None
            ),
            boltz2_dg_r_squared=(
                r_squared(b2_dgs, exp_dgs_b2) if len(b2_dgs) >= 2 else None
            ),
            flexaidds_roc_auc=(
                roc_auc(fa_scores_enrich, labels_fa) if labels_fa else None
            ),
            boltz2_roc_auc=(
                roc_auc(b2_scores_enrich, labels_b2) if labels_b2 else None
            ),
            flexaidds_enrichment_factor_1pct=(
                enrichment_factor(fa_scores_enrich, labels_fa, 0.01) if labels_fa else None
            ),
            boltz2_enrichment_factor_1pct=(
                enrichment_factor(b2_scores_enrich, labels_b2, 0.01) if labels_b2 else None
            ),
            flexaidds_mean_time_seconds=_mean(fa_times),
            boltz2_mean_time_seconds=_mean(b2_times),
        )

    def to_records(self) -> List[Dict[str, Any]]:
        """One record per system with both methods' results side-by-side."""
        records = []
        for sr in self.systems:
            rec: Dict[str, Any] = {
                "system_id": sr.system.system_id,
                "experimental_dg_kcal_mol": sr.system.experimental_dg_kcal_mol,
                "is_active": sr.system.is_active,
            }
            for prefix, mr in [("flexaidds", sr.flexaidds_result),
                               ("boltz2", sr.boltz2_result)]:
                if mr is not None:
                    rec[f"{prefix}_rmsd"] = mr.best_pose_rmsd_angstrom
                    rec[f"{prefix}_dg"] = mr.predicted_dg_kcal_mol
                    rec[f"{prefix}_score"] = mr.predicted_score
                    rec[f"{prefix}_n_poses"] = mr.n_poses
                    rec[f"{prefix}_time"] = mr.wall_time_seconds
                else:
                    rec[f"{prefix}_rmsd"] = None
                    rec[f"{prefix}_dg"] = None
                    rec[f"{prefix}_score"] = None
                    rec[f"{prefix}_n_poses"] = None
                    rec[f"{prefix}_time"] = None
            records.append(rec)
        return records

    def to_json(self, path: Union[str, Path, None] = None, **kwargs) -> Optional[str]:
        """Serialize benchmark results to JSON."""
        summary = self.summary()
        payload = {
            "n_systems": self.n_systems,
            "rmsd_threshold_angstrom": self.rmsd_threshold_angstrom,
            "metadata": self.metadata,
            "summary": {
                k: v for k, v in summary.__dict__.items() if v is not None
            },
            "systems": self.to_records(),
        }
        kwargs.setdefault("indent", 2)
        text = json.dumps(payload, **kwargs)
        if path is None:
            return text
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
            fh.write("\n")
        return None

    def to_csv(self, path: Union[str, Path, None] = None) -> Optional[str]:
        """Write per-system results to CSV."""
        import csv

        records = self.to_records()
        if not records:
            return "" if path is None else None

        fieldnames = list(records[0].keys())
        if path is None:
            buf = _io.StringIO()
            writer = csv.DictWriter(buf, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
            return buf.getvalue()

        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)
        return None

    def to_dataframe(self):
        """Convert to pandas DataFrame (one row per system)."""
        try:
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "pandas is required for to_dataframe(); use to_records() instead."
            ) from exc
        return pd.DataFrame(self.to_records())

    @classmethod
    def from_json(cls, source: Union[str, Path]) -> "BenchmarkResult":
        """Load from JSON produced by :meth:`to_json`."""
        source_path = Path(source)
        if source_path.is_file():
            text = source_path.read_text(encoding="utf-8")
        else:
            text = str(source)

        payload = json.loads(text)
        systems = []
        for rec in payload.get("systems", []):
            sys_obj = BenchmarkSystem(
                system_id=rec["system_id"],
                protein_pdb_path=Path("."),
                protein_sequence="",
                ligand_mol2_path=Path("."),
                ligand_smiles="",
                reference_pose_pdb_path=Path("."),
                experimental_dg_kcal_mol=rec.get("experimental_dg_kcal_mol"),
                is_active=rec.get("is_active"),
            )
            fa_result = None
            if rec.get("flexaidds_rmsd") is not None or rec.get("flexaidds_dg") is not None:
                fa_result = MethodResult(
                    method="flexaidds",
                    system_id=rec["system_id"],
                    best_pose_rmsd_angstrom=rec.get("flexaidds_rmsd"),
                    predicted_dg_kcal_mol=rec.get("flexaidds_dg"),
                    predicted_score=rec.get("flexaidds_score"),
                    n_poses=rec.get("flexaidds_n_poses", 0),
                    wall_time_seconds=rec.get("flexaidds_time", 0.0),
                )
            b2_result = None
            if rec.get("boltz2_rmsd") is not None or rec.get("boltz2_dg") is not None:
                b2_result = MethodResult(
                    method="boltz2",
                    system_id=rec["system_id"],
                    best_pose_rmsd_angstrom=rec.get("boltz2_rmsd"),
                    predicted_dg_kcal_mol=rec.get("boltz2_dg"),
                    predicted_score=rec.get("boltz2_score"),
                    n_poses=rec.get("boltz2_n_poses", 0),
                    wall_time_seconds=rec.get("boltz2_time", 0.0),
                )
            systems.append(SystemBenchmarkResult(
                system=sys_obj,
                flexaidds_result=fa_result,
                boltz2_result=b2_result,
            ))

        return cls(
            systems=tuple(systems),
            rmsd_threshold_angstrom=payload.get("rmsd_threshold_angstrom", 2.0),
            metadata=payload.get("metadata", {}),
        )


# ---------------------------------------------------------------------------
# Dataset I/O
# ---------------------------------------------------------------------------


def load_benchmark_dataset(path: Path) -> List[BenchmarkSystem]:
    """Load a benchmark dataset from a JSON file.

    Paths in the JSON are resolved relative to the file's parent directory.
    """
    path = Path(path)
    base_dir = path.parent
    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)
    return [BenchmarkSystem.from_dict(s, base_dir) for s in data["systems"]]


def save_benchmark_dataset(systems: List[BenchmarkSystem], path: Path) -> None:
    """Save a benchmark dataset to JSON."""
    payload = {"systems": [s.to_dict() for s in systems]}
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2)
        fh.write("\n")


# ---------------------------------------------------------------------------
# Runner functions
# ---------------------------------------------------------------------------


def run_flexaidds(
    system: BenchmarkSystem,
    *,
    binary: Optional[str] = None,
    timeout: int = 3600,
    temperature: float = 300.0,
) -> MethodResult:
    """Run FlexAIDdS on one system and extract benchmark metrics.

    Uses :func:`flexaidds.dock` to invoke the C++ engine, then extracts
    the top binding mode's free energy and best pose RMSD vs reference.
    """
    from . import dock

    t0 = time.monotonic()
    population = dock(
        receptor=str(system.protein_pdb_path),
        ligand=str(system.ligand_mol2_path),
        temperature=temperature,
        binary=binary,
        timeout=timeout,
    )
    elapsed = time.monotonic() - t0

    # Extract top mode
    ranked = population.rank_by_free_energy()
    if not ranked:
        return MethodResult(
            method="flexaidds", system_id=system.system_id,
            wall_time_seconds=elapsed,
        )

    top = ranked[0]
    best_cf = min((p.energy for p in top._poses), default=None)

    # RMSD: use the best-energy pose coordinates vs reference
    rmsd_val = None
    try:
        ref_coords = extract_ligand_coords_from_pdb(system.reference_pose_pdb_path)
        # Find best pose RMSD among top mode poses
        best_rmsd = float("inf")
        for pose in top._poses:
            if pose.coordinates is not None:
                pred_coords = pose.coordinates
                if pred_coords.shape == ref_coords.shape:
                    r = compute_rmsd(pred_coords, ref_coords)
                    if r < best_rmsd:
                        best_rmsd = r
        if best_rmsd < float("inf"):
            rmsd_val = best_rmsd
    except (ValueError, FileNotFoundError, OSError):
        pass

    return MethodResult(
        method="flexaidds",
        system_id=system.system_id,
        best_pose_rmsd_angstrom=rmsd_val,
        predicted_dg_kcal_mol=top.free_energy,
        predicted_score=best_cf,
        n_poses=top.n_poses,
        wall_time_seconds=elapsed,
    )


def run_boltz2(
    system: BenchmarkSystem,
    *,
    client: Optional[Any] = None,
    predict_affinity: bool = True,
    diffusion_samples: int = 5,
) -> MethodResult:
    """Run Boltz-2 on one system and extract benchmark metrics.

    Uses :class:`Boltz2Client` to call the NIM API, parses the mmCIF
    output for RMSD, and extracts affinity predictions.
    """
    from .boltz2 import Boltz2Client

    if client is None:
        client = Boltz2Client()

    t0 = time.monotonic()
    result = client.predict_protein_ligand(
        protein_sequence=system.protein_sequence,
        ligand_smiles=system.ligand_smiles,
        predict_affinity=predict_affinity,
        pocket_residues=list(system.pocket_residues) if system.pocket_residues else None,
        diffusion_samples=diffusion_samples,
    )
    elapsed = time.monotonic() - t0

    # Extract RMSD from best structure sample
    rmsd_val = None
    try:
        ref_coords = extract_ligand_coords_from_pdb(system.reference_pose_pdb_path)
        best_rmsd = float("inf")
        for struct in result.structures:
            try:
                pred_coords = extract_ligand_coords_from_mmcif(struct)
                if pred_coords.shape == ref_coords.shape:
                    r = compute_rmsd(pred_coords, ref_coords)
                    if r < best_rmsd:
                        best_rmsd = r
            except (ValueError, IndexError):
                continue
        if best_rmsd < float("inf"):
            rmsd_val = best_rmsd
    except (ValueError, FileNotFoundError, OSError):
        pass

    # Extract affinity
    predicted_dg = None
    predicted_score = None
    for aff in result.affinities.values():
        if aff.pic50:
            mean_pic50 = sum(aff.pic50) / len(aff.pic50)
            predicted_score = mean_pic50
            predicted_dg = pic50_to_dg(mean_pic50)
            break

    return MethodResult(
        method="boltz2",
        system_id=system.system_id,
        best_pose_rmsd_angstrom=rmsd_val,
        predicted_dg_kcal_mol=predicted_dg,
        predicted_score=predicted_score,
        n_poses=len(result.structures),
        wall_time_seconds=elapsed,
    )


# ---------------------------------------------------------------------------
# Cloud RAID storage helpers
# ---------------------------------------------------------------------------


def _detect_cloud_drives() -> List[Path]:
    """Detect available iCloud and Google Drive paths on macOS.

    Returns a list of writable cloud drive directories (up to 2: iCloud, Google
    Drive).  Falls back to an empty list on non-macOS or when drives are not
    mounted.
    """
    import glob

    candidates: List[Path] = []

    # iCloud Drive
    icloud = Path.home() / "Library" / "Mobile Documents" / "com~apple~CloudDocs"
    if icloud.is_dir():
        candidates.append(icloud)

    # Google Drive (any account)
    gd_pattern = str(Path.home() / "Library" / "CloudStorage" / "GoogleDrive-*" / "My Drive")
    for p in sorted(glob.glob(gd_pattern)):
        gp = Path(p)
        if gp.is_dir():
            candidates.append(gp)
            break  # use first account

    return candidates


def _mirror_write(data: str, filename: str, paths: List[Path]) -> None:
    """Write *data* to *filename* in every directory in *paths* (RAID-1 mirror).

    Uses atomic temp-file + rename for crash safety.
    """
    for base in paths:
        base.mkdir(parents=True, exist_ok=True)
        target = base / filename
        fd, tmp = tempfile.mkstemp(dir=str(base), suffix=".tmp")
        try:
            os.write(fd, data.encode("utf-8"))
            os.close(fd)
            os.replace(tmp, str(target))
        except OSError:
            try:
                os.unlink(tmp)
            except OSError:
                pass


# ---------------------------------------------------------------------------
# Checkpoint helpers
# ---------------------------------------------------------------------------


def _save_checkpoint(
    checkpoint_path: Path,
    results_by_idx: Dict[int, SystemBenchmarkResult],
    systems: List[BenchmarkSystem],
    mirror_dirs: Optional[List[Path]] = None,
) -> None:
    """Atomically write partial benchmark results for crash recovery.

    If *mirror_dirs* is provided the checkpoint is also mirrored to those
    directories (iCloud + Google Drive RAID-1).
    """
    completed = [results_by_idx[i] for i in sorted(results_by_idx)]
    partial = BenchmarkResult(systems=tuple(completed))
    text = partial.to_json() or ""

    # Write primary checkpoint atomically
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(
        dir=str(checkpoint_path.parent), suffix=".tmp",
    )
    try:
        os.write(fd, text.encode("utf-8"))
        os.close(fd)
        os.replace(tmp, str(checkpoint_path))
    except OSError:
        try:
            os.unlink(tmp)
        except OSError:
            pass

    # Mirror to cloud drives
    if mirror_dirs:
        _mirror_write(text, checkpoint_path.name, mirror_dirs)


def _load_checkpoint(
    checkpoint_path: Path,
    systems: List[BenchmarkSystem],
) -> Dict[int, SystemBenchmarkResult]:
    """Load previously completed systems from a checkpoint file.

    Returns a mapping of system index → result for systems already finished.
    """
    if not checkpoint_path.is_file():
        return {}

    try:
        prev = BenchmarkResult.from_json(checkpoint_path)
    except (json.JSONDecodeError, KeyError, ValueError):
        return {}

    done_ids = {sr.system.system_id for sr in prev.systems}
    idx_map: Dict[int, SystemBenchmarkResult] = {}

    for idx, sys in enumerate(systems):
        if sys.system_id in done_ids:
            for sr in prev.systems:
                if sr.system.system_id == sys.system_id:
                    idx_map[idx] = SystemBenchmarkResult(
                        system=sys,
                        flexaidds_result=sr.flexaidds_result,
                        boltz2_result=sr.boltz2_result,
                    )
                    break
    return idx_map


# ---------------------------------------------------------------------------
# Parallel worker function (must be top-level for pickling)
# ---------------------------------------------------------------------------


def _run_single_system(
    system: BenchmarkSystem,
    methods: Tuple[str, ...],
    flexaidds_binary: Optional[str],
    timeout_per_system: int,
    boltz2_predict_affinity: bool,
) -> SystemBenchmarkResult:
    """Execute all requested methods for one benchmark system.

    This function is the unit of work submitted to ``ProcessPoolExecutor``.
    It must remain a module-level function so that it is picklable.
    """
    fa_result = None
    b2_result = None

    if "flexaidds" in methods:
        try:
            fa_result = run_flexaidds(
                system, binary=flexaidds_binary, timeout=timeout_per_system,
            )
        except Exception:
            pass

    if "boltz2" in methods:
        try:
            b2_result = run_boltz2(
                system, predict_affinity=boltz2_predict_affinity,
            )
        except Exception:
            pass

    return SystemBenchmarkResult(
        system=system,
        flexaidds_result=fa_result,
        boltz2_result=b2_result,
    )


# ---------------------------------------------------------------------------
# Main benchmark runner
# ---------------------------------------------------------------------------


def run_benchmark(
    systems: List[BenchmarkSystem],
    *,
    methods: Tuple[str, ...] = ("flexaidds", "boltz2"),
    flexaidds_binary: Optional[str] = None,
    boltz2_client: Optional[Any] = None,
    boltz2_predict_affinity: bool = True,
    timeout_per_system: int = 3600,
    on_error: str = "skip",
    progress_callback: Optional[Callable[[str, int, int], None]] = None,
    max_workers: Optional[int] = None,
    checkpoint_path: Optional[Union[str, Path]] = None,
    use_cloud_raid: bool = False,
) -> BenchmarkResult:
    """Run the full benchmark across all systems.

    Args:
        systems: List of benchmark systems.
        methods: Which methods to run (``"flexaidds"``, ``"boltz2"``, or both).
        flexaidds_binary: Path to FlexAID executable.
        boltz2_client: Pre-configured Boltz2Client instance (ignored when
            *max_workers* > 1; each worker creates its own client).
        boltz2_predict_affinity: Request affinity from Boltz-2.
        timeout_per_system: Wall-clock timeout per system for FlexAIDdS.
        on_error: ``"skip"`` to continue on failure, ``"raise"`` to propagate.
        progress_callback: Called with ``(system_id, completed_count, total)``.
        max_workers: Parallelism level.  ``None`` or ``1`` for sequential
            execution (backward-compatible default).  ``0`` to auto-detect
            via :func:`auto_workers`.  ``>=2`` for that exact worker count.
        checkpoint_path: If set, partial results are saved here after each
            system completes so that a crashed run can be resumed.
        use_cloud_raid: If ``True``, checkpoint files are mirrored to iCloud
            Drive and Google Drive (RAID-1) to avoid local SSD usage.

    Returns:
        :class:`BenchmarkResult` with all per-system results.
    """
    total = len(systems)

    # Resolve checkpoint and cloud RAID paths
    ckpt: Optional[Path] = Path(checkpoint_path) if checkpoint_path else None
    mirror_dirs: Optional[List[Path]] = None
    if use_cloud_raid:
        cloud_drives = _detect_cloud_drives()
        if cloud_drives:
            mirror_dirs = [d / "FlexAIDdS" / "benchmarks" for d in cloud_drives]

    # Load any previous checkpoint
    results_by_idx: Dict[int, SystemBenchmarkResult] = {}
    if ckpt is not None:
        results_by_idx = _load_checkpoint(ckpt, systems)
        if results_by_idx:
            n_resumed = len(results_by_idx)
            if progress_callback:
                for idx in sorted(results_by_idx):
                    progress_callback(systems[idx].system_id, n_resumed, total)

    # Determine effective worker count
    effective_workers = 1
    if max_workers is not None:
        effective_workers = auto_workers() if max_workers == 0 else max_workers

    # --- Parallel path ---
    if effective_workers > 1:
        pending = {
            idx: sys for idx, sys in enumerate(systems)
            if idx not in results_by_idx
        }

        with ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                executor.submit(
                    _run_single_system,
                    sys,
                    methods,
                    flexaidds_binary,
                    timeout_per_system,
                    boltz2_predict_affinity,
                ): idx
                for idx, sys in pending.items()
            }

            completed_count = len(results_by_idx)
            for future in as_completed(futures):
                idx = futures[future]
                try:
                    result = future.result()
                except Exception:
                    if on_error == "raise":
                        raise
                    result = SystemBenchmarkResult(system=systems[idx])
                results_by_idx[idx] = result
                completed_count += 1
                if progress_callback:
                    progress_callback(
                        systems[idx].system_id, completed_count, total,
                    )
                if ckpt is not None:
                    _save_checkpoint(ckpt, results_by_idx, systems, mirror_dirs)

    # --- Sequential path (backward-compatible) ---
    else:
        for idx, system in enumerate(systems):
            if idx in results_by_idx:
                continue

            if progress_callback:
                progress_callback(system.system_id, len(results_by_idx), total)

            fa_result = None
            b2_result = None

            if "flexaidds" in methods:
                try:
                    fa_result = run_flexaidds(
                        system, binary=flexaidds_binary,
                        timeout=timeout_per_system,
                    )
                except Exception:
                    if on_error == "raise":
                        raise

            if "boltz2" in methods:
                try:
                    b2_result = run_boltz2(
                        system, client=boltz2_client,
                        predict_affinity=boltz2_predict_affinity,
                    )
                except Exception:
                    if on_error == "raise":
                        raise

            results_by_idx[idx] = SystemBenchmarkResult(
                system=system,
                flexaidds_result=fa_result,
                boltz2_result=b2_result,
            )
            if ckpt is not None:
                _save_checkpoint(ckpt, results_by_idx, systems, mirror_dirs)

    # Reassemble in original order
    ordered = [results_by_idx[i] for i in range(total) if i in results_by_idx]
    return BenchmarkResult(systems=tuple(ordered))
