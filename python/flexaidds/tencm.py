"""Torsional Elastic Network Model (TorsionalENM) for backbone flexibility.

Provides Pythonic wrappers around the C++ TorsionalENM engine, plus a
pure-Python fallback for environments without the compiled extension.

Reference:
    Delarue & Sanejouand (2002) J. Mol. Biol. 320:1011-24.
    Yang, Song & Cui (2009) Biophys. J. 97:2327-37.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

try:
    from . import _core
    _HAS_CORE = _core is not None
except ImportError:
    _core = None
    _HAS_CORE = False


@dataclass
class TorsionalNormalMode:
    """A single torsional normal mode.

    Attributes:
        eigenvalue:  Stiffness in kcal mol⁻¹ rad⁻².
        eigenvector: Displacement vector over torsion DOFs.
    """
    eigenvalue: float = 0.0
    eigenvector: List[float] = field(default_factory=list)

    def __repr__(self) -> str:
        return f"<TorsionalNormalMode λ={self.eigenvalue:.6g}>"


@dataclass
class Conformer:
    """A perturbed backbone conformation from TorsionalENM sampling.

    Attributes:
        delta_theta:    Torsion perturbations in radians.
        ca_positions:   Perturbed Cα positions as (N, 3) list.
        strain_energy:  ½ δθᵀ H δθ in kcal/mol.
    """
    delta_theta: List[float] = field(default_factory=list)
    ca_positions: List[List[float]] = field(default_factory=list)
    strain_energy: float = 0.0


@dataclass
class FullThermoResult:
    """Result from the ShannonThermoStack pipeline.

    Attributes:
        deltaG:                Base ΔG plus calibrated entropy terms (kcal/mol).
        shannonEntropy:        Shannon configurational entropy (nats).
        torsionalVibEntropy:   Torsional entropy; heuristic unless calibrated.
        entropyContribution:   Applied −T·S term; excludes heuristic S_vib.
        report:                Human-readable summary.
    """
    deltaG: float = 0.0
    shannonEntropy: float = 0.0
    torsionalVibEntropy: float = 0.0
    entropyContribution: float = 0.0
    report: str = ""

    def __repr__(self) -> str:
        return (
            f"<FullThermoResult ΔG={self.deltaG:.4f} "
            f"H_shannon={self.shannonEntropy:.4f} nats "
            f"S_vib_heuristic={self.torsionalVibEntropy:.6f} kcal/(mol·K)>"
        )


# ── Constants ────────────────────────────────────────────────────────────────
kB_kcal = 0.001987206  # kcal mol⁻¹ K⁻¹
kB_SI = 1.380649e-23   # J K⁻¹
hbar_SI = 1.054571817e-34  # J·s
DEFAULT_CUTOFF = 9.0    # Å
DEFAULT_K0 = 1.0        # kcal mol⁻¹ Å⁻²


def _read_ca_coords(pdb_path: str) -> np.ndarray:
    """Read Cα coordinates from a PDB file. Returns (N, 3) array."""
    coords = []
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        raise ValueError(f"No Cα atoms found in {pdb_path}")
    return np.array(coords, dtype=np.float64)


def _build_torsional_hessian(
    coords: np.ndarray,
    cutoff: float = DEFAULT_CUTOFF,
    k0: float = DEFAULT_K0,
) -> tuple:
    """Build torsional Hessian from Cα coordinates (pure-Python fallback).

    Returns (bonds, hessian, contacts) where bonds is a list of bond axes,
    hessian is the (n_bonds, n_bonds) matrix, and contacts is a list of
    (i, j, k_ij, r0) tuples.
    """
    n = len(coords)
    if n < 4:
        raise ValueError(f"Need at least 4 Cα atoms, got {n}")

    # Build pseudo-bonds (Cα_k → Cα_{k+1})
    bonds = []
    for k in range(n - 1):
        diff = coords[k + 1] - coords[k]
        length = np.linalg.norm(diff)
        if length < 1e-8:
            continue
        axis = diff / length
        pivot = (coords[k] + coords[k + 1]) / 2.0
        bonds.append((k, axis, pivot))

    n_bonds = len(bonds)
    if n_bonds == 0:
        return bonds, np.zeros((0, 0)), []

    # Build contacts
    contacts = []
    for i in range(n):
        for j in range(i + 1, n):
            r = np.linalg.norm(coords[j] - coords[i])
            if r <= cutoff:
                k_ij = k0 * (cutoff / r) ** 6
                contacts.append((i, j, k_ij, r))

    # Jacobian: ∂r_atom / ∂θ_bond for torsional DOF
    def jac(bond_idx: int, atom_idx: int) -> np.ndarray:
        bk, axis, pivot = bonds[bond_idx]
        # Atom is downstream of bond if atom_idx > bk+1
        if atom_idx <= bk:
            return np.zeros(3)
        r = coords[atom_idx] - pivot
        return np.cross(axis, r)

    # Assemble Hessian
    H = np.zeros((n_bonds, n_bonds))
    for ci, cj, k_ij, r0 in contacts:
        for bk in range(n_bonds):
            J_k_i = jac(bk, ci)
            J_k_j = jac(bk, cj)
            dJ_k = J_k_j - J_k_i
            for bl in range(bk, n_bonds):
                J_l_i = jac(bl, ci)
                J_l_j = jac(bl, cj)
                dJ_l = J_l_j - J_l_i
                val = k_ij * np.dot(dJ_k, dJ_l)
                H[bk, bl] += val
                if bk != bl:
                    H[bl, bk] += val

    return bonds, H, contacts


class TorsionalENM:
    """Torsional elastic network model for backbone flexibility.

    Mirrors the C++ TorsionalENM API. When the C++ extension is available,
    delegates to it for performance; otherwise uses a pure-Python fallback.

    Example:
        >>> tenm = TorsionalENM()
        >>> tenm.build_from_pdb("receptor.pdb")
        >>> print(f"Modes: {tenm.n_modes}")
        >>> for mode in tenm.modes[:5]:
        ...     print(f"  λ={mode.eigenvalue:.6g}")
    """

    def __init__(self) -> None:
        self._built = False
        self._modes: List[TorsionalNormalMode] = []
        self._n_residues = 0
        self._n_bonds = 0
        self._cpp_engine = None

    def build_from_pdb(
        self,
        pdb_path: str,
        cutoff: float = DEFAULT_CUTOFF,
        k0: float = DEFAULT_K0,
    ) -> None:
        """Build the torsional ENM from a PDB file.

        Args:
            pdb_path: Path to a PDB file with Cα atoms.
            cutoff:   Contact cutoff in Angstroms (default 9.0).
            k0:       Base spring constant (default 1.0).
        """
        if _HAS_CORE:
            try:
                self._cpp_engine = _core.TorsionalENM()
                self._cpp_engine.build_from_pdb(pdb_path, cutoff, k0)
                if self._cpp_engine.is_built:
                    self._built = True
                    self._n_residues = self._cpp_engine.n_residues
                    self._n_bonds = self._cpp_engine.n_bonds
                    self._modes = [
                        TorsionalNormalMode(
                            eigenvalue=m.eigenvalue,
                            eigenvector=list(m.eigenvector),
                        )
                        for m in self._cpp_engine.modes
                    ]
                    return
            except (AttributeError, RuntimeError):
                self._cpp_engine = None

        # Pure-Python fallback
        coords = _read_ca_coords(pdb_path)
        self._n_residues = len(coords)
        bonds, H, contacts = _build_torsional_hessian(coords, cutoff, k0)
        self._n_bonds = len(bonds)

        if self._n_bonds == 0:
            self._built = False
            return

        eigenvalues, eigenvectors = np.linalg.eigh(H)

        # Keep non-trivial modes (positive eigenvalues)
        self._modes = []
        for idx in range(len(eigenvalues)):
            lam = float(eigenvalues[idx])
            if lam > 1e-8:
                self._modes.append(TorsionalNormalMode(
                    eigenvalue=lam,
                    eigenvector=eigenvectors[:, idx].tolist(),
                ))
        self._built = True

    @property
    def is_built(self) -> bool:
        return self._built

    @property
    def n_residues(self) -> int:
        return self._n_residues

    @property
    def n_bonds(self) -> int:
        return self._n_bonds

    @property
    def modes(self) -> List[TorsionalNormalMode]:
        return self._modes

    @property
    def n_modes(self) -> int:
        return len(self._modes)

    def __repr__(self) -> str:
        status = "built" if self._built else "not built"
        return f"<TorsionalENM {status} residues={self._n_residues} modes={self.n_modes}>"


def compute_shannon_entropy(
    values: List[float],
    num_bins: int = 20,
) -> float:
    """Compute Shannon entropy of a continuous distribution.

    Bins the values into a histogram and computes H = -Σ p_i ln(p_i).

    Args:
        values:    List of continuous values.
        num_bins:  Number of histogram bins (default 20).

    Returns:
        Shannon entropy in nats (natural log).
    """
    if _HAS_CORE:
        try:
            return _core.compute_shannon_entropy(values, num_bins)
        except (AttributeError, RuntimeError):
            pass

    if not values:
        return 0.0

    arr = np.array(values, dtype=np.float64)
    counts, _ = np.histogram(arr, bins=num_bins)
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    return -float(np.sum(probs * np.log(probs)))


def compute_torsional_vibrational_entropy(
    modes: List[TorsionalNormalMode],
    temperature_K: float = 298.15,
    eigenvalue_to_omega: float = 1.0,
) -> float:
    """Compute torsional vibrational entropy from normal modes.

    Sums classical harmonic oscillator entropy for each torsional mode:
        S = k_B × [1 + ln(k_B T / (ħ ω))]
        ω = eigenvalue_to_omega × sqrt(λ)

    The default eigenvalue_to_omega=1.0 is a model-scale heuristic, not an
    absolute frequency calibration.

    Args:
        modes:          List of TorsionalNormalMode objects.
        temperature_K:  Temperature in Kelvin.
        eigenvalue_to_omega: rad/s per sqrt(model eigenvalue).

    Returns:
        Torsional vibrational entropy in kcal/(mol·K); model-scale heuristic
        unless eigenvalue_to_omega is physically calibrated.
    """
    if not math.isfinite(eigenvalue_to_omega) or eigenvalue_to_omega <= 0.0:
        raise ValueError("eigenvalue_to_omega must be finite and positive")

    if _HAS_CORE and eigenvalue_to_omega == 1.0:
        try:
            cpp_modes = []
            for m in modes:
                cm = _core.TorsionalNormalMode()
                cm.eigenvalue = m.eigenvalue
                cm.eigenvector = m.eigenvector
                cpp_modes.append(cm)
            return _core.compute_torsional_vibrational_entropy(cpp_modes, temperature_K)
        except (AttributeError, RuntimeError):
            pass

    # Pure-Python fallback
    if not modes:
        return 0.0

    S_total = 0.0
    for mode in modes[6:]:
        if mode.eigenvalue <= 1e-8:
            continue
        omega = eigenvalue_to_omega * math.sqrt(mode.eigenvalue)
        ratio = (kB_SI * temperature_K) / (hbar_SI * omega)
        if ratio > 0.0 and math.isfinite(ratio):
            S_total += kB_kcal * (1.0 + math.log(ratio))

    return S_total


def run_shannon_thermo_stack(
    energies: List[float],
    tencm_model: Optional[TorsionalENM] = None,
    base_deltaG: float = 0.0,
    temperature_K: float = 298.15,
    use_super_cluster: bool = False,
) -> FullThermoResult:
    """Run the full ShannonThermoStack pipeline.

    Combines Shannon configurational entropy from the GA ensemble with
    torsional vibrational entropy from the ENCoM backbone model.

    Formula:
        S_conf = k_B * H_nats           (nats → physical units)
        S_total = S_conf + S_vib         (additive for independent DOFs)
        ΔG = base_ΔG - T * S_total

    Args:
        energies:          List of pose energies from the GA ensemble.
        tencm_model:       Built TorsionalENM (optional; if None, vibrational
                           contribution is zero).
        base_deltaG:       Base enthalpy-dominated ΔG from scoring (kcal/mol).
        temperature_K:     Simulation temperature in Kelvin.
        use_super_cluster: When True, pre-filter energies through super-cluster
                           extraction for ~40%% faster Shannon collapse.

    Returns:
        FullThermoResult with combined thermodynamic quantities.
    """
    # Optional super-cluster pre-filtering
    sc_info = ""
    if use_super_cluster and len(energies) > 4:
        from .supercluster import SuperCluster
        sc = SuperCluster(energies)
        filtered = sc.filter_energies()
        sc_info = f"  SuperCluster pre-filter   = {sc.n_selected}/{sc.n_total} poses\n"
        energies = filtered

    # Shannon entropy
    H_shannon = compute_shannon_entropy(energies) if energies else 0.0
    S_conf_phys = H_shannon * kB_kcal

    # Torsional vibrational entropy
    S_vib = 0.0
    if tencm_model is not None and tencm_model.is_built:
        S_vib = compute_torsional_vibrational_entropy(
            tencm_model.modes, temperature_K)

    # The torsional ENCoM path is model-scale only unless a calibrated frequency
    # path is supplied. Keep S_vib as a relative diagnostic, but exclude it from
    # kcal/mol free-energy arithmetic in this stack.
    total_S = S_conf_phys
    entropy_contribution = -temperature_K * total_S
    deltaG = base_deltaG + entropy_contribution

    report = (
        f"ShannonThermoStack (T={temperature_K:.1f} K)\n"
        f"{sc_info}"
        f"  Shannon conf entropy    = {H_shannon:.4f} nats\n"
        f"  Torsional vib entropy   = {S_vib:.6f} kcal/(mol·K) "
        f"(model-scale heuristic; excluded from dG)\n"
        f"  Entropy contribution    = {entropy_contribution:.4f} kcal/mol (-TΔS)\n"
        f"  Total ΔG (F + applied entropy) = {deltaG:.4f} kcal/mol\n"
    )

    return FullThermoResult(
        deltaG=deltaG,
        shannonEntropy=H_shannon,
        torsionalVibEntropy=S_vib,
        entropyContribution=entropy_contribution,
        report=report,
    )
