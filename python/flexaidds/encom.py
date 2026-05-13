"""ENCoM normal-mode analysis and vibrational entropy for FlexAID∆S.

Provides Pythonic wrappers around the C++ ENCoMEngine, plus a pure-Python
quasi-harmonic fallback for environments where the compiled extension is not
available.

Reference:
    Frappier et al. (2015). Proteins 83(11):2073-82.
    DOI: 10.1002/prot.24922
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

import numpy as np

try:
    from . import _core
    _HAS_CORE = _core is not None
except ImportError:
    _core = None
    _HAS_CORE = False

# Physical constants (used in pure-Python fallback)
kB_kcal  = 0.001987206      # kcal mol⁻¹ K⁻¹
kB_SI    = 1.380649e-23     # J K⁻¹
hbar_SI  = 1.054571817e-34  # J·s
NA       = 6.02214076e23    # mol⁻¹


@dataclass
class NormalMode:
    """A single normal mode from an elastic network calculation.

    Attributes:
        index:      1-based mode index.
        eigenvalue: λ_i in ENCoM arbitrary units.
        frequency:  sqrt(λ_i) in model scale unless calibrated to rad/s.
        eigenvector: Displacement vector with 3N components.
    """
    index: int = 0
    eigenvalue: float = 0.0
    frequency: float = 0.0
    eigenvector: List[float] = field(default_factory=list)

    def __repr__(self) -> str:
        return f"<NormalMode {self.index} λ={self.eigenvalue:.6g}>"


@dataclass
class FrequencyCalibration:
    """Eigenvalue-to-angular-frequency calibration metadata.

    ENCoM/tENCoM eigenvalues are model stiffness quantities. A calibrated
    scale is required before reporting absolute SI-frequency or absolute
    vibrational entropy claims.
    """
    eigenvalue_to_omega: float = 1.0
    calibrated: bool = False
    label: str = "model-scale"
    provenance: str = (
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied"
    )

    @classmethod
    def model_scale(cls) -> "FrequencyCalibration":
        return cls()

    @classmethod
    def calibrated_scale(
        cls,
        eigenvalue_to_omega: float,
        label: str,
        provenance: str = "User-supplied eigenvalue-to-frequency calibration",
    ) -> "FrequencyCalibration":
        if not math.isfinite(eigenvalue_to_omega) or eigenvalue_to_omega <= 0.0:
            raise ValueError("eigenvalue_to_omega must be finite and positive")
        return cls(
            eigenvalue_to_omega=eigenvalue_to_omega,
            calibrated=True,
            label=label or "calibrated",
            provenance=provenance,
        )

    @property
    def status(self) -> str:
        return "calibrated" if self.calibrated else "model_scale_heuristic"


@dataclass
class VibrationalEntropy:
    """Quasi-harmonic vibrational entropy from ENCoM normal modes.

    Attributes:
        S_vib_kcal_mol_K: Vibrational entropy in kcal mol⁻¹ K⁻¹.
        S_vib_J_mol_K:    Vibrational entropy in J mol⁻¹ K⁻¹.
        omega_eff:         Effective angular frequency ω_eff (rad/s).
        n_modes:           Number of non-trivial normal modes (3N − 6).
        temperature:       Temperature in Kelvin.
        eigenvalue_to_omega: Calibration scale used for omega_eff.
        calibrated:         True only when the scale has physical/empirical provenance.
    """
    S_vib_kcal_mol_K: float = 0.0
    S_vib_J_mol_K: float = 0.0
    omega_eff: float = 0.0
    n_modes: int = 0
    temperature: float = 300.0
    eigenvalue_to_omega: float = 1.0
    calibrated: bool = False
    calibration_label: str = "model-scale"
    calibration_provenance: str = (
        "No mass/inertia or empirical eigenvalue-to-frequency calibration supplied"
    )

    @property
    def free_energy_correction(self) -> float:
        """−T·S_vib vibrational free energy correction (kcal/mol)."""
        return -self.temperature * self.S_vib_kcal_mol_K

    @property
    def absolute_claim_allowed(self) -> bool:
        return self.calibrated

    @property
    def calibration_status(self) -> str:
        return "calibrated" if self.calibrated else "model_scale_heuristic"

    def __repr__(self) -> str:
        return (
            f"<VibrationalEntropy n_modes={self.n_modes} "
            f"S_vib={self.S_vib_kcal_mol_K:.6f} kcal/(mol·K) "
            f"T={self.temperature:.1f}K status={self.calibration_status}>"
        )


def _python_compute_vibrational_entropy(
    modes: List[NormalMode],
    temperature_K: float = 300.0,
    eigenvalue_cutoff: float = 1e-6,
    frequency_calibration: Optional[FrequencyCalibration] = None,
) -> VibrationalEntropy:
    """Pure-Python classical harmonic-oscillator S_vib approximation.

    Used when the C++ extension is not available.  Matches the C++ engine
    formula:
        ω_eff = scale × geometric_mean(sqrt(λ_i) for non-trivial modes)
        S_vib = n × k_B × [1 + ln(k_B T / (ħ ω_eff))]
    """
    calibration = frequency_calibration or FrequencyCalibration.model_scale()
    if not math.isfinite(calibration.eigenvalue_to_omega) or calibration.eigenvalue_to_omega <= 0.0:
        raise ValueError("eigenvalue_to_omega must be finite and positive")
    non_trivial = [m for m in modes if m.eigenvalue > eigenvalue_cutoff]
    if not non_trivial:
        return VibrationalEntropy(
            temperature=temperature_K,
            eigenvalue_to_omega=calibration.eigenvalue_to_omega,
            calibrated=calibration.calibrated,
            calibration_label=calibration.label,
            calibration_provenance=calibration.provenance,
        )

    log_sum = sum(0.5 * math.log(m.eigenvalue) for m in non_trivial)
    omega_eff = calibration.eigenvalue_to_omega * math.exp(log_sum / len(non_trivial))

    kT = kB_SI * temperature_K
    x = kT / (hbar_SI * omega_eff)
    if x <= 0.0 or not math.isfinite(x):
        return VibrationalEntropy(
            temperature=temperature_K,
            eigenvalue_to_omega=calibration.eigenvalue_to_omega,
            calibrated=calibration.calibrated,
            calibration_label=calibration.label,
            calibration_provenance=calibration.provenance,
        )

    n = len(non_trivial)
    S_vib_J = n * kB_SI * (1.0 + math.log(x))
    S_vib_J_mol = S_vib_J * NA
    S_vib_kcal_mol = S_vib_J_mol / 4184.0

    return VibrationalEntropy(
        S_vib_kcal_mol_K=S_vib_kcal_mol,
        S_vib_J_mol_K=S_vib_J_mol,
        omega_eff=omega_eff,
        n_modes=n,
        temperature=temperature_K,
        eigenvalue_to_omega=calibration.eigenvalue_to_omega,
        calibrated=calibration.calibrated,
        calibration_label=calibration.label,
        calibration_provenance=calibration.provenance,
    )


def _read_ca_coords(pdb_path: str) -> np.ndarray:
    """Read Cα coordinates from a PDB file.

    Returns:
        (N, 3) NumPy array of Cα positions in Angstroms.
    """
    coords = []
    with open(pdb_path) as fh:
        for line in fh:
            if (line.startswith("ATOM") and
                    line[12:16].strip() == "CA"):
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    if not coords:
        raise ValueError(f"No Cα atoms found in {pdb_path}")
    return np.array(coords, dtype=np.float64)


def _build_enm_modes(
    coords: np.ndarray,
    cutoff: float = 9.0,
) -> List[NormalMode]:
    """Build an elastic network model and return its normal modes.

    Uses the Anisotropic Network Model (ANM) with uniform spring constant.
    The Hessian is 3N×3N; eigenvalues below a trivial threshold (first 6
    modes = translations/rotations) are excluded.

    Args:
        coords:  (N, 3) array of Cα coordinates.
        cutoff:  Distance cutoff for spring contacts in Angstroms.

    Returns:
        List of NormalMode objects (excluding 6 trivial modes).
    """
    n = len(coords)
    hessian = np.zeros((3 * n, 3 * n), dtype=np.float64)

    for i in range(n):
        for j in range(i + 1, n):
            diff = coords[j] - coords[i]
            dist = np.linalg.norm(diff)
            if dist > cutoff:
                continue
            # Spring constant scaled by 1/dist² (ENCoM-like)
            k = 1.0 / (dist * dist)
            outer = np.outer(diff, diff) * (k / (dist * dist))
            ii, jj = 3 * i, 3 * j
            hessian[ii:ii+3, ii:ii+3] += outer
            hessian[jj:jj+3, jj:jj+3] += outer
            hessian[ii:ii+3, jj:jj+3] -= outer
            hessian[jj:jj+3, ii:ii+3] -= outer

    eigenvalues, eigenvectors = np.linalg.eigh(hessian)

    # Skip first 6 trivial modes (translation + rotation)
    modes: List[NormalMode] = []
    for idx in range(6, len(eigenvalues)):
        lam = float(eigenvalues[idx])
        if lam < 0:
            lam = 0.0
        modes.append(NormalMode(
            index=idx - 5,
            eigenvalue=lam,
            frequency=math.sqrt(lam),
            eigenvector=eigenvectors[:, idx].tolist(),
        ))
    return modes


class ENCoMEngine:
    """ENCoM quasi-harmonic entropy calculator.

    All methods are static; the class is a namespace mirroring the C++ API
    so that calling code is identical whether the compiled extension is
    available or not.

    Example:
        >>> modes = ENCoMEngine.load_modes("eigenvalues.txt", "eigenvectors.txt")
        >>> vs = ENCoMEngine.compute_vibrational_entropy(modes, temperature_K=300.0)
        >>> print(f"S_vib = {vs.S_vib_kcal_mol_K:.6f} kcal/(mol·K)")
        >>> print(f"F correction = {vs.free_energy_correction:.3f} kcal/mol")
    """

    @staticmethod
    def load_modes(
        eigenvalue_file: str,
        eigenvector_file: str,
    ) -> List[NormalMode]:
        """Load normal modes from ENCoM output files.

        Args:
            eigenvalue_file:  Plain text file with one eigenvalue per line.
            eigenvector_file: Plain text file with one mode per row,
                              space-separated displacement components.

        Returns:
            List of NormalMode objects sorted by mode index.
        """
        if _HAS_CORE:
            try:
                cpp_modes = _core.ENCoMEngine.load_modes(eigenvalue_file, eigenvector_file)
                return [
                    NormalMode(
                        index=m.index,
                        eigenvalue=m.eigenvalue,
                        frequency=m.frequency,
                        eigenvector=list(m.eigenvector),
                    )
                    for m in cpp_modes
                ]
            except (AttributeError, RuntimeError):
                # Fall through to pure-Python implementation if C++ unavailable
                pass

        # Pure-Python fallback
        eigenvalues: List[float] = []
        with open(eigenvalue_file) as fh:
            for line in fh:
                stripped = line.strip()
                if stripped and not stripped.startswith("#"):
                    eigenvalues.append(float(stripped.split()[0]))

        modes: List[NormalMode] = []
        with open(eigenvector_file) as fh:
            for i, (line, lam) in enumerate(
                zip(fh, eigenvalues), start=1
            ):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                evec = [float(x) for x in stripped.split()]
                freq = math.sqrt(max(lam, 0.0))
                modes.append(NormalMode(index=i, eigenvalue=lam, frequency=freq, eigenvector=evec))

        return modes

    @staticmethod
    def compute_vibrational_entropy(
        modes: List[NormalMode],
        temperature_K: float = 300.0,
        eigenvalue_cutoff: float = 1e-6,
        frequency_calibration: Optional[FrequencyCalibration] = None,
    ) -> VibrationalEntropy:
        """Schlitter quasi-harmonic S_vib from a list of NormalMode objects.

        Args:
            modes:             List of NormalMode objects.
            temperature_K:     Temperature in Kelvin (default 300 K).
            eigenvalue_cutoff: Modes with λ ≤ cutoff are treated as trivial
                               (translations/rotations) and excluded.

        Returns:
            VibrationalEntropy with S_vib in both kcal/(mol·K) and J/(mol·K).
            Values are model-scale heuristic unless frequency_calibration is
            calibrated.
        """
        calibration = frequency_calibration or FrequencyCalibration.model_scale()
        if _HAS_CORE:
            try:
                cpp_modes = []
                for m in modes:
                    cm = _core.NormalMode()
                    cm.index = m.index
                    cm.eigenvalue = m.eigenvalue
                    cm.frequency = m.frequency
                    cm.eigenvector = m.eigenvector
                    cpp_modes.append(cm)
                eng = _core.ENCoMEngine(
                    eigenvalue_cutoff=eigenvalue_cutoff,
                    eigenvalue_to_omega=calibration.eigenvalue_to_omega,
                    calibrated=calibration.calibrated,
                    calibration_label=calibration.label,
                    calibration_provenance=calibration.provenance,
                )
                vs_cpp = eng.compute_vibrational_entropy(
                    cpp_modes, temperature_K
                )
                return VibrationalEntropy(
                    S_vib_kcal_mol_K=vs_cpp.S_vib_kcal_mol_K,
                    S_vib_J_mol_K=vs_cpp.S_vib_J_mol_K,
                    omega_eff=vs_cpp.omega_eff,
                    n_modes=vs_cpp.n_modes,
                    temperature=vs_cpp.temperature,
                    eigenvalue_to_omega=getattr(vs_cpp, "eigenvalue_to_omega", calibration.eigenvalue_to_omega),
                    calibrated=getattr(vs_cpp, "calibrated", calibration.calibrated),
                    calibration_label=getattr(vs_cpp, "calibration_label", calibration.label),
                    calibration_provenance=getattr(vs_cpp, "calibration_provenance", calibration.provenance),
                )
            except (AttributeError, RuntimeError, TypeError):
                # Fall through to pure-Python implementation if C++ unavailable
                pass

        return _python_compute_vibrational_entropy(
            modes,
            temperature_K,
            eigenvalue_cutoff,
            calibration,
        )

    @staticmethod
    def total_entropy(
        S_conf_kcal_mol_K: float,
        S_vib_kcal_mol_K: float,
    ) -> float:
        """S_total = S_conf + S_vib  (kcal mol⁻¹ K⁻¹).

        Args:
            S_conf_kcal_mol_K: Configurational entropy from StatMechEngine.
            S_vib_kcal_mol_K:  Vibrational entropy from compute_vibrational_entropy.

        Returns:
            Total entropy in kcal mol⁻¹ K⁻¹.
        """
        if _HAS_CORE:
            try:
                return _core.ENCoMEngine.total_entropy(S_conf_kcal_mol_K, S_vib_kcal_mol_K)
            except (AttributeError, RuntimeError):
                pass
        return S_conf_kcal_mol_K + S_vib_kcal_mol_K

    @staticmethod
    def compute_delta_s(
        apo_pdb: str,
        holo_pdb: str,
        eigenvalue_cutoff: float = 1e-6,
        temperature_K: float = 300.0,
    ) -> float:
        """Compute ΔS_vib between apo and holo receptor conformations.

        Builds coarse-grained elastic network models for both structures,
        computes normal modes via eigenvalue decomposition, then returns
        the difference in quasi-harmonic vibrational entropy:
            ΔS_vib = S_vib(holo) − S_vib(apo)

        Args:
            apo_pdb:            Path to apo (unbound) receptor PDB.
            holo_pdb:           Path to holo (ligand-bound) receptor PDB.
            eigenvalue_cutoff:  Modes with λ ≤ cutoff treated as trivial.
            temperature_K:      Temperature in Kelvin (default 300 K).

        Returns:
            ΔS_vib in kcal mol⁻¹ K⁻¹ (positive = holo more flexible).
        """
        coords_apo = _read_ca_coords(apo_pdb)
        coords_holo = _read_ca_coords(holo_pdb)

        modes_apo = _build_enm_modes(coords_apo)
        modes_holo = _build_enm_modes(coords_holo)

        vs_apo = ENCoMEngine.compute_vibrational_entropy(
            modes_apo, temperature_K, eigenvalue_cutoff)
        vs_holo = ENCoMEngine.compute_vibrational_entropy(
            modes_holo, temperature_K, eigenvalue_cutoff)

        return vs_holo.S_vib_kcal_mol_K - vs_apo.S_vib_kcal_mol_K

    @staticmethod
    def free_energy_with_vibrations(
        F_electronic: float,
        S_vib_kcal_mol_K: float,
        temperature_K: float,
    ) -> float:
        """F_total = F_elec − T·S_vib  (kcal/mol).

        Args:
            F_electronic:      Electronic free energy from StatMechEngine (kcal/mol).
            S_vib_kcal_mol_K:  Vibrational entropy (kcal mol⁻¹ K⁻¹).
            temperature_K:     Temperature in Kelvin.

        Returns:
            Total free energy including vibrational correction (kcal/mol).
        """
        if _HAS_CORE:
            try:
                return _core.ENCoMEngine.free_energy_with_vibrations(
                    F_electronic, S_vib_kcal_mol_K, temperature_K
                )
            except (AttributeError, RuntimeError):
                pass
        return F_electronic - temperature_K * S_vib_kcal_mol_K
