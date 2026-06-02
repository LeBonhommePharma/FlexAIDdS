"""Statistical mechanics and thermodynamic analysis for FlexAID∆S.

Provides Pythonic wrappers around C++ StatMechEngine with NumPy integration.
"""

from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple
from dataclasses import dataclass, field

try:
    import numpy as np
except ImportError:  # pragma: no cover
    np = None  # type: ignore[assignment]

try:
    from . import _core
except ImportError:
    _core = None

# Physical constants
if _core is not None:
    kB_kcal: float = _core.kB_kcal  # kcal mol⁻¹ K⁻¹
    kB_SI: float = _core.kB_SI      # J K⁻¹
else:
    kB_kcal = 0.001987206
    kB_SI = 1.380649e-23


@dataclass
class Thermodynamics:
    """Complete thermodynamic properties of a conformational ensemble.
    
    Attributes:
        temperature: Temperature in Kelvin
        log_Z: Natural logarithm of partition function
        free_energy: Helmholtz free energy F = -kT ln Z (kcal/mol)
        mean_energy: Boltzmann-weighted average energy ⟨E⟩ (kcal/mol)
        mean_energy_sq: ⟨E²⟩ for variance calculation
        heat_capacity: Cv = (⟨E²⟩ - ⟨E⟩²) / (kT²) (kcal mol⁻¹ K⁻²)
        entropy: Configurational entropy S = (⟨E⟩ - F) / T (kcal mol⁻¹ K⁻¹)
        std_energy: Standard deviation of energy σ_E (kcal/mol)
    """
    temperature: float
    log_Z: float
    free_energy: float
    mean_energy: float
    mean_energy_sq: float
    heat_capacity: float
    entropy: float
    std_energy: float
    
    @property
    def binding_free_energy(self) -> float:
        """Alias for free_energy (common in docking context)."""
        return self.free_energy
    
    @property
    def entropy_term(self) -> float:
        """Entropic contribution to free energy: TΔS (kcal/mol)."""
        return self.temperature * self.entropy
    
    def __repr__(self) -> str:
        return (
            f"<Thermodynamics T={self.temperature:.1f}K "
            f"F={self.free_energy:.3f} S={self.entropy:.6f} "
            f"Cv={self.heat_capacity:.6f}>"
        )

    def to_dict(self) -> dict:
        """Convert to dictionary for serialization (includes all fields for roundtrip)."""
        return {
            'temperature_K': self.temperature,
            'log_Z': self.log_Z,
            'free_energy_kcal_mol': self.free_energy,
            'enthalpy_kcal_mol': self.mean_energy,
            'mean_energy_sq': self.mean_energy_sq,
            'entropy_kcal_mol_K': self.entropy,
            'heat_capacity_kcal_mol_K2': self.heat_capacity,
            'std_energy_kcal_mol': self.std_energy,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Thermodynamics":
        """Construct from either legacy raw keys or unit-suffixed JSON keys."""
        """Construct a Thermodynamics instance from a dictionary.

        Accepts the key format produced by :meth:`to_dict` (suffixed keys such
        as ``temperature_K``, ``free_energy_kcal_mol``, …) as well as the raw
        attribute names (``temperature``, ``free_energy``, …).  Suffixed keys
        take priority when both forms are present.
        """
        def _get(suffixed: str, raw: str) -> float:
            if suffixed in data:
                return float(data[suffixed])
            if raw in data:
                return float(data[raw])
            raise KeyError(
                f"Missing required key: expected '{suffixed}' or '{raw}'"
            )

        # mean_energy_sq may be absent in legacy dicts; default safe 0.0
        def _get_opt(suffixed: str, raw: str, default: float = 0.0) -> float:
            if suffixed in data:
                return float(data[suffixed])
            if raw in data:
                return float(data[raw])
            return default

        return cls(
            temperature=_get("temperature_K", "temperature"),
            log_Z=_get("log_Z", "log_Z"),
            free_energy=_get("free_energy_kcal_mol", "free_energy"),
            mean_energy=_get("enthalpy_kcal_mol", "mean_energy"),
            mean_energy_sq=_get_opt("mean_energy_sq", "mean_energy_sq"),
            heat_capacity=_get("heat_capacity_kcal_mol_K2", "heat_capacity"),
            entropy=_get("entropy_kcal_mol_K", "entropy"),
            std_energy=_get("std_energy_kcal_mol", "std_energy"),
        )


# ─── ThermodynamicBreakdown (Task 2 Python exposure + parity with C++) ──────
# Mirrors the C++ statmech::ThermodynamicBreakdown exactly for round-trip and
# C++/Python parity tests. All field names and units match the JSON contract
# in the roadmap. This is the pure-Python path; when C++ _core is available
# the C++ version (once bound) will be preferred for compute, but this
# dataclass remains the canonical serialisation shape.
@dataclass
class ThermodynamicBreakdown:
    """Auditable thermodynamic ledger for a binding mode or ensemble.

    All quantities follow the invariants in docs/dev/thermo_invariants.md.
    G_total = G_config + G_vib + G_natural + G_other (always defined).

    This is the shape emitted under "thermodynamics" in JSON output (Task 2+).
    """
    temperature_K: float = 300.0

    logZ_config: float = 0.0
    G_config_kcal_mol: float = 0.0
    H_eff_kcal_mol: float = 0.0
    S_config_kcal_mol_K: float = 0.0
    minus_T_S_config_kcal_mol: float = 0.0
    Cv_kcal_mol_K: float = 0.0
    sigma_E_kcal_mol: float = 0.0

    G_vib_kcal_mol: float = 0.0
    G_natural_kcal_mol: float = 0.0
    G_other_kcal_mol: float = 0.0
    G_total_kcal_mol: float = 0.0

    has_vib: bool = False
    has_natural: bool = False
    has_other: bool = False

    # Task 3: component-wise Boltzmann averages (when available)
    component_means: Dict[str, float] = field(default_factory=dict)
    component_sum_kcal_mol: float = 0.0
    components_complete: bool = False

    # Task 6: Standard-state affinity calibration (safe / experimental)
    affinity: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        """Exact JSON shape required by roadmap (no legacy aliases here)."""
        return {
            "temperature_K": self.temperature_K,
            "logZ_config": self.logZ_config,
            "G_config_kcal_mol": self.G_config_kcal_mol,
            "H_eff_kcal_mol": self.H_eff_kcal_mol,
            "S_config_kcal_mol_K": self.S_config_kcal_mol_K,
            "minus_T_S_config_kcal_mol": self.minus_T_S_config_kcal_mol,
            "Cv_kcal_mol_K": self.Cv_kcal_mol_K,
            "sigma_E_kcal_mol": self.sigma_E_kcal_mol,
            "G_vib_kcal_mol": self.G_vib_kcal_mol,
            "G_natural_kcal_mol": self.G_natural_kcal_mol,
            "G_other_kcal_mol": self.G_other_kcal_mol,
            "G_total_kcal_mol": self.G_total_kcal_mol,
            "component_sum_kcal_mol": self.component_sum_kcal_mol,
            "components_complete": self.components_complete,
            "component_means": self.component_means,
            "affinity": self.affinity,
        }

    @classmethod
    def from_thermodynamics(cls, thermo: "Thermodynamics",
                            G_vib: float = 0.0, has_vib: bool = False,
                            G_natural: float = 0.0, has_natural: bool = False,
                            G_other: float = 0.0, has_other: bool = False) -> "ThermodynamicBreakdown":
        """Factory mirroring C++ make_breakdown() for pure-Python parity."""
        b = cls(
            temperature_K=thermo.temperature,
            logZ_config=thermo.log_Z,
            G_config_kcal_mol=thermo.free_energy,
            H_eff_kcal_mol=thermo.mean_energy,
            S_config_kcal_mol_K=thermo.entropy,
            minus_T_S_config_kcal_mol=thermo.free_energy - thermo.mean_energy,
            Cv_kcal_mol_K=thermo.heat_capacity,
            sigma_E_kcal_mol=thermo.std_energy,
            G_vib_kcal_mol=G_vib,
            G_natural_kcal_mol=G_natural,
            G_other_kcal_mol=G_other,
            G_total_kcal_mol=thermo.free_energy + G_vib + G_natural + G_other,
            has_vib=has_vib,
            has_natural=has_natural,
            has_other=has_other,
            component_means={},
            component_sum_kcal_mol=0.0,
            components_complete=False,
        )
        return b


# ─── Diagnostic-only enthalpy–entropy metrics (Task 4) ───────────────────────
# These are NEVER to be used for ranking, pose selection, or optimization.
# They are purely for analysis and must be labelled as diagnostics in all
# output and documentation. compensation_score == 1 when H and -TS perfectly
# cancel (G≈0); low when one term dominates.
EPS = 1e-12

def entropy_fraction(H_eff: float, minus_T_S: float) -> float:
    """| -TΔS | / (|H_eff| + |-TΔS| + eps)  — diagnostic only."""
    return abs(minus_T_S) / (abs(H_eff) + abs(minus_T_S) + EPS)

def enthalpy_fraction(H_eff: float, minus_T_S: float) -> float:
    """|H_eff| / (|H_eff| + |-TΔS| + eps) — diagnostic only."""
    return abs(H_eff) / (abs(H_eff) + abs(minus_T_S) + EPS)

def compensation_score(G_config: float, H_eff: float, minus_T_S: float) -> float:
    """1 - |G| / (|H| + |-T S| + eps) clamped to [0,1].

    High when enthalpy and entropy compensate (G small relative to parts).
    FORBIDDEN for ranking or affinity claims.
    """
    denom = abs(H_eff) + abs(minus_T_S) + EPS
    score = 1.0 - (abs(G_config) / denom)
    return max(0.0, min(1.0, score))  # clamp numerical noise


# ─── Task 6: Pure-Python affinity calibration (parity with C++) ──────────────
def deltaG_standard_to_Kd_M(deltaG_kcal_mol: float, T_K: float, c0_M: float = 1.0) -> float:
    """Safe conversion. Raises on invalid inputs."""
    if T_K <= 0:
        raise ValueError("Temperature must be > 0 K")
    if c0_M <= 0:
        raise ValueError("c0_M must be > 0")
    RT = kB_kcal * T_K
    return c0_M * math.exp(deltaG_kcal_mol / RT)

def Kd_M_to_deltaG_standard(Kd_M: float, T_K: float, c0_M: float = 1.0) -> float:
    """Safe conversion. Raises on invalid inputs."""
    if T_K <= 0:
        raise ValueError("Temperature must be > 0 K")
    if Kd_M <= 0:
        raise ValueError("Kd must be > 0 M")
    if c0_M <= 0:
        raise ValueError("c0_M must be > 0")
    RT = kB_kcal * T_K
    return RT * math.log(Kd_M / c0_M)


class _PyStatMechEngine:
    """Pure-Python canonical-ensemble engine (fallback when C++ _core is absent).

    Uses the log-sum-exp trick for numerical stability.
    """

    def __init__(self, temperature_K: float) -> None:
        self._T = float(temperature_K)
        if self._T <= 0.0:
            raise ValueError("StatMechEngine: temperature must be > 0")
        self._beta = 1.0 / (kB_kcal * self._T)
        self._states: List[Tuple[float, float]] = []

    # ------------------------------------------------------------------
    # sample accumulation
    # ------------------------------------------------------------------

    def add_sample(self, energy: float, multiplicity: float = 1.0) -> None:
        count = float(multiplicity)
        if count < 0.0:
            raise ValueError("StatMechEngine: multiplicity must be non-negative")
        self._states.append((float(energy), count))

    def clear(self) -> None:
        self._states.clear()

    @property
    def size(self) -> int:
        return len(self._states)

    @property
    def temperature(self) -> float:
        return self._T

    @property
    def beta(self) -> float:
        return self._beta

    # ------------------------------------------------------------------
    # thermodynamic computation
    # ------------------------------------------------------------------

    def compute(self) -> Thermodynamics:
        if not self._states:
            raise RuntimeError("No samples added to StatMechEngine before compute()")

        energies = [energy for energy, _count in self._states]
        log_w = [
            (math.log(count) - self._beta * energy) if count > 0.0 else -math.inf
            for energy, count in self._states
        ]
        max_log_w = max(log_w)
        if not math.isfinite(max_log_w):
            raise RuntimeError("StatMechEngine: ensemble has zero total multiplicity")

        log_Z = max_log_w + math.log(sum(math.exp(w - max_log_w) for w in log_w))

        # Boltzmann weights
        weights = [math.exp(w - log_Z) if math.isfinite(w) else 0.0 for w in log_w]

        mean_e = sum(wi * ei for wi, ei in zip(weights, energies))
        mean_e2 = sum(wi * ei * ei for wi, ei in zip(weights, energies))
        var_e = mean_e2 - mean_e ** 2
        std_e = math.sqrt(max(0.0, var_e))

        free_energy = -kB_kcal * self._T * log_Z
        heat_capacity = var_e / (kB_kcal * self._T ** 2)
        entropy = (mean_e - free_energy) / self._T

        return Thermodynamics(
            temperature=self._T,
            log_Z=log_Z,
            free_energy=free_energy,
            mean_energy=mean_e,
            mean_energy_sq=mean_e2,
            heat_capacity=heat_capacity,
            entropy=entropy,
            std_energy=std_e,
        )

    def compute_breakdown(
        self,
        G_vib_kcal_mol: float = 0.0,
        G_natural_kcal_mol: float = 0.0,
        G_other_kcal_mol: float = 0.0,
        has_vib: bool = False,
        has_natural: bool = False,
        has_other: bool = False,
    ) -> ThermodynamicBreakdown:
        thermo = self.compute()
        return ThermodynamicBreakdown(
            temperature_K=thermo.temperature,
            logZ_config=thermo.log_Z,
            G_config_kcal_mol=thermo.free_energy,
            H_eff_kcal_mol=thermo.mean_energy,
            S_config_kcal_mol_K=thermo.entropy,
            minus_T_S_config_kcal_mol=thermo.free_energy - thermo.mean_energy,
            Cv_kcal_mol_K=thermo.heat_capacity,
            sigma_E_kcal_mol=thermo.std_energy,
            G_vib_kcal_mol=G_vib_kcal_mol,
            G_natural_kcal_mol=G_natural_kcal_mol,
            G_other_kcal_mol=G_other_kcal_mol,
            G_total_kcal_mol=thermo.free_energy + G_vib_kcal_mol + G_natural_kcal_mol + G_other_kcal_mol,
            has_vib=has_vib,
            has_natural=has_natural,
            has_other=has_other,
        )

    def component_averages(self, components: List[EnergyComponents]) -> ComponentAverages:
        if len(components) != len(self._states):
            raise ValueError("component count must match ensemble size")
        if not components:
            raise ValueError("component list must not be empty")

        weights = self.boltzmann_weights()
        avg = ComponentAverages(component_completeness_flag=True)
        for weight, comp in zip(weights, components):
            avg.mean_CF_kcal_mol += weight * comp.cf
            avg.mean_receptor_strain_kcal_mol += weight * comp.receptor_strain
            avg.mean_ligand_internal_kcal_mol += weight * comp.ligand_internal
            avg.mean_hbond_kcal_mol += weight * comp.hbond
            avg.mean_gist_kcal_mol += weight * comp.gist
            avg.mean_metal_kcal_mol += weight * comp.metal
            avg.mean_water_kcal_mol += weight * comp.water
            avg.mean_other_kcal_mol += weight * comp.other
            avg.component_completeness_flag = avg.component_completeness_flag and comp.complete

        avg.component_sum_kcal_mol = (
            avg.mean_CF_kcal_mol
            + avg.mean_receptor_strain_kcal_mol
            + avg.mean_ligand_internal_kcal_mol
            + avg.mean_hbond_kcal_mol
            + avg.mean_gist_kcal_mol
            + avg.mean_metal_kcal_mol
            + avg.mean_water_kcal_mol
            + avg.mean_other_kcal_mol
        )
        avg.component_status = "available" if avg.component_completeness_flag else "included_in_other"
        return avg

    def boltzmann_weights(self) -> List[float]:
        if not self._states:
            return []
        thermo = self.compute()
        log_Z = thermo.log_Z
        return [
            count * math.exp(-self._beta * energy - log_Z) if count > 0.0 else 0.0
            for energy, count in self._states
        ]

    def delta_G(self, other: "_PyStatMechEngine") -> float:
        return self.compute().free_energy - other.compute().free_energy


class StatMechEngine:
    """Statistical mechanics engine for conformational ensembles.

    Computes partition functions, free energies, entropies, and heat capacities
    from sampled configurations using canonical ensemble formalism.

    Example:
        >>> engine = StatMechEngine(temperature_K=300.0)
        >>> engine.add_samples([-10.5, -9.8, -10.2, -11.0])  # energies in kcal/mol
        >>> thermo = engine.compute()
        >>> print(f"Free energy: {thermo.free_energy:.2f} kcal/mol")
        >>> print(f"Entropy: {thermo.entropy:.5f} kcal/(mol·K)")
    """

    def __init__(self, temperature_K: float = 300.0):
        """Initialize engine at specified temperature.

        Args:
            temperature_K: Simulation temperature in Kelvin (default 300K)
        """
        if _core is not None:
            self._engine = _core.StatMechEngine(temperature_K)
        else:
            self._engine = _PyStatMechEngine(temperature_K)
    
    def add_sample(self, energy: float, multiplicity: int = 1) -> None:
        """Add a single sampled configuration.
        
        Args:
            energy: Configuration energy in kcal/mol (negative = favorable)
            multiplicity: Degeneracy/sampling count (default 1)
        """
        self._engine.add_sample(energy, multiplicity)
    
    def add_samples(self, energies) -> None:
        """Add multiple configurations from a sequence or NumPy array.

        Args:
            energies: Iterable of configuration energies (kcal/mol)
        """
        if np is not None:
            energies = np.asarray(energies, dtype=np.float64)
        for e in energies:
            self._engine.add_sample(float(e))
    
    def compute(self) -> Thermodynamics:
        """Compute full thermodynamics from current ensemble.
        
        Returns:
            Thermodynamics object with F, S, H, Cv, etc.
        """
        thermo_cpp = self._engine.compute()
        return Thermodynamics(
            temperature=thermo_cpp.temperature,
            log_Z=thermo_cpp.log_Z,
            free_energy=thermo_cpp.free_energy,
            mean_energy=thermo_cpp.mean_energy,
            mean_energy_sq=thermo_cpp.mean_energy_sq,
            heat_capacity=thermo_cpp.heat_capacity,
            entropy=thermo_cpp.entropy,
            std_energy=thermo_cpp.std_energy,
        )

    def compute_breakdown(
        self,
        G_vib_kcal_mol: float = 0.0,
        G_natural_kcal_mol: float = 0.0,
        G_other_kcal_mol: float = 0.0,
        has_vib: bool = False,
        has_natural: bool = False,
        has_other: bool = False,
    ) -> ThermodynamicBreakdown:
        """Compute the explicit thermodynamic ledger for the ensemble."""
        if hasattr(self._engine, "compute_breakdown"):
            result = self._engine.compute_breakdown(
                G_vib_kcal_mol,
                G_natural_kcal_mol,
                G_other_kcal_mol,
                has_vib,
                has_natural,
                has_other,
            )
            return ThermodynamicBreakdown(
                temperature_K=result.temperature_K,
                logZ_config=result.logZ_config,
                G_config_kcal_mol=result.G_config_kcal_mol,
                H_eff_kcal_mol=result.H_eff_kcal_mol,
                S_config_kcal_mol_K=result.S_config_kcal_mol_K,
                minus_T_S_config_kcal_mol=result.minus_T_S_config_kcal_mol,
                Cv_kcal_mol_K=result.Cv_kcal_mol_K,
                sigma_E_kcal_mol=result.sigma_E_kcal_mol,
                G_vib_kcal_mol=result.G_vib_kcal_mol,
                G_natural_kcal_mol=result.G_natural_kcal_mol,
                G_other_kcal_mol=result.G_other_kcal_mol,
                G_total_kcal_mol=result.G_total_kcal_mol,
                has_vib=result.has_vib,
                has_natural=result.has_natural,
                has_other=result.has_other,
                components=(
                    _component_averages_from_cpp(result.components)
                    if getattr(result, "has_components", False)
                    else None
                ),
                has_components=getattr(result, "has_components", False),
            )

        thermo = self.compute()
        return ThermodynamicBreakdown(
            temperature_K=thermo.temperature,
            logZ_config=thermo.log_Z,
            G_config_kcal_mol=thermo.free_energy,
            H_eff_kcal_mol=thermo.mean_energy,
            S_config_kcal_mol_K=thermo.entropy,
            minus_T_S_config_kcal_mol=thermo.free_energy - thermo.mean_energy,
            Cv_kcal_mol_K=thermo.heat_capacity,
            sigma_E_kcal_mol=thermo.std_energy,
            G_vib_kcal_mol=G_vib_kcal_mol,
            G_natural_kcal_mol=G_natural_kcal_mol,
            G_other_kcal_mol=G_other_kcal_mol,
            G_total_kcal_mol=thermo.free_energy + G_vib_kcal_mol + G_natural_kcal_mol + G_other_kcal_mol,
            has_vib=has_vib,
            has_natural=has_natural,
            has_other=has_other,
        )

    def component_averages(self, components: List[EnergyComponents]) -> ComponentAverages:
        """Boltzmann-weight component diagnostics over the current ensemble."""
        if hasattr(self._engine, "component_averages"):
            cpp_components = []
            for comp in components:
                cpp_comp = _core.EnergyComponents() if _core is not None else None
                if cpp_comp is None:
                    break
                cpp_comp.total = comp.total
                cpp_comp.cf = comp.cf
                cpp_comp.receptor_strain = comp.receptor_strain
                cpp_comp.ligand_internal = comp.ligand_internal
                cpp_comp.hbond = comp.hbond
                cpp_comp.gist = comp.gist
                cpp_comp.metal = comp.metal
                cpp_comp.water = comp.water
                cpp_comp.other = comp.other
                cpp_comp.complete = comp.complete
                cpp_components.append(cpp_comp)
            if len(cpp_components) == len(components):
                return _component_averages_from_cpp(self._engine.component_averages(cpp_components))

        return self._engine.component_averages(components)
    
    def boltzmann_weights(self):
        """Get Boltzmann weights for all samples.

        Returns:
            NumPy array of normalized weights (sum to 1.0), or a plain list
            when NumPy is not available.
        """
        weights = self._engine.boltzmann_weights()
        if np is not None:
            return np.array(weights)
        return list(weights)
    
    def delta_G(self, reference: 'StatMechEngine') -> float:
        """Compute relative free energy to another ensemble.
        
        Args:
            reference: Reference StatMechEngine
        
        Returns:
            ΔG = F_this - F_reference (kcal/mol)
        """
        return self._engine.delta_G(reference._engine)
    
    def clear(self) -> None:
        """Remove all samples from ensemble."""
        self._engine.clear()
    
    @property
    def temperature(self) -> float:
        """Temperature in Kelvin."""
        return self._engine.temperature
    
    @property
    def beta(self) -> float:
        """Thermodynamic beta = 1/(kT) in (kcal/mol)⁻¹."""
        return self._engine.beta
    
    @property
    def n_samples(self) -> int:
        """Number of configurations in ensemble."""
        return self._engine.size
    
    def __len__(self) -> int:
        return self.n_samples
    
    def __repr__(self) -> str:
        return f"<StatMechEngine T={self.temperature:.1f}K n_samples={self.n_samples}>"

    # Task 7/9: temperature scan + ΔCp (C++ path when bindings present; pure fallback for fit)
    def temperature_scan(self, temperatures: list[float]) -> list[dict]:
        """Recompute G/H/S/Cv at multiple temperatures using fixed ensemble energies (C++ only).

        Returns list of dicts with keys T_K, logZ, G_kcal_mol, H_kcal_mol, S_kcal_mol_K, Cv_kcal_mol_K.
        All new thermodynamic terms carry explicit units per roadmap invariants.
        """
        if _core is None or not hasattr(self._engine, "temperature_scan"):
            raise NotImplementedError(
                "temperature_scan requires the C++ extension. "
                "Run `pip install -e .` with BUILD_PYTHON_BINDINGS=ON."
            )
        points = self._engine.temperature_scan(temperatures)
        return [
            {
                "T_K": p.T_K,
                "logZ": p.logZ,
                "G_kcal_mol": p.G_kcal_mol,
                "H_kcal_mol": p.H_kcal_mol,
                "S_kcal_mol_K": p.S_kcal_mol_K,
                "Cv_kcal_mol_K": p.Cv_kcal_mol_K,
            }
            for p in points
        ]

    @staticmethod
    def fit_delta_Cp(scan_points: list[dict], T_ref_K: float) -> dict:
        """Linear regression ΔCp fit (requires >=4 points). Pure-Python (no C++ needed).

        This is model-derived / experimental-diagnostic output (see roadmap Task 7).
        Always labelled experimental=true, model_derived=true. Never used for ranking.
        Units: delta_Cp_kcal_mol_K, rmse_kcal_mol, T_ref_K.
        """
        if len(scan_points) < 4:
            raise ValueError("ΔCp fit requires at least 4 temperature points")
        n = len(scan_points)
        sum_x = sum_y = sum_xx = sum_xy = 0.0
        h0 = scan_points[0].get("H_kcal_mol", 0.0)
        for p in scan_points:
            x = p["T_K"] - T_ref_K
            y = p.get("H_kcal_mol", 0.0) - h0
            sum_x += x
            sum_y += y
            sum_xx += x * x
            sum_xy += x * y
        denom = n * sum_xx - sum_x * sum_x
        delta_Cp = (n * sum_xy - sum_x * sum_y) / denom if abs(denom) > 1e-12 else 0.0
        # RMSE
        sse = 0.0
        for p in scan_points:
            x = p["T_K"] - T_ref_K
            pred = h0 + delta_Cp * x
            sse += (p.get("H_kcal_mol", 0.0) - pred) ** 2
        rmse = (sse / n) ** 0.5
        return {
            "delta_Cp_kcal_mol_K": delta_Cp,
            "T_ref_K": T_ref_K,
            "rmse_kcal_mol": rmse,
            "model_derived": True,
            "experimental": True,
        }


class BoltzmannLUT:
    """Pre-tabulated Boltzmann factors for fast inner-loop evaluation.
    
    Provides O(1) lookup for exp(-βE) over a specified energy range.
    
    Example:
        >>> lut = BoltzmannLUT(beta=1.688, e_min=-20.0, e_max=5.0, n_bins=10000)
        >>> weight = lut(-12.5)  # exp(-β × -12.5)
    """
    
    def __init__(self, beta: float, e_min: float, e_max: float, n_bins: int = 10000):
        """Initialize lookup table.
        
        Args:
            beta: 1/(kT) in (kcal/mol)⁻¹
            e_min: Minimum energy (kcal/mol)
            e_max: Maximum energy (kcal/mol)
            n_bins: Number of table entries (default 10000)
        """
        if _core is None:
            raise RuntimeError("C++ bindings not available")
        self._lut = _core.BoltzmannLUT(beta, e_min, e_max, n_bins)
    
    def __call__(self, energy: float) -> float:
        """Look up Boltzmann factor for given energy.
        
        Args:
            energy: Energy in kcal/mol
        
        Returns:
            exp(-βE)
        """
        return self._lut(energy)


def helmholtz_from_energies(energies, temperature: float = 300.0) -> float:
    """Convenience function: Helmholtz free energy from energy array.

    Works with or without C++ bindings (uses pure-Python engine as fallback).

    Args:
        energies: Iterable of configuration energies (kcal/mol)
        temperature: Temperature in Kelvin

    Returns:
        Helmholtz free energy F (kcal/mol)
    """
    engine = StatMechEngine(temperature)
    engine.add_samples(energies)
    return engine.compute().free_energy
