"""Statistical mechanics and thermodynamic analysis for FlexAID∆S.

Provides Pythonic wrappers around C++ StatMechEngine with NumPy integration.
"""

import math
from typing import List, Optional, Tuple
from dataclasses import dataclass

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
        """Convert to dictionary for serialization."""
        return {
            'temperature_K': self.temperature,
            'log_Z': self.log_Z,
            'free_energy_kcal_mol': self.free_energy,
            'enthalpy_kcal_mol': self.mean_energy,
            'entropy_kcal_mol_K': self.entropy,
            'heat_capacity_kcal_mol_K2': self.heat_capacity,
            'std_energy_kcal_mol': self.std_energy,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Thermodynamics":
        """Construct a Thermodynamics instance from a dictionary.

        Accepts the key format produced by :meth:`to_dict` (suffixed keys such
        as ``temperature_K``, ``free_energy_kcal_mol``, …) as well as the raw
        attribute names (``temperature``, ``free_energy``, …).  Suffixed keys
        take priority when both forms are present.

        Args:
            data: Dictionary with thermodynamic quantities.

        Returns:
            A new :class:`Thermodynamics` instance.

        Raises:
            KeyError: If a required field is missing under both key forms.
        """
        def _get(suffixed: str, raw: str) -> float:
            if suffixed in data:
                return float(data[suffixed])
            if raw in data:
                return float(data[raw])
            raise KeyError(
                f"Missing required key: expected '{suffixed}' or '{raw}'"
            )

        return cls(
            temperature=_get("temperature_K", "temperature"),
            log_Z=_get("log_Z", "log_Z"),
            free_energy=_get("free_energy_kcal_mol", "free_energy"),
            mean_energy=_get("enthalpy_kcal_mol", "mean_energy"),
            mean_energy_sq=_get("mean_energy_sq", "mean_energy_sq"),
            heat_capacity=_get("heat_capacity_kcal_mol_K2", "heat_capacity"),
            entropy=_get("entropy_kcal_mol_K", "entropy"),
            std_energy=_get("std_energy_kcal_mol", "std_energy"),
        )


@dataclass
class EnergyComponents:
    """Per-microstate energy components in kcal/mol.

    ``complete`` means the listed terms partition the microstate total. When
    false, component means are diagnostics and ``component_sum_kcal_mol`` must
    not be interpreted as the ensemble ``H_eff``.
    """

    total: float = 0.0
    cf: float = 0.0
    receptor_strain: float = 0.0
    ligand_internal: float = 0.0
    hbond: float = 0.0
    gist: float = 0.0
    metal: float = 0.0
    water: float = 0.0
    other: float = 0.0
    complete: bool = False


@dataclass
class ComponentAverages:
    """Boltzmann-weighted energy component means in kcal/mol."""

    mean_CF_kcal_mol: float = 0.0
    mean_receptor_strain_kcal_mol: float = 0.0
    mean_ligand_internal_kcal_mol: float = 0.0
    mean_hbond_kcal_mol: float = 0.0
    mean_gist_kcal_mol: float = 0.0
    mean_metal_kcal_mol: float = 0.0
    mean_water_kcal_mol: float = 0.0
    mean_other_kcal_mol: float = 0.0
    component_sum_kcal_mol: float = 0.0
    component_completeness_flag: bool = False
    component_status: str = "not_computed"

    def to_dict(self) -> dict:
        return {
            "mean_CF_kcal_mol": self.mean_CF_kcal_mol,
            "mean_receptor_strain_kcal_mol": self.mean_receptor_strain_kcal_mol,
            "mean_ligand_internal_kcal_mol": self.mean_ligand_internal_kcal_mol,
            "mean_hbond_kcal_mol": self.mean_hbond_kcal_mol,
            "mean_gist_kcal_mol": self.mean_gist_kcal_mol,
            "mean_metal_kcal_mol": self.mean_metal_kcal_mol,
            "mean_water_kcal_mol": self.mean_water_kcal_mol,
            "mean_other_kcal_mol": self.mean_other_kcal_mol,
            "component_sum_kcal_mol": self.component_sum_kcal_mol,
            "component_completeness_flag": self.component_completeness_flag,
            "component_status": self.component_status,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ComponentAverages":
        return cls(
            mean_CF_kcal_mol=float(data.get("mean_CF_kcal_mol", 0.0)),
            mean_receptor_strain_kcal_mol=float(data.get("mean_receptor_strain_kcal_mol", 0.0)),
            mean_ligand_internal_kcal_mol=float(data.get("mean_ligand_internal_kcal_mol", 0.0)),
            mean_hbond_kcal_mol=float(data.get("mean_hbond_kcal_mol", 0.0)),
            mean_gist_kcal_mol=float(data.get("mean_gist_kcal_mol", 0.0)),
            mean_metal_kcal_mol=float(data.get("mean_metal_kcal_mol", 0.0)),
            mean_water_kcal_mol=float(data.get("mean_water_kcal_mol", 0.0)),
            mean_other_kcal_mol=float(data.get("mean_other_kcal_mol", 0.0)),
            component_sum_kcal_mol=float(data.get("component_sum_kcal_mol", 0.0)),
            component_completeness_flag=bool(data.get("component_completeness_flag", False)),
            component_status=str(data.get("component_status", "not_computed")),
        )


@dataclass
class ThermodynamicBreakdown:
    """Auditable thermodynamic ledger with explicit units and corrections.

    ``G_config_kcal_mol`` is the canonical configurational free energy from
    the sampled scoring-energy ensemble.  ``G_total_kcal_mol`` is the sum of
    that configurational term plus explicitly flagged correction terms.  These
    fields are not calibrated affinity estimates.
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
    components: Optional[ComponentAverages] = None
    has_components: bool = False

    def to_dict(self) -> dict:
        data = {
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
            "has_vib": self.has_vib,
            "has_natural": self.has_natural,
            "has_other": self.has_other,
            "has_components": self.has_components,
        }
        if self.components is not None:
            data["components"] = self.components.to_dict()
        return data

    @classmethod
    def from_dict(cls, data: dict) -> "ThermodynamicBreakdown":
        components_data = data.get("components")
        return cls(
            temperature_K=float(data.get("temperature_K", 300.0)),
            logZ_config=float(data.get("logZ_config", 0.0)),
            G_config_kcal_mol=float(data.get("G_config_kcal_mol", 0.0)),
            H_eff_kcal_mol=float(data.get("H_eff_kcal_mol", 0.0)),
            S_config_kcal_mol_K=float(data.get("S_config_kcal_mol_K", 0.0)),
            minus_T_S_config_kcal_mol=float(data.get("minus_T_S_config_kcal_mol", 0.0)),
            Cv_kcal_mol_K=float(data.get("Cv_kcal_mol_K", 0.0)),
            sigma_E_kcal_mol=float(data.get("sigma_E_kcal_mol", 0.0)),
            G_vib_kcal_mol=float(data.get("G_vib_kcal_mol", 0.0)),
            G_natural_kcal_mol=float(data.get("G_natural_kcal_mol", 0.0)),
            G_other_kcal_mol=float(data.get("G_other_kcal_mol", 0.0)),
            G_total_kcal_mol=float(data.get("G_total_kcal_mol", 0.0)),
            has_vib=bool(data.get("has_vib", False)),
            has_natural=bool(data.get("has_natural", False)),
            has_other=bool(data.get("has_other", False)),
            components=(
                ComponentAverages.from_dict(components_data)
                if isinstance(components_data, dict)
                else None
            ),
            has_components=bool(data.get("has_components", isinstance(components_data, dict))),
        )


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
