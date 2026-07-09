"""Grand-canonical (μVT) competitive-binding helpers for FlexAIDdS.

This module does **not** reimplement the C++ GrandPartitionFunction.
It provides:

1. A pure-Python `CompetitiveSite` mirror of the single-site Ξ math
   (for NRGsuite / notebooks without compiling bindings).
2. `set_concentration(...)` — the user-facing multi-ligand API.

When C++ `_core` eventually exposes GrandPartitionFunction, prefer that
path; this module stays as a fallback and documentation surface.

Thermodynamics (single site, ideal solution, c° = 1 M):

    Ξ = 1 + Σ_i (c_i/c°) Z_i
    p_i = z_i Z_i / Ξ
    ⟨N⟩ = 1 − p_empty

See docs/theory.md.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union

# Match C++ statmech::kB_kcal
kB_kcal: float = 0.001987206
C_STANDARD_M: float = 1.0


@dataclass
class LigandEntry:
    name: str
    log_Z: float
    concentration_M: float = 1.0

    @property
    def log_c(self) -> float:
        if self.concentration_M <= 0.0:
            raise ValueError("concentration_M must be > 0")
        return math.log(self.concentration_M / C_STANDARD_M)

    @property
    def log_zZ(self) -> float:
        return self.log_c + self.log_Z


@dataclass
class OccupancyPoint:
    concentration_M: float
    p_bound: float
    p_species: float
    mean_N: float


@dataclass
class CompetitiveSite:
    """Pure-Python single-site grand partition function (μVT).

    Mirrors `target::GrandPartitionFunction` semantics for offline analysis.
    """

    temperature_K: float = 300.0
    ligands: Dict[str, LigandEntry] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.temperature_K <= 0.0:
            raise ValueError("temperature_K must be > 0")

    @property
    def kT(self) -> float:
        return kB_kcal * self.temperature_K

    def add(self, name: str, log_Z: float, c_M: float = 1.0) -> None:
        if name in self.ligands:
            raise ValueError(f"Ligand '{name}' already registered")
        if c_M <= 0.0:
            raise ValueError("c_M must be > 0")
        if c_M > 1e3:
            raise ValueError("c_M > 1000 M — convert µM/nM to M first")
        self.ligands[name] = LigandEntry(name, log_Z, c_M)

    def set_concentration(self, name: str, c_M: float) -> None:
        if name not in self.ligands:
            raise KeyError(f"Ligand '{name}' not found")
        if c_M <= 0.0:
            raise ValueError("c_M must be > 0")
        if c_M > 1e3:
            raise ValueError("c_M > 1000 M")
        self.ligands[name].concentration_M = c_M

    def log_Xi(self) -> float:
        # log_sum_exp(0, log_zZ_i...)
        terms = [0.0] + [L.log_zZ for L in self.ligands.values()]
        m = max(terms)
        s = sum(math.exp(t - m) for t in terms)
        return m + math.log(s)

    def binding_probability(self, name: str) -> float:
        if name not in self.ligands:
            raise KeyError(f"Ligand '{name}' not found")
        return math.exp(self.ligands[name].log_zZ - self.log_Xi())

    def empty_probability(self) -> float:
        return math.exp(-self.log_Xi())

    def mean_N(self) -> float:
        return 1.0 - self.empty_probability()

    def mixing_entropy(self) -> float:
        lxi = self.log_Xi()
        probs = [math.exp(-lxi)]
        probs += [math.exp(L.log_zZ - lxi) for L in self.ligands.values()]
        s = 0.0
        for p in probs:
            if p > 1e-300:
                s -= p * math.log(p)
        return kB_kcal * s

    def ligand_entropy_collapse(self) -> float:
        M = len(self.ligands)
        if M <= 1:
            return 0.0
        p_bound = self.mean_N()
        if p_bound < 1e-15:
            return 0.0
        lxi = self.log_Xi()
        s = 0.0
        for L in self.ligands.values():
            pt = math.exp(L.log_zZ - lxi) / p_bound
            if pt > 1e-300:
                s -= pt * math.log(pt)
        s_max = math.log(M)
        collapse = 1.0 - s / s_max
        return max(0.0, min(1.0, collapse))

    def occupancy_vs_concentration(
        self, titrate: str, concentrations_M: Sequence[float]
    ) -> List[OccupancyPoint]:
        if titrate not in self.ligands:
            raise KeyError(f"Titrate '{titrate}' not found")
        # snapshot fixed concentrations
        fixed = {n: L.concentration_M for n, L in self.ligands.items()}
        log_Z = {n: L.log_Z for n, L in self.ligands.items()}
        curve: List[OccupancyPoint] = []
        for c in concentrations_M:
            if c <= 0.0:
                continue
            tmp = CompetitiveSite(self.temperature_K)
            for n, lz in log_Z.items():
                tmp.add(n, lz, c if n == titrate else fixed[n])
            curve.append(
                OccupancyPoint(
                    concentration_M=c,
                    p_bound=tmp.mean_N(),
                    p_species=tmp.binding_probability(titrate),
                    mean_N=tmp.mean_N(),
                )
            )
        return curve


def set_concentration(
    site: CompetitiveSite,
    concentrations: Union[Mapping[str, float], Sequence[float]],
    names: Optional[Sequence[str]] = None,
) -> CompetitiveSite:
    """NRGsuite-facing API: ``set_concentration([L1, L2, ...])``.

    Parameters
    ----------
    site
        CompetitiveSite (or any object with set_concentration(name, c_M)).
    concentrations
        Either a dict ``{name: c_M}`` or a sequence of concentrations
        parallel to ``names`` (or to ``sorted(site.ligands)`` if names is None).
    names
        Required when ``concentrations`` is a sequence without a dict.

    Returns
    -------
    The same site instance (mutated) for chaining.
    """
    if isinstance(concentrations, Mapping):
        for name, c in concentrations.items():
            site.set_concentration(name, float(c))
        return site

    # sequence form
    conc_list = [float(x) for x in concentrations]
    if names is None:
        names = list(site.ligands.keys())
    if len(names) != len(conc_list):
        raise ValueError("names and concentrations length mismatch")
    for n, c in zip(names, conc_list):
        site.set_concentration(n, c)
    return site


def plot_occupancy_curve(
    points: Sequence[OccupancyPoint],
    *,
    title: str = "Fractional occupancy vs concentration",
    species_label: str = "p_species",
) -> Optional[object]:
    """Optional matplotlib visualization of occupancy_vs_concentration results.

    Returns the Axes if matplotlib is available, else None.
    """
    try:
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        return None

    xs = [p.concentration_M for p in points]
    ys = [p.p_species for p in points]
    yb = [p.p_bound for p in points]
    fig, ax = plt.subplots()
    ax.semilogx(xs, ys, "o-", label=species_label)
    ax.semilogx(xs, yb, "s--", label="p_bound (any)")
    ax.set_xlabel("concentration (M)")
    ax.set_ylabel("fractional occupancy")
    ax.set_title(title)
    ax.set_ylim(-0.05, 1.05)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    return ax
