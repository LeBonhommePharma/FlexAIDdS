"""Pure-Python fallback for Grand Canonical Partition Function (Ξ).

This module provides _PyGrandPartitionFunction (and public alias when no C++ grand
bindings) that exactly mirrors the C++ target::GrandPartitionFunction API
and numerical behavior (log-space, log-sum-exp with empty-site anchor=0,
intrinsic vs apparent selectivity, LigandRank, etc.).

It is additive and does not affect single-ligand canonical paths.
Use HAS_GRAND_BINDINGS (set when C++ bindings for grand are present) to select.

Intended for:
- compute_grand_partition( ligand_logZs or StatMech results, concs )
- augmenting DockingResult / competitive manifests
- roundtrips and tests before/after C++ pybind11 wiring (P2 parallel track)

See GPF_IMPLEMENTATION_PLAN.md and LIB/GrandPartitionFunction.{h,cpp}
for contract, math, and edge cases.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

from .thermodynamics import (
    StatMechEngine,
    Thermodynamics,
    kB_kcal,
)
from .models import LigandSpec  # canonical definition (also re-exported at package level)

# Standard reference concentration c° = 1 M (IUPAC)
C_STANDARD = 1.0


@dataclass(frozen=True)
class LigandRank:
    """One entry from GrandPartitionFunction.rank().

    Mirrors C++ GrandPartitionFunction::LigandRank.
    Sorted by dG ascending (most favorable first).
    """
    name: str
    log_Z: float          # intrinsic ln(Z_i)
    dG: float             # -kT * log_Z   (F_bound)
    p_bound: float        # z_i Z_i / Ξ  (p_bind_like; CF-proxy occupancy, not calibrated ΔG)


def _logsumexp_with_anchor(log_zZ_values: List[float]) -> float:
    """Replicate C++ compute_log_Xi_fresh exactly.

    Ξ = 1 + Σ zZ_i
    lnΞ = logsumexp( 0 , log_zZ_0, log_zZ_1, ... )
    Anchor 0 for empty (apo) site is *always* considered for max and sum.
    """
    if not log_zZ_values:
        return 0.0
    max_val = 0.0
    for v in log_zZ_values:
        if v > max_val:
            max_val = v
    # empty contribution
    s = math.exp(0.0 - max_val)
    for v in log_zZ_values:
        s += math.exp(v - max_val)
    return max_val + math.log(s)


class _PyGrandPartitionFunction:
    """Pure-Python implementation of GrandPartitionFunction.

    All operations in log-space. Matches C++ numerical results within fp tolerance.
    Not thread-safe (no mutex); use one per binding site in Python usage.
    """

    def __init__(self, temperature_K: float = 300.0):
        if temperature_K <= 0.0:
            raise ValueError(f"Temperature must be positive (got {temperature_K})")
        self._T: float = temperature_K
        self._beta: float = 1.0 / (kB_kcal * temperature_K)
        # name -> {'log_Z': float, 'log_c': float, 'log_zZ': float}
        self._ligands: Dict[str, Dict[str, float]] = {}
        self._cached_log_xi: Optional[float] = None

    @property
    def temperature(self) -> float:
        return self._T

    def _invalidate(self) -> None:
        self._cached_log_xi = None

    # ── Registration (mirror signatures) ─────────────────────────────────

    def add_ligand(
        self,
        name: str,
        log_Z_or_engine: Union[float, "StatMechEngine"],
        concentration_M: float = 1.0,
    ) -> None:
        """Register ligand.

        Accepts either log_Z (float) or StatMechEngine (extracts .compute().log_Z).
        concentration_M in molar (1.0 = standard state).
        """
        if concentration_M <= 0.0:
            raise ValueError(
                f"concentration_M must be > 0 (got {concentration_M} M)"
            )
        if concentration_M > 1000.0:
            raise ValueError(
                "Concentration > 1000 M — did you pass µM or nM without conversion to M?"
            )

        if isinstance(log_Z_or_engine, (int, float)):
            log_Z = float(log_Z_or_engine)
        elif hasattr(log_Z_or_engine, "compute"):
            # StatMechEngine-like
            thermo = log_Z_or_engine.compute()
            log_Z = float(thermo.log_Z)
        else:
            # allow objects with .log_Z or thermodynamics attr
            if hasattr(log_Z_or_engine, "log_Z"):
                log_Z = float(log_Z_or_engine.log_Z)
            else:
                thermo = getattr(log_Z_or_engine, "get_thermodynamics", lambda: None)()
                if thermo is None:
                    raise TypeError("add_ligand expects float log_Z or StatMechEngine-like object")
                log_Z = float(thermo.log_Z)

        log_c = math.log(concentration_M / C_STANDARD)
        log_zZ = log_c + log_Z

        if name in self._ligands:
            raise ValueError(f"Ligand '{name}' already registered")
        self._ligands[name] = {"log_Z": log_Z, "log_c": log_c, "log_zZ": log_zZ}
        self._invalidate()

    def add_or_overwrite(
        self,
        name: str,
        log_Z: float,
        concentration_M: float = 1.0,
    ) -> None:
        if concentration_M <= 0.0:
            raise ValueError(f"concentration_M must be > 0 (got {concentration_M} M)")
        if concentration_M > 1000.0:
            raise ValueError(
                "Concentration > 1000 M — did you pass µM or nM without conversion to M?"
            )
        log_c = math.log(concentration_M / C_STANDARD)
        log_zZ = log_c + log_Z
        if name in self._ligands:
            self._ligands[name]["log_Z"] = log_Z
            self._ligands[name]["log_c"] = log_c
            self._ligands[name]["log_zZ"] = log_zZ
        else:
            self._ligands[name] = {"log_Z": log_Z, "log_c": log_c, "log_zZ": log_zZ}
        self._invalidate()

    def overwrite_ligand(self, name: str, new_log_Z: float) -> None:
        if name not in self._ligands:
            raise ValueError(f"Ligand '{name}' not found")
        entry = self._ligands[name]
        entry["log_Z"] = new_log_Z
        entry["log_zZ"] = entry["log_c"] + new_log_Z
        self._invalidate()

    def merge_ligand(self, name: str, new_log_Z: float) -> None:
        if name not in self._ligands:
            raise ValueError(f"Ligand '{name}' not found")
        entry = self._ligands[name]
        log_c = entry["log_c"]
        a = entry["log_zZ"]
        b = log_c + new_log_Z
        max_val = max(a, b)
        merged_zZ = max_val + math.log(math.exp(a - max_val) + math.exp(b - max_val))
        entry["log_zZ"] = merged_zZ
        entry["log_Z"] = merged_zZ - log_c
        self._invalidate()

    def remove_ligand(self, name: str) -> None:
        if name not in self._ligands:
            raise ValueError(f"Ligand '{name}' not found")
        del self._ligands[name]
        self._invalidate()

    # ── Queries (log-space, match C++ contract) ──────────────────────────

    def log_Xi(self) -> float:
        if self._cached_log_xi is None:
            self._cached_log_xi = self._compute_log_Xi_fresh()
        return self._cached_log_xi

    def _compute_log_Xi_fresh(self) -> float:
        if not self._ligands:
            return 0.0
        log_zZs = [e["log_zZ"] for e in self._ligands.values()]
        return _logsumexp_with_anchor(log_zZs)

    def binding_probability(self, name: str) -> float:
        if name not in self._ligands:
            raise ValueError(f"Ligand '{name}' not found")
        log_xi = self.log_Xi()
        log_zZ = self._ligands[name]["log_zZ"]
        # guard extreme
        if log_xi > 700 and (log_zZ - log_xi) < -700:
            return 0.0
        return math.exp(log_zZ - log_xi)

    def empty_probability(self) -> float:
        log_xi = self.log_Xi()
        return math.exp(-log_xi)

    def mean_occupancy(self) -> float:
        return 1.0 - self.empty_probability()

    def occupancy_variance(self) -> float:
        mu = self.mean_occupancy()
        return mu * (1.0 - mu)

    def F_bound(self, name: str) -> float:
        if name not in self._ligands:
            raise ValueError(f"Ligand '{name}' not found")
        log_Z = self._ligands[name]["log_Z"]
        kT = 1.0 / self._beta
        return -kT * log_Z

    def delta_G_bind(self, name: str, F_ref: float = 0.0) -> float:
        return self.F_bound(name) - F_ref

    def selectivity(self, a: str, b: str) -> float:
        diff = self.log_selectivity(a, b)
        if diff > 700.0:
            return float("inf")  # but C++ uses DBL_MAX sentinel; for py inf ok in most contexts, or sys.float_info.max
        if diff < -700.0:
            return 0.0
        return math.exp(diff)

    def log_selectivity(self, a: str, b: str) -> float:
        if a not in self._ligands:
            raise ValueError(f"Ligand '{a}' not found")
        if b not in self._ligands:
            raise ValueError(f"Ligand '{b}' not found")
        return self._ligands[a]["log_zZ"] - self._ligands[b]["log_zZ"]

    def log_intrinsic_selectivity(self, a: str, b: str) -> float:
        # MUST be log_Z only, independent of conc (see C++ comments + tests)
        if a not in self._ligands:
            raise ValueError(f"Ligand '{a}' not found")
        if b not in self._ligands:
            raise ValueError(f"Ligand '{b}' not found")
        return self._ligands[a]["log_Z"] - self._ligands[b]["log_Z"]

    def rank(self) -> List[LigandRank]:
        log_xi = self.log_Xi()
        kT = 1.0 / self._beta
        ranks: List[LigandRank] = []
        for name, entry in self._ligands.items():
            p = math.exp(entry["log_zZ"] - log_xi) if log_xi != float("-inf") else 0.0
            ranks.append(
                LigandRank(
                    name=name,
                    log_Z=entry["log_Z"],
                    dG=-kT * entry["log_Z"],
                    p_bound=p,
                )
            )
        ranks.sort(key=lambda r: r.dG)
        return ranks

    def num_ligands(self) -> int:
        return len(self._ligands)

    def has_ligand(self, name: str) -> bool:
        return name in self._ligands

    def all_log_Z(self) -> List[Tuple[str, float]]:
        return [(name, e["log_Z"]) for name, e in self._ligands.items()]

    def all_log_zZ(self) -> List[Tuple[str, float]]:
        return [(name, e["log_zZ"]) for name, e in self._ligands.items()]

    # Convenience: return dict of current p(bound) for all + p(empty)
    def probabilities(self) -> Dict[str, float]:
        """Return {ligand_name: p_bound, '__empty__': p_empty}."""
        res: Dict[str, float] = {"__empty__": self.empty_probability()}
        for name in self._ligands:
            res[name] = self.binding_probability(name)
        return res


# Public alias used by __init__.py conditional
GrandPartitionFunction = _PyGrandPartitionFunction


def compute_grand_partition(
    ligands: List[Union[LigandSpec, Tuple[str, float, float]]],
    *,
    temperature_K: float = 300.0,
    # alternative: pass precomputed log_Z map
    log_Z_map: Optional[Dict[str, float]] = None,
) -> _PyGrandPartitionFunction:
    """High-level helper: build a GPF from specs or logZ map.

    ligands: list of LigandSpec or (name, log_Z, conc_M) tuples.
    If log_Z_map provided, use it for log_Z lookup by name (conc from spec).

    Returns a populated _PyGrandPartitionFunction (or C++ when wired).
    """
    gpf = _PyGrandPartitionFunction(temperature_K=temperature_K)
    for item in ligands:
        if isinstance(item, LigandSpec):
            name = item.name
            conc = item.concentration_M
            if log_Z_map and name in log_Z_map:
                log_Z = log_Z_map[name]
            else:
                # caller must supply via map or use add after
                raise ValueError(
                    f"No log_Z for {name}; pass log_Z_map or pre-populate"
                )
            gpf.add_ligand(name, log_Z, conc)
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            name = item[0]
            log_Z = item[1]
            conc = item[2] if len(item) > 2 else 1.0
            gpf.add_ligand(name, float(log_Z), float(conc))
        else:
            raise TypeError("ligands items must be LigandSpec or (name, log_Z, [conc])")
    return gpf


# Note: when C++ grand bindings land (future), we will do:
# try:
#     from ._core import GrandPartitionFunction as _CppGrand...
#     HAS_GRAND_BINDINGS = True
# except:
#     ...
#     HAS_GRAND_BINDINGS = False
# and GrandPartitionFunction = _Cpp... if available else _Py...
# For P2 we expose the pure always under the name, flag starts False.
