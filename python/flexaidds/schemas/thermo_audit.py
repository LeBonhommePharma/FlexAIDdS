"""entropy.help / A1.1 thermodynamic audit schema.

This module defines the auditable TotalSampledPartitionFunction and
ThermodynamicOutput types used for public thermodynamic validation.

It is deliberately separate from the existing ThermodynamicBreakdown
(per-mode + corrections ledger).
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, TypedDict

# Physical constant (mirrors statmech / thermodynamics.py)
kB_kcal = 0.001987206


# ──────────────────────────────────────────────────────────────────────────────
# TypedDicts (for JSON / strict schema validation)
# ──────────────────────────────────────────────────────────────────────────────

class TotalSampledPartitionFunction(TypedDict):
    """Raw-ensemble configurational thermodynamics (pre-clustering)."""
    logZ_total_sampled: float
    F_config_kcal_mol: float
    H_eff_kcal_mol: float
    S_config_kcal_mol_K: float


class Provenance(TypedDict):
    """Self-describing audit metadata required for credibility."""
    temperature_K: float
    n_samples: int
    git_sha: str
    timestamp: str
    gate_results: Dict[str, Any]
    seed: str | int | None
    runner_info: str | None
    engine_version: str | None


class ThermodynamicOutput(TypedDict):
    """Top-level auditable container for entropy.help / public ledger."""
    total_sampled: TotalSampledPartitionFunction
    temperature_K: float
    n_samples_raw: int
    provenance: Provenance
    raw_ensemble_digest: str | None


# ──────────────────────────────────────────────────────────────────────────────
# Dataclasses (ergonomic runtime usage, matches project conventions)
# ──────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class TotalSampledPartitionFunctionDC:
    logZ_total_sampled: float
    F_config_kcal_mol: float
    H_eff_kcal_mol: float
    S_config_kcal_mol_K: float

    def to_dict(self) -> TotalSampledPartitionFunction:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "TotalSampledPartitionFunctionDC":
        return cls(
            logZ_total_sampled=float(data["logZ_total_sampled"]),
            F_config_kcal_mol=float(data["F_config_kcal_mol"]),
            H_eff_kcal_mol=float(data["H_eff_kcal_mol"]),
            S_config_kcal_mol_K=float(data["S_config_kcal_mol_K"]),
        )


@dataclass(frozen=True)
class ProvenanceDC:
    temperature_K: float
    n_samples: int
    git_sha: str
    timestamp: str
    gate_results: Dict[str, Any] = field(default_factory=dict)
    seed: str | int | None = None
    runner_info: str | None = None
    engine_version: str | None = None

    def to_dict(self) -> Provenance:
        d = asdict(self)
        # Ensure gate_results is always a dict
        if d.get("gate_results") is None:
            d["gate_results"] = {}
        return d  # type: ignore[return-value]

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ProvenanceDC":
        return cls(
            temperature_K=float(data["temperature_K"]),
            n_samples=int(data["n_samples"]),
            git_sha=str(data["git_sha"]),
            timestamp=str(data["timestamp"]),
            gate_results=dict(data.get("gate_results", {})),
            seed=data.get("seed"),
            runner_info=data.get("runner_info"),
            engine_version=data.get("engine_version"),
        )


@dataclass(frozen=True)
class ThermodynamicOutputDC:
    """Primary public API type for entropy.help audits."""

    total_sampled: TotalSampledPartitionFunctionDC
    temperature_K: float
    n_samples_raw: int
    provenance: ProvenanceDC
    raw_ensemble_digest: str | None = None

    def to_dict(self) -> ThermodynamicOutput:
        return {
            "total_sampled": self.total_sampled.to_dict(),
            "temperature_K": self.temperature_K,
            "n_samples_raw": self.n_samples_raw,
            "provenance": self.provenance.to_dict(),
            "raw_ensemble_digest": self.raw_ensemble_digest,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ThermodynamicOutputDC":
        return cls(
            total_sampled=TotalSampledPartitionFunctionDC.from_dict(
                data["total_sampled"]
            ),
            temperature_K=float(data["temperature_K"]),
            n_samples_raw=int(data["n_samples_raw"]),
            provenance=ProvenanceDC.from_dict(data["provenance"]),
            raw_ensemble_digest=data.get("raw_ensemble_digest"),
        )

    # ─── Validation (enforces rules from THERMODYNAMIC_OUTPUT_SCHEMA.md) ────
    def validate(self) -> None:
        """Raise on any violation of the audit-grade schema contract."""
        ts = self.total_sampled
        prov = self.provenance

        if self.temperature_K <= 0:
            raise ValueError("temperature_K must be > 0")
        if self.n_samples_raw < 1:
            raise ValueError("n_samples_raw must be >= 1")

        # Consistency: F = -kT * logZ (within floating-point tolerance)
        expected_F = -kB_kcal * self.temperature_K * ts.logZ_total_sampled
        if abs(ts.F_config_kcal_mol - expected_F) > 1e-8:
            raise ValueError(
                f"F_config inconsistency: got {ts.F_config_kcal_mol}, "
                f"expected {expected_F}"
            )

        # S derivation check
        derived_S = (ts.H_eff_kcal_mol - ts.F_config_kcal_mol) / self.temperature_K
        if abs(ts.S_config_kcal_mol_K - derived_S) > 1e-8:
            raise ValueError("S_config derivation does not match (H - F)/T")

        # Gate results presence (minimum for audit credibility)
        if not isinstance(prov.gate_results, dict):
            raise ValueError("provenance.gate_results must be a dict")

        # All values must be finite
        for name, val in [
            ("logZ_total_sampled", ts.logZ_total_sampled),
            ("F_config_kcal_mol", ts.F_config_kcal_mol),
            ("H_eff_kcal_mol", ts.H_eff_kcal_mol),
            ("S_config_kcal_mol_K", ts.S_config_kcal_mol_K),
        ]:
            if not (val == val and abs(val) != float("inf")):  # NaN / Inf check
                raise ValueError(f"{name} must be finite")


# ──────────────────────────────────────────────────────────────────────────────
# Convenience constructors (used by runners / StatMechEngine wrappers)
# ──────────────────────────────────────────────────────────────────────────────

def make_total_sampled_output(
    *,
    logZ: float,
    mean_energy: float,
    temperature_K: float,
    n_samples: int,
    git_sha: str,
    timestamp: str,
    gate_results: Dict[str, Any],
    raw_ensemble_digest: str | None = None,
    **provenance_extras: Any,
) -> ThermodynamicOutputDC:
    """Factory used by the future raw-ensemble path in A2.1 / runners."""
    F = -kB_kcal * temperature_K * logZ
    S = (mean_energy - F) / temperature_K

    total = TotalSampledPartitionFunctionDC(
        logZ_total_sampled=logZ,
        F_config_kcal_mol=F,
        H_eff_kcal_mol=mean_energy,
        S_config_kcal_mol_K=S,
    )

    prov = ProvenanceDC(
        temperature_K=temperature_K,
        n_samples=n_samples,
        git_sha=git_sha,
        timestamp=timestamp,
        gate_results=gate_results,
        **provenance_extras,
    )

    out = ThermodynamicOutputDC(
        total_sampled=total,
        temperature_K=temperature_K,
        n_samples_raw=n_samples,
        provenance=prov,
        raw_ensemble_digest=raw_ensemble_digest,
    )
    out.validate()
    return out
