# Thermodynamic Output Schema (TotalSampledPartitionFunction)

**Version**: 1.0  
**Status**: Draft for implementation (A1.1)  
**Last Updated**: 2026-05-27  
**Related**: [GitHub #219](https://github.com/LeBonhommePharma/FlexAIDdS/issues/219), A2.1, B2.1 (public audit ledger), B4.1 (report template)

---

## 1. Purpose

This schema defines the **auditable, engine-agnostic** representation of configurational thermodynamics computed over the *raw* (pre-clustering) conformational ensemble sampled by a docking run.

It is deliberately distinct from:
- `ThermodynamicBreakdown` (existing per-mode + correction ledger)
- `GrandPartitionFunction` (competitive binding / grand canonical; see GPF_IMPLEMENTATION_PLAN.md P3/P5, competition_example.yaml, grand_calibrate.py, *_grand_summary.csv emission, --conc CLI, DatasetRunner grand_summary for P3+ outputs + verified 3-ligand cases)

The new construct is the **TotalSampledPartitionFunction** — the finite-sample canonical partition function over every pose and multiplicity generated during the GA run.

## 2. Core Concepts

- **Z_sampled** (Total Sampled Partition Function): `Σ n_i * exp(−β E_i)` over the complete raw ensemble.
- **F_config**: −kT ln Z_sampled (configurational Helmholtz free energy)
- **S_config**: (⟨E⟩ − F_config) / T   (or equivalently −kB Σ p_i ln p_i)
- **Provenance + Gate Results**: Mandatory self-describing metadata for public audit credibility (Gates 5 & 6 from the 7-day launch plan).

All quantities are **configurational only** until explicit corrections (vib, natural, etc.) are added in higher layers.

## 3. Units & Precision

- Energies: kcal/mol
- Entropy: kcal mol⁻¹ K⁻¹
- Temperature: K
- All floating-point values must be finite (no NaN/Inf).
- Recommended output precision: 6–8 decimal places in JSON for audit reproducibility.

## 4. Schema Definitions

### 4.1 TotalSampledPartitionFunction

```typescript
interface TotalSampledPartitionFunction {
  logZ_total_sampled: number;      // ln(Z_sampled) — numerically stable storage
  F_config_kcal_mol: number;       // −kT * logZ_total_sampled
  H_eff_kcal_mol: number;          // Boltzmann-weighted mean energy ⟨E⟩
  S_config_kcal_mol_K: number;     // (H_eff − F_config) / T
}
```

### 4.2 Provenance

```typescript
interface Provenance {
  temperature_K: number;
  n_samples: number;               // raw pose count (with multiplicity)
  git_sha: string;                 // full preferred
  timestamp: string;               // ISO-8601
  gate_results: Record<string, any>; // structured results for Gates 5/6
  seed?: number | string;
  runner_info?: string;
  engine_version?: string;
}
```

### 4.3 ThermodynamicOutput (Top Level)

```typescript
interface ThermodynamicOutput {
  total_sampled: TotalSampledPartitionFunction;
  temperature_K: number;
  n_samples_raw: number;
  provenance: Provenance;
  raw_ensemble_digest?: string;    // optional sha256 of sorted energies for extra repro
  // Future: per_mode_refs, component_corrections, etc.
}
```

## 5. Validation Rules

- `temperature_K > 0`
- `n_samples_raw >= 1`
- `|F_config_kcal_mol + (kB_kcal * temperature_K * logZ_total_sampled)| < 1e-9` (consistency)
- `S_config_kcal_mol_K` derivation must match within 1e-10
- `gate_results` must contain at minimum keys for the active Gates (e.g. `gate5_convergence`, `gate6_crosscheck`)
- All energies must be finite

## 6. Serialization Examples

### JSON (for ledger / signed sidecar)

See `audit-report-example.json` in this directory for a realistic populated instance.

### Python (TypedDict + dataclass)

See `python/flexaidds/schemas/thermo_audit.py` (companion implementation).

### C++ (header sketch for A2.1)

```cpp
// In LIB/ or new audit/ namespace — sketch only
namespace statmech {

struct TotalSampledPartitionFunction {
    double logZ_total_sampled = 0.0;
    double F_config_kcal_mol = 0.0;
    double H_eff_kcal_mol = 0.0;
    double S_config_kcal_mol_K = 0.0;
};

struct Provenance {
    double temperature_K = 300.0;
    std::size_t n_samples = 0;
    std::string git_sha;
    std::string timestamp;
    // gate_results stored as nlohmann::json or simple map<string, variant>
};

struct ThermodynamicOutput {
    TotalSampledPartitionFunction total;
    Provenance provenance;
    std::string raw_ensemble_digest;   // optional

    std::string to_json_string() const; // implementation in .cpp
    bool validate() const;
};

} // namespace statmech
```

Implementation note for A2.1: Source `logZ_total_sampled` from a new raw-ensemble accumulator path (distinct from `BindingMode::get_thermodynamic_breakdown` and `get_global_ensemble`).

## 7. Integration Points

- **Runner / DatasetRunner**: After a docking job, call `compute_total_sampled_output()` on the raw pose list + multiplicities before clustering.
- **Python results loader**: Attach `ThermodynamicOutput` (or its dict form) to `DockingResult` / `BindingPopulation` metadata.
- **Public Ledger (B2)**: Each published audit stores the full `ThermodynamicOutput` as the canonical signed artifact.
- **Report Template (B4.1)**: The human-readable section "4. Thermodynamic Ledger (TotalSampled)" is populated directly from this schema.

## 8. Relationship to Existing Types

| Existing Type            | Scope                  | When to use vs. new schema |
|--------------------------|------------------------|----------------------------|
| `Thermodynamics`         | Single ensemble        | Internal StatMechEngine    |
| `ThermodynamicBreakdown` | Mode + corrections     | Per-mode reporting + vib/natural |
| `TotalSampled...` (new)  | Raw full GA ensemble   | Audit, reproducibility, external verification |

Never conflate the raw `TotalSampled` Z with any clustered / binding-mode Z.

## 9. Future Evolution

- v1.1: Add optional `raw_energies_histogram` for visual audit.
- Component corrections (vib, natural) will be added as sibling fields under a `corrections` object, not mutating the core `total_sampled`.

---

**Implementation owner**: A2.1 (C++ core) + Python/TS bindings.  
**Reviewers**: Thermodynamics maintainers + entropy.help auditors.

This document is the single source of truth for the schema until superseded by a formal JSON Schema file in a later revision.