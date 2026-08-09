# Thermodynamic Output Schema (TotalSampledPartitionFunction)

**Version**: 1.0  
**Status**: Draft for implementation (A1.1)  
**Last Updated**: 2026-05-27  
**Related**: [GitHub #219](https://github.com/LeBonhommePharma/FlexAIDdS/issues/219), A2.1, B2.1 (public audit ledger), B4.1 (report template)

---

> **Non-claiming design document.** This schema is not evidence that the proposed
> raw-ensemble accumulator exists in a validated production path. In particular,
> FlexAIDdS CF/contact-function scores are ranking proxies, not automatically
> kcal/mol energies. Fields bearing `*_kcal_mol` names may be populated only when
> the input-energy calibration and reference-state provenance justify those units;
> otherwise an implementation must label the output as a score-space proxy.

## 1. Purpose

This schema proposes an **auditable, engine-agnostic** representation of finite-sample ensemble calculations over the *raw* (pre-clustering) states retained by a docking run.

It is deliberately distinct from:
- `ThermodynamicBreakdown` (existing per-mode + correction ledger)
- `GrandPartitionFunction` (competitive binding / grand canonical; see its implementation plan and fixtures, whose validation status must be established from their own artifacts)

The proposed construct is the **TotalSampledPartitionFunction** — a finite-sample partition calculation over every deposited pose and declared multiplicity. Whether it is canonical thermodynamics or a score-space diagnostic depends on the provenance and units of its inputs.

## 2. Core Concepts

- **Z_sampled** (Total Sampled Partition Function): `Σ n_i * exp(−β E_i)` over the deposited raw ensemble under a declared energy/score convention.
- **F_config**: −kT ln Z_sampled, only a configurational Helmholtz free energy when `E_i` and β have validated physical units.
- **S_config**: (⟨E⟩ − F_config) / T (or equivalently −kB Σ p_i ln p_i), subject to the same unit and sampling assumptions.
- **Provenance + Gate Results**: Mandatory self-describing metadata for public audit credibility (Gates 5 & 6 from the 7-day launch plan).

These quantities cover the declared finite sample only. They exclude vibrational, solvation, concentration, standard-state, and bound/unbound reference terms unless those appear as separately sourced corrections.

### 2.1 Scientific provenance and claim validity (schema v2)

Every runtime thermodynamic record carries a nested `scientific_provenance`
object with these source fields:

```json
{
  "schema_version": 2,
  "energy_domain": "unclassified | cf_arbitrary_units | calibrated_kcal_per_mol | model_scale",
  "ensemble_measure": "unclassified | optimizer_samples | enumerated_microstates | weighted_quadrature",
  "reference_state": "none | bound_only | matched_association_cycle",
  "energy_provenance": "sha256:<64 hex>",
  "measure_provenance": "sha256:<64 hex>",
  "reference_provenance": "sha256:<64 hex>",
  "claim_validity": "proxy_only | canonical_physical | binding_physical"
}
```

`claim_validity` is derived output and is never trusted on input. Canonical
physical status requires calibrated kcal/mol energies, an enumerated or
weighted-quadrature measure, and nontrivial SHA-256 identities for both source
artifacts. Binding status additionally requires a matched association cycle and
its artifact identity. All-zero, low-diversity, malformed, and the historical
example filler digest are rejected. The digest check is only a syntactic gate;
publication must also deposit and verify the referenced artifacts. A result
marked unavailable can never authorize a physical claim.

## 3. Units & Precision

- Physical energy fields: kcal/mol, permitted only with explicit unit provenance
- Physical entropy fields: kcal mol⁻¹ K⁻¹, permitted only with explicit unit provenance
- Score-space calculations: use an explicit proxy label and do not serialize values into `*_kcal_mol` fields
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

See `audit-report-example.json` for a synthetic, null-valued shape example. It is intentionally unsigned and carries no result or gate verdict.

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

- **Runner / DatasetRunner (planned)**: A future implementation would compute the output from the deposited raw pose list + multiplicities before clustering.
- **Python results loader (planned)**: A future implementation may attach `ThermodynamicOutput` (or its dict form) to result metadata.
- **Public Registry (planned)**: Any future published audit must store the full output plus separate provenance and human-readable artifacts.
- **Report Template**: The human-readable section "4. Thermodynamic Ledger (TotalSampled)" is a placeholder for values derived from a real artifact.

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

This document is a draft design note until superseded by an implemented, tested formal JSON Schema. It is not a completion receipt.
