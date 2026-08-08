# Entropy Audit Report

**Status**: TEMPLATE — NOT AN AUDIT RESULT
All values and verdict choices below are placeholders. A populated copy must not
be described as complete, published, validated, verified, or reproducible until
its linked JSON report and separate provenance record pass
`scripts/validate_thermo_claims.py`.

**Audit ID**: AUD-YYYY-NNN  
**Target**: [Receptor] + [Ligand]  
**Engine**: FlexAIDdS / [Other]  
**Report Version**: 1.0  
**Date**: YYYY-MM-DD  
**Auditor**: entropy.help (automated + manual review)

---

## 1. Summary

| Quantity                  | Value                  | Notes |
|---------------------------|------------------------|-------|
| Total Sampled Partition Function (log Z) | `logZ_total_sampled` | Raw ensemble |
| F_config (kcal/mol)       | `F_config_kcal_mol`    | −kT ln Z_sampled |
| S_config (kcal mol⁻¹ K⁻¹) | `S_config_kcal_mol_K`  | (⟨E⟩ − F_config) / T |
| H_eff (kcal/mol)          | `H_eff_kcal_mol`       | Boltzmann-weighted mean energy |
| Temperature (K)           | `temperature_K`        | — |
| Raw poses sampled         | `n_samples_raw`        | Before clustering |

**Gate 5 (Partition Convergence)**: [PASSED / FAILED] — delta_logZ = X.XXe-YY  
**Gate 6 (F/S Cross-check)**: [PASSED / FAILED] — max deviation = X.XX kcal/mol

**Overall Verdict**: [UNVERIFIED / MINOR DISCREPANCY / SIGNIFICANT DEVIATION / REPRODUCIBLE only with deposited evidence]

---

## 2. Provenance

- **Git SHA**: `full_or_short_sha`
- **Timestamp (ISO-8601)**: `...`
- **Runner / Campaign**: `...`
- **Random seed(s)**: `...`
- **Host / Hardware**: `...`
- **Docking binary / version**: `...`

**Raw ensemble digest** (placeholder; a published quantitative audit requires a provenance receipt here): `[sha256 of the deposited ensemble]`

---

## 3. Input Data

- Receptor file(s): `path(s)`
- Ligand file(s): `path(s)`
- Full docking output directory: `path`
- Reference experimental data (if ITC or literature): `...`

---

## 4. Thermodynamic Ledger (TotalSampled)

```json
{
  "total_sampled": {
    "logZ_total_sampled": null,
    "F_config_kcal_mol": null,
    "H_eff_kcal_mol": null,
    "S_config_kcal_mol_K": null
  },
  "temperature_K": null,
  "n_samples_raw": null,
  "provenance": { ... }
}
```

(A completed report requires a provenance receipt: link the real JSON sidecar and the separate provenance JSON here.)

---

## 5. Comparison vs Original Engine Output

| Metric                    | Original Report | This Audit | Δ |
|---------------------------|-----------------|------------|---|
| Top pose free energy      |                 |            |   |
| Selected binding mode     |                 |            |   |
| Rank of crystallographic pose |             |            |   |
| S_config contribution     |                 |            |   |

**Explanation of discrepancies** (if any):

---

## 6. Methodology Notes

- Partition function computed via log-sum-exp over the complete raw pose ensemble (pre-clustering).
- Multiplicity handling: each GA chromosome contributed with its sampling weight.
- Vibrational / other corrections: [none / tENCoM applied / ...]
- Concentration / reference state handling: [N/A for this audit / ...]

---

## 7. Signature

**Content hash (SHA-256 of canonical JSON)**: `[computed 64-hex digest; never a placeholder]`

**Auditor signature** (optional GPG or equivalent): attach a real detached
signature, or omit this field. Never publish placeholder signature armor.

**Report URL** (permanent, only after artifact publication): `[repository-backed URL]`

---

*This template is version-controlled alongside the entropy.help tooling.  
Suggestions and improvements are welcome via the public coordination issue:  
https://github.com/LeBonhommePharma/FlexAIDdS/issues/219*
