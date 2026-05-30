# Entropy Audit Report

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

**Overall Verdict**: [REPRODUCIBLE / MINOR DISCREPANCY / SIGNIFICANT DEVIATION]

---

## 2. Provenance

- **Git SHA**: `full_or_short_sha`
- **Timestamp (ISO-8601)**: `...`
- **Runner / Campaign**: `...`
- **Random seed(s)**: `...`
- **Host / Hardware**: `...`
- **Docking binary / version**: `...`

**Raw ensemble digest** (optional, for high-reproducibility audits): `sha256:...`

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
    "logZ_total_sampled": -45.237,
    "F_config_kcal_mol": -13.42,
    "H_eff_kcal_mol": -16.81,
    "S_config_kcal_mol_K": 0.0113
  },
  "temperature_K": 300.0,
  "n_samples_raw": 12480,
  "provenance": { ... }
}
```

(Full signed JSON sidecar attached as `audit-report-*.json`.)

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

**Content hash (SHA-256 of canonical JSON)**: `...`

**Auditor signature** (optional GPG or equivalent):  
```
-----BEGIN PGP SIGNATURE-----
...
-----END PGP SIGNATURE-----
```

**Report URL** (permanent): `https://entropy.help/audits/AUD-YYYY-NNN`

---

*This template is version-controlled alongside the entropy.help tooling.  
Suggestions and improvements are welcome via the public coordination issue:  
https://github.com/LeBonhommePharma/FlexAIDdS/issues/219*
