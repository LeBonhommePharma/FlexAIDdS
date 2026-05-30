# FlexAID∆S Thermodynamics

This document describes the thermodynamic and ensemble analysis capabilities in FlexAID∆S (post Phase 1–3 + roadmap Tasks 0–9).

## Core Ensemble Thermodynamics

FlexAID∆S computes the canonical ensemble over the genetic algorithm conformational ensemble using the StatMechEngine.

### Key Quantities (with units)

- **logZ** — Natural log of the partition function (dimensionless)
- **G_config** (or F) — Configurational Helmholtz free energy: `G_config = -kB T logZ` (kcal/mol)
- **H_eff** (⟨E⟩) — Boltzmann-weighted mean energy (kcal/mol). This is an *effective scoring-energy proxy*, **not** true calorimetric enthalpy unless calibration is performed.
- **S_config** — Configurational entropy: `S = (H_eff - G_config) / T` (kcal mol⁻¹ K⁻¹)
- **-T S_config** — Entropic contribution to free energy (kcal/mol)
- **Cv** — Heat capacity at constant volume (kcal mol⁻¹ K⁻¹)
- **σ_E** — Standard deviation of energy (kcal/mol)

All calculations use numerically stable log-sum-exp.

## Additive Corrections (in ThermodynamicBreakdown)

The full free energy is decomposed as:

```
G_total = G_config + G_vib + G_natural + G_other
```

- **G_vib**: Vibrational entropy correction (ENCoM / tENCoM)
- **G_natural**: NATURaL co-translational / receptor strain correction
- **G_other**: Future terms (e.g., explicit GIST)

Presence is tracked with `has_vib`, `has_natural`, `has_other` flags.

## Component-wise Boltzmann Averages (Task 3)

For any per-pose energy decomposition, ensemble averages are computed as:

```
<X> = Σ_i p_i X_i
```

Where `p_i` are the Boltzmann weights from the canonical ensemble.

See `EnergyComponents` struct and `component_means` in `ThermodynamicBreakdown`.

If components are incomplete, `components_complete = false` and `component_sum` may legitimately differ from `H_eff`.

## Diagnostic Metrics (Task 4) — Use With Caution

These are **strictly diagnostic** and must never be used for ranking, pose selection, or affinity claims:

- `entropy_fraction`
- `enthalpy_fraction`
- `compensation_score = 1 - |G_config| / (|H_eff| + |-T S_config| + ε)` (clamped [0,1])

**Forbidden**: Using `(ΔH + TΔS)/ΔG` or any compensation metric as a score.

## Joint Receptor–Ligand Ensemble (Task 5 — Experimental)

`JointEnsembleResult` provides:

- Joint entropy S_joint
- Marginal receptor and ligand entropies
- Mutual information I(R;L)

Fallback mode (when receptor conformer IDs are unavailable): S_receptor = 0, I(R;L) = 0.

## Standard-State Affinity Calibration (Task 6 — Safe / Experimental)

Utilities:

- `deltaG_standard_to_Kd_M(ΔG°, T, c0)`
- `Kd_M_to_deltaG_standard(Kd, T, c0)`

**Hard rules**:
- T > 0 and Kd > 0 required (otherwise exception)
- `predicted_Kd_M` is only meaningful when `calibrated = true`
- Relative free energies may be emitted but must be labelled "relative"

No uncalibrated docking scores are ever presented as real Kd/Ki.

## Temperature Scan & Model-Derived ΔCp (Task 7 — Experimental)

- `temperature_scan(T_grid)` recomputes all quantities at multiple temperatures on **fixed** microstate energies.
- `fit_delta_Cp(...)` performs linear regression (requires ≥4 points) and reports RMSE.
- Output is always labelled `model_derived = true` and `experimental = true`.
- **Single-temperature Cv is not experimental binding heat capacity change.**

## Cleft Annotation & Flexible Residue Selection (Task 8 — Preprocessing / Experimental)

`CleftAnnotation` + `select_flexible_residues()` convert geometry into:

- Classification (orthosteric / allosteric / unknown)
- Recommended flexible residues (respecting fixed/forced lists, Gly/Ala rules, distance shells)

This is **preprocessing only** — it does not affect scoring.

External annotations must come from user-supplied files.

## Reporting & Visualization (Task 9)

All reporting is driven exclusively from `DockingResult` / JSON:

- `generate_pymol_script()` — Labels G_total, H_eff, -T*S with units and experimental warnings
- `generate_markdown_report()`
- `generate_temperature_scan_plot()`
- `write_all_reports()`

No new energy calculations that could affect ranking are performed.

## Support Classification (as of end of roadmap)

| Feature                        | Classification          |
|--------------------------------|-------------------------|
| ThermodynamicBreakdown (G, H, S, Cv, σ) | Core |
| Component-wise averages        | Core diagnostic (when data available) |
| Diagnostic compensation metrics| Diagnostic only |
| Joint receptor–ligand ensemble | Experimental |
| AffinityCalibration / Kd       | Experimental (safe only) |
| Temperature scan + ΔCp fit     | Experimental / model-derived |
| CleftAnnotation + selector     | Experimental preprocessing |
| PyMOL / Markdown / plot reports| Visualization / reporting only |

See also:
- `docs/VALIDATED_CAPABILITIES.md`
- `docs/EXPERIMENTAL_CAPABILITIES.md`
- `docs/KNOWN_LIMITATIONS.md`
- `docs/dev/thermo_invariants.md`
- `docs/dev/benchmark_acceptance_checklist.md`

## Backwards Compatibility

All legacy fields (`free_energy`, `enthalpy`, `entropy`, etc.) on `BindingModeResult` and REMARKs are preserved. New fields are additive.

---

**End of document.** All claims in this file are supported by unit tests and the implementation PRs (#212–#218).