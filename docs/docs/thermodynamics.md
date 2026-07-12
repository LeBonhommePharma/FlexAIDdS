# Thermodynamics (MkDocs summary)

**Last actualized**: 2026-07-12  

Full reference: [`docs/thermodynamics.md`](../thermodynamics.md).

## Pipeline placement

1. **GA search** ranks candidates with the **CF / Voronoi contact-function scoring proxy**.
2. **Binding modes** cluster poses from the ensemble.
3. **StatMechEngine** computes the canonical partition function over sampled energies.
4. Optional **vibrational** (ENCoM / tENCoM) and experimental corrections enter `ThermodynamicBreakdown`.

## Core identities

\[
G_{\mathrm{config}} = -k_B T \log Z,\quad
S_{\mathrm{config}} = (H_{\mathrm{eff}} - G_{\mathrm{config}}) / T
\]

- \(H_{\mathrm{eff}} = \langle E \rangle\) is a Boltzmann-weighted **effective energy**, not automatically calorimetric ΔH.
- Numerics use log-sum-exp for stability.
- Single-temperature \(C_v\) is **not** experimental binding ΔCp.

## Support classification (short)

| Feature | Class |
|---------|-------|
| ThermodynamicBreakdown (G, H, S, Cv, σ) | Core analysis |
| Component averages | Core diagnostic when data present |
| Compensation / enthalpy–entropy fractions | Diagnostic only — never rank |
| Joint receptor–ligand ensemble | Experimental |
| Affinity / Kd converters | Experimental (safe when calibrated) |
| Temperature scan + ΔCp fit | Experimental / model-derived |

## Audit schema

`make_total_sampled_output` / `ThermodynamicOutputDC` enforce consistency gates used by entropy.help packaging. Validated by `python/tests/test_thermo_schema.py` and C++ `test_thermo_ledger` / `test_statmech`.

## Language checklist

| Prefer | Avoid |
|--------|--------|
| CF / contact-function scoring proxy | “true ΔG from CF alone” |
| ensemble-derived free energy estimate | “experimental binding free energy” without ITC calibration |
| thermodynamic ledger (F, H, −TS, Cv) | conflating score rank with thermo rank |
