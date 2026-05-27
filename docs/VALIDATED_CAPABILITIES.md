# Validated Capabilities

This file lists the capability surface that the repository is willing to present as validated for **Core 1.0**.

A capability belongs here only if it is both documented and exercised by automation or release validation.

## Core execution surface

The following are the validated execution surfaces for Core 1.0:

- command-line docking workflows through `FlexAIDdS`
- legacy-compatible command-line workflows through `FlexAID`
- command-line vibrational entropy workflows through `tENCoM`
- Python package workflows through `flexaidds`
- documented JSON configuration workflows that map to supported CLI and Python use cases

## Repository-level guarantees

The following repository-level guarantees are part of the validated surface:

- supported-platform installation documentation
- an explicit support matrix
- an explicit security policy
- an explicit reproducibility policy for benchmark claims
- explicit known-limitations documentation

## Benchmark-facing guarantees

Only benchmark-facing outputs backed by repository artifacts should be treated as validated from the repository itself.

That means:

- the benchmark has a reproducibility bundle under `benchmarks/`
- commands and expected outputs are documented
- the metric computation path is discoverable in the repository

## Exclusions

Anything not listed here should be treated as either planned, provisional, or experimental.

## Thermodynamic Capabilities (post-roadmap)

The following thermodynamic features have completed implementation + unit testing and are considered validated for reporting/analysis use (see `docs/thermodynamics.md` for full details):

- Canonical ensemble quantities via `StatMechEngine` and `ThermodynamicBreakdown` (G_config, H_eff, S_config, Cv, σ_E)
- Additive corrections (G_vib, G_natural) with presence flags
- Component-wise Boltzmann averages (`EnergyComponents` + `component_means`)
- Pure-Python parity for all core quantities
- JSON/CSV round-tripping of thermodynamic results

These are safe for use in analysis, visualization, and publication when properly labelled.
