# Validated Capabilities

**Last actualized**: 2026-07-12

This file lists the capability surface that the repository is willing to present as validated for **Core 1.0**.

A capability belongs here only if it is both **documented** and **exercised by automation or release validation**.

## Core execution surface

The following are the validated execution surfaces for Core 1.0:

- command-line docking workflows through `FlexAIDdS`
- legacy-compatible command-line workflows through `FlexAID`
- command-line vibrational entropy workflows through `tENCoM`
- Python package workflows through `flexaidds`
- documented JSON configuration workflows that map to supported CLI and Python use cases

## Repository-level guarantees

- supported-platform installation documentation (`docs/INSTALLATION.md`)
- explicit support matrix (`docs/SUPPORT_MATRIX.md`)
- explicit security policy and buffer-safety regressions (`tests/test_buffer_safety.cpp`)
- explicit reproducibility policy for benchmark claims (`docs/REPRODUCIBILITY.md`)
- explicit known-limitations documentation (`docs/KNOWN_LIMITATIONS.md`)
- testing guide and coverage inventory (`docs/TESTING.md`, `docs/TEST_COVERAGE_ANALYSIS.md`)

## Thermodynamic capabilities

The following thermodynamic features have completed implementation + unit testing and are considered validated for **reporting/analysis** use (see `docs/thermodynamics.md`):

| Capability | Primary automation |
|------------|--------------------|
| Canonical ensemble via `StatMechEngine` | `tests/test_statmech.cpp`, Python thermo tests |
| Thermodynamic ledger (`ThermodynamicBreakdown`) | `tests/test_thermo_ledger.cpp`, `python/tests/test_thermo_breakdown.py` |
| Additive corrections (G_vib, G_natural) with flags | thermo ledger tests |
| Component-wise Boltzmann averages | C++/Python thermo suites |
| Pure-Python parity for core quantities | `test_py_statmech`, `test_thermodynamics*` |
| JSON/CSV round-tripping of results | `test_models*`, `test_results*` |
| A1.1 audit schema identities | `python/tests/test_thermo_schema.py` |

These are safe for use in analysis, visualization, and publication when properly labelled (ensemble-derived; CF is a scoring proxy).

## Benchmark-facing guarantees

Only benchmark-facing outputs backed by repository artifacts should be treated as validated from the repository itself.

That means:

- the benchmark has a reproducibility bundle under `benchmarks/`
- commands and expected outputs are documented
- the metric computation path is discoverable in the repository
- PoseBusters + RMSD criteria apply when the benchmarking skill claims docking success (see `.agents/skills/flexaidds-benchmarking/SKILL.md`)

## GA and scoring automation (supporting Core)

Validated *as software contracts* (not as a universal accuracy claim on all targets):

| Area | Primary tests |
|------|----------------|
| GA population / operators / validation | `test_ga_*`, `test_gaboom`, `test_production_blockers` |
| BindingMode clustering + I/O | `test_binding_mode_*` |
| MOL2/SDF readers | `test_mol2_sdf_reader` |
| Hardware dispatch reporting | `test_hardware_dispatch`, `test_unified_dispatch` |
| Dataset adapter affinity normalization | `python/tests/test_dataset_adapters.py` |
| Ranking-bias diagnostics (post-hoc CF) | `python/tests/test_ranking_bias.py` |

## Exclusions

Anything not listed here should be treated as either planned, provisional, or experimental.

See also:

- `docs/EXPERIMENTAL_CAPABILITIES.md`
- `docs/KNOWN_LIMITATIONS.md`
- `docs/SUPPORT_MATRIX.md`
