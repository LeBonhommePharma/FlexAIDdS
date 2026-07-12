# Architecture Overview

**Last actualized**: 2026-07-12

## Execution layers (order)

```text
Inputs (PDB/CIF/MOL2/SDF/SMILES)
        │
        ▼
┌───────────────────┐
│  Site / cavity    │  CleftDetector, CavityDetect (optional / config)
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ Genetic algorithm │  gaboom / GAContext — explores conformational space
│ CF scoring proxy  │  Vcontacts / Voronoi CF (+ optional H-bond, GIST, …)
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ Pose clustering   │  BindingMode / FAST-OPTICS / density peaks
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ StatMechEngine    │  partition function, G/F, H_eff, S, Cv
│ Thermo ledger     │  ThermodynamicBreakdown + optional vib / NATURaL terms
└─────────┬─────────┘
          ▼
┌───────────────────┐
│ Outputs           │  PDBs + REMARKs, JSON, Python loaders, reports, figures
└───────────────────┘
```

## Repository map (simplified)

| Path | Role |
|------|------|
| `LIB/` | C++26 core engine (GA, scoring, thermo, I/O, HW dispatch) |
| `LIB/tENCoM/` | Torsional ENCoM vibrational module |
| `LIB/ShannonThermoStack/` | Shannon entropy + GPU bridges |
| `LIB/LigandRingFlex/`, `ChiralCenter/`, `NATURaL/`, `CavityDetect/` | Specialized modules |
| `python/flexaidds/` | Analysis package, DatasetRunner, schemas, diagnostics |
| `python/bindings/` | pybind11 `_core` |
| `tests/` | GoogleTest + skill Python tests |
| `benchmarks/` | Datasets, runners, reproducibility bundles |
| `.grok/skills/flexaidds/` | Agent skill (validation, data ensure, campaign helpers) |
| `docs/` | Product + MkDocs documentation |

## Ranking vs thermodynamics

- **During search**: individuals are compared primarily by CF (and configured score terms). Changing this requires explicit user request + tests + feature flags (`AGENTS.md`).
- **After ensemble collection**: BindingMode / StatMech produce thermodynamic quantities used for analysis, reporting, and best-mode extraction protocols that sort by free energy **when configured**.

## Hardware dispatch

`UnifiedHardwareDispatch` / Shannon stack select backends at runtime when built:

CUDA → Metal → AVX-512 → AVX2 → OpenMP → scalar

Support level depends on platform and CI matrix — see Support Matrix.

## Python boundary

- Pure Python can load results, compute pure-Python StatMech, run adapters, diagnostics, and DatasetRunner orchestration.
- Full docking still requires the native `FlexAIDdS` (or compatible) binary and runtime data files.
