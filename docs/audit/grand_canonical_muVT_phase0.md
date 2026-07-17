# Phase 0 Audit — Canonical (NVT) vs Grand-Canonical (μVT)

**Branch:** `feature/grand-canonical-muVT-v1`  
**Date:** 2026-07-17  
**Scope:** Document existing Z / Ξ code and the PartitionFunctionBase refactor.

## Inventory of partition-function code

| Component | Path | Ensemble | Role |
|-----------|------|----------|------|
| `StatMechEngine` | `LIB/statmech.h/.cpp` | **NVT** | Canonical Z from GA pose energies; F, S, C_v |
| `log_sum_exp` | `StatMechEngine` → `flexaids::log_sum_exp_dispatch` | — | Shannon-stabilized LSE (AVX-512 / Metal / scalar) |
| `UnifiedHardwareDispatch` | `LIB/UnifiedHardwareDispatch.*` | — | Hardware LSE + entropy kernels |
| `GrandPartitionFunction` | `LIB/GrandPartitionFunction.*` | **μVT** (single site) | Ξ = 1 + Σ z_i Z_i competitive binding |
| `MultiSiteGPF` | `LIB/MultiSiteGPF.*` | **μVT** multi-site | Per-site GPF + cooperativity |
| `TargetServer` | `LIB/TargetServer.*` | owns live Ξ | Session registry + NRGsuite concentrations |
| `GrandCanonicalEngine` | `LIB/GrandCanonicalEngine.*` | **μVT** | Multi-N channel + GPF composition + `PartitionFunctionBase` |
| `PartitionFunctionBase` | `LIB/PartitionFunctionBase.h` | abstract | Shared surface: `log_partition()`, `free_energy()`, T |
| `CanonicalPartitionAdapter` | same header | NVT adapter | Wraps log_Z or StatMechEngine result |
| `BindingMode` | `LIB/BindingMode.*` | NVT per mode | Optional concentration metadata for Ξ feed |
| `VoronoiCFBatch` | `LIB/VoronoiCFBatch.h` | CF proxy | **Does not** fold fugacity into CF |
| Python mirror | `python/flexaidds/grand_canonical.py` | μVT helper | `CompetitiveSite`, `set_concentration` |

## Refactor decisions

1. **Do not fork Ξ math** — multi-species competitive binding remains in
   `GrandPartitionFunction`. `GrandCanonicalEngine` *composes* it.
2. **Add multi-N channel** — classical Σ_N λ^N Z_N for occupancy numbers,
   OpenMP outer sum + Shannon LSE.
3. **Extract `PartitionFunctionBase`** — polymorphic ledger surface for NVT
   adapters and μVT engines without forcing `StatMechEngine` ABI breakage
   (adapter pattern preferred over invasive inheritance).
4. **GPF uses dispatch LSE** — `compute_log_Xi_fresh()` calls
   `flexaids::log_sum_exp_dispatch` (same kernel as NVT Z).

## Test gates (Phase 0/1)

- `tests/test_grand_partition.cpp` — GPF analytic + MOR/naloxone NVT vs μVT
- `tests/test_grand_canonical_engine.cpp` — engine + multi-N + 5-HT2A toy
- `python/tests/test_grand_canonical.py` — pure-Python mirror

## Ranking guardrail

GA / CF ranking paths are **unchanged**. Concentration and fugacity only
affect post-hoc Ξ, p_i, ⟨N⟩, and diagnostic entropy metrics.
