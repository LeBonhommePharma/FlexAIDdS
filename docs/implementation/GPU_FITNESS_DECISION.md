# GPU FITNESS_EVAL decision (Chunk 6)

**Date:** 2026-08-13  
**Status:** FITNESS_EVAL stays on CPU. Not a default flip.

## Decision

`UnifiedHardwareDispatch::best_backend(KernelType::FITNESS_EVAL)` continues to
return `select_cpu_backend()`, including when a GPU override is set.

## Why (measured-gate, not a ceiling claim)

Chunks 1–5 of the acceleration stack (rigid fast path, two-stage screen, RNG
stream fix, cheap LUT/LSE, niche/ca_rec/pProp) have **unit** coverage only in
this worktree. There is no Astex-85 A/B of GPU vs CPU fitness, and the existing
GPU kernels still:

1. Decode only the first 3 genes as translation (no `buildcc` torsion geometry).
2. Use a drift-tolerant C0 contact-area model instead of Voronoi CF.
3. Are not a CI-run fitness path.

Until those close **and** an ON-arm A/B is measured, `best_backend(FITNESS_EVAL)`
must stay CPU. Shannon/cavity kernels may still use Metal/CUDA.

## How to revisit

Port full `buildcc` gene decode, close channel gaps, FP64 thermo reductions,
then Astex-85 A/B vs CPU. Only then consider flipping FITNESS_EVAL.

Pinned by `tests/test_unified_dispatch.cpp` (`FitnessEvalDefaultsToCpuBackend`,
`GpuOverrideDoesNotApplyToFitnessEval`).
