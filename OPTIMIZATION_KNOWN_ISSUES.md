# Optimization branch — known issues

## FLEXAIDDS_PARALLEL_REPRODUCE (Opt1, commit df8c87feb) — DEFECTIVE, DO NOT ENABLE

The GA offspring-CF-eval parallelization is **gated OFF by default** and must stay off.

**Measured defect (1G9V, FLEXAID_SEED=12345):**
- Serial / flag OFF: CF -51.9293 (bit-identical to main, 10/10 poses).
- Flag ON @1 thread: CF -5.2330, 0/10 poses match serial — WRONG result even without
  any thread concurrency, so the patch altered algorithm semantics, not just parallelism.
- Flag ON @4 threads: non-reproducible (run1 -5.2330, run2 -12.1093).

**Status:** the offspring loop does not preserve serial GA semantics (RNG draw order and/or
population-write ordering). Needs a correctness fix (per-offspring deterministic RNG stream
keyed by offspring index, ordered reduction) before the flag can be enabled. Until then the
headline GA speedup is NOT available.

All other branch optimizations (contacts memset->epoch, hoist rigid index, precompute
PoseBust vdW radius) are verified bit-identical to main with default flags.
