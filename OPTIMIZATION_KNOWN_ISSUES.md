# Optimization branch — known issues

## FLEXAIDDS_PARALLEL_REPRODUCE (Opt1) — partially fixed, still GATED OFF

The GA offspring-CF-eval parallelization remains **OFF by default** and must stay off
until the two remaining issues below are resolved.

### Fixed (branch opt1/stale-status-fix, commit 9b514ca82)
The deferred path left offspring `status` untouched, assuming it was `' '` ("needs eval").
But `chrom[num_chrom+i]` is REUSED memory from the previous generation (status typically
`'n'`), and `calculate_fitness`'s eval loop SKIPS `status=='n'` — so deferred offspring kept
the prior occupant's stale CF. This produced wrong results even single-threaded
(CF -5.23 vs serial -51.93). Fix: explicitly set `status=' '` on deferred offspring.

**Verified after fix (1G9V, FLEXAID_SEED=12345):**
- parallel @1 thread is now run-to-run DETERMINISTIC (10/10 identical) and close to serial
  (-33.61 vs -36.07). The stale-CF defect is gone.

### Remaining (why the flag is still OFF)
1. **parallel @1 thread is not yet BIT-identical to serial** (-33.61 vs -36.07). The deferred
   batch-eval changes the order in which offspring reach `calculate_fitness`'s `QuickSort`,
   altering tie-ordering and thus which chromosomes survive BOOM selection. Deterministic,
   but not equal to the serial reference — so "bit-identical speedup" is not yet proven.
2. **parallel @4 threads is non-reproducible** (run1 -60.9 vs run2 -45.8). A shared-state
   write inside the `eval_chromosome`/`vcfunction`/Vcontacts hot path (NOT
   `ccbm_inject_strain`, which is multi-model-only and read-only on FA here) races under
   `schedule(dynamic)`. Needs the CF-eval shared-write audit (the same Vcontacts static
   buffers the perf swarm flagged) before the flag can be enabled.

**Bottom line:** the headline GA speedup is NOT yet realized. One real bug fixed; two
remaining before FLEXAIDDS_PARALLEL_REPRODUCE can be turned on.

All other merged optimizations (contacts memset->epoch, hoist rigid index, precompute
PoseBust vdW radius) are verified bit-identical to main with default flags (10/10 parity).
