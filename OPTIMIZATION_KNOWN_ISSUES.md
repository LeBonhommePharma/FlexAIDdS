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
   batch-eval uses the per-thread `tl_fa[0]` copy and defers `ring_load_chrom_to_fa`, a
   different but valid GA trajectory. Deterministic, but not equal to the serial reference —
   so "bit-identical speedup" is not yet proven.
2. **parallel @>1 thread is non-reproducible — but this is PRE-EXISTING, not an Opt1 defect.**
   Verified 2026-07-18 (1G9V, FLEXAID_SEED=12345): with the flag **OFF** (pure main behaviour),
   a 4-thread dock already gives -30.23 vs -40.55 run-to-run. Opt1 only routes more work
   through the parallel `calculate_fitness` path that main already ships and that is already
   non-deterministic across threads.

   **ThreadSanitizer (Debug, Metal OFF) result:** races appear ONLY in
   `CleftDetector.cpp:87,129` (startup probe-merge); the GA CF-eval hot path
   (`gaboom.cpp`/`vcfunction.cpp`/`Vcontacts.cpp`/`ic2cf.cpp`) is **race-clean**. The 4-thread
   non-determinism is therefore **floating-point reduction/accumulation order** under
   `schedule(dynamic)`, not a memory race. True >1-thread bit-reproducibility needs an
   ordered/deterministic FP reduction (per-thread partials combined in fixed thread-index
   order) in the parallel eval — the real remaining engineering task.

   Ruled out this pass (all reverted): Voronoi degeneracy jitter reseed (the `edgenum>=200`
   failsafe fires 0x on 1G9V), `schedule(static)` (did not restore reproducibility),
   CleftDetector deterministic probe sort (helps cleft determinism but does not fix full-dock
   4-thread reproducibility, and changes cleft output so it needs its own Astex-85 A/B).

   Full write-up: `OPT1_RACE_INVESTIGATION.md` (artifact).

**Bottom line:** the headline GA speedup is NOT yet realized. The stale-CF bug is fixed;
1-thread bit-parity and a deterministic FP reduction for >1 thread remain before
FLEXAIDDS_PARALLEL_REPRODUCE can be turned on. The >1-thread non-determinism is a
pre-existing engine property, exposed — not caused — by Opt1.

All other merged optimizations (contacts memset->epoch, hoist rigid index, precompute
PoseBust vdW radius) are verified bit-identical to main with default flags (10/10 parity).
