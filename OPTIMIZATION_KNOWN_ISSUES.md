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
2. **parallel @>1 thread is non-reproducible — narrowed to ≈0.2% of chromosomes.**
   Gen-0 population CF checksum (order-independent sum over all chromosomes), 1G9V,
   FLEXAID_SEED=12345:
   - flag OFF, 4 threads: −302083.242271 — **reproducible**
   - flag ON, 1 thread: −302096.627183 — **reproducible** (stale-status fix)
   - flag ON, 4 threads: −302130.80 / −302302.87 — **DIVERGES**

   So at the population-eval level the divergence is **flag-ON specific** (a real Opt1
   property), correcting an earlier note that called it purely pre-existing. Per-chromosome
   decomposition: **998/1000 chromosomes carry identical CF; 2/1000 genuinely differ**
   (one with a large `wal` swing — a rare clash/penalty branch). Offspring genes are created
   **serially** in `reproduce()`, so the population contents are deterministic; only
   `calculate_fitness` is parallel. The residual is therefore a **rare, branch/boundary-
   sensitive numeric nondeterminism** in the multi-thread batch-eval copy path.

   **ThreadSanitizer (Debug, Metal OFF) result:** races appear ONLY in
   `CleftDetector.cpp:87,129` (startup probe-merge); the GA CF-eval hot path
   (`gaboom.cpp`/`vcfunction.cpp`/`Vcontacts.cpp`/`ic2cf.cpp`) is **race-clean** — so the
   2-chromosome divergence is an ordering / rare-uninitialised-read subtlety in the per-thread
   `tl_fa`/`tl_vc`/`tl_atoms` copies, not a coarse race. Pinning the exact field/branch is a
   multi-hour instrumented bisection with diminishing certainty — the real remaining task.
   (A separate, milder full-dock 4-thread divergence exists even with the flag OFF, arising
   later in the GA; that one is pre-existing.)

   Ruled out this pass (all reverted): Voronoi degeneracy jitter reseed (the `edgenum>=200`
   failsafe fires 0x on 1G9V), `schedule(static)` (did not restore reproducibility),
   CleftDetector deterministic probe sort (helps cleft determinism but does not fix full-dock
   4-thread reproducibility, and changes cleft output so it needs its own Astex-85 A/B).

   Full write-up: `OPT1_RACE_INVESTIGATION.md` (artifact).

### Latent (separate) concurrency bugs found while investigating — not Opt1-specific
- **CleftDetector probe-merge race** (`CleftDetector.cpp:87,129`): TSan-confirmed; per-thread
  `local` vectors merged in thread-arrival order under `omp critical`, and `cluster_probes()`
  is order-sensitive single-linkage. Candidate fix (`cleftdetector_deterministic_order.patch`,
  not committed): canonical geometric sort of probes; changes cleft output, needs Astex-85 A/B.
- **FlexDEE shared linked list** (`ic2cf.cpp:461,526`): when `FA->psFlexDEENode != NULL`
  (flexible-sidechain / DEE targets), the parallel eval walks and **splices/grows** a shared
  linked list (`FlexDEE_Nodes++`, `->next`/`->prev` writes) via each thread's `tl_fa` pointer
  copy that still aliases shared nodes — a genuine race. Inactive on 1G9V (DEE off), so not the
  1G9V cause, but it WILL corrupt flexible-sidechain docks under Opt1. Must be per-thread or
  serialized before Opt1 can be enabled for flexible-residue runs.

**Bottom line:** the headline GA speedup is NOT yet realized. The stale-CF bug is fixed and makes
Opt1 **1-thread** deterministic; **>1-thread** still has a narrow (≈0.2% of chromosomes) numeric
nondeterminism in the batch-eval copy path, plus the latent FlexDEE race for flexible-residue
targets. Keep `FLEXAIDDS_PARALLEL_REPRODUCE` OFF until both are resolved.

All other merged optimizations (contacts memset->epoch, hoist rigid index, precompute
PoseBust vdW radius) are verified bit-identical to main with default flags (10/10 parity).
