# TASK 1 — Land CleftDetector deterministic probe ordering
Role: software engineer (implementation). Model: **Sonnet 5** (mechanical + test wiring; low scientific-judgment content).
Repo: /Users/lp.more/Projects/FlexAIDdS   Branch off: main (currently 356d77b0c)

## Problem (established by OPS, do not re-litigate)
`LIB/CleftDetector.cpp::generate_probes()` fills a per-thread `local` probe vector under
`#pragma omp for schedule(dynamic,64) nowait`, then merges into the shared `probes` vector
under `#pragma omp critical` in THREAD-ARRIVAL order. Downstream `cluster_probes()`
(single-linkage) is order-sensitive, so the cleft grid (`.rrg`) is NOT reproducible across
thread counts. OPS measurement (FLEXAID_SEED=12345, 4 threads, 2 runs):
  baseline .rrg: DIFFERS on 4/5 targets (1G9V,1GM8,1GPK,1HP0), identical only 1HNN.

## Fix (already proven for reproducibility — implement cleanly)
Sort `probes` by a canonical geometric key (center x,y,z then radius) immediately before
`return probes;`. OPS verified: 5/5 targets byte-identical .rrg after the sort, and grid-point
COUNT is unchanged (1027/15224/3972/1003) — sort changes order/tie-break, not cleft volume.
Reference patch: artifact cleftdetector_deterministic_order.patch (version 5e3c3d7b).

## Deliverables
1. Apply the sort in generate_probes().
2. Add a ctest determinism test: dock a fixed seed on ≥3 targets at OMP_NUM_THREADS=1 and =4,
   assert the emitted .rrg is byte-identical across both thread counts AND run-to-run.
3. Keep the change behavior-neutral in serial (must not change 1-thread output vs its own prior
   1-thread output beyond the deterministic reordering).

## DO NOT
- Do not merge to main. Push a branch `fix/cleft-deterministic-order` and stop.
- Landing is gated on the OPS Astex-85 autonomous A/B (accuracy). OPS runs it, not you.

## OPS acceptance gate (I verify before merge)
- ctest determinism test green.
- Astex-85 autonomous top-1 RMSD (spyrmsd, 2 Å): sorted engine within noise of baseline
  (no target flips success→fail attributable to the sort).
