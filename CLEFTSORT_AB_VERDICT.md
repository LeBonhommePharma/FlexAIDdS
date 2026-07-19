# CleftDetector Probe-Sort — A/B Verdict

## Status: reproducibility PROVEN; accuracy A/B inconclusive (infrastructure); **LANDED flag-gated, OFF by default** (`FLEXAIDDS_CLEFT_SORT`)

## IMPLEMENTED
The sort is wired behind `FLEXAIDDS_CLEFT_SORT` (OFF by default) in `LIB/CleftDetector.cpp`.
Validated on 1G9V (engine md5 `891fcc25f25c996dc6a640d52f3ea08b`):
- **Gate OFF, 1-thread grid == true-baseline (`7f1f10a0`) grid** — byte-identical (md5 `37aec8fe...`).
  Default path reproduces prior behavior exactly; all prior single-thread numbers stay valid.
- **Gate OFF, 4-thread:** DIFFERS run-to-run (baseline nondeterminism, unchanged).
- **Gate ON, 4-thread:** IDENTICAL run-to-run (deterministic cleft grid).
- **Volume preserved:** OFF and ON produce the same grid-point count (939), differing only in order.

## 1. What the sort does
`CleftDetector::generate_probes` merges per-thread probe lists under `#pragma omp critical`
in **thread-arrival order**, then `cluster_probes` (single-linkage) is order-sensitive. Under
`--omp-threads > 1` this makes the SURFNET cleft grid (`.rrg`) **non-deterministic** run-to-run.
The candidate patch adds a canonical geometric sort of probes before return, fixing the order.

## 2. Reproducibility result — PROVEN (prior measurement, 5 targets, 4-thread)
- **Baseline engine:** `.rrg` cleft grid DIFFERS across 4-thread runs on **4/5** targets
  (1G9V, 1GM8, 1GPK, 1HP0 divergent; 1HNN identical). 1G9V byte-375 diff, md5 df14ee3e vs 469aaa7f.
- **Sorted engine:** `.rrg` grid **IDENTICAL** across 4-thread runs on **5/5** targets.
- **Grid-point count unchanged** (1027/15224/3972/1003 identical baseline vs sorted): the sort
  changes probe *order / tie-break*, not cleft *volume*. Volume-preserving.

## 3. Accuracy A/B — INCONCLUSIVE (infrastructure)
Attempted 5 targets × 2 engines, autonomous/blind, 2000-gen, 1-thread. The parallel run hit
machine oversubscription (load ~68) compounded by an orphaned process from a prior run that the
first cleanup missed and that could not be reaped from the analysis sandbox (`ps`/`kill` EPERM).
Only **1/5 pairs completed** (1G9V) before the run stalled.

**1G9V (rank-0 elected pose, symmetry-corrected heavy-atom RMSD, 1-thread):**
- baseline: **5.13 Å** (fail @ 2.0 Å)
- sorted:   **7.54 Å** (fail @ 2.0 Å)
- No success flip (both miss — 1G9V is a hard blind target: pocket localization, not scoring).
- **But the sort changed the pose** (5.13 → 7.54 Å) even at 1 thread — the probe reorder is
  NOT inert; it perturbs the blind search trajectory.

One hard-target datapoint where both arms fail is too noisy to conclude regression or safety.

## 4. Decision rationale — why NOT default-merge
1. **The sort targets a mode the canonical protocol does not use.** The reproducible benchmark
   protocol runs **single-thread per restart** (`--omp-threads 1`, `FLEXAIDDS_RESTARTS=10`,
   `FLEXAID_SEED`), where cleft detection is **already deterministic**. The sort only matters for
   `--omp-threads > 1`.
2. **The sort changes single-thread results** (1G9V: 5.13 → 7.54 Å) because it reorders probe
   tie-breaks even in serial execution. Merging it as default would **change / invalidate every
   prior single-thread benchmark number** — to fix a non-canonical mode.
3. This trades a known, controlled baseline for an unvalidated one. Not acceptable as a default.

## 5. Recommendation
**Land the sort behind an OFF-by-default env gate** (e.g. `FLEXAIDDS_CLEFT_SORT=1`), consistent
with the established CF-fix pattern (all fixes opt-in, defaults reproduce prior behavior):
- Default OFF → canonical single-thread results **unchanged**, all prior numbers stay valid.
- Set ON → deterministic multi-thread cleft detection for anyone running `--omp-threads > 1`.
- Revisit default-ON only after a **clean-machine, cognate-site A/B** (known pocket → most
  targets succeed → a success flip is real signal, not blind-search noise) shows no regression
  across the full 85 × 10-restart protocol.

## 6. Provenance
- Baseline engine md5 `7f1f10a0f10b682b33a76622a40f1a60`; sorted md5 `b2dd05e347e5e859cd818632027942d1`.
- Candidate patch: `cleftdetector_deterministic_order.patch` (artifact 80eceadd), NOT committed.
- 1G9V scoring: `cleftsort_ab_1G9V.csv` (this run).
