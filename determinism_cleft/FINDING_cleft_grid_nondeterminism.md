# The Astex benchmark variance is in cleft detection, not the GA

**Date:** 2026-08-09
**Engine:** `build_tests/FlexAIDdS` (working tree at `c5395779` + uncommitted audit instrumentation)
**Status:** root cause reproduced on 2 targets; a working fix already exists in-tree, gated OFF

---

## Summary

Run-to-run variance in multi-threaded Astex campaigns originates in **SURFNET cleft
detection**, upstream of the genetic algorithm. `CleftDetector::generate_probes` merges
per-thread probe lists in *thread-arrival* order; that order permutes the cleft grid, and
the cleft grid index **is** GA gene 0. The same seed therefore addresses different 3D
anchors from run to run.

The GA is not the source. Held at a fixed thread count with the grid proven identical, the
GA is bit-reproducible.

Two consequences:

1. **Run-to-run nondeterminism** appears only when the machine is loaded enough to perturb
   thread arrival order. `R>1` restarts supply exactly that load, which is why the defect
   hides on a quiet single-restart probe and surfaces in an 84-target campaign.
2. **Results at different `--omp-threads` are not comparable at all** — the binding-site
   grid itself differs, so `omp=1` and `omp=4` are searching differently-indexed spaces.

---

## Measurements

All runs: 1GPK / 1G9V from `benchmarks/astex_diverse/astex_diverse/`, oracle mode
(`FLEXAIDDS_ORACLE_SITE`), `FLEXAID_SEED=12345`, 200 chromosomes × 30 generations,
`FLEXAIDDS_PARALLEL_RESTARTS=0`. Grid = md5 of the emitted `.rrg`. No grid cache was in
play — every run computed and saved its own grid (verified in the logs).

### The probe set is stable; only its order is not

`CleftDetector: 540 gap-spheres survived shrinking` in every 1GPK run, at every thread
count. The probe **multiset** is thread-invariant; the merge **order** is not. The
per-cluster sphere dump (`FLEXAIDDS_CLEFT_DUMP`) differs between two 4-thread runs, and
union-find cluster labels permute (6/40 in one run, 14/16 in the next) while cluster
membership, centroids and `ligandable_score` values are unchanged.

### Thread count changes the grid (phenomenon A)

| threads | replicates | grid md5 | poses md5 |
|---|---|---|---|
| 4 | ×2 | `2ca980fc` | `37852730` |
| 1 | ×2 | `aa1a261f` | `ca0860fe` |

Reproducible within a thread count, **different across thread counts**. This is the
mechanism behind the reported `omp=1 → 5.41 Å`, `omp=2 → 8.96 Å`, `omp=3 → 9.92 Å`.

### Load makes it diverge run-to-run (phenomenon B)

Six identical 4-thread runs, gate OFF, with 8 background CPU hogs:

| target | distinct grids | breakdown |
|---|---|---|
| 1GPK | **2 of 6** | 4× `2ca980fc`, 2× `5c769048` |
| 1G9V | **4 of 6** | 3× `9a9d74e8`, and `1c3ee1dc`, `2fe3e261`, `4275ed04` |

Same binary, same seed, same inputs, same thread count. Without the background load, the
1GPK 4-thread pair agreed — which is precisely why quiet-box probes look clean.

### The existing gate fixes both (`FLEXAIDDS_CLEFT_SORT=1`)

| condition | distinct grids | distinct pose sets |
|---|---|---|
| 1GPK, 4 threads under load, ×6 | **1** | **1** |
| 1GPK, {1,4} threads × 2 reps | **1** (`4d55826e`) | **1** (`cbed58cc`) |
| 1G9V, 4 threads under load, ×3 | **1** | — |

With the gate on, 1-thread and 4-thread runs produce **byte-identical grids and poses**.
That is thread-count invariance, not merely reproducibility.

---

## Mechanism

1. `LIB/CleftDetector.cpp:87-130` — `#pragma omp for schedule(dynamic,64) nowait` fills a
   thread-private `std::vector<Probe> local`; `#pragma omp critical` appends each thread's
   block to the shared `probes`. Block order = thread completion order.
2. `LIB/CleftDetector.cpp:388` — the sphere linked list is built by walking `probes[]` and
   **prepending**, so list order is the reverse of the merged probe order.
3. `LIB/generate_grid.cpp:33-81` — `generate_grid` walks that list and appends each newly
   seen lattice point at `cleftgrid[FA->num_grd++]`. The point **set** is order-independent;
   the **index assignment** is not.
4. GA gene 0 indexes `cleftgrid`. Permuting the index assignment silently redefines what
   every gene value means — the same seed explores a different space.

Cluster membership is permutation-invariant (union-find connected components), so this is
purely an ordering effect. It does not change *which* pocket is found in oracle mode; it
changes the **addressing** of the search space.

## Why the existing determinism gates do not help

Both live in the GA, downstream of the defect:

- `FLEXAID_DETERMINISTIC` (`gaboom.cpp:3199-3220`) computes `eval_threads` but only the
  pragma at `gaboom.cpp:3221` carries `num_threads(eval_threads)`; loops at 3338, 3449 and
  4022 are untouched. Either way it cannot affect a grid built before the GA starts.
- `FLEXAIDDS_PARALLEL_REPRODUCE` (`gaboom.cpp:2152-2170`) pins the CF-eval reduction path
  only.

## Refuted — do not spend time here

- **Static scheduling on the four GA pragmas** (the original fix request). The loop bodies
  write per-index (`chrom[ii].evalue`, `pshare_out[pi]`), so the result is order-independent.
  Confirmed empirically: with the grid held constant, the GA is bit-reproducible at a fixed
  thread count.
- **Nested OpenMP across restarts.** Restarts are separate OS processes
  (`DatasetRunner.cpp:6437` `fork_exec`); there is no nesting. `R>1` matters as a source of
  CPU *contention*, not nesting.
- **Restart-timeout truncation of the pose pool.** The `remaining` budget at
  `DatasetRunner.cpp:6476-6479` was misread; the drain does not silently SIGKILL restarts.
  Differing pose counts follow from the differing grid.
- **Unstable MIF energy sort** (`MIFGrid.h:117-122`) as an independent cause.

## Open

- `FLEXAIDDS_CLEFT_SORT=1` changed the 1GPK pose set even at 1 thread
  (`ca0860fe` → `cbed58cc`), consistent with `CLEFTSORT_AB_VERDICT.md` (1G9V 5.13 → 7.54 Å).
  The sort reorders tie-breaks in serial too, so **enabling it invalidates prior
  single-thread numbers**. That is a protocol decision, not a bug fix.
- The stated rationale in `CLEFTSORT_AB_VERDICT.md` §4.1 for keeping it off — "the canonical
  protocol runs single-thread, where cleft detection is already deterministic" — does not
  cover the campaigns actually being run at `omp=3`/`omp=4`.
- Untested here: whether blind/autonomous mode (no oracle pre-filter, many clusters) can also
  flip *which* cleft is elected, not just the addressing. Cluster-size ties break on
  union-find labels, which are probe indices.

## Reproduce

    determinism_cleft/repro_cleft_determinism.sh 1GPK

Prints distinct grid/pose checksum counts for gate OFF vs ON, under load, at 1 and 4 threads.
