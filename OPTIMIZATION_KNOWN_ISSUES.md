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
   - flag OFF, 4 threads: −302083.242271 — **reproducible** ← SCOPE-LIMITED, see
     "Default path (flag OFF) is NOT reproducible on the flexible arm" below.
     This measurement stands as taken (1G9V, gen-0, rigid), but it does **not**
     generalise: flag-OFF multi-thread diverges on the flexible-sidechain arm.
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
- **FlexDEE shared linked list** (`ic2cf.cpp:456–535`): initially flagged as a latent race, but
  on closer read it is **already mitigated** — the whole mutation block is guarded by
  `FA->useflexdee > 0 && rclash && !omp_in_parallel()` (ic2cf.cpp:416), so it is **skipped under
  Opt1's parallel eval**. Further, the DEE-list *consumer* in `gaboom.cpp:1889` is commented out
  (dead code) and the parallel eval path never calls `check_clash`/consumes the list, so the
  linked list is effectively write-only in the current engine. **No live race; no fix needed.**
  (`check_clash.cpp:54` uses only the scalar `FA->dee_clash` threshold, not the list.) If DEE
  pruning is ever re-enabled inside the parallel region, the list must be made per-thread first.

**Bottom line:** the headline GA speedup is NOT yet realized. The stale-CF bug is fixed and makes
Opt1 **1-thread** deterministic; **>1-thread** still has a narrow (≈0.2% of chromosomes) numeric
nondeterminism in the batch-eval copy path, plus the latent FlexDEE race for flexible-residue
targets. Keep `FLEXAIDDS_PARALLEL_REPRODUCE` OFF until both are resolved.

All other merged optimizations (contacts memset->epoch, hoist rigid index, precompute
PoseBust vdW radius) are verified bit-identical to main with default flags (10/10 parity).

> ⚠️ **Read "with default flags" literally.** For `FLEXAIDDS_CONTACTS_EPOCH` the default is OFF,
> so the 10/10 parity run exercised the legacy memset path and said **nothing** about the epoch
> path. With the flag ON the optimization was in fact **incorrect**: the epoch counter lived in
> `FA_Global`, which the threaded GA re-snapshots from the master every generation
> (`gaboom.cpp`, `tl_fa[t] = *FA;`) while the stamp buffer it points at stays resident across
> generations. The counter therefore rewound at every generation boundary against stale
> high-water stamps, contacts were silently skipped, and CF came out wrong.
>
> Fixed by moving the counter **inside** the buffer it stamps (`CONTACTS_EPOCH_SLOT`, see
> `flexaid.h`), so a struct copy can no longer separate the two. Regression coverage:
> `tests/test_contacts_epoch.cpp` plus the paired `WILL_FAIL` target that re-creates the
> pre-fix layout. The flag remains **default OFF** pending the ON-vs-OFF parity evidence
> required by `METHODOLOGY.md` §1.
>
> **Lesson for this table: a parity run with a flag OFF is not evidence about the flag ON.**
> Any future entry here must state which arm was actually exercised.

---

## Default path (flag OFF) is NOT reproducible on the flexible arm

*Added 2026-09-07. Everything above this line concerns `FLEXAIDDS_PARALLEL_REPRODUCE`
(Opt1), which is gated OFF. This section concerns the **default** configuration — Opt1
unset — and is a separate defect with a separate cause.*

**Why this section exists.** The `FLEXAID_DETERMINISTIC` comment in `gaboom.cpp` used to
describe the default multi-thread path as accepting "the ~0.2% chromosome-level numeric
drift documented in OPTIMIZATION_KNOWN_ISSUES.md." That citation was misattributed: the
0.2% figure above is called **flag-ON specific** in its own section, and the same table
measures flag-OFF at 4 threads as reproducible. The number never described the default
path. The source comment has been corrected.

**What is measured.** 1P2Y, `defined-cleft-redock`, flexible sidechains, gen=1000,
`FLEXAID_SEED=12345`, `FLEXAIDDS_PARALLEL_REPRODUCE` unset:

| omp threads | per-restart pose-set jaccard | verdict |
|---|---|---|
| 3 | 0.976 / 0.000 / 0.138 | **NOT reproducible** |
| 1 | 1.000 / 1.000 / 1.000 | reproducible |

Best CF across four runs of identical inputs: **−254.86 / −230.86 / −231.03 / −231.41**.
Initial poses were byte-identical in every run, so setup and seeding are deterministic and
the divergence begins inside the search. This is a **different GA trajectory**, not a
rounding difference — restart 1 shares *zero* pose-set members between runs.

**Cause.** The Voronoi non-convergence failsafe (`Vcontacts.cpp`, `voronoi_poly2`,
`edgenum >= 200`) perturbs an atom's coordinates by ±0.005 Å drawn from a **thread_local**
RNG. The CF-eval loop is `schedule(dynamic)`, so which thread — and therefore which draw —
a chromosome receives is not reproducible run-to-run. `recalc == 1` on the GA path, so the
hull is recomputed **with** the perturbed coordinate and the returned CF depends on that
draw. Restarts share one process and one RNG stream, so a single differing cluster in
restart 0 amplifies to complete disagreement by restart 1.

**Target specificity** tracks the failsafe firing rate, not flexibility: 1P2Y fires 3,575
times against a **median of 42** across the 84-target set — a 15× outlier with no target
in the gap.

**Three further defects at the same site**, measured with `FLEXAIDDS_VCT_FAILSAFE_DIAG`
(18,771 firings across 10,762 calls in one 1P2Y run):

1. `origcoor` is re-saved from the already-perturbed coordinates on each pass, so a call
   taking N passes restores as of pass N and leaves **N−1 perturbations applied**. Measured
   maximum: 1,730 passes → 1,729 surviving perturbations (RMS 0.21 Å/axis). 99.665% of
   calls take one pass, but 26 calls take ≥10 and produce 41% of all firings.
2. `origcoor` is zeroed at the top of `RESTART`, so a pass reaching the restore without
   re-entering the failsafe writes **(0,0,0)** and teleports the atom to the origin —
   **79 occurrences** in that run.
3. The `return -1` exit restores nothing at all.

**None of this self-heals**, because the GA's per-chromosome reset is *selective*:
`gaboom.cpp` builds `dirty_atm` from `mov[]` and `map_par[].atm` — "atom indices modified
by **ic2cf**" — while this failsafe modifies whichever atom's hull failed. Measured
intersection on 1P2Y: 11 ligand atoms (90001–90011) are covered; **receptor atoms 743 and
2293 are not**. `use_selective` is therefore an optimisation resting on an unstated
invariant that another translation unit violates.

**Impact is not noise.** With the coordinate guard on, 1P2Y goes from **rmsd 4.1597 Å
(failure)** to **0.8939 Å (success)**, while best CF moves from −254.86 to −231.47. The
corrupted-geometry search was finding poses that score *better* against a structure with
displaced atoms — a better number for a worse answer, which is why no score-distribution
check ever flagged it. Independently, the keyed-jitter arm lands at −231.31: two different
fixes converging on the same basin.

**Fixes, all gated, defaults preserved and proven inert by pose-hash identity:**

| flag | effect | default |
|---|---|---|
| `FLEXAIDDS_VCT_COORD_GUARD` | arm-once RAII guard; restores the pristine coordinate at scope exit, so it wins on every exit path including `return -1` | OFF |
| `FLEXAIDDS_VORONOI_KEYED_JITTER` | derives the displacement from atom identity + seed instead of a thread_local stream | OFF |
| `FLEXAIDDS_SELECTIVE_INVARIANT_CHECK` | per-generation check that no atom outside `dirty_atm` has moved in `tl_atoms` | OFF |
| `FLEXAIDDS_VCT_FAILSAFE_DIAG` | logs each firing with its pass number, branch and displacement | OFF |
| `FLEXAIDDS_DIRTY_SET_DUMP` | one-shot dump of the dirty set (both index and atom number) | OFF |

**`FLEXAID_DETERMINISTIC` suppresses the divergence but does not fix it** — it removes the
scheduling variable by serialising evaluation (`eval_threads = 1`), paying for determinism
in throughput. Measured cost on 1P2Y: +6.3% wall. That target has `n_flex_bonds = 1`;
on a target where OMP actually scales the cost is the whole parallel speedup. The
coordinate guard addresses the cause and keeps the threads.

**Open at time of writing:** whether the leak on receptor atoms 743/2293 materialises in
`tl_atoms` across chromosomes (the invariant checker's full-budget run decides it), and
the per-target effect of the guard across all 84 targets. The guard has been measured on
**one** target — the 15× outlier — so its default stays OFF until that table exists.
