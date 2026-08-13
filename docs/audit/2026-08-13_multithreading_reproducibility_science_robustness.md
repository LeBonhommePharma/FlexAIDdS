# Deep code audit — multithreading, reproducibility & science robustness

**Date:** 2026-08-13
**Scope:** `LIB/` C++ engine (GA, scoring, statistical mechanics, entropy stacks, clustering,
parallel orchestration), the CMake build flags, and the Python/DatasetRunner reproducibility
contract.
**Base commit:** `fe25507e` (detached from `main`).
**Method:** first-hand source inspection with real-execution grounding — a full test build
(`clang-18`, `-DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release`) with `ctest` **89/89 passing**, the
Python suite **1253 passed / 74 skipped / 1 environmental failure**, and a standalone compiled
reproduction of the headline RNG defect (see F1).

This document is descriptive/diagnostic. It changes no engine behaviour, ranking, or
thermodynamics. Per `AGENTS.md`/`METHODOLOGY.md`, any fix that could move a score, a pose, or a
ledger value must land behind an opt-in flag that defaults OFF and pass the §1 parity gate; the
recommendations below are written with that in mind.

---

## 1. Executive summary

The engine's *published* reproducibility model is fundamentally sound and, in most respects,
carefully engineered:

- **Restarts are process-isolated**, each seeded deterministically from `(target_id, restart,
  seed_base)` (`LIB/DatasetRunner.cpp:177`). This is the strongest determinism guarantee in the
  system and it sidesteps every in-process threading hazard for the canonical Astex path.
- **GA reproduction (crossover/mutation/selection) runs serially** on the master thread; only
  fitness evaluation is parallelised, and each chromosome's contact-function (CF) score is computed
  independently on per-thread scratch, so CF values are thread-count-invariant.
- **Cleft-grid probe generation already merges deterministically** (ascending-`ii` buckets,
  `LIB/CleftDetector.cpp:85-137`), closing a previously-known multi-thread nondeterminism.
- **Numerics are defensively written** where it matters most: log-sum-exp with max-shift
  (`LIB/statmech.cpp`, `LIB/UnifiedHardwareDispatch.cpp`), `p>0` guards before `p·log p`
  (`LIB/BindingMode.cpp`), a relative eigenvalue floor before `kT/λ` in tENCoM
  (`LIB/tENCoM/tencm.cpp:1020-1031`), and a fail-closed scientific-claim firewall
  (`LIB/statmech.cpp:126-157`).

Against that backdrop, the audit found **one high-severity latent RNG-correctness bug**, several
medium-severity reproducibility/robustness gaps, and a set of low-severity hardening items. None of
them invalidate the published Astex-85 self-docking numbers (which do not exercise the affected
feature paths), but they do affect: the ring-flexibility feature, the experimental grid-decomposed
`ParallelDock` mode, the thermodynamic ledger's bit-reproducibility across thread counts, and the
engine's resilience to a non-finite score.

### Findings at a glance

| ID | Severity | Area | One-line |
|----|----------|------|----------|
| **F1** | **High** | RNG correctness / reproducibility | `lazy_thread_rng` multiplexes one generator across all "streams" and re-seeds on every stream switch → interleaved streams on one thread produce repeating draws (reproduced). |
| **F2** | Medium | Reproducibility (MT) | Voronoi hull-failure jitter draws RNG *inside* the `schedule(dynamic)` eval loop; value depends on thread scheduling/count and (via F1) resets the GA's main stream. |
| **F3** | Medium | Reproducibility (MT) | Boltzmann/log-sum-exp OpenMP reductions sum in thread-count-dependent order → F/S/Cv/weights not bit-identical across thread counts for large ensembles. |
| **F4** | Medium | Reproducibility (FP) | Engine built with `-ffast-math` (reassociation + FMA contraction) plus runtime SIMD-width dispatch → not bit-reproducible across binaries/architectures. |
| **F5** | Medium | Science robustness | `QuickSort` comparator has no finite-check; a NaN/Inf CF sorts unpredictably and can be elected rank-0. |
| **F6** | Low–Med | Reproducibility (MT) | `ParallelDock::run_region` ignores its `rng_seed` argument; region RNG keyed on `omp_get_thread_num()` under `schedule(dynamic,1)` → grid-decomposed mode non-reproducible run-to-run. |
| **F7** | Low | Science robustness | `Cv = var/(kT²)` uses an unclamped variance (can be slightly negative); all-zero multiplicities yield NaN moments. |
| **F8** | Low | Multithreading | `UnifiedHardwareDispatch::detect()` guarded by a non-atomic `detected_` flag; a data race if the first dispatch call ever happens concurrently. |
| **F9** | Low | Test coverage | The only `lazy_thread_rng` test repeats a *single* stream, so it cannot catch F1; it gives false confidence in stream independence. |

An automated breadth sweep (three independent read-only passes) corroborated F1–F9 and added
**F10–F24** — see §8. The highest-impact additions are **F10** (Critical, but scoped to the
experimental `ParallelDock` concurrent-GA mode: a shared function-`static` eval workspace and shallow
region copies) and **F13** (an inverted mutual-information sign in the experimental joint-ensemble
API). Positive/no-action-needed observations are collected in §4.

---

## 2. What was inspected & verified

- **Parallelism surfaces:** every `#pragma omp` region in `gaboom.cpp` (fitness eval at
  `:3252`, `:4054`; PSHARE/SMFREE fitness at `:3369`, `:3480`), the histogram/reduction kernels in
  `UnifiedHardwareDispatch.cpp` and `statmech.cpp`, `ParallelDock.cpp`, `CleftDetector.cpp`, and the
  `static thread_local` scratch in `Vcontacts.cpp`.
- **RNG:** `LIB/RngSeed.h` and all call sites (`RandomDouble`, `Vcontacts.cpp:416`,
  `ChiralCenterGene.cpp:22`, `SugarPucker.cpp:119`, `RingConformerLibrary.cpp:69/74`,
  `FOPTICS.cpp:15`, `shannon_ga.cpp:32`) plus the DatasetRunner per-restart seeding.
- **Numerics:** `statmech.cpp` (partition function, moments, Cv, joint ensemble, affinity),
  `UnifiedHardwareDispatch.cpp` (LSE/Boltzmann/Shannon dispatch), `BindingMode.cpp` (clustering,
  conformational/mutual-information entropy), `tENCoM/tencm.cpp` (vibrational modes).
- **Build:** `CMakeLists.txt` and `cmake/FlexAIDOptions.cmake` FP/optimisation flags.
- **Execution:** clean configure + build; `ctest` 89/89; Python `pytest` 1253 passed; standalone
  RNG reproduction compiled against `LIB/RngSeed.h`.

---

## 3. Findings

### F1 — `lazy_thread_rng` collapses all logical streams into one generator (High)

**Location:** `LIB/RngSeed.h:118-131`; live trigger `LIB/gaboom.cpp:122-130`.

`lazy_thread_rng(stream)` is written as if each `stream` constant had its own generator, but it keeps
a *single* `thread_local mt19937` plus a cached stream id and **re-seeds the generator whenever the
requested `stream` differs from the cached one**:

```118:131:LIB/RngSeed.h
inline std::mt19937& lazy_thread_rng(std::uint64_t stream)
{
    thread_local std::uint64_t cached_stream = ~0ULL;
    thread_local std::uint64_t cached_epoch  = ~0ULL;
    thread_local std::mt19937 rng = make_thread_rng(stream);

    const std::uint64_t epoch = g_seed_epoch.load(std::memory_order_acquire);
    if (cached_stream != stream || cached_epoch != epoch) {
        rng = make_thread_rng(stream);
        cached_stream = stream;
        cached_epoch = epoch;
    }
    return rng;
}
```

Because `make_thread_rng(stream)` is deterministic (same seed for a given stream/thread/epoch), each
switch back to a stream **restarts that stream from its seed**. Two call sites with different stream
ids, interleaved on the same thread, therefore never advance — they each keep returning their first
draw.

At least seven distinct stream constants exist (`0x9A800D` GA-main, `0x0C0A11` Voronoi jitter,
`0xC417A1` chiral, `0x5A6A9` sugar pucker, `0x516`/`0x515` ring library, `0xF0701C5` FOPTICS). The
**live** interleave is in `ring_mutate_chrom`, which alternates the GA-main stream and the
sugar-pucker stream within a single function on a single thread:

```122:130:LIB/gaboom.cpp
	for (int i = 0; i < ns && i < MAX_RING_FLEX; ++i)
		if (RandomDouble() < pucker_mut_prob)                 // stream 0x9A800D
			c->ring_phases[i] = sugar_pucker::mutate_phase(c->ring_phases[i]);   // stream 0x5A6A9 -> re-seed
	for (int i = 0; i < n6 && i < MAX_RING_FLEX; ++i)
		if (RandomDouble() < ring_mut_prob)                   // stream 0x9A800D -> re-seed back to its SEED
			c->ring_six[i] = static_cast<uint8_t>(RandomDouble() * lib.n_six());
```

**Reproduction (compiled against the real header, run on this VM):**

```
baseline  A stream : 0.773906 0.047997 0.376984 0.445347 0.912045
interleaved A draws: a0=0.773906 a1=0.773906 a2=0.773906
interleaved B draws: b0=0.409877 b1=0.409877
RESULT: A_stuck=YES(bug)  B_stuck=YES(bug)  A_never_advances=YES(bug)
```

The baseline advances normally; the moment a second stream is touched between draws, both streams
freeze at their first value.

**Impact.**
- **Sugar-pucker / ring-flex docking:** every pucker mutation draws the *same* phase, and the GA's
  main RNG is reset mid-generation, collapsing the diversity of standard-gene mutation for any ligand
  with a flexible furanose ring. This is a scientific-correctness defect for that feature (the GA is
  far less stochastic than intended and biased toward a fixed pucker).
- **Cross-thread reproducibility (with F2):** when the Voronoi jitter (`0x0C0A11`) fires on the
  master thread inside the eval region, the next serial `RandomDouble()` (`0x9A800D`) re-seeds the
  GA-main stream to its start. Whether that happens depends on which thread scored the degenerate
  chromosome — i.e. on `schedule(dynamic)` timing — so it becomes a source of run-to-run and
  thread-count divergence.
- **Latent:** `chiral::ChiralCenterGene::mutate/crossover` and
  `RingConformerLibrary::random_six/five_index` are defined but *not* currently called from the live
  GA path (verified: no callers under `LIB/**/*.cpp`), so their streams are dormant — but they are
  landmines for anyone who wires them in.

**Why the published Astex-85 numbers are unaffected:** standard self-docking without ring flex only
ever touches `0x9A800D` on the master thread (plus the rare jitter path), so no *sustained*
interleave occurs. The bug is real but its blast radius is the ring-flex feature and MT determinism,
not the headline benchmark.

**Recommendation.** Give each `stream` its own generator instead of multiplexing one. E.g. a
`thread_local std::unordered_map<uint64_t, std::mt19937>` (or a small fixed array keyed by a stream
enum) seeded lazily per stream, re-seeded per stream on epoch change. This is behaviour-preserving
for any path that only ever uses one stream per thread (the published path), so it can pass §1 parity
with the fix ON, while fixing ring-flex sampling and removing the cross-boundary reset. Add a
regression test that interleaves two streams and asserts each advances independently (see F9).

---

### F2 — Voronoi hull-failure jitter uses a per-thread RNG inside the dynamic eval loop (Medium)

**Location:** `LIB/Vcontacts.cpp:405-436` (drawn from `lazy_thread_rng(0x0C0A11)`), reached from the
`schedule(dynamic)` fitness loops at `LIB/gaboom.cpp:3252` and `:4054`.

When a Voronoi cell fails to converge (≥200 hull edges — a degenerate/clashing geometry), the code
perturbs the atom coordinates by ±0.005 Å and recomputes:

```414:419:LIB/Vcontacts.cpp
				// perturb atom coordinates (deterministic when ga.seed/FLEXAID_SEED set)
				thread_local std::uniform_real_distribution<float> vc_dist(-0.005f, 0.005f);
				auto& vc_rng = flexaids_rng::lazy_thread_rng(0x0C0A11ULL);
				VC->Calc[atomzero].atom->coor[0] += vc_dist(vc_rng);
				VC->Calc[atomzero].atom->coor[1] += vc_dist(vc_rng);
				VC->Calc[atomzero].atom->coor[2] += vc_dist(vc_rng);
```

The comment claims determinism, but the draw happens on whichever worker thread the chromosome was
scheduled to (`schedule(dynamic)`), and the per-thread stream is salted by `omp_get_thread_num()`
(`RngSeed.h:68-76`). So on the degenerate path the perturbation — and therefore the resulting CF and
coordinates for that pose — depends on the thread that scored it, which is timing-dependent. This
violates the byte-identical-across-thread-counts expectation of `METHODOLOGY.md §2`, and it is the
concrete mechanism that (via F1) can reset the GA-main stream when it fires on the master thread.

**Impact.** Only the degenerate-hull path is affected, but such hulls are common early in a GA run
when random poses clash, so this can fire in practice. Effect: non-reproducible CF for the perturbed
pose across thread counts / runs; amplified by F1.

**Recommendation.** Make the jitter independent of the scheduling thread: derive the perturbation
from a stream keyed on stable pose identity (e.g. chromosome index and atom index) rather than the
thread-salted generator, or apply a fixed deterministic offset. Guard behind the parity gate.

---

### F3 — Thread-count-dependent floating-point reductions in the dispatch layer (Medium)

**Location:** `LIB/UnifiedHardwareDispatch.cpp:582-595` (`lse_openmp`), `:724-748`
(`boltzmann_openmp`), `:751-826` (`boltzmann_avx512` with the OpenMP branch).

These kernels accumulate a partition-function sum with an OpenMP reduction:

```588:591:LIB/UnifiedHardwareDispatch.cpp
    #pragma omp parallel for reduction(+:sum) schedule(static)
    for (int i = 0; i < n; ++i)
        sum += std::exp(values[i] - x_max);
    return x_max + std::log(sum);
```

Floating-point addition is not associative, and the order in which per-thread partial sums are
combined by an OpenMP `reduction` is unspecified and varies with thread count. So `log_Z`, and every
quantity derived from it (Helmholtz F, entropy S, heat capacity Cv, Boltzmann weights), is **not
bit-identical across thread counts** for ensembles large enough to take these paths
(`N ≥ 4096` for the OpenMP branch; the AVX-512 branch always chunks per thread). The same file also
selects among scalar/Eigen/AVX-512/Metal/OpenMP paths *by runtime detection and by size thresholds*
(`N ≥ 16/256/1024/4096`), so a single machine can switch reduction/rounding strategy at a
data-dependent boundary.

**Impact.** The GA's pose *ranking* is CF-based and computed elsewhere, so the *elected pose* is not
affected by this. The affected outputs are the **thermodynamic ledger values** (F, S, Cv, weights)
reported in receipts/REMARKs; two runs at different `--omp-threads` can print ledger numbers that
differ in the low-order digits. Given docking ensembles are usually small (retained poses 10–50,
population 1000), the `N≥4096` OpenMP paths rarely trigger in the canonical flow, which is why this
is Medium rather than High.

**Recommendation.** For the thermodynamic layer, prefer a deterministic reduction: accumulate into
per-thread partials and combine them in a fixed index order (as `boltzmann_pmf` and the Shannon
kernels already do with per-thread histogram buckets), or use a compensated/tree reduction with a
fixed shape. Document that the thermo ledger is only bit-stable at a fixed backend + thread count
otherwise.

---

### F4 — `-ffast-math` + runtime SIMD dispatch breaks cross-binary bit-reproducibility (Medium)

**Location:** `CMakeLists.txt:286, 567, 609, 669, 784, 869, 992, 1040, 1065, 1093, 1150`;
`cmake/FlexAIDOptions.cmake:70-115, 247-270`.

Every engine target is compiled with:

```286:286:CMakeLists.txt
    target_compile_options(FlexAID PRIVATE -Wall -O3 -ffast-math -fno-finite-math-only)
```

`-ffast-math` enables FP reassociation (`-fassociative-math`), FMA contraction, reciprocal
approximations, and `-fno-signed-zeros`. Combined with the runtime SIMD-width dispatch (scalar vs
AVX2 vs AVX-512 vs NEON tails, all with different associativity) and the native-arch tuning — which is **on by default for the shipped binary**: `BUILD_FLEXAIDDS_FAST`
defaults **ON** (`CMakeLists.txt:821`) and adds `-flto` + `-march=native`/`-mcpu=native` + `-DNDEBUG`
to the `FlexAIDdS` target on any host (`:869-873`), and `FLEXAIDS_MCPU_NATIVE` defaults **ON** for
`flexaid_core` (`FlexAIDOptions.cmake:81`) — the exact CF score and its tie-breaks can differ across compilers,
optimisation levels, SIMD paths, and machines. This is consistent with `REPRODUCIBILITY.md`'s
"agree within floating-point rounding" caveat, but it is worth stating precisely because CF ties
decide pose order.

**Positive note (deliberately correct):** the targets append `-fno-finite-math-only`, which overrides
the `-ffinite-math-only` implied by `-ffast-math` and **keeps `isnan`/`isinf` working**. That is the
right call and it is what makes an F5-style finite-check meaningful rather than optimised away. Order
matters and the flags are appended in the correct order.

**Impact.** Same-binary/same-machine/same-thread-count runs remain deterministic (the compiler's
reassociation is fixed at compile time). The concern is cross-binary and cross-architecture
reproducibility, and SIMD-width-dependent tie-breaks in ranking.

**Recommendation.** No behaviour change proposed. Consider (a) documenting that the parity gate is
per-toolchain, and (b) for the *thermodynamic* translation units only, evaluating
`-ffp-contract=off` (or a non-fast-math subset) so the audited ledger identities in
`thermo_invariants.md` hold to IEEE precision. Keep `-fno-finite-math-only` everywhere.

---

### F5 — Ranking sort has no NaN/finite guard on the CF score (Medium)

**Location:** `LIB/gaboom.cpp:4693-4740` (`QuickSort`) with comparator macros
`LIB/gaboom.h:43-45`.

The GA orders chromosomes by `evalue` with:

```43:45:LIB/gaboom.h
#define QS_TYPE double
#define QS_ASC(a,b) ((a)-(b))
#define QS_DSC(a,b) ((b)-(a))
```

If any `evalue` is NaN, both scan predicates `(NaN - piv) <= 0` and `(NaN - piv) > 0` evaluate
false, so a NaN key partitions arbitrarily. The bounds checks (`l<=r`) prevent an out-of-bounds
crash, but the resulting order is undefined — a NaN-scoring pose can land at rank 0 and be **elected
as the top pose**. `-fno-finite-math-only` (F4) means NaNs are *not* assumed away, so this path is
reachable if scoring ever produces a non-finite CF (degenerate geometry, a division/`log`/`sqrt`
domain slip in the contact function, an Inf strain term).

**Impact.** Silent mis-ranking / election of an invalid pose, with no diagnostic. Low probability but
high consequence when it happens.

**Recommendation.** Before ranking (or when writing `evalue`), map non-finite CF to a defined
worst-case sentinel (e.g. `+HUGE_VAL` for the ascending "lower is better" convention) and emit a
one-line warning. This is a pure robustness guard with no effect on finite scores, so it is
parity-safe.

---

### F6 — `ParallelDock::run_region` ignores its seed argument; region RNG is thread-assignment-dependent (Low–Medium)

**Location:** `LIB/ParallelDock.cpp:114-131, 152-250`.

The grid-decomposed docking mode computes a per-region seed and passes it to `run_region`:

```117:122:LIB/ParallelDock.cpp
    for (int r = 0; r < n_regions; r++) {
        unsigned int seed;
        #pragma omp critical
        { seed = seed_gen() + r; }

        results_[r] = run_region(regions_[r], seed, target);
```

but `run_region(const GridRegion&, unsigned int rng_seed, ...)` never applies `rng_seed` — it runs
`GA(&ws.fa, ...)` on a shallow copy whose `ws.gb.seed` is inherited unchanged from the parent for
every region. Meanwhile the GA's per-thread RNG is salted by `omp_get_thread_num()` and the loop is
`schedule(dynamic, 1)`, so which region gets which stream depends on thread arrival → the
grid-decomposed results are **non-reproducible run-to-run**. Additionally, each worker's
`ga()` calls the *global* `flexaids_rng::set_master_seed()` (`gaboom.cpp:361`), so concurrent region
GAs race on the global seed/epoch atomics (same value, so numerically benign, but it forces
cross-thread epoch re-seeds).

**Impact.** Confined to the experimental octree grid-decomposition path (not the published restart
path). Effect: non-deterministic pose search in that mode and a dead determinism knob.

**Recommendation.** Apply `rng_seed` to `ws.gb.seed` in `run_region` and drive each region GA from a
per-region master seed rather than the process-global `set_master_seed`, or run regions as
subprocesses like the DatasetRunner restart loop. Guard behind parity.

---

### F7 — StatMech edge cases: unclamped Cv variance and all-zero multiplicities (Low)

**Location:** `LIB/statmech.cpp:300-310` and `:216-233`.

`std_energy` is clamped to `≥0` but the variance fed to heat capacity is not:

```299:310:LIB/statmech.cpp
    double kT  = kB_kcal * T_;
    double var = E2_avg - E_avg * E_avg;
    ...
    th.heat_capacity  = var / (kB_kcal * T_ * T_);
    ...
    th.std_energy     = std::sqrt(std::max(0.0, var));
```

Floating-point cancellation in `E2_avg − E_avg²` can make `var` slightly negative for a
near-degenerate ensemble, yielding a tiny negative `Cv`. Separately, `compute()` builds log-weights
as `counts.log() − β·energies`; if all multiplicities are 0 (or negative), `lnZ = −∞` and the
moments become `exp(−∞ − (−∞)) = exp(NaN) = NaN`, producing NaN free energy/entropy rather than a
clean error. `add_sample` does not validate multiplicity.

**Impact.** Minor: a small negative Cv is cosmetic; the all-zero-multiplicity case is an unusual input
but would emit NaN thermodynamics instead of failing closed.

**Recommendation.** Clamp `var` with `std::max(0.0, var)` before dividing for Cv, and validate
`multiplicity > 0` in `add_sample` (or fail closed in `compute()` when `Σ counts ≤ 0`).

---

### F8 — Unsynchronised lazy hardware detection in the dispatch singleton (Low)

**Location:** `LIB/UnifiedHardwareDispatch.cpp:71-77` and every `if (!detected_) detect();` entry
(`:531, :623, :653, :873, :906, :1034, :1107`).

`detect()` is guarded by a plain `bool detected_` (`check-then-act`) with no atomicity or lock. The
`instance()` function-local `static` is thread-safe to construct, but if the *first* dispatch call
ever occurs concurrently from two threads, both can run `detect()` and race on `info_`/`detected_`.

**Impact.** In the current flow the dispatch layer is first touched serially (post-GA thermodynamics),
so the race is not currently triggered; `detect()` is also idempotent (writes the same values), so
worst case is redundant work / a torn read of an `info_` field. Latent, not active.

**Recommendation.** Call `detect()` once explicitly during startup, or make `detected_` a
`std::atomic<bool>` / use `std::call_once`.

---

### F9 — The `lazy_thread_rng` test cannot catch F1 (Low)

**Location:** `tests/test_ga_context.cpp:104-107`.

The only test exercising the RNG draws the **same** stream three times:

```104:107:tests/test_ga_context.cpp
    const double a = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
    ...
    const double b = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
    const double c = dist(flexaids_rng::lazy_thread_rng(0x9A800DULL));
```

Single-stream use is exactly the case F1 does *not* break, so this test passes while the defect is
live. It gives false confidence in stream independence.

**Recommendation.** Add a test that interleaves two distinct streams and asserts each advances
independently and reproducibly (the reproduction in F1 is a ready-made assertion).

---

## 4. Things that are correct (no action needed)

- **Deterministic cleft-probe merge.** `CleftDetector::generate_probes` fills per-`ii` buckets in
  parallel and concatenates them in ascending `ii` order, reproducing the serial order bit-exactly at
  any thread count (`LIB/CleftDetector.cpp:85-137`). Cluster tally/selection iterates `std::map`/
  `std::set` in key order (`:257-317`). Deterministic.
- **GA fitness eval scratch.** The `schedule(dynamic)` eval loops write only `chrom[ii]` and
  per-thread `tl_*[tid]`/`p_*[tid]` scratch, reset per iteration; no shared writes
  (`gaboom.cpp:3260-3318, 4060-4112`). CF values are thread-count-invariant.
- **Log-sum-exp everywhere it matters.** `statmech.cpp:203-207`, `BindingMode.cpp:1127-1132`, and the
  `E_min`-shifted Boltzmann batch (`UnifiedHardwareDispatch.cpp:828-844`) all avoid overflow/underflow.
- **Entropy guards.** Every `p·log p` sum guards `p>0` (`BindingMode.cpp:1188, 1217, 1252, 1296,
  1311`, each preceded by an `if (p > 0.0)`; `statmech.cpp:975-978` `safe_entropy`).
- **tENCoM eigenvalues.** Negative/near-zero modes are dropped via a *relative* floor
  `max(1e-8, λ_max·1e-9)` before `kT/λ` (`tencm.cpp:1020-1031`, guard also at `:906`).
- **Scientific-claim firewall.** `StatMechEngine` fails closed to proxy-only provenance unless a
  calibrated-energy + physical-measure + matched-reference witness with a real artifact SHA256 is
  present (`statmech.cpp:53-157`), and `merge`/`merge_samples` downgrade on witness mismatch
  (`:778-813`). This is exactly the "separate the scoring proxy from thermodynamics" guardrail.
- **Affinity/temperature guards.** `T ≤ 0` throws in the engine ctor and in the affinity converters
  (`statmech.cpp:191-192, 1084-1106`); the joint-ensemble path sorts map ids for determinism
  (`:1027-1044`).
- **`Vcontacts` scratch** uses `static thread_local` (`:24, :1817`) and once-only `static const`
  env caches (`:1659, :2050`) — thread-safe.
- **`SharedPosePool` is fully synchronised** — all mutating/reading methods take `mtx_`, move-assign
  locks both mutexes in address order to avoid deadlock, and `publish` **rejects non-finite energies**
  (`SharedPosePool.cpp:56`) with a comment noting it uses `DBL_MAX` rather than `infinity`
  "(UB under `-ffast-math`/`-ffinite-math-only`)" (`:106`). This is precisely the finite-check pattern
  that F5's `QuickSort` path lacks and should adopt.
- **`TargetServer` / grand partition function concurrency is correct.** `create_session` uses an
  atomic id counter; `register_result` protects `knowledge_` under `knowledge_mtx_`
  (`TargetServer.cpp:60-68`) and delegates to `GrandPartitionFunction`, which has its own
  `mutable std::mutex mtx_` (`GrandPartitionFunction.h:192`). No unsynchronised shared writes found.
- **Published restart determinism.** Per-restart process isolation with distinct deterministic seeds
  (`DatasetRunner.cpp:177-181`, restart subprocesses at `:5844-5860`) is the model reproducibility
  pattern here.

---

## 5. Coverage vs `METHODOLOGY.md §1–§2`

- **§1 parity gate** (single-thread, `FLEXAID_SEED`, byte-identical elected poses) does its job for
  the *default* flag configuration, and would **not** flag F1/F2 because it runs single-thread and
  without ring flex.
- **§2 multi-thread determinism** (cleft `.rrg` byte-identical across 1 vs 4 threads) is satisfied by
  the deterministic probe merge (§4), **but** the GA-population clause ("all 10 elected poses
  byte-identical across two 4-thread runs") is exactly what F2 (and F1 via F2) can violate on the
  degenerate-hull path, and what F1 violates outright once ring flex is enabled. The gate should be
  run with ring flex ON and with a target that forces Voronoi hull failures to exercise these.

---

## 6. Prioritised recommendations

1. **F1 (High):** rework `lazy_thread_rng` to keep an independent generator per stream; add the
   interleave regression test (F9). Parity-safe because single-stream paths are unchanged.
2. **F5 (Medium):** clamp non-finite CF to a worst-case sentinel before `QuickSort`; warn once.
   Parity-safe (no effect on finite scores).
3. **F2 (Medium):** make the Voronoi jitter independent of the scheduling thread (key on pose/atom
   identity). Behind the parity gate.
4. **F3 (Medium):** deterministic fixed-order reduction for the thermodynamic partition sums; or
   document thermo-ledger bit-stability as fixed-backend/fixed-thread only.
5. **F6 (Low–Med):** apply the per-region seed in `ParallelDock::run_region` (or subprocess-isolate).
6. **F7/F8 (Low):** clamp Cv variance, validate multiplicities, and make hardware detection
   `call_once`/atomic.
7. **F4 (Low, doc):** evaluate `-ffp-contract=off` for the thermodynamic TUs; document the
   per-toolchain nature of the parity gate. Keep `-fno-finite-math-only`.

---

## 7. Verification performed in this audit

| Check | Result |
|-------|--------|
| Configure + build (`clang-18`, Release, `BUILD_TESTING=ON`, GPU off) | Success, no errors |
| `ctest --output-on-failure` | **89/89 passed** (4.9 s) |
| `pytest python/tests/` | **1253 passed, 74 skipped, 1 failed** |
| The 1 Python failure | Environmental only — Ubuntu `/bin/sh`→`dash` symlink breaks `test_bare_name_on_PATH_still_resolves`'s `.name=="sh"` assertion; not a code defect |
| F1 reproduction (compiled vs `LIB/RngSeed.h`) | Confirmed: interleaved streams freeze at their first draw |

No engine sources were modified; this audit is analysis only.

---

## 8. Addendum — automated breadth sweep (corroboration + new findings)

Three read-only exploration passes independently swept the tree for the same three axes. They
**corroborated F1–F9** (all three flagged the thread-salted RNG + Vcontacts jitter, the OpenMP
partition-sum reductions, `-ffast-math`, and the NaN-in-ranking hazard), and surfaced additional
concrete issues. The items below were **verified first-hand** against the cited lines before being
recorded here.

### F10 — Concurrent-GA (`ParallelDock`) path has structural data races (Critical for that mode)

The single-GA OpenMP eval is well-engineered, but `ParallelDock` runs *many* full `GA()` instances
concurrently (`ParallelDock.cpp:115` `#pragma omp parallel for`), and several process-global pieces
of GA state are not isolated per instance:

- **`static ParEvalWS ws`** (`gaboom.cpp:3129`) is a function-static evaluation workspace, explicitly
  documented as "resident across generations … GA parallelism lives strictly inside this loop." Under
  concurrent `GA()` calls each instance sees a different `ws.fa`/`ws.atoms`, fails the `ws_valid`
  check, and executes `ws = ParEvalWS{}` (full reallocation, `:3140-3142`) while another thread is
  indexing into it → heap corruption / wrong CF / crash.
- **Shallow FA/VC copy** in `create_workspace` (`ParallelDock.cpp:63-72`) leaves `optres`,
  `contacts`, `Calc`, DEE lists, etc. pointing at shared master arrays.
- **`omp_set_num_threads(2)`** (`gaboom.cpp:218-220`) and the global `srand`/`set_master_seed`
  (`:360-361`) are process-global side effects invoked from every concurrent GA.
- **DEE update is gated on `!omp_in_parallel()`** (`ic2cf.cpp`), so under `ParallelDock`'s outer
  parallel region FlexDEE pruning silently never runs — a different search than the serial dock.
- **`atoms_copy` sized `atm_cnt`** rather than `atm_cnt+1` (`ParallelDock.cpp:69`) against the
  engine's 1-based atom convention (reported off-by-one; medium confidence).

**Impact.** Confined to the experimental octree grid-decomposition mode (not the published restart
path, which is process-isolated). In that mode: crashes, torn scores, non-reproducible search.
**Recommendation.** Make `ParEvalWS` non-static (thread/instance-local), deep-copy the region
workspace, and drive per-region seeds/threads locally — or run regions as subprocesses like the
DatasetRunner restart loop.

### F11 — TOCTOU race in classic clustering (High)

`cluster.cpp:145-167` checks `Clus_GAPOP[i] == -1` *outside* the `#pragma omp critical`, then assigns
inside it; two threads can both pass the check and double-assign the same index, double-decrementing
`n_unclus` and inflating `Clus_FRE`. **Impact:** wrong/non-reproducible cluster sizes when the classic
clustering path runs multi-threaded.

### F12 — Population-init eval omits the `metal_coord` reset (Medium, correctness)

The steady-state eval clears `cf.metal_coord` in its per-thread reset (`gaboom.cpp:3290-3303`), but
the population-init eval's reset block (`:4087-4099`) does **not**. A stale metal-coordination term
can leak into the initial-population CF for metal-containing ligands. **Recommendation:** add the
`metal_coord = 0.0` reset to the populate path for parity with the steady-state path.

### F13 — `compute_joint_ensemble` mutual-information sign is inverted (High, experimental API)

```1065:1066:LIB/statmech.cpp
    // Step 5: Mutual information I(R;L) = S_joint - S_receptor - S_ligand  (in nats, dimensionless after scaling)
    result.mutual_information_dimensionless = S_joint - S_receptor - S_ligand;
```

Mutual information is `I(R;L) = S_R + S_L − S_joint ≥ 0`; the code computes its negation.
`BindingMode.cpp:1256-1257` uses the correct sign, so the two APIs disagree. **Impact:** the
experimental joint-ensemble MI is negated (and would print negative "information"). **Recommendation:**
flip to `S_receptor + S_ligand − S_joint`; add a test asserting `I ≥ 0`.

### F14 — `BindingMode` energy sort has no tie-break → nondeterministic elected mode (Medium)

`EnergyComparator` (`BindingMode.h:250-260`) compares cached energy with a strict `<` and no secondary
key. `BindingModes` is ordered with `std::sort` (`BindingMode.cpp:100`), which is not stable, so when
two modes share a cached free energy the "top" (elected) mode is implementation-/order-dependent.
**Impact:** which pose is emitted rank-0 can flip on ties across std::sort implementations or input
order. **Recommendation:** add a deterministic secondary key (e.g. representative CF then a stable id)
or use `std::stable_sort`.

### F15 — `geometry.cpp` bond angle: unclamped `acos`, unguarded zero denominator (Medium, NaN source)

```127:128:LIB/geometry.cpp
  cosa /= sqrt(absu*absv);
  cosa = (float)(acos(cosa)*180.0/PI);
```

If two bond vectors are (near-)collinear/degenerate, `absu*absv` can be 0 (→ Inf) and FP error can push
`cosa` outside `[-1,1]` (→ `acos` returns NaN). `Vcontacts.cpp:1607-1612` and `hbond_potential.h:40-41`
clamp; this shared helper does not. This is a concrete route to a non-finite geometric term feeding the
CF — i.e. an upstream producer of the NaN that F5's `QuickSort` cannot then order. **Recommendation:**
guard the denominator and `std::clamp(cosa, -1.0, 1.0)` before `acos`.

### F16 — Two β conventions on the same scale; physical `kB` must not touch CF units (Medium)

The GA/cluster/soft-β layers use `β = 1/T` over the unitless CF landscape (`read_input.cpp:253-258`,
deliberately — folding in `kB_kcal` over-sharpens weights ~503× at 300 K), while `StatMechEngine`
defaults to physical `β = 1/(kB·T)` (`statmech.cpp:186-188`) and also exposes `selection_weights()`
at `β = 1/T`. The hazard is calling the physical-`β` `compute()` on CF-unit energies and printing the
result as a kcal/mol ledger. This is a labelling/consistency risk rather than an arithmetic bug — the
claim firewall marks such output proxy-only — but the two "temperatures" (`thermo_T_eff = 0.596`
CF-units vs Kelvin `TEMPER`) are easy to cross. **Recommendation:** assert at the call site which β a
given ensemble was built for, and keep CF-scored ensembles on `selection_weights()`/`β=1/T`.

### F17 — `boltzmann_pmf` and other thermo entry points: empty-bin and `T≤0` gaps (Medium)

`boltzmann_pmf` (the deprecated-name "WHAM", `statmech.cpp:626-732`) assigns empty bins `F=0` (so an
unvisited bin looks like a ground state), never uses `f_old` in the estimator (the "self-consistency"
iteration only measures convergence), and does not guard `T≤0` (`:640`). `helmholtz` (`:573-576`),
`init_replicas` (`:591-595`) and `compute_joint_ensemble` (`:994`) likewise compute `β=1/(kB·T)`
without a `T>0` check, unlike the engine ctor. **Impact:** misleading PMF on sparse coordinates; Inf β
/ NaN thermodynamics at `T=0`. **Recommendation:** guard `T>0` uniformly; mark empty PMF bins as
undefined rather than 0.

### F18 — Grand-canonical sums are hash-order- and overflow-sensitive (Medium)

`GrandPartitionFunction` accumulates `log_Xi` and per-ligand sums by iterating `std::unordered_map`
(`GrandPartitionFunction.cpp:130-136, 287-298`); non-associative `sum += exp(...)` over hash order
makes `log_Xi`/`p_bind` values (not the final `dG`-sorted rank) depend on insertion/hash order.
Separately, `MultiSiteGPF.cpp:167-174` leaves log-space to form `exp(log_xi)−1`, which overflows to
Inf once `log_Xi ≳ 700` (strong binding). **Impact:** competitive-binding ledger noise; Inf multi-site
cooperativity. The single-ligand canonical Astex path is unaffected (GPF unused). **Recommendation:**
accumulate in a fixed key order and keep the `Ξ−1` correction in log-space (`log1p(-exp(-log_xi))`
form).

### F19 — Python result aggregation depends on filesystem order (Medium)

`python/flexaidds/docking.py:742-747` discovers pose PDBs ordered by `st_mtime`, and
`python/flexaidds/figures.py:195` consumes an **unsorted** `rglob`. Which pose is treated as primary,
or which audit JSON "wins", therefore depends on filesystem timestamps/enumeration rather than
rank/name. **Recommendation:** sort by rank/name explicitly at these sites (the C++ batch/`results.py`
paths already sort).

### F20 — PoseBusters temp path collides under same-process concurrency (Medium, claim gate)

`LIB/PoseBust/ChecksChemistry.cpp:483-486` builds an InChI temp file from `getpid()`, so two PoseBust
checks in the *same* process (concurrent docks) can clobber one SDF and flip the PB pass/fail — which
is a benchmark **claim gate**, not the CF search. **Recommendation:** include a thread id / unique
counter (or `mkstemp`) in the temp name.

### F21 — tENCoM/ENCoM absolute vibrational entropy is model-scale, and `tencom_diff` abs-folds negative modes (Low, scientific caveat)

Absolute `S_vib` in `ShannonThermoStack.cpp:332-359` combines a model-scale `ω = √λ` with SI `ħ`, so
the absolute magnitude is heuristic (~0.06 kcal/mol/K per mode offset); only *differences* cancel the
scale. `tENCoM/tencom_diff.cpp:25-27` sets `frequency = √|λ|`, folding indefinite/negative Hessian
eigenvalues into real frequencies rather than filtering them. There is also a floor inconsistency:
`sample()` uses an absolute `λ<1e-8` cut (`tencm.cpp:906`) while `bfactors()` uses a relative floor
(`:1020-1021`). **Recommendation:** treat absolute `S_vib` as proxy unless calibrated; filter negative
λ (don't `abs`); unify the soft-mode floor.

### F22 — FOPTICS clustering draws the (F1-affected) thread RNG (Low)

`FOPTICS.cpp:15,846` uses `lazy_thread_rng`/`RandomDouble` for random splits, so cluster membership —
and therefore the OPTICS representative selected as the elected pose (`BindingMode.cpp:640-650`) — is
subject to the F1 multiplexing defect and the thread-count salt when run under OpenMP.

### F23 — `simd_distance.h` Boltzmann helper lacks the max-shift (Low, latent)

`simd_distance.h:485-486` computes `exp(-β·E)` with no `E_min` shift (unlike the production
`compute_boltzmann_batch`). Harmless where currently used, but a latent overflow/underflow footgun if
reused on raw CF energies.

### F24 — Default (unseeded) master seed is wall-clock `time(0)` (Low, expected)

With neither `ga.seed` nor `FLEXAID_SEED` set, the GA seeds from `time(0)` (`gaboom.cpp:349-361`).
Expected behaviour, but worth stating: reproducibility requires an explicit seed, and the effective
seed is recorded in the pose REMARK. (`srand` there is vestigial — no `rand()` consumers exist in
`LIB/`.)

### Corroboration & scope notes

- All three sweeps independently confirmed there are **no unsynchronised shared writes in the
  single-GA OpenMP eval** (per-`tid` scratch + `default(none)`), and that `SharedPosePool`,
  `TargetServer`, `GrandPartitionFunction`, `InStreamClustering`, and the thread pools are
  mutex/atomic-correct — matching §4.
- The net reproducibility contract is: **bit-identical docking holds only under `FLEXAID_SEED` +
  `OMP_NUM_THREADS=1` + `FLEXAIDDS_PARALLEL_REPRODUCE` off + identical binary/flags/host.**
  `FLEXAID_DETERMINISTIC` pins only the main CF eval loop to one thread; it does not cover the niche
  loops, clustering, FOPTICS RNG, the thermo OpenMP reductions, or host math flags.
- Prioritisation update: **F10 (ParallelDock) and F13 (MI sign)** join **F1/F5** as the items with the
  clearest correctness impact; F10 is Critical but scoped to the experimental grid-decomposition mode.

