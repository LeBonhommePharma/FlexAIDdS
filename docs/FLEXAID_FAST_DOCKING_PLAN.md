# FlexAID-Fast — Ultra-Optimized Classic FlexAID Docking Plan

**Branch:** `claude/flexaid-optimization-docking-t2h83d`
**Goal:** Take the *classic FlexAID* docking engine (Voronoi complementarity function + genetic algorithm) — behavior kept "almost as is" — and accelerate it to the same level of performance engineering that the rest of FlexAIDdS already carries. **Target: raw docking throughput.**

## Design decisions (locked)

These framing choices drive every priority below:

| Decision | Choice | Consequence |
|---|---|---|
| **Fidelity** | *Small numerical drift allowed everywhere* | May use `float` geometry, reordered parallel reductions, fast/approx transcendentals, and GPU. Ranking and poses must stay essentially equivalent; scores need not match the legacy engine bit-for-bit. |
| **Hardware** | *CPU **and** GPU (CUDA/Metal)* | Do the CPU work **and** bring the existing GPU batch kernels up to full CF fidelity so GPU docking becomes usable, then auto-dispatch per machine. |
| **Scope** | *Keep features, just speed the core* | ΔS entropy machinery (statmech / Shannon / tENCoM / ENCoM) stays available but **off-by-default and zero-cost when unused**; accelerate the shared GA + Voronoi core. |

Because drift is allowed, the historical blocker on `FLEXAIDDS_PARALLEL_REPRODUCE` (a ~0.2% chromosome-level numeric non-determinism under multi-thread offspring eval) is **no longer a gate** — it can be turned on, with a separate deterministic-reduction switch kept only for regression benchmarking.

---

## 1. Where the time actually goes (measured hot path)

```
main()                         LIB/top.cpp:405
 └─ GA(...)                    LIB/top.cpp:2915          scoring fn = ic2cf
     └─ for gen in 0..max_gen  LIB/gaboom.cpp:691        SERIAL generation loop (GA data dependency)
         └─ calculate_fitness  LIB/gaboom.cpp:2307
             └─ #pragma omp parallel for over population   gaboom.cpp:2773-2781   ← only fine-grained parallelism
                 └─ eval_chromosome  gaboom.cpp:3183
                     └─ ic2cf → buildcc (gene→Cartesian) → vcfunction → Vcontacts
```

- **Loop budget:** default `num_chrom = 1000`, `max_generations = 2000` (`top.cpp:2812-2813`); production configs reach ~1750×2000. Each generation scores the whole population → ~10⁹+ per-atom Voronoi cell constructions for one full dock.
- **Dominant cost — Voronoi geometry:** `Vcontacts()` (`Vcontacts.cpp:14`), specifically `voronoi_poly2()` (`:172`) and `calc_areas()` (`:871`). `calc_areas` does **5× `sqrt` + an `atan(sqrt(...))` per face triangle**, all in `double`.
- **Second cost — contact-list scoring walk:** `vcfunction.cpp:392` chases a `ca_rec[...].prev` **linked list** per atom (pointer-chasing, cache-hostile), with a strided `energy_matrix` lookup (`get_yval`) per contact and optional `exp`/H-bond/metal transcendentals.
- **Per-pose heap churn:** `index_protein` `malloc`+`memset` of the cubed bounding box **every** `Vcontacts` call (`Vcontacts.cpp:1738, 1773`); `get_contlist4` `malloc(100000*int)` per call (`:1268`); `save_areas` grows `ca_rec` by realloc.
- **Precision:** geometry is `double` although atom coordinates originate as `float` (`atom.coor`).
- **Parallelism today:** OpenMP across the population *within* a generation, made thread-safe by per-thread deep copies of nearly every mutable structure (`gaboom.cpp:2703-2734`) — the cloning itself costs memory bandwidth. Offspring deferred-eval (`FLEXAIDDS_PARALLEL_REPRODUCE`) is currently **OFF**.

## 2. What optimization machinery already exists (reuse, don't reinvent)

| Asset | State | Action |
|---|---|---|
| `flexaids_configure_simd()` + AVX2/512/NEON (`cmake/FlexAIDOptions.cmake`) | Applied to all targets | Reuse as-is |
| `simd_distance.h` (920 lines, mature AVX-512/AVX2/NEON kernels) | Used only in clustering/RMSD/tENCoM — **not** CF scoring | Extend into scoring core |
| SoA Voronoi path (`FLEXAIDS_USE_SOA_DISTANCES`, `VoronoiCFBatch_SoA.h`, `AtomSoA.h`) | Built but **default OFF**, experimental (`Vcontacts.cpp:1972`) | Finish + default ON |
| `VoronoiCFBatch.h` span batch | Standalone, stub-CF benchmark only, dormant | Wire real `ic2cf` batch into GA |
| CUDA/Metal batch eval (`cuda_eval.cu`, `metal_eval.mm`, `GPUContextPool`) | Dispatch-connected but **reduced fidelity** (com/wal/sas only; zeros clash/con/hbond/gist — `gaboom.cpp:2428-2432`) | Bring to full CF parity |
| OpenMP GA population parallelism | **Real, default ON** (`gaboom.cpp:2773`) | Reduce per-thread copy cost |
| LTO/IPO + `-march/-mcpu=native` + `-DNDEBUG` + `-s` | On `FlexAIDdS`/`tENCoM` **only**; classic `FlexAID` is `-O3 -ffast-math` (`CMakeLists.txt:225-235` vs `766-793`) | Apply to fast target |
| `FLEXAIDDS_PARALLEL_REPRODUCE` (offspring deferred eval) | Gated OFF for reproducibility | Now unblocked by drift allowance |

**Bottom line:** the reusable, production-grade pieces are the OpenMP per-thread-workspace scheme, the SIMD CMake helper, the `simd_distance.h` kernel library, and the `FlexAIDdS` LTO/native flag block. The batch/GPU/SoA layers are already built and tested but deliberately kept off the default CF path pending fidelity parity — which the drift allowance now grants us.

---

## 3. Phased plan (prioritized)

Each phase is independently shippable, each guarded by the same acceptance gate (§4). Rough speedup estimates are *hypotheses to be measured*, not promises.

### P0 — Build-level free wins  *(risk: none · est. 1.1–1.3×)*
1. Give the classic `FlexAID` target the `FlexAIDdS` optimization flag block: `-flto`, CMake `INTERPROCEDURAL_OPTIMIZATION`, `-march=native` / `-mcpu=native`, `-DNDEBUG`, link `-s`. This is copying the block at `CMakeLists.txt:766-793` onto the `FlexAID` target (`:225-235`).
2. Add **Profile-Guided Optimization (PGO)** as a CMake option: instrumented build → run a representative dock (e.g. 1G9V, fixed seed) → rebuild with the profile. A single tight hot loop is the ideal PGO case; typically 10–20%.
3. Verify LTO doesn't break the `flexaid_core` OBJECT-library linkage caveat noted at `CMakeLists.txt:786-790` (link `.o` without `-flto`, or LTO the whole fast target as a monolith).

### P1 — Eliminate per-pose allocation & exploit the rigid receptor  *(risk: low · est. 1.3–2×)*
1. **Kill the per-call bounding-box `malloc`+`memset`** (`Vcontacts.cpp:1738, 1773`). The receptor index box is invariant across all poses (rigid receptor); allocate once and reuse via the existing `FA->vindex` / `hoist_receptor_index` caching path. Replace `memset`-per-call with **epoch stamping** (the repo already did `contacts memset→epoch`; extend the same trick to the box).
2. **Per-thread scratch, allocated once**: `get_contlist4`'s `malloc(100000*int)` (`Vcontacts.cpp:1268`) and `ca_rec` growth become preallocated per-thread buffers.
3. **Cache receptor-only Voronoi cells.** For a rigid receptor, receptor↔receptor contacts are pose-invariant; only ligand↔receptor and ligand↔ligand (+ flexible-residue) cells change. Recompute only moved atoms' cells; reuse the static receptor partition. (Ties into P5.)
4. Hoist invariants out of the inner loop: `4πr²` sphere areas, `energy_matrix` row base pointers, radius+`Rw` sums.

### P2 — Precision, data layout, and SIMD in the scoring core  *(risk: medium, drift allowed · est. 1.5–3×)*
1. **`double` → `float` in Voronoi geometry** (`voronoi_poly2`, `calc_areas`). Halves memory traffic and doubles SIMD lane count. Allowed under the drift decision.
2. **Flatten the contact linked list into a contiguous per-atom SoA array** (leveraging `AtomSoA.h` / `VoronoiCFBatch_SoA.h`). Stops pointer-chasing and makes the `com/wal/elec/hbond` accumulation vectorizable.
3. **Vectorize `calc_areas`**: batch face triangles; replace `sqrt` with SIMD `rsqrt` (+1 Newton step) and `atan` with a polynomial approximation (drift allowed). Add an area kernel alongside the existing `simd_distance.h` kernels.
4. **Default-enable and finish `FLEXAIDS_USE_SOA_DISTANCES`** (`Vcontacts.cpp:1972`), removing its experimental status.
5. **Batch scoring across the population**: evaluate a generation's poses as a SoA batch (real `ic2cf`, not the `VoronoiCFBatch.h` stub) so the static receptor stays hot in cache and wide SIMD applies — the `cpu_eval.cpp` / batch scaffolding already exists but is GPU-gated.

### P3 — Parallelism unlocked by the drift allowance  *(risk: medium · est. near-linear in cores)*
1. **Turn on `FLEXAIDDS_PARALLEL_REPRODUCE`** (offspring deferred eval). The ~0.2% chromosome drift documented in `OPTIMIZATION_KNOWN_ISSUES.md` is now acceptable. Keep a `FLEXAID_DETERMINISTIC` reduction mode for regression benchmarking only.
2. **Shrink the per-thread deep copy** (`gaboom.cpp:2703-2734`). The receptor state is read-only during scoring → share it across threads instead of cloning; clone only the small mutable ligand/pose/flex state. This is the biggest single lever on the existing OpenMP path's overhead.
3. Confirm load balance (`schedule(dynamic)` already present) and NUMA-friendliness for many-core hosts.

### P4 — Full-fidelity GPU CF + auto-dispatch  *(risk: high, largest ceiling · est. 5–20× on batch/screening)*
1. **Close the GPU fidelity gap.** Current kernels compute only `com/wal/sas` and zero `pb_clash / con / hbond / gist` (`gaboom.cpp:2428-2432`), which is why the GPU path is off the default. Add the missing terms to `cuda_eval.cu` / `metal_eval.mm`.
2. **Voronoi on GPU is branchy and unfriendly** — evaluate two routes and pick by measured parity:
   - (a) direct port of the analytic Voronoi (accurate, poor GPU occupancy), vs
   - (b) a **drift-tolerant contact-area model** (grid/neighbor-based or a per-pair area LUT) that reproduces the CF within tolerance and maps cleanly to GPU. The drift allowance makes (b) viable and is the recommended first attempt.
3. **Keep receptor + energy matrix resident on device**; batch the whole population per generation (the `pack_genes_batch` → `cuda_eval_batch` scaffolding at `gaboom.cpp:2411, 2537` already exists). Use Metal's `metal_eval_batch_multi` for many-ligand single-dispatch screening.
4. **Flip on auto-dispatch** in `UnifiedHardwareDispatch` once the Astex A/B parity gate (§4) passes: GPU when a fidelity-matched kernel exists, optimized CPU otherwise. `FLEXAIDDS_FORCE_CPU` stays as the escape hatch.

### P5 — Throughput/algorithmic wins (biggest for screening)  *(risk: medium · est. 10–100× VS throughput)*
1. **Two-stage screening**: use the in-repo NRGRank rigid `CoarseScreen` / `TwoStageScreen` (tests at `CMakeLists.txt:2799`) to prune a library before the expensive GA. Order-of-magnitude for virtual screening, no change to per-dock physics.
2. **Adaptive generations / early convergence**: detect plateau in best-CF and stop the GA early instead of always running `max_generations`.
3. **Incremental Voronoi** (with P1.3): recompute only moved-atom cells across generations where the ligand barely moves.

### Cross-cutting — "keep features, zero-cost when off"
Audit the GA hot path so entropy/thermodynamics code is **skipped**, not merely multiplied by zero, when disabled: confirm `statmech` / `ShannonThermoStack` / `tENCoM` / `encom` calls sit behind runtime flags that short-circuit before any per-pose work. This is what makes fast FlexAID genuinely "classic FlexAID" speed while retaining FlexAIDdS features on demand.

---

## 4. Acceptance gate (applies to every phase)

Because drift is allowed, "correct" means **rank-and-pose-equivalent**, not bit-identical. Every optimization must clear:

1. **Speedup:** wall-clock improvement on the reference dock (1G9V, fixed `FLEXAID_SEED`) and on an Astex-85 subset.
2. **Parity:** on Astex-85 — success rate (top-pose RMSD < 2 Å) **non-regressing**, top-pose RMSD distribution unchanged within noise, best-CF within a stated tolerance band. Reuse `scripts/reproduce_astex85.sh`, `benchmarks/m3pro/`, and `.github/workflows/benchmark-tier1.yml` / `benchmark-tier2.yml` as the regression harness.
3. **Determinism (benchmark mode only):** `FLEXAID_DETERMINISTIC` build reproduces run-to-run for CI A/B; production fast mode need not.
4. **Build health:** fresh CMake build clean on Linux GCC/Clang + macOS Clang; `ctest` green; `-DNDEBUG` and LTO builds link.

No phase merges without showing *both* a speedup number and a passing parity table — per the repo's "verify with actual execution" rule.

## 5. Deliverable shape

Since drift is allowed globally, the cleanest packaging is to **optimize the `FlexAID` target directly** (rather than fork a separate "fast" binary), exposing:
- a `FLEXAIDDS_FORCE_CPU` / backend-select flag (already present) for CPU-vs-GPU,
- a `FLEXAID_DETERMINISTIC` switch that pins serial-equivalent reductions for regression testing,
- entropy features behind their existing runtime flags, off by default.

## 6. Suggested sequencing & rough budget

| Order | Phase | Effort | Est. incremental | Cumulative single-dock (CPU) |
|---|---|---|---|---|
| 1 | P0 build flags + PGO | days | 1.1–1.3× | ~1.2× |
| 2 | P1 alloc removal + receptor cache | ~1 wk | 1.3–2× | ~2–2.5× |
| 3 | P2 float + SoA + SIMD scoring | 2–3 wk | 1.5–3× | ~3–6× |
| 4 | P3 parallel offspring + lean copies | ~1 wk | cores-bound | scales with cores |
| 5 | P4 GPU full-fidelity CF | 3–5 wk | 5–20× (batch) | screening ceiling |
| 6 | P5 two-stage + adaptive | 1–2 wk | 10–100× (VS) | throughput ceiling |

**Honest expectation:** a realistic **~3–6× single-dock CPU speedup** by end of P3, and an **order-of-magnitude+ screening throughput** gain once P4/P5 land — all subject to the §4 parity gate. Estimates must be replaced with measured numbers as each phase completes.

## 7. First concrete step

P0 is a self-contained, low-risk PR: lift the `FlexAIDdS` LTO/native/NDEBUG flag block onto the `FlexAID` target, add the PGO CMake option, and stand up the Astex-85 speedup+parity harness as the CI gate that all later phases report against.
