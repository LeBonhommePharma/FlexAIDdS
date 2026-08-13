# Test coverage quality — real gates vs false confidence

**Date:** 2026-08-13
**Branch:** `cursor/audit-mt-reproducibility-science-1c2f`
**Companion:** [`2026-08-13_multithreading_reproducibility_science_robustness.md`](2026-08-13_multithreading_reproducibility_science_robustness.md)
**Question:** after a green `ctest` / `pytest`, would this suite actually catch the audit bugs, or is the pass count inflated by tautologies and files that never run?

This note is diagnostic. It changes no engine behaviour. Adding the high-value tests listed in §7
would currently **fail** on F1, F5, F13 (and would not exercise F10 without a concurrent-`GA()`
harness). They must not land as red CI; they belong with the corresponding parity-gated fixes.

---

## 1. Verdict

The suite is **not empty theater**. StatMech identities, the claim firewall, several numeric
regressions (fractional multiplicity → `log(0)`, Shannon outlier-bin clamp, SoftBeta NaN skip), and
the BindingMode CCBM mutual-information path are real gates that would fail if those implementations
broke.

It **systematically misses the audit's highest-impact bugs** (F1 stream collapse, F5 NaN ranking,
F10 concurrent `GA()`, F13 inverted `compute_joint_ensemble` MI, F14 energy-sort ties). In three
places the existing tests **give false confidence**: they execute the buggy function on a happy-path
input that cannot expose the defect, then pass.

`ctest` **89/89** is a count of **binaries / discovered wrappers**, not of GoogleTest cases. Inside
those binaries there are **1722** `TEST()`/`TEST_F()` cases in 83 registered files. Separately,
**94** `TEST()` cases live in three `.cpp` files that CMake never builds (`build_sources.ignore`).
Counting those files as "the suite" is coverage theater.

Line-coverage (`gcov`) was **not** used as the quality metric. It would overstate protection: the F1
RNG test *does* call `lazy_thread_rng`, and the QuickSort tests *do* execute `QS_ASC`, without being
able to fail when the bugs are present. The right question is: **would this test go red if the bug
existed?** For F1, F5, F10, F13, F14 the answer is no.

---

## 2. What was actually run (this session)

| Suite | Result | Notes |
|-------|--------|--------|
| C++ `ctest --test-dir build --output-on-failure` | **89/89 passed**, ~3.9 s | `clang-18` Release, `BUILD_TESTING=ON`, Metal/CUDA off. Re-run 2026-08-13. |
| Python subset | 14 passed, 17 skipped, **1 failed** | Reproduced `python/tests/test_runner_absolute_paths.py::test_bare_name_on_PATH_still_resolves`: Ubuntu `/bin/sh` → `dash`, so `Path(binary).name == "sh"` is false. Environment assertion, not an engine defect. |
| Prior full `pytest python/tests/` (same branch, earlier this day) | **1253 passed, 74 skipped, 1 env fail** | Not re-run in full this pass; the env failure is unchanged. |

CMake registers most `tests/test_*.cpp` files as one `add_test()` per binary (GoogleTest then runs
every `TEST()` inside). `gtest_discover_tests` is used only for CMA-ES. That is why `ctest` prints
89 items while `test_statmech.cpp` alone contains 107 cases.

---

## 3. Inventory

| Bucket | Files | `TEST()` / `def test_` |
|--------|-------|-------------------------|
| C++ on disk (`tests/test_*.cpp`) | 86 | 1816 |
| C++ **registered in `CMakeLists.txt` and therefore in CI** | 83 | **1722** |
| C++ listed in `build_sources.ignore` — **never compiled, never run** | 3 | **94** |
| Python `python/tests/test_*.py` | 64 | 1257 `def test_` |

Ignored (dead) C++ files:

| File | Cases | What they would have tested |
|------|-------|------------------------------|
| `tests/test_ga_population.cpp` | 61 | `fitness_stats`, `adapt_prob`, `roullete_wheel`, `crossover`, `mutate`, `remove_dups` |
| `tests/test_binding_mode_io.cpp` | 29 | BindingMode I/O |
| `tests/test_production_blockers.cpp` | 4 | Roulette bounds, converged `adapt_prob` finiteness, GPU pool singleton |

`test_production_blockers.cpp` is the most misleading: its name claims production-blocker
regressions, but it is not in any CMake target. The GPU-pool singleton case is duplicated by
`tests/test_ga_context.cpp` (`GPUContextPoolTest.SingletonInstanceNoGPU`), which *does* run and is
still tautological (`&pool1 == &pool2`).

---

## 4. Would the suite catch the audit findings?

| ID | Bug | Would current tests catch it? | Why |
|----|-----|-------------------------------|-----|
| **F1** | `lazy_thread_rng` stream collapse | **No — false confidence** | `RngSeedTest.LazyThreadRngRespectsMasterSeedEpoch` (`test_ga_context.cpp:101-109`) draws **one** stream (`0x9A800D`) three times. After `set_master_seed(11)` it asserts `a==b` then `a!=c` on that same stream. That is the behaviour of a single `mt19937`. F1 is interleaved **distinct** streams on one thread. No test mentions `0x5A6A9` (sugar-pucker) or two stream ids. |
| **F2** | Voronoi hull jitter inside `schedule(dynamic)` | **No** | No test of `0x0C0A11` jitter under varying `OMP_NUM_THREADS`. `MasterSeedOverridesRandomDevice` only checks `seed_from_env_or_random(0x0C0A11)` equality, not hull perturbation. |
| **F3** | OpenMP `reduction(+:sum)` order in LSE/Boltzmann | **No** | No thread-count bit-identity test on `UnifiedHardwareDispatch` LSE. Hardware-dispatch tests cover Shannon identities and a real outlier-bin regression (`RobustBinning.FarOutlierLandsInTheTopBinNotOnTopOfTheBulk`), not reduction associativity. |
| **F4** | `-ffast-math` + native arch | **N/A (build flag)** | Tests run against whatever binary CI built. They cannot detect cross-host bit drift. |
| **F5** | `QS_ASC` / `QuickSort` + NaN CF → possible rank-0 | **No** | `QuickSortTest` (`test_gaboom.cpp:283-360`) and duplicate `QuickSortGA` (`test_ga_core.cpp:41-128`) use **finite** energies only. Zero `isnan` / `quiet_NaN` in `test_gaboom.cpp`. Duplicate happy-path sorts inflate the case count. Contrast: `SharedPosePool.RejectsNaN` and `SoftBetaIdentity.RejectsNaNInfMembers` **are** real NaN gates — on other APIs. |
| **F6** | `ParallelDock::run_region` ignores `rng_seed` | **No** | `test_parallel_dock.cpp` tests `GridDecomposer`, `SharedPosePool`, StatMech merge. It never constructs `ParallelDock` or calls `run_region`. |
| **F7** | Unclamped `Cv`; multiplicity ≤ 0 → NaN | **Partial** | `FractionalMultiplicityNoNaN` is a **real** regression (int truncation of 0.5 → `log(0)`). `HeatCapacityNonNegative` uses well-separated energies and will not see FP-cancellation negative variance. No `add_sample(E, 0)` / all-zero counts case. |
| **F8** | Unsynchronized `detected_` | **No** | No concurrent first-call `detect()` test. |
| **F9** | RNG test cannot catch F1 | **Confirmed** | See F1. This finding *is* about the test suite. |
| **F10** | Concurrent `GA()` + `static ParEvalWS` | **No** | Parallel-dock tests do not run two `GA()` calls. No mention of `ParEvalWS`. |
| **F11** | `cluster.cpp` TOCTOU | **No** | No concurrent clustering test of `Clus_GAPOP`. |
| **F12** | Populate-eval omits `metal_coord` | **No** | Metal-coordination tests exist (`test_metal_coordination.cpp`) but do not compare the two eval-reset sites in `gaboom.cpp`. |
| **F13** | `compute_joint_ensemble` MI = `S_joint − S_R − S_L` | **No — false confidence** | `JointEnsemble_SingleReceptorFallback` forces `receptor_conformer_id = -1`, which **zeroes MI in the fallback branch** (`statmech.cpp:1073-1076`) and never hits the inverted formula. `JointEnsemble_ProbabilitiesSumToOne` only checks `Σp = 1`. `CCBMTest.MutualInformationNonNegative` **does** assert `MI ≥ 0` — on `BindingMode::ligand_receptor_mutual_information()`, which uses the **correct** sign (`BindingMode.cpp:1256-1257`, clamped). The production path is tested; the experimental API with the bug is not. |
| **F14** | `EnergyComparator` strict `<`, unstable `std::sort` | **No** | No BindingMode tie-break / equal-energy election test. |
| **F15** | Unclamped `acos` / unguarded `sqrt` in `geometry.cpp` | **No** | No geometry-angle NaN test. |
| **F16** | Dual β (CF-unit `1/T` vs `1/(kB·T)`) | **Partial** | Claim-firewall tests correctly refuse to treat CF scores as physical ΔG. They do not catch mixing the two β conventions inside one formula. |
| **F17** | `boltzmann_pmf` empty-bin F=0 / `T≤0` | **Partial** | Engine ctor throws `T≤0` (`InvalidTemperatureThrows`). `helmholtz` / `init_replicas` / `compute_joint_ensemble` do not; those entry points are not tested for `T≤0`. Affinity conversion **does** throw (`AffinityCalibrationTest.RejectsInvalidTemperature`). |
| **F18** | GPF `unordered_map` FP-sum order | **No** | Grand-partition tests check probabilities sum to 1 and some bit-identity for **identical numeric inputs**, not hash-order accumulation across insertion orders. |
| **F19** | `MultiSiteGPF` `exp(log_xi)-1` overflow | **No** | No huge-`log_xi` case found. |
| **F20** | Python `st_mtime` / unsorted `rglob` | **No** | Results I/O tests do not pin aggregation order. |
| **F21** | Model-scale absolute `S_vib`; `sqrt(|λ|)` | **Partial** | tENCoM tests exist (`test_tencom_diff.cpp`, `test_tencom_entropy_diff.cpp`) including a **real** NaN-overlap case for mismatched dimensionality. Absolute vs relative `S_vib` claim language is not gated. |
| **F22** | FOPTICS uses F1 RNG | **No** | `test_fast_optics.cpp` does not interleave RNG streams. |
| **F23** | `simd_distance.h` Boltzmann without max-shift | **No** | Unhit helper. |
| **F24** | Unseeded `time(0)` | **N/A** | Expected; `FLEXAID_SEED` tests exist and are real (`FlexaidSeedMakesStreamSeedsRepeatable`). |

---

## 5. What is real (keep, do not dismiss)

These would fail if the corresponding production code broke. They are the opposite of fluff.

- **Two-state analytical StatMech** (`test_statmech.cpp:549-575`): F, ⟨E⟩, Cv against a closed form.
  The comments record that a previous hard-coded `expected_Cv` **masked a formula bug**. That is
  how a useful test looks.
- **Single-state identities:** F = E, S = 0, Cv = 0 (`:512-525`); mirrored in
  `test_thermo_ledger.cpp` (`G = −kT log Z`, `−TS = G − H`).
- **Fractional multiplicity NaN regression** (`FractionalMultiplicityNoNaN`).
- **LSE / extreme-T finiteness** (`ExtremelyLowTemperatureFinite`, `ExtremeEnergySpreadLogsumexpStable`).
- **Claim firewall** (C++ `ScientificProvenanceTest.*` in `test_statmech.cpp`; Python
  `test_results_io.py` hostile-PDB cannot self-authorize physical claims;
  `test_thermodynamics_dataclass.py` / `test_models_deserialization.py` round-trips). Fail-closed
  on missing/malformed evidence is the point of the feature.
- **Thermo-impossibility gate** (`test_thermo_gate.cpp`): four quadrants + strict zero boundaries +
  the Shannon-H-cannot-trip-the-gate pin. Small file, high signal.
- **Affinity ΔG ↔ Kd round-trip** and invalid T/Kd throws.
- **Shannon robust binning** (`test_hardware_dispatch.cpp` `FarOutlierLandsInTheTopBinNotOnTopOfTheBulk`):
  documents a real INT_MAX/INT_MIN clamp bug that differed x86 vs ARM.
- **SoftBeta / classic entropy ranking NaN skip** (`test_classic_entropy_ranking.cpp`).
- **SharedPosePool NaN rejection** (`test_knowledge_pool.cpp`).
- **BindingMode CCBM MI ≥ 0** on the **correct** formula (does not cover F13).
- **Cleft-grid deterministic merge** — previously a real MT bug; tests remain around cleft/cavity.

---

## 6. What is fluff, duplicate, or false confidence

### False confidence (worse than missing)

1. **F1 RNG test** — executes `lazy_thread_rng` and cannot fail F1. A reviewer who greps for
   `lazy_thread_rng` will conclude streams are tested.
2. **F13 joint-ensemble tests** — named after the experimental API, then avoid the inverted-sign
   line by using the single-receptor fallback or by asserting only `Σp = 1`.
3. **QuickSort duplicates** — two files, finite energies only. NaN ranking (F5) looks "covered"
   because `QuickSort` appears in the test list.

### Tautologies that always pass

- `GPUContextPoolTest.SingletonInstanceNoGPU`: two references to a Meyers singleton.
- `GAContextTest.NotCopyable` / `IsMovable`: `std::is_copy_constructible_v` on a type with
  `unique_ptr`. Compile-time traits, not runtime GA re-entrancy.
- `GAContextTest.DefaultConstruction` / `IndependentCounters`: field zeros on two stack objects.
  Useful as a smoke for the struct layout; not a ParallelDock race test (the file's header comment
  claims that mission).

### Dead files counted as coverage if you `ls tests/`

- 61 + 29 + 4 GoogleTest cases never linked. `test_ga_population.cpp` in particular overlaps
  `test_ga_operators.cpp` (which **is** registered: `calc_poss`, `set_bins`, `validate_dups`).

### Python property-based "examples" that never call FlexAIDdS

`python/tests/test_property_based_example.py` is the clearest fluff. Hypothesis generates inputs,
then the assertions test **Python itself**:

```python
# test_boltzmann_weight_normalized — docstring claims Boltzmann weights
weight = 1.0  # Single state always has weight 1.0
assert weight == 1.0
```

```python
# test_energy_calculation_finite — never imports flexaidds
calculated = energy_val * 2.0
assert not (calculated != calculated)
```

```python
# test_config_value_lengths_bounded — never calls the C++/Python config parser
safe_value = config_line[:buffer_size-1]
assert len(safe_value) < buffer_size
```

`test_partition_function_positive` computes `sum(2.718**(-(e/0.5)))` in pure Python and asserts
`z > 0`. That cannot fail for finite negative energies. `test_addition_associativity` asserts
`a+b == b+a` for IEEE floats.

If Hypothesis is installed these still **collect and pass**, padding the 1253 count. They are
examples of a testing style, not gates on this repository.

### Python env-strictness (real intent, brittle assertion)

`test_bare_name_on_PATH_still_resolves` correctly wants PATH resolution of `sh`. On Ubuntu
`Path(resolved).name == "sh"` fails because `/bin/sh` is `dash`. The behaviour under test is
sound; the basename check is host-dependent. That is a test bug, not engine fluff.

---

## 7. High-value tests that do not exist (parity-safe to add *with* the fixes)

Do **not** merge these as failing tests. Pair each with the corresponding default-OFF / methodology
fix from the companion audit §6.

| Gap | Minimal test that would have caught it |
|-----|----------------------------------------|
| F1 | Same thread: draw stream A, stream B, stream A again. Assert the third A draw is **not** equal to the first A draw (independent generators) **or**, if the API is specified as "one RNG, stream is a seed key", document that and stop advertising independent streams. Today's test is the latter behaviour accidentally. |
| F5 | `QuickSort` a 3-chromosome array with `evalue = {NaN, -10, -5}`. Assert rank-0 is finite (or that NaN is rejected before sort, matching `SharedPosePool`). |
| F13 | Two receptor ids × two ligand ids, non-uniform joint `p(r,i)`. Assert `I ≈ S_R + S_L − S_joint` and `I ≥ −ε`. The fallback test can stay; it is a different branch. |
| F7 | `add_sample(-1.0, 0.0)` and an all-zero-count ensemble: finite Cv or a defined throw. |
| F10 | Two threads calling `GA()` (or a unit that `ParEvalWS` is not `function-static`). ThreadSanitizer, not just GoogleTest. |
| F14 | Two BindingModes with identical cached energy; assert election is stable (index or insertion order). |
| F6 | `run_region(..., rng_seed=S)` twice with `OMP_NUM_THREADS>1`; seeds must be honoured. |

Wiring `test_ga_population.cpp` / `test_production_blockers.cpp` into CMake is **not** a substitute
for the table above. Those files test other functions; they still would not catch F1/F5/F13.

---

## 8. How to read a green `ctest` after this

- **89/89** means every **registered binary** exited 0. It does not mean F1–F24 are absent.
- **1722 GoogleTest cases** is a large, mixed-quality corpus: genuine analytic StatMech + claim
  firewall + a handful of regressions, plus duplicates, trait checks, and experimental-API tests
  that miss the experimental-API bug.
- **94 additional cases** in `build_sources.ignore` are not CI coverage.
- **Python 1253** includes a real claim-firewall suite and a file of Hypothesis demos that do not
  import the engine. Subtract those mentally when using the number as a health metric.
- The honest coverage statement for the audit bugs: **the published single-GA, single-stream,
  finite-CF, BindingMode-MI path is well tested; the interleaved-RNG, NaN-rank, concurrent-GA, and
  `compute_joint_ensemble` paths are not.**
