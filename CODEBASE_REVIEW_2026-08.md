# FlexAIDdS — In-Depth Code & Scientific Review (2026-08)

Scope: full-tree review of the current `main`/`claude/codebase-review-gs804p` state —
~107 K lines C++, ~91 K lines Python, plus CUDA/Metal/HIP GPU code. Seven parallel
review passes (statmech core, scoring/entropy, GA & memory safety, Python package,
build/CI, tests & validation methodology, GPU/concurrency), each finding cross-checked
against the source. No files were modified.

---

## 0. Headline: the accuracy regression you're chasing

**Symptom (yours):** success rate regressed; last-good = "before #372"; macOS Metal build.

**Key structural facts established during review:**

1. **#372 and #370 landed 2 minutes apart** (2026-08-01, 21:10 and 21:12). "Before #372"
   is therefore *before both*, so both are in the regression window:
   - **#372** (`4974ac4`, mif-default) edited **only** `LIB/config_defaults.h` — it flips
     `seeding.mif_enabled` to `true` on the **unconfigured / CI-gate path**.
   - **#370** (`6fcaf09`, FlexAID-Fast) rewrote the CPU hot path: `Vcontacts.cpp` (+117),
     `gaboom.cpp` (+488), and is **default-ON** (`BUILD_FLEXAIDDS_FAST=ON`,
     `CMakeLists.txt:811`).

2. **On your macOS build, the GPU GA-fitness path is unreachable.**
   `UnifiedHardwareDispatch::best_backend(FITNESS_EVAL)` always returns a CPU backend
   (`UnifiedHardwareDispatch.cpp:207-222`). So #370's giant `metal_eval.mm` rewrite (+932)
   is **dead code for GA scoring** — it cannot be your regressor. The live parts of #370
   are the **CPU** changes to `Vcontacts.cpp`/`gaboom.cpp` and `FLEXAIDDS_PARALLEL_REPRODUCE`
   (now default-ON, "drift allowed").

3. **The mif-seeding mechanism is net-negative-prone.** `gaboom.cpp:3731-3742` replaces the
   anchor gene (`genes[0]`) of every seeded initial chromosome with a MIF-pocket-weighted
   inverse-CDF sample. The whole initial population is biased toward the *predicted* pocket:
   a big win when the MIF field is right (#372's only validation, 1MQ6: 7.3→2.7 Å), a loss
   when it's wrong. **#372 validated on n=1** — no 85-set evidence. On a diverse set this is
   exactly how a change improves individual targets yet regresses the aggregate.

4. **Attribution caveat.** If you run the **campaign path** (`DatasetRunner` /
   `reproduce_astex85.sh`), `DatasetRunner` already emitted `"mif_enabled": true` *before*
   #372 (#372 didn't touch `DatasetRunner.cpp`). On that path #372 is a **no-op** and the
   real change in the window is **#370's CPU rewrite**. If you run the **unconfigured
   engine / tier-1 gate**, #372 genuinely turned mif seeding on.

### The isolation experiment (single-factor, run on your Mac — GCC≥14)

Both suspects toggle **without bisecting**, so a 3-run matrix on the subset of targets you
saw regress disambiguates cleanly. Same `FLEXAID_SEED=12345`, same restarts, in-place
Hungarian RMSD (name the instrument, METHODOLOGY §0):

| run | `seeding.mif_enabled` | `BUILD_FLEXAIDDS_FAST` | isolates |
|---|---|---|---|
| baseline | `false` | `OFF` | both reverted |
| A | `true` | `OFF` | **#372 only** |
| B | `false` | `ON` | **#370 only** |

- If **A** regresses vs baseline → mif-seeding (#372) is the cause. Ship it opt-in/off by
  default, or improve the MIF field before biasing the anchor.
- If **B** regresses vs baseline → FlexAID-Fast (#370) is the cause; then narrow with
  `FLEXAIDDS_PARALLEL_REPRODUCE=0` and `FLEXAIDS_USE_SOA_DISTANCES=OFF` to find which sub-change.
- Use the repo's own METHODOLOGY §1 parity harness (byte-identical elected poses when the
  flag is OFF) as the instrument.

**Fastest first probe (no rebuild):** run your regressed subset with a config that sets
`seeding.mif_enabled: false`. If accuracy returns, you've confirmed #372 and can stop.

---

## 1. Critical / High findings (cross-verified)

### Validation integrity

- **[CRITICAL] Flagship "94.1 % Astex" is 90 %-native-seeded; the repo's own contract scores
  it 0/85.** `REPRODUCIBILITY.md:26` publishes 80/85 with `FLEXAIDDS_NATIVE_SEED_FRAC=0.90` +
  `SEED_ELITISM=1` (`reproduce_astex85.sh:420,462`); multiple targets read exactly 0.00 Å
  (crystal ligand surviving as elected pose). This contradicts `METHODOLOGY.md:38` ("NEVER
  report seed-elitism RMSD as the result") and `P0_CLAIM_CONTRACT_REPORT.md` (same run →
  0/85 STRICT). *You've acknowledged this — it's the single most important number to restate
  honestly.* Report the unseeded rate as the headline; seeded is an oracle ceiling.
  `reproduce_astex85.sh:564` also "verifies" against a ±3 % band around 94.1 rather than
  measuring independently.

- **[HIGH] Success gate threshold (0.70) sits above the method's own best figure and is
  unsourced.** `astex_diverse.yaml expected_baselines.docking_power_top1: 0.70` gates a
  regression below 0.665, while the YAML's own comment admits 0.70 is above FlexAID (0.66)
  and FlexAID_dS (0.69) "with no cited source." A genuine ~0.66 run is flagged as regression;
  only a seeded run clears it. The disclaimer exists only in `benchmarks/datasets/…` — it's
  **absent from the package copy** `python/flexaidds/dataset_runner/datasets/astex_diverse.yaml`
  that the runtime actually loads.

- **[HIGH] Silent-pass tests give false CI confidence.** `tests/test_aggregate_oracle_ceiling.py`
  and `test_posebust_upstream_parity.py` do a bare `return` (not `pytest.skip`) when private
  result dirs are absent → recorded as **PASSED**. The "92 % ceiling" they assert is an
  **oracle** metric (best pose selected by knowing true RMSD) sitting next to the 94.1 %
  headline — high risk of being read as achieved docking power.

- **[HIGH] License gate is a no-op on modern scancode.** `scripts/check_licenses.py:46` reads
  the pre-v32 `licenses[]` schema; `license-scan.yml:22` installs `scancode-toolkit` **unpinned**
  (v32+ dropped that array → every file yields `[]` → zero violations). Compounded by
  `scancode … || true` and exit-0 on parse error. `ALLOWED_LICENSES` is dead code. The
  no-GPL clean-room policy has **no working automated enforcement**.

- **[HIGH] Only end-to-end integration test never runs in CI.** `test_docking_pipeline`
  (ProcessLigand→StatMech→BindingMode) and `test_ligand_ring_flex` are defined inside
  `if(BUILD_SWIFT_BRIDGE)` (default OFF, `CMakeLists.txt:2992`); no CI job sets it ON.

### Memory safety (C++ core)

- **[HIGH] Unchecked stack overflow on `num_genes > MAX_NUM_GENES` (100).**
  `gaboom.cpp:2057` (`chrop1_gen[MAX_NUM_GENES]`), `:3519` (`icv[MAX_NUM_GENES]` in
  `eval_chromosome`), `coarse_init.cpp:198`. `FA->npar` grows unbounded via `realloc_par`
  with no cap; `reproduce()` `memcpy`s `num_genes*sizeof(gene)` into the fixed array. A large
  flexible ligand + many `FLEXSC` residues (>100 DOF) overflows every generation. Production
  path has no guard (only an `assert` compiled out under `NDEBUG`).

- **[HIGH] STEADY reproduction aliases `genes` pointers + leaks.** `gaboom.cpp:2247` shallow-
  copies `chromosome` structs (`chrom[a]=chrom[b]`), so two slots share one heap `genes`
  buffer and the original is orphaned; the next `reproduce()` `memcpy` corrupts the aliased
  live individual. Default `rep_model=BOOM` is safe; triggers only under STEADY. One-line fix:
  `swap_chrom(...)`.

- **[HIGH] `MPITransport::gather_results` loses/garbles results.** `MPITransport.cpp:83-126`
  broadcasts each region from `owner = r % world` — correct only for `ParallelDock`'s exact
  round-robin; any other decomposition broadcasts empty data. Non-owner ranks also reset
  `multiplicities` to all-1, corrupting the partition function. Reachable only under MPI.

### Scientific / GPU divergence

- **[HIGH] Metal Shannon/Boltzmann/LSE compute in FP32; thermodynamic outputs depend on which
  backend ran.** `shannon_metal.metal` + `ShannonMetalBridge.mm` do min/max/binning and
  `exp`/sum in `float`; the CPU reference uses `double`. `compute_boltzmann_batch` dispatches
  to Metal at `size ≥ 256` (`UnifiedHardwareDispatch.cpp:865`) and `ParallelDock.cpp:142`
  registers that `log_Z` as the reported ΔG. **An Apple run and a Linux run of the same input
  report different ΔG/entropy** — and, because election can depend on thermodynamic ranking,
  potentially different elected poses. This is *reachable* (unlike the GPU fitness path) and is
  directly relevant to reproducibility of your macOS results. Keep these reductions in FP64 on
  the host, or gate Metal off for values feeding thermodynamic outputs.

- **[HIGH, legacy/dead] `spfunction.cpp:172` wall term is dimensionally wrong** (r⁻⁶ − cr⁻¹²,
  non-vanishing at contact) and **`cffunction.cpp:96`** reads `distances[-1]` (OOB) / **:80**
  unbounded 100-element stack buffers. All in the dead `cffunction`/`spfunction` legacy path —
  latent, but should be deleted to prevent reactivation.

---

## 2. Medium findings

- **[MED, active] GIST off-by-one drops the last ligand atom.** `GISTEvaluator.cpp:203` loops
  `for(i=0;i<atm_cnt_real;++i)` on a 1-based `atoms[]` array — reads unused `atoms[0]`, skips
  the atom at `atm_cnt_real`. Disagrees with the correct GIST loop at `vcfunction.cpp:893`.

- **[MED, science] tENCoM/Shannon drop a fixed 6 "rigid-body" modes in *internal* coordinates.**
  `ShannonThermoStack.cpp:294` (`for m=6`) and `tencm.cpp:845,950` skip the first 6 modes as
  rigid-body, but the torsional Hessian is built in dihedral space (`tencm.cpp:128`), which has
  **no** 6-fold zero-mode manifold — so this removes the 6 *softest real* modes (largest
  entropy weight). The ENCoM path (`encom.cpp:140`) selects by eigenvalue cutoff instead, so
  the two entropy implementations use different mode sets on the same spectrum. Worth author
  confirmation — likely a systematic under-estimate of torsional S.

- **[MED, science] Physical-kB thermodynamics emitted with physical labels on arbitrary-unit CF
  scores.** `BindingMode::get_thermodynamics()` runs `β=1/(kB_kcal·T)` on CF (not kcal/mol) and
  `output_BindingMode` writes `REMARK free_energy/entropy/heat_capacity` as if kcal/mol
  (`BindingMode.cpp:745-762`) — the repo's own `SoftBetaFreeEnergy.h:11` / `read_input.cpp:253`
  say this β convention is invalid on CF units. Downstream `load_results()` reads them as
  physical. Label as CF-unit diagnostics or use the soft-β convention.

- **[MED, numerics] Heat capacity can go negative.** `statmech.cpp:165,173` uses uncentered
  `var = ⟨E²⟩−⟨E⟩²` with no clamp (while `std_energy` right below *is* clamped). With
  `⟨E⟩~−100`, cancellation can drive `var<0` → unphysical negative Cv in the ledger. Use
  centered variance + `max(0,·)`.

- **[MED, latent T=0] Constructing any `BindingMode` under a `Temperature==0` population throws.**
  `StatMechEngine` rejects `T≤0` (`statmech.cpp:56`) but `BindingMode` eagerly builds
  `engine_(Temperature)` (`:282`) while the rest of the stack treats T=0 as a valid "entropy
  off" mode (`read_input.cpp:248`, `BindingMode.cpp:1001`). Shielded only by convention.

- **[MED, parsers] Large-ligand OOB writes into `num_atm[]` (100000).** `Mol2Reader.cpp:296`
  and `SdfReader.cpp` V3000 path (`:152-209`) have no atom-count cap; ligand numbering
  `90001+ai` overflows at ≥9999 atoms. `read_lig.cpp:330` `altfdih` overflows a 3-int alloc
  when a FLEDIH record lists >3 atoms.

- **[MED, correctness] `abs()` on a double truncates to int** (compiler-confirmed).
  `gaboom.cpp:4608` (`remove_dups`) and `:4015` bind unqualified `abs` to `int abs(int)`, so
  the intended 0.1-IC dedup tolerance is effectively ~1.0 — duplicate collapse is far too
  aggressive, silently reducing GA diversity. Use `std::fabs`. *(This one is worth checking as a
  secondary contributor to accuracy — over-aggressive dedup thins the population.)*

- **[MED, Python] `_success` counts negative sentinel RMSDs as successes.** `benchmark.py:823`
  lacks the `0.0 <= r` lower bound that its sibling `metrics.docking_power:334` has; the
  runner's `-1.0`/`999.0` miss sentinels loaded via external records would score a miss as a
  success, inflating the rate.

- **[MED, Python] CSV thermodynamics not round-trippable / "frozen" dataclasses aren't hashable.**
  `models.py:437` writes the dataclass repr into CSV (JSON path is correct); `PoseResult`/
  `BindingModeResult` raise `TypeError` on `hash()` despite the docstring claiming they're usable
  as dict keys.

- **[MED, build] `cmake/MetalAcceleration.cmake` is dead & broken** (targets non-existent
  `flexaid_lib`, sources a missing file) yet CLAUDE.md advertises it as the Metal helper.
  `LIB/rocm_detect.cpp` is orphaned and the "strict source validator" misses it via a
  basename match against a dead patch file. 13 test targets re-list an identical ~35-file
  source closure by copy-paste.

- **[MED, GPU] `UnifiedHardwareDispatch::detect()` races** (plain-bool guard, no sync) despite
  a "thread-safe" doc; called lazily inside the OpenMP region. Use `call_once`/magic-static.

---

## 3. Low / notable

- **Determinism claim vs FP reductions.** `lse_openmp`/`boltzmann_openmp`/`rmsd_openmp` use
  `reduction(+:)` and AVX-512 block sums whose bit-result depends on thread count/vectorization
  (`UnifiedHardwareDispatch.cpp:576,727,1007,909`). `ParallelDock.cpp:117` seeds regions in
  thread-arrival order. These undercut the "byte-identical across thread counts" claim for the
  float partition-function paths (integer histograms are fine).
- **Two independent Hungarian-RMSD implementations, still no cross-check** between the pose-PDB
  `REMARK RMSD` and the `result.csv` value (METHODOLOGY §0 flags this; #371 added a test — verify
  it actually compares 4↔5).
- **Default build is non-reproducible** (`-march=native -flto -ffast-math` ON by default), so the
  git-provenance embedding can't be paired with a stable binary hash unless disabled.
- **TSan runs with OpenMP OFF** (`tsan.yml:33`) — the real threading model isn't exercised.
- **Claude workflows use mutable action tags** (`claude.yml`, `claude-code-review.yml`) while every
  other workflow SHA-pins; SHA-pin them.
- Classical (not quantum) vibrational-entropy oscillator in `encom.cpp:187` (documented heuristic,
  acceptable); Metal batch-completion synchronization keys on reused buffer pointers
  (`GPUContextPool.h:207`, screening path only); `SharedPosePool::serialize` ships uninitialized
  padding bytes.

---

## 4. What is genuinely strong (keep it)

- **Verification culture is a real asset.** `METHODOLOGY.md`/`VALIDATION.md` self-flag their own
  integrity risks (🔴), and the git history shows those gaps being closed: RMSD symmetry
  correction (element-blocked Hungarian, `benchmark.py:_symmetry_permutation`), pocket-field
  default (#372, with a before/after measurement), timeout→liveness accounting (`runner.py:1415`),
  liveness/productivity/completeness gates (#326). Several doc 🔴s are now *stale in the safe
  direction* (describing already-fixed issues as open).
- **Load-bearing statmech physics is correct** (independently re-derived): log-sum-exp with
  max-subtraction, parallel-tempering acceptance `Δ=(β_a−β_b)(E_a−E_b)` (detailed balance holds),
  grand partition function anchoring `Ξ≥1`, `S=(⟨E⟩−F)/T ≥ 0`, ΔG↔Kd standard-state.
- **Pure-Python StatMech fallback mirrors the C++ formula-by-formula** — no silent divergence in
  the highest-risk parity area.
- **C++ physics tests assert real properties** (`test_statmech.cpp` ~235 asserts: `S=k_B ln g`,
  Σw=1, Cv analytics, ΔG antisymmetry); `test_rmsd_symmetry.py` is exemplary (control test that a
  real 3 Å shift is *not* rescued, forbids cross-element swaps).
- **P0 claim contract is the right instrument** (frozen 85 denominator, missing=failure, real
  PoseBusters discrimination). Dataset self/cross-docking semantics gate is fail-closed and correct.
- OpenMP GA eval isolates per-thread buffers correctly; `roullete_wheel` selection is unbiased;
  GA-level elitism deep-copies; `cmaes_search.cpp`, `CifReader.cpp`, `SharedPosePool` locking are clean.

---

## 5. Suggested priority order

1. **Regression:** run the §0 3-run matrix on your regressed subset (start with the no-rebuild
   `mif_enabled:false` probe). Decide #372 vs #370.
2. **Restate the headline** unseeded (CRITICAL); fix the ±3 % self-confirmation and the silent-pass
   oracle tests (HIGH).
3. **Fix the license gate** (schema + pin scancode + drop `|| true`) — it's currently blind.
4. **Guard `num_genes ≤ MAX_NUM_GENES`** and the STEADY `swap_chrom` fix (memory safety).
5. **FP64 for Metal thermodynamic reductions** (macOS reproducibility of your own ΔG/entropy).
6. **`std::fabs` in the GA dedup** (diversity), GIST off-by-one, negative-Cv clamp.
7. Un-nest the integration test; wire it into CI.
