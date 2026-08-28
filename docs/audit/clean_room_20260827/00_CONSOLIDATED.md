# FlexAIDdS clean-room audit — consolidated, 2026-08-27

## Status of this audit

**Partial but source-verified.** Eight parallel review slices were dispatched.
The `reviewer` subagent provider (`opencode-go/qwen3.8-max`) returned
`401 Insufficient balance`; six slices died mid-read. The surviving leads were
then **re-verified by the lead agent against the source**. Every finding below
quotes code read this pass. Nothing is inherited from a prior audit document.

`SecurityHygieneAudit` was still running at write time; its results are not
folded in.

---

## Summary (highest severity first)

1. **HIGH / determinism** — Default RNG stream interleaves three consumers on
   one generator and re-seeds on every stream switch, collapsing each stream
   to its first draw. The correct per-stream generators exist and are **opt-in**.
2. **HIGH / determinism** — `-ffast-math` is on `flexaid_core` by default.
3. **HIGH / determinism** — `-mcpu=native` defaults ON, making the engine SHA
   and numerics host-dependent.
4. **HIGH / determinism** — Fixed-order log-sum-exp exists and is **opt-in**.
5. **HIGH / race** — `UnifiedHardwareDispatch::detected_` is an unsynchronized
   `bool` guarding lazy init from seven public entry points.
6. **HIGH / claim surface** — CSV column `success` is RMSD-only, sitting next
   to the contract-correct `success_pb` / `claim_ready`.
7. **HIGH / I/O** — CIF `get_float` returns `0.0f` for missing/`?`/`.` and for
   `atof` failure, placing atoms at the origin with no warning.
8. **MEDIUM / science** — Metal Boltzmann/LSE path is FP32; GPU ≠ CPU bits.
9. **MEDIUM / science** — `kB_kcal` is defined at least four times with three
   distinct values; CMA-ES uses a different constant from `statmech`.
10. **MEDIUM / science** — `BinarySnapshotRecord.score` is documented as
    `kcal/mol` for a CF proxy.

**Pattern:** in F-01, F-02, F-03, F-08 the *correct, reproducible mechanism
exists and the default selects the fast, irreproducible one*. The comments
know this. The defaults do not.

By contrast the **claim-integrity layer is genuinely well built** (see
Non-findings).

---

## Findings

### F-01 — `-ffast-math` is applied to the science-critical core library by default
- **Severity**: HIGH
- **Category**: determinism / science
- **Location**: `LIB/CMakeLists.txt:230-235`; also root `CMakeLists.txt:286, 567, 609, 669, 869, 992, 1040, 1065, 1093, 1150`
- **Evidence**:
```cmake
    # -fno-finite-math-only: keep -ffast-math optimizations but preserve
    # infinity()/NaN semantics used as sentinels in fast_optics, BindingMode,
    # SharedPosePool, GrandPartitionFunction, UnifiedHardwareDispatch, etc.
    target_compile_options(flexaid_core PRIVATE
        -Wall -O3 -ffast-math -fno-finite-math-only)
```
- **Why it matters**: `flexaid_core` carries scoring and thermodynamics.
  `-ffast-math` licenses reassociation and FMA contraction, so accumulation
  results become a function of compiler version, optimization decisions, and
  target CPU. AGENTS.md requires pose ranking to be bit-identical across runs.
  `-fno-finite-math-only` correctly rescues inf/NaN sentinels — the
  reassociation hazard remains.
- **Fix**: drop `-ffast-math` from `flexaid_core` and any scoring/thermo
  target. If the performance is required, replace it with the narrow subset
  that does not permit reassociation, and add a test asserting a fixed input
  produces a byte-identical score vector. Changes scoring numerics →
  **REQUIRES a feature flag + new tests per AGENTS.md.**
- **Confidence**: verified-by-reading

### F-02 — `-mcpu=native` defaults ON, making the binary and its numerics host-dependent
- **Severity**: HIGH
- **Category**: determinism
- **Location**: `cmake/FlexAIDOptions.cmake:71-81`, `106-118`; `LIB/CMakeLists.txt:237-243`
- **Evidence**:
```cmake
option(FLEXAIDS_MCPU_NATIVE "Compile flexaid_core with -mcpu=native on Apple/Clang arm64 (perf; not portable across machines, may change FP codegen under -ffast-math)" ON)
```
  The same file correctly defaults `BUILD_FLEXAID_FAST` **OFF** and documents
  that `-march=native` makes two correct builds of identical source produce
  different md5s, so METHODOLOGY §1 cannot compare them across machines.
- **Why it matters**: stacked on F-01 the *numerical output* moves too. The
  strict option is defaulted safe; the equally hazardous one is not.
- **Fix**: default `FLEXAIDS_MCPU_NATIVE` to `OFF`. Record effective ISA flags
  in the engine-SHA pin so a pinned result cannot be silently compared across
  CPU tunings.
- **Confidence**: verified-by-reading

### F-03 — Deterministic fixed-order log-sum-exp is opt-in
- **Severity**: HIGH
- **Category**: determinism / science
- **Location**: `LIB/log_sum_exp.h:18-21`; consumer `LIB/UnifiedHardwareDispatch.cpp:908-914`
- **Evidence**:
```cpp
inline bool fixed_order_lse_enabled() noexcept
{
    return env_bool("FLEXAIDDS_FIXED_ORDER_LSE", false);
}
```
```cpp
    if (flexaids::fixed_order_lse_enabled())
        return flexaids::log_sum_exp_fixed_order(values);

    if (!detected_) detect();
```
- **Why it matters**: log-sum-exp is how the partition function is accumulated.
  A correct fixed-order implementation exists — it runs only when an env var is
  set. Default falls through to SIMD/OpenMP/Metal, whose summation order
  depends on ISA and thread count. AGENTS.md requires bit-identity across
  thread counts.
- **Fix**: invert the default. Record the env var in the run receipt. Add a
  cross-thread-count bit-identity test. Changes thermodynamic numerics →
  **REQUIRES a feature flag + new tests per AGENTS.md** (the flag already
  exists; this is a default flip).
- **Confidence**: verified-by-reading

### F-04 — `UnifiedHardwareDispatch::detected_` is an unsynchronized `bool`
- **Severity**: HIGH
- **Category**: correctness (data race)
- **Location**: `LIB/UnifiedHardwareDispatch.h:216`; `LIB/UnifiedHardwareDispatch.cpp:72-78`; racing call sites at `534, 626, 656, 876, 912, 1040, 1113`
- **Evidence**:
```cpp
    bool         detected_  = false;
```
```cpp
void UnifiedHardwareDispatch::detect() {
    if (detected_) return;
    detect_cpu();
    detect_gpu();
    detect_libraries();
    detected_ = true;
}
```
- **Why it matters**: two threads entering any two of those seven entry points
  concurrently both observe `detected_ == false` and both write `info_`.
  `detected_ = true` is a plain store with no release ordering, so a second
  thread can observe `detected_ == true` while `info_` writes are not yet
  visible, then dispatch on a partially-initialized `HardwareInfo`. Failure
  mode is a silently wrong backend, not a crash.
- **Fix**: `std::call_once` / `std::once_flag`, or a function-local `static`
  `HardwareInfo` (C++11 thread-safe init).
- **Confidence**: verified-by-reading

### F-05 — CSV column `success` is RMSD-only
- **Severity**: HIGH
- **Category**: science (claim surface)
- **Location**: `LIB/DatasetRunner.cpp:7529-7534`, emitted at `7887-7888` / `7923-7925` and `8375-8376` / `8413-8415`
- **Evidence**:
```cpp
                docking_completed && rmsd_report >= 0.0f &&
                rmsd_report <= 2.0f && !result.seed_echo &&
                !result.pose_sha256.empty() &&
                result.rmsd_pose_sha256 == result.pose_sha256;
            // Legacy column remains exactly the same-pose ordered-RMSD gate.
            result.success = result.success_rmsd;
```
```cpp
                           "wall_time_s,success,success_rmsd,pb_pass,success_pb,claim_ready,"
```
- **Why it matters**: AGENTS.md states a pose is successful **only** when
  RMSD ≤ 2.0 Å **and** PoseBusters passes. The gating logic itself is sound
  (see Non-findings). The emitted CSV still carries a column named plain
  `success` whose value is the RMSD-only verdict, two columns away from
  `success_pb`. Any downstream consumer that does the obvious thing and reads
  `success` obtains an RMSD-only number. The comment marks it "legacy" — a
  known-debt signal, not a mitigation. The markdown summary *does* headline
  `claim_ready` correctly (`8230-8239`); the CSV does not.
- **Fix**: rename the column to `success_rmsd_legacy` (or drop it — `success_rmsd`
  already carries the value). Add a validator assertion that rejects any
  summary computed from a bare `success` column. Touches the claim path →
  **REQUIRES new tests per AGENTS.md**; migrate every reader (clean cutover).
- **Confidence**: verified-by-reading

### F-06 — Metal GPU thermodynamics run in FP32
- **Severity**: MEDIUM
- **Category**: science / determinism
- **Location**: `LIB/ShannonThermoStack/shannon_metal.metal:9-11, 57-70`; host bridge `LIB/ShannonThermoStack/ShannonMetalBridge.mm:6-8, 121-128`
- **Evidence**:
```cpp
// Note: Metal Shading Language does not support double; all FP is float (FP32).
// The host bridge (ShannonMetalBridge.mm) converts double→float before upload.
```
```metal
    weights[gid] = exp(neg_beta * (energies[gid] - E_min));
```
- **Why it matters**: Boltzmann weights, the summation reduction, and
  log-sum-exp all execute in single precision on the GPU, while the CPU path
  is `double`. Combined with F-03, the same input can yield different
  thermodynamic ledger values depending on which backend was selected —
  automatically. The limitation is honestly documented; nothing prevents a
  GPU-computed value from being reported as equivalent to a CPU one. The
  `E_min` shift is correctly applied, so the GPU path is overflow-safe.
- **Fix**: record the backend that produced each thermodynamic quantity in the
  run receipt. Restrict the Metal path to non-claim work, or add a
  tolerance-bounded CPU/GPU agreement test. Do not mix backends within one
  campaign.
- **Confidence**: verified-by-reading

### F-07 — Default RNG stream interleaves three consumers and re-seeds on switch
- **Severity**: HIGH
- **Category**: determinism
- **Location**: `LIB/RngSeed.h:121-182`; consumers `LIB/gaboom.cpp:5030-5033` (stream `0x9A800D`), `LIB/FOPTICS.cpp:846` (`RandomDouble()`)
- **Evidence**:
```cpp
// FLEXAIDDS_RNG_STREAM_FIX — DEFAULT OFF.
//
// OFF reproduces the historical single-generator behaviour bit-for-bit. That
// behaviour is defective: three streams interleave on one thread (GA 0x9A800D,
// Vcontacts 0x0C0A11 inside chromosome evaluation, FOPTICS 0xF0701C5) and the
// generator re-seeds on every stream switch, so each stream collapses to its
// first draw forever. It is nonetheless the behaviour every frozen reference
// number in this repo was produced under, so it remains the default until a
// baseline re-run retires it.
```
```cpp
    if (!rng_stream_fix_enabled()) {
        // ---- LEGACY PATH (default) — byte-identical to the pre-fix code. ----
        thread_local std::mt19937 rng = make_thread_rng(stream);
        if (cached_stream != stream || cached_epoch != epoch) {
            rng = make_thread_rng(stream);   // RE-SEED
            cached_stream = stream;
            cached_epoch = epoch;
        }
        return rng;
    }
```
  Roulette selection (`gaboom.cpp:2217-2218`) and FOPTICS random projection
  splits (`FOPTICS.cpp:846`) both call `RandomDouble()`, which uses stream
  `0x9A800D` via `lazy_thread_rng`.
- **Why it matters**: the comment is an unusually honest self-diagnosis and it
  is accurate. On the default path, every time evaluation (Vcontacts) and
  search (GA) or clustering (FOPTICS) interleave, the generator is reset, so
  subsequent draws in a stream are not independent — they are the first draw
  of that stream, forever. Clustering output therefore depends on how many
  times the stream was switched during the preceding GA, i.e. on evaluation
  count and thread interleaving. The correct per-stream `std::map` of
  generators exists (`RngSeed.h:185-201`) and is **opt-in**. Enabling it
  "changes the draw sequence for a given `FLEXAID_SEED` and therefore
  invalidates comparison against those numbers" — which is why it is off, and
  which is also why every frozen reference number in the repo was produced
  under a defective RNG.
- **Fix**: enable `FLEXAIDDS_RNG_STREAM_FIX` by default **after** a new
  baseline campaign is frozen under the fixed generator, with the flag
  recorded in every receipt. Until then, document on every published number
  that it was produced under the defective stream. Changes search trajectory
  → **REQUIRES a feature flag + new tests + a new baseline per AGENTS.md.**
- **Confidence**: verified-by-reading

### F-08 — CIF coordinates silently become the origin on missing or unparsable values
- **Severity**: HIGH
- **Category**: correctness / api-safety
- **Location**: `LIB/CifReader.cpp:134-138, 333-341`
- **Evidence**:
```cpp
// Get float value, return 0.0 if missing
static float get_float(const std::vector<std::string>& toks, int col_idx) {
    std::string v = get_val(toks, col_idx);
    if (v == "?" || v == ".") return 0.0f;
    return static_cast<float>(std::atof(v.c_str()));
}
```
```cpp
        a.x = get_float(toks, cols.Cartn_x);
        a.y = get_float(toks, cols.Cartn_y);
        a.z = get_float(toks, cols.Cartn_z);
        ...
        cif_atoms.push_back(a);
```
- **Why it matters**: `atof` returns `0.0` on a non-numeric string with no
  error signal. Combined with the explicit `?`/`.` → `0.0f` branch, a CIF
  atom with a missing or garbage coordinate is accepted as `(0,0,0)` with no
  warning. An atom at the origin corrupts Voronoi contacts, CF scores, and
  RMSD. The file is accepted; the molecule is wrong. Silent wrongness is
  worse than a hard error. The same helper is used for occupancy and B-factor
  (`336-337`), which is less severe but the same pattern.
- **Fix**: reject the atom (or the file) if `Cartn_*` is `?`/`.` or if
  `strtod` leaves a non-empty remainder. Do not default coordinates to origin.
- **Confidence**: verified-by-reading

### F-09 — `kB_kcal` is defined at least four times with three distinct values
- **Severity**: MEDIUM
- **Category**: science
- **Location**:
  - `LIB/statmech.h:28` — `0.001987206` (canonical)
  - `LIB/DiFT/DiFT.h:50` — `0.001987206` ("matches statmech")
  - `LIB/ShannonThermoStack/ShannonThermoStack.h:24` — `0.001987206`
  - `LIB/ReferenceEntropy.h:53` — `0.001987206`
  - `LIB/cmaes_search.cpp:31` — `0.001987204258` (**different**)
  - `LIB/MIFGrid.h:132` — `0.001987` (truncated)
  - `LIB/NATURaL/disco_natural.hpp:206, 230` — `0.001987`
  - `LIB/CavityDetect/CavityDetect.metal:144` — `kT = 0.592f` at 300 K
    (implies kB ≈ 0.001973, a fourth value)
- **Evidence**:
```cpp
// LIB/statmech.h:28
inline constexpr double kB_kcal = 0.001987206;   // kcal mol⁻¹ K⁻¹

// LIB/cmaes_search.cpp:31
constexpr double kB_kcal = 0.001987204258;
```
```metal
    float kT = 0.592f;  // k_B * T (kcal/mol at 300K)
```
- **Why it matters**: CMA-ES uses kB to convert temperature into an energy
  scale for its entropy trace (`cmaes_search.cpp:503-505, 657-659`). A 1 ppm
  difference from `statmech` is not going to reorder poses, but it is a
  second source of truth for a physical constant, and the Metal cavity kernel
  is ~0.7% off (0.592 vs 0.596). `config_parser.cpp:178-181` documents a
  *real* historical bug in which folding `kB_kcal` into GA β made it ~503×
  too large and collapsed Boltzmann weights — the comment is the reason this
  constant must live in one place.
- **Fix**: one `inline constexpr` in `statmech.h`; every other file includes
  it. Delete the locals. The Metal `0.592f` should be `kB_kcal * T` passed in
  as a buffer, not a magic number.
- **Confidence**: verified-by-reading

### F-10 — CF proxy labelled `kcal/mol` on the snapshot record
- **Severity**: MEDIUM
- **Category**: science (overclaim)
- **Location**: `LIB/BinarySnapshot.h:56`
- **Evidence**:
```cpp
    float    score;          ///< energy / CF score (kcal/mol)
```
- **Why it matters**: AGENTS.md forbids presenting the CF/contact-function
  scoring proxy as a free energy. `DatasetRunner.h:133-142` gets this right
  (`best_score ≡ CF … NOT ΔG`; `predicted_dG` is a historical name). The
  on-disk snapshot header does not. Anyone reading a `.bin` snapshot will
  treat the field as kcal/mol. Nearby, `MultiModelDock.cpp:213-269` correctly
  stamps `claim_validity: proxy_only` / `energy_domain: cf_arbitrary_units`
  — the snapshot header is the outlier.
- **Fix**: change the comment to `CF/contact-function scoring proxy (a.u.)`.
  Do not change the on-disk layout.
- **Confidence**: verified-by-reading

### F-11 — Unseeded `time(0)` fallback when `GB->seed == 0`
- **Severity**: MEDIUM
- **Category**: determinism
- **Location**: `LIB/gaboom.cpp:349-362`
- **Evidence**:
```cpp
	unsigned int tt;
	if (GB->seed == 0) {
		if (flexaids_rng::has_master_seed()) {
			tt = static_cast<unsigned int>(flexaids_rng::master_seed());
		} else {
			tt = static_cast<unsigned int>(time(0));
		}
	} else {
		tt = GB->seed;
	}
	printf("srand=%u\n", tt);
	srand(tt);
	flexaids_rng::set_master_seed(static_cast<std::uint64_t>(tt));
```
- **Why it matters**: DatasetRunner *does* wire a deterministic per-target
  seed (`DatasetRunner.cpp:178-181, 6198-6199`), so the claim path is
  protected. The classic `./FlexAID cfg.inp ga.inp` path, and any caller that
  leaves `GB->seed == 0` without setting `FLEXAID_SEED`, silently time-seeds.
  `srand(tt)` is also still called; `rand()` itself has no remaining call
  sites in `LIB/` (verified this pass), so `srand` is currently a no-op
  consumer — but it is a trap for the next person who reaches for `rand()`.
  The seed *is* printed (`srand=%u`), so an unseeded run is at least
  reconstructable from the log.
- **Fix**: refuse to start when neither `GB->seed` nor `FLEXAID_SEED` nor a
  master seed is set, instead of `time(0)`. Delete the `srand` call.
- **Confidence**: verified-by-reading

---

## Non-findings / verified-good

1. **The PoseBusters conjunction is properly enforced.**
   `LIB/DatasetRunner.cpp:7658-7663`:
   ```cpp
   // success_pb := success_rmsd && pb_pass; never imply pass without pb_ran.
   {
       result.pb_pass = false;
       result.pb_ran  = false;
       result.success_pb = false;
   ```
   All three flags **default to false**. Absent ≠ pass. `pb_ran` is tracked
   separately from `pb_pass`.

2. **The denominator does not shrink on input failure.**
   `DatasetRunner.cpp:5593-5597` records `success = false`,
   `rmsd_to_crystal = -1.0f`, `rmsd_fail_reason = "input_missing"` into
   `report.results[idx]`. N is preserved.

3. **The success gate is anti-gaming beyond the stated contract.**
   `7529-7532` additionally requires `!result.seed_echo`, a non-empty
   `pose_sha256`, and `rmsd_pose_sha256 == pose_sha256` (RMSD measured on the
   *same* pose that was elected). Elected pose is persisted with SHA256
   (`7442-7444`). Stronger than AGENTS.md demands.

4. **The 2.0 Å comparison is on an unrounded float** (`7529-7530`).

5. **`claim_ready` is the actual headline in the markdown summary**
   (`8230-8255`): RMSD ∧ PoseBusters ∧ tENCoM/Eigen ∧ hash receipts ∧
   protocol eligibility. Rates are divided by `report.total_systems`, not by
   surviving rows (`8165-8167`).

6. **Budget semantics match AGENTS.md on the DatasetRunner path.**
   `DatasetRunner.cpp:5811-5860`: default `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1`
   scales **population** by `max(1, n_flex_bonds/4)` and keeps generations
   fixed. Mode 0 (gen-scale) is labelled "NOT for claim runs". Mode −1 is
   labelled oracle-ceiling only. `FLEXAIDDS_BUDGET_SCALE` (default ON) is an
   extra **population** multiplier for `n_genes >= 14`, not a generation
   multiplier. This is the contract.

7. **CMA-ES and GA agree on fitness sign.** CMA-ES initializes
   `best_cf = +infinity` and updates on `fit[k] < result->best_cf`
   (`cmaes_search.cpp:552, 610-611`) — minimize CF. GA uses `evalue` as CF
   (lower is better) and a separate `fitnes` field for roulette. No sign
   inversion.

8. **`-fno-finite-math-only` is paired with every `-ffast-math`**
   (`LIB/CMakeLists.txt:231-235`). NaN-safety was anticipated.

9. **x86 SIMD flags are force-disabled on Apple Silicon**
   (`cmake/FlexAIDOptions.cmake:161-165`).

10. **A correct fixed-order log-sum-exp implementation exists**
    (`LIB/log_sum_exp.h:23-27`). Defect is the default, not the algorithm.

11. **Metal FP32 limitation is documented at both kernel and bridge**, and
    `E_min` is applied before `exp()`.

12. **`simd_distance.h` AVX-512 / AVX2 nesting is correct.**
    `#if FLEXAIDS_HAS_AVX512` (line 87) / `#elif FLEXAIDS_HAS_AVX2` (line 235).
    The earlier "unreachable AVX-512" lead is **dropped**. Tail loops are
    scalar (`153-155`, `288-291`); no over-read of the allocation.

13. **No `std::unordered_map` / `unordered_set` in clustering sources**
    (`cluster.cpp`, `FOPTICS.cpp`, `DensityPeak_Cluster.cpp`). FOPTICS
    random-projection *is* stochastic (F-07), but not via unordered
    iteration.

14. **No `pull_request_target` in `.github/`.**

15. **`rand()` has no remaining call sites in `LIB/`.** `srand` is still
    called (F-11) but currently has no consumer.

16. **`predicted_dG` is documented as a historical name, not experimental ΔG**
    (`DatasetRunner.h:133-142`, `DatasetRunner.cpp:6900-6905`). Fallback to
    CF when the ledger is absent is explicit.

---

## Units ledger (partial)

| Constant | Location | Value | Units | Status |
|---|---|---|---|---|
| `statmech::kB_kcal` | `LIB/statmech.h:28` | `0.001987206` | kcal mol⁻¹ K⁻¹ | canonical |
| `DiFT::kB_kcal` | `LIB/DiFT/DiFT.h:50` | `0.001987206` | same | PASS (duplicate) |
| `ShannonThermoStack::kB_kcal` | `LIB/ShannonThermoStack/ShannonThermoStack.h:24` | `0.001987206` | same | PASS (duplicate) |
| `reference_entropy::kB_kcal` | `LIB/ReferenceEntropy.h:53` | `0.001987206` | same | PASS (duplicate) |
| CMA-ES `kB_kcal` | `LIB/cmaes_search.cpp:31` | `0.001987204258` | same | **FAIL** (F-09) |
| MIFGrid `kBT` | `LIB/MIFGrid.h:132` | `0.001987 * T` | kcal/mol | **FAIL** (truncated) |
| NATURaL `kB` | `LIB/NATURaL/disco_natural.hpp:206` | `0.001987` | same | **FAIL** (truncated) |
| Cavity Metal `kT` | `LIB/CavityDetect/CavityDetect.metal:144` | `0.592f` at 300 K | kcal/mol | **FAIL** (implies 0.001973) |
| GA β | `LIB/config_parser.cpp:178-181` | `1/T` (T in Kelvin, CF a.u.) | dimensionless / K⁻¹ | PASS (deliberately not kB) |
| GPF β | `LIB/GrandPartitionFunction.cpp:16` | `1/(kB_kcal * T)` | mol/kcal | PASS |
| BindingMode β | `LIB/BindingMode.cpp:1180` | `1/(kB_kcal * T)` | mol/kcal | PASS |

GA β being `1/T` rather than `1/(kB T)` is **intentional and documented**
(`config_parser.cpp:178-181`): CF is in arbitrary units, not kcal/mol.
Folding kB in historically collapsed the ensemble. Do not "fix" it.

---

## Claim-path trace (DatasetRunner)

| Site | What it emits | Verdict |
|---|---|---|
| `5593-5597`, `5604-5608` | `success = false` on missing input | failure recorded; N preserved |
| `7529-7534` | `success = success_rmsd` (RMSD ≤ 2 ∧ same-pose hash ∧ ¬seed_echo) | **RMSD-ONLY** (legacy column) |
| `7658-7690` | `pb_pass`/`pb_ran`/`success_pb` default false; then `validate_elected_pose` | **ENFORCES-AND** |
| `7790-7793` | `claim_ready = success_pb ∧ protocol_claim_eligible ∧ score_pose_consistent ∧ …` | **ENFORCES-AND** (STRICT) |
| `7887-7926` | CSV row: `success`, `success_rmsd`, `pb_pass`, `success_pb`, `claim_ready` | mixed (F-05) |
| `8230-8255` | markdown headline is `claim_ready` | **ENFORCES-AND** |
| `8165-8167` | rates ÷ `total_systems` | denominator intact |

Answers:
- (a) Can success be reported without PoseBusters passing? **The column named
  `success` can (F-05). `success_pb` and `claim_ready` cannot.**
- (b) Is the 2.0 Å comparison on an unrounded value? **Yes.**
- (c) Can the denominator shrink on failure? **Not at the sites checked.**
- (d) Is ground truth used to SELECT the reported pose? **Not at the sites
  checked.** Election is CF/Softβ; RMSD is applied after, and the SHA check
  binds RMSD to the elected pose. (The bulk of `run()` was not read; this is
  not a proof of absence.)
- (e) Symmetry-aware RMSD exhaustive or approximate? **Not examined**
  (`calc_rmsd.cpp` unread this pass).

---

## Compiler flag matrix (science-critical targets)

| Target | Flags | Reproducibility risk |
|---|---|---|
| `flexaid_core` | `-O3 -ffast-math -fno-finite-math-only`; `-mcpu=native` if `FLEXAIDS_MCPU_NATIVE` (default ON, Apple/Clang arm64) | **HIGH** — F-01 + F-02 |
| `FlexAIDdS` | `-O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG`; `-march=native` (non-Apple) / `-mcpu=native` (Apple arm64) | **HIGH** |
| `FlexAID` | `-O3 -ffast-math -fno-finite-math-only`; extra LTO/native only if `BUILD_FLEXAID_FAST` (default OFF) | medium (fast-math only) |
| `tENCoM`, `tencom_entropy_diff` | `-O3 -ffast-math -fno-finite-math-only -flto -DNDEBUG` + `-march=native` (non-Apple) | **HIGH** |
| `benchmark_datasets`, `cavity_detect_cli`, others | `-O3 -ffast-math -fno-finite-math-only` | medium |

---

## Determinism verdict

**Results are not bit-identical across thread counts or machines in the
default configuration.** Independently sufficient reasons:

1. `-ffast-math` on `flexaid_core` (F-01)
2. `-mcpu=native` default ON (F-02)
3. Order-dependent log-sum-exp unless `FLEXAIDDS_FIXED_ORDER_LSE=1` (F-03)
4. Defective RNG stream interleaving unless `FLEXAIDDS_RNG_STREAM_FIX=1` (F-07)
5. Metal FP32 vs CPU FP64 if the GPU backend is selected (F-06)
6. Unsynchronized hardware detection (F-04) can pick different backends
   under concurrent first use

The DatasetRunner claim path *does* pin a per-target seed
(`deterministic_ga_seed`). That is necessary and not sufficient.

---

## Coverage

**Read and verified this pass (lead agent):**
- `LIB/CMakeLists.txt:230-243`; root `CMakeLists.txt` fast-math/native sites
- `cmake/FlexAIDOptions.cmake:71-118, 161-165, 207-273`
- `LIB/log_sum_exp.h` (full)
- `LIB/UnifiedHardwareDispatch.h:215-219`; `.cpp:72-78` + seven `detected_` sites
- `LIB/DatasetRunner.cpp` claim/budget/seed windows: `178-182`, `5593-5608`,
  `5810-5863`, `6198-6199`, `7439-7444`, `7529-7534`, `7658-7718`,
  `7790-7793`, `7887-7926`, `8165-8167`, `8230-8255`, `8375-8416`
- `LIB/DatasetRunner.h:133-166`
- `LIB/ShannonThermoStack/shannon_metal.metal`, `ShannonMetalBridge.mm`
- `LIB/simd_distance.h:1-300` (AVX-512 and AVX2 blocks + tails)
- `LIB/gaboom.cpp:349-366, 2186-2238, 2514-2538, 5026-5034`
- `LIB/RngSeed.h` (full)
- `LIB/cmaes_search.cpp:1-32, 408-422, 551-614, 748-811`
- `LIB/statmech.h:27-29`; `config_parser.cpp:178-181`
- `LIB/CifReader.cpp:134-145, 247-252, 333-341`
- `LIB/FOPTICS.cpp:821-860`
- `LIB/BinarySnapshot.h:55-57`
- `LIB/CavityDetect/CavityDetect.metal:134-159`
- `LIB/GrandPartitionFunction.cpp:15-19`
- kB definitions listed in the units ledger
- `.github/` for `pull_request_target` (none)
- `LIB/` for `rand(` (none remaining)

**Not examined — no findings should be inferred:**
`Vcontacts.cpp` / `vcfunction.cpp` (Voronoi degeneracies, NaN),
`calc_rmsd.cpp` (symmetry permutation), bulk of `DatasetRunner.cpp` (I/O
safety, `system()`/`popen()`), `top.cpp`, `BindingMode.cpp` beyond the β
sites, `statmech.cpp` overflow/Cv estimator, GPF vs single-ligand isolation,
`DensityPeak_Cluster.cpp`, `CleftDetector.cpp`, `build_rotamers.cpp`,
`TurboQuant.h`, `ParallelDock.cpp`, `GPUContextPool.h`, Mol2/SDF/PDB
readers beyond CIF `get_float`, `python/`, `scripts/` metric aggregation,
CI token permissions, dependency licenses, committed secrets.

**Unverified leads (do not treat as findings):**
- Scoring-proxy slice mentioned a "potential stale-plane issue" in the
  `Vcontacts.cpp` hull code before dying.
- `SecurityHygieneAudit` was still running at write time.

---

## Recommended next actions (priority order)

1. **Do not flip `FLEXAIDDS_RNG_STREAM_FIX` on a live campaign.** Freeze a
   new baseline under the fixed generator, then invert the default (F-07).
2. Default `FLEXAIDS_MCPU_NATIVE` OFF; drop `-ffast-math` from
   `flexaid_core` or replace it with a non-reassociating subset (F-01, F-02).
3. Default `FLEXAIDDS_FIXED_ORDER_LSE` ON (F-03).
4. `std::call_once` the hardware detector (F-04).
5. Rename CSV `success` → `success_rmsd_legacy`; reject summaries that read
   the bare column (F-05).
6. Reject unparsable CIF coordinates instead of origin-filling (F-08).
7. One `kB_kcal` (F-09); fix the snapshot comment (F-10); refuse unseeded
   classic-path starts (F-11).
8. Restore subagent credits and re-run the killed slices, especially
   Vcontacts degeneracies and `calc_rmsd.cpp`.
