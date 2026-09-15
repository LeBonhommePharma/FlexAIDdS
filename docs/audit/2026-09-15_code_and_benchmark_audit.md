# FlexAIDdS — Code and Benchmark Audit (2026-09-15)

**Tip audited:** `959dd8eb60b0ed035e2076426f6eafc32c80d62c` (`main`, merge of #495)
**Auditor:** Cursor Grok 4.6 (cloud agent), first-hand source + CI + tests this session
**Mode:** Diagnostic. This PR does not change engine, ranking, thermodynamics, or election.
**Machine:** Linux cloud VM. **Mac live result trees are not here.**
**Identity:** GitHub `LeBonhommePharma`

Prior audits this document supersedes *as a snapshot*, not as a rewrite of their evidence:

- `docs/audit/2026-08-16_science_and_code_audit.md` (tip then `5c891e6e`) — several user-facing rate findings are **closed** by #437
- `docs/audit/2026-08-18_past_week_science_and_code_audit.md` (tip then `bde7908c`)
- `docs/audit/26h-swarm/SCIENCE_AUDIT.md`

Every factual statement below was re-read from files, `git log`, GitHub Actions, or a test run in this session. **No Astex dock was executed.** No docking-power percentage is reported as current product performance.

---

## 1. Executive verdict

The engine on this tip is a **CF-search / layered-ranking** docking system with a real claim firewall on the pages a reviewer actually opens (`README.md`, `docs/BENCHMARK.md`, `REPRODUCIBILITY.md`): this repository **publishes no receipted Astex-85 rate**. That honesty is the most important product fact. The scientific *idea* (search on the Voronoi CF proxy; optional density clustering; optional soft-β \(\tilde G\) on CF; post-election PoseBusters; torsional tENCoM as a ledger diagnostic) is coherent and much of the July–August P0 science debt is closed in code. What is **not** coherent is the **operator surface**: three ranking layers that docs still collapse into one word (“entropy”); a DatasetRunner default that emits ACF-ordered heads at `T=300` while electing min-CF behind a frequency gate; two energy-matrix MD5s still named as “the” pin; campaign-status files that still tell agents to write live GA to iCloud and to use TEMPER 298 / gen 6000; and a Python public API that still prints \(\Delta G\) in kcal/mol. **ΔS_vib is not claim-ready.** Shannon on the live path is **not** pose-space density. Cartesian ligand ENCoM is still the HVIB default. Flexed-receptor GA is **off** on DatasetRunner (`autoflex_max=0`). CI on `main` at this SHA is green (FlexAID CI 44m55s success after #495), with Eigen **5.0.1 hash-pinned** — that is real engineering. It does not mint a docking-power number. Treat every numeric campaign rate that is not a literature anchor or an explicitly withdrawn figure as **Mac-local, other-SHA, or comment-cited**. Reproduce from a local-first `result.csv` + `RUN_RECEIPT.json` on *this* binary, or refuse the claim.

---

## 2. What this VM could and could not see

| Path / source | Status on 2026-09-15 |
|---------------|----------------------|
| Git `main` @ `959dd8eb` | Present, clean |
| `MC_st0r5.2_6.dat` (repo root and `WRK/`) | Present; **md5 `9dc93717dfed0698006d88dd6a9627bc`** (the “9dc9” pin) |
| `~/flexaidds_results` | **MISSING** |
| iCloud / `Mobile Documents/…/FlexAIDdS_benchmarks` | **MISSING** |
| `/workspace/results` | **MISSING** (gitignored) |
| `benchmarks/astex_repro/full/` pose tree | **MISSING** (CSVs that point at it are orphans) |
| GitHub Actions on `main` after #495 | **FlexAID CI, sanitizers, tsan, coverage, CodeQL, license, tier-2 dry-run: success** |
| Live `result.csv` / campaign `RUN_RECEIPT.json` | **None tracked in git** |

Consequence: §6 truth table is **repo + CI + published/withdrawn text + in-code comments that name Mac arms**. It is not a re-score of those arms.

---

## 3. Architecture (do not conflate these layers)

| Layer | What it optimizes / reports | Default on DatasetRunner JSON path | Elects rank-0? |
|-------|------------------------------|-------------------------------------|----------------|
| **L1 Search** | Voronoi **CF** proxy (`vcfunction` / `ic2cf`) | Always CF | No |
| **L2 Cluster** | CF / FO / DP modes | **CF** (`clustering_algorithm{"CF"}`) | Groups |
| **L3a Engine emission** | Soft-β ACF order when `FA->temperature>0` && `!force_cf_rank_emission` | **ON** (`T=300`, `classic_entropy_ranking: true`, `force_cf_rank_emission: false`) | Yes for `_N.pdb` order |
| **L3b DatasetRunner S1** | Re-elect among emitted heads (Fix B freq-gate + min-CF, or opt-in Softβ) | Softβ **OFF**; Fix B **ON** | Yes for `elected_pose.pdb` / `result.csv` |
| **L3c Classic arm B** | `TEMPER 21` + `CLUSTA FO` + engine ACF | **Not** this path | Yes inside classic FlexAID |
| **L4 Ledger** | StatMech F/H/S/Cv; tENCoM λ; ShannonThermoStack `deltaG` | Off / diagnostic | **No** (vib correction fail-closed 0) |
| **L5 Validators** | RMSD + PoseBusters `bust_cli` | Post-election; missing `bust` → `pb_pass=false` | Does not re-rank |

**Softβ S1 ≠ FO@TEMPER21 ≠ engine ACF@300.** Citing `_0.pdb` as “the elected pose” on a DatasetRunner run is wrong: the harness copies `elected_pose.pdb` and says so (`LIB/DatasetRunner.cpp` ~8746–8747).

Default DatasetRunner JSON still writes:

```text
temperature: <DockingConfig, default 300>
clustering_algorithm: CF
classic_entropy_ranking: true
force_cf_rank_emission: false
mif_enabled: true
coarse_init: ON
seed_fraction: 0.0
autoflex_max: 0 (unless FLEXAIDDS_AUTOFLEX_MAX)
```

(`LIB/DatasetRunner.h` temperature default; `LIB/DatasetRunner.cpp` ~7244–7320.)

C0 claim script forces Softβ S1 off and uses `TEMP=298` (`scripts/run_C0_claim_clean.sh`) — engine ACF **still on** because `T>0`.

---

## 4. Critical findings (P0 / P1)

### P0-1. Three ranking products, one marketing word

**Evidence.** Engine ACF when `T>0`: `LIB/cluster.cpp` 377–394, 528–543. DatasetRunner Softβ S1 default OFF: `LIB/ProtocolConfig.cpp` 246–266. Classic arm B: `docs/implementation/3dsig_red_pair_protocol.md` §2.

README “What’s New” still says FlexAID∆S ranking **is** cluster-local \(\tilde G\) and that a two-gate spread guard is a product feature (`README.md` 88–96). Spread guard default is **off** (`cluster_spread_max{0.0f}`, `LIB/ProtocolConfig.h` 108–114). DatasetRunner S1 is min finite head CF when Softβ is off (`LIB/DatasetRunner.cpp` 1524–1527).

A reader who quotes README’s table has described **engine emission at T>0**, not the claim-harness elector, and has described a spread guard that is **inert**.

**Action.** Split the table into “engine emission / DatasetRunner S1 / classic arm B / spread guard (opt-in)”. Do not “fix” ranking to match the table.

### P0-2. Matrix pin fork is still live in the instruction surface

Repo files hash to **9dc9** (`9dc93717dfed0698006d88dd6a9627bc`). Claim generators and the 3Dsig protocol pin 9dc9 (`scripts/generate_flexaid_inp.py` `MATRIX_MD5_PIN`; `docs/implementation/3dsig_red_pair_protocol.md` line 20; `docs/implementation/arm_pins.json`).

Still pinning **72d7** (`72d7c7396702331d96ff12d18f831796`):

| Surface | Pin |
|---------|-----|
| `.grok/skills/flexaidds/SKILL.md` | 72d7 |
| `.grok/skills/flexaidds-dataset-runner/SKILL.md` | 72d7 |
| `.agents/skills/flexaidds-benchmarking/LIVE_QUEUE.md` | 72d7 |
| `METHODOLOGY.md` §3 blind republish | 72d7 |
| `docs/implementation/BLIND_ASTEX85_RECEIPT_PROTOCOL.md` | 72d7 |
| `scripts/blind_astex85_receipt_protocol.py` | 72d7 |
| `benchmarks/protocols/three_engine_entropy_comparison.md` (default / § notes) | 72d7 *or* 9dc9 |

Those are **different matrices**. Mixing them is mixing physics. `docs/implementation/arm_pins.json` calls 72d7 the “packing-sweetened fork” and forbids it for comparative arms.

**Action.** One sentence in `METHODOLOGY.md` and both skills: claim/comparative = 9dc9; 72d7 is a named historical fork. Update `LIVE_QUEUE.md` or mark it **superseded**. Do not re-fit either matrix.

### P0-3. `LIVE_QUEUE.md` is an ops landmine

`.agents/skills/flexaidds-benchmarking/LIVE_QUEUE.md` is dated **2026-07-15**. It still says:

- Storage: **iCloud only**
- Search: pop 1000 · **gen 6000** · restarts 5 · **T 298**
- Matrix: **72d7**
- C0 full85: **LIVE**
- Primary KPIs: S1/S2 (not 3Dsig S_top10)

That contradicts `AGENTS.md` (local-first, thin iCloud), the 3Dsig protocol (1000×2000, TEMPER **21** for arm B, S_top10, 9dc9), and `docs/implementation/CAMPAIGN_STATUS_2026-07-25.md` (no live science dock; C0 suspended). An agent that “reads the benchmarking skill then LIVE_QUEUE” will launch the wrong experiment on the wrong filesystem.

**Action.** Banner at top of `LIVE_QUEUE.md`: superseded by `3dsig_red_pair_protocol.md` + `CAMPAIGN_STATUS_2026-07-25.md` + `AGENTS.md` storage rules. Do not delete the historical queue; quarantine it.

### P0-4. Softβ policy doc disagrees with the binary

`docs/implementation/softbeta_election_policy.md` §4 still documents **`FLEXAIDDS_ACF_STRICT`** (lines 91–93). No reader for that name exists in `LIB/`. The live gate is **`FLEXAIDDS_ELECT_LEGACY_ACF`**, default **strict ON** (`LIB/cluster.cpp` 380–394; `tests/test_wave_gates.py` 38–42). Setting the documented name is a no-op.

§6 still says the Softβ-OFF path “keeps legacy ZH composite” (lines 121–133). The binary logs and elects **min finite head CF**; ZH is removed (`LIB/DatasetRunner.cpp` 1524–1527).

DatasetRunner Softβ ON uses **non-strict** `free_energy()` (`LIB/DatasetRunner.cpp` ~1585–1607). Engine emission uses `free_energy_strict(UniqueGeometry)` unless `FLEXAIDDS_ELECT_LEGACY_ACF=1`. Enabling Softβ S1 therefore elects a **different** \(\tilde G\) than FO@TEMPER21.

**Action.** Rewrite §4 and §6 to match the binary. Do not change the binary to match the doc.

### P0-5. ΔS_vib / Shannon are not the quantities the thesis slogan needs

**Vibrational.** `BindingMode::compute_vibrational_correction()` always caches **0.0** — `atom::eigen` holds eigenvectors, not eigenvalues (`LIB/BindingMode.cpp` 1507–1551). Real λ is ledger-only behind `FLEXAIDDS_LEDGER_TENCOM_LAMBDA` (default OFF), tagged inert on election (`LIB/tencom_ledger.h` 1–8). Search-path ligand ANM is skipped unless `FLEXAIDDS_TENCOM_WEIGHT>0` (default 0; `LIB/ic2cf.cpp` 36–47, 917–919). HVIB default ON is **Cartesian** unless `FLEXAIDDS_TENCOM_BASIS=torsional` (missing `_ligtopo.json` is a hard error — good). `encom.cpp` still filters eigenvalues by sign-like cutoff; tENCoM now classifies zeros by **magnitude** (`LIB/tENCoM/tencm.cpp` `classify_spectrum`, 590–617). Two ENMs, two zero-mode rules.

**Conformational Shannon.** Live Shannon is **not** pose-space density:

1. GA `H_final`: 20-bin histogram of CF `evalue` (`LIB/gaboom.cpp` ~971–975)
2. Softβ \(\tilde S\): \(-\sum p\ln p\) over Boltzmann weights of **CF** inside a cluster (score-space occupancy)
3. `BindingPopulation::get_shannon_entropy()`: same occupancy Shannon **× `kB_kcal`** (`LIB/BindingMode.cpp` ~290–298) — unearned kcal
4. `ThermodynamicEngine` gene-space 256-bin Shannon: `FLEXAIDDS_THERMO` only
5. FastOPTICS **is** pose-space density clustering — that is L2, not Shannon

Thesis constraint for this audit’s recommendations: **conformational entropy = pose-space density; vibrational = torsional-consistent tENCoM.** Neither is the default ranking term on this SHA.

**Action.** Keep ranking fail-closed. Next science work is FO density modes (arm B) and torsional tENCoM as a **validator**, not a CF weight.

### P0-6. Public Python example still prints kcal/mol ΔG

`python/flexaidds/__init__.py` 219–221:

```python
print(f"Mode: ΔG={mode.free_energy:.2f} kcal/mol")
```

`docs/SCORING.md` and README physics correctly call this proxy. The import docstring does not. That is the quote that will land in a slide.

### P0-7. No receipted Astex-85 docking-power number exists at this tip

User-facing pages are fail-closed (`README.md` badge + prose; `docs/BENCHMARK.md` withdrawal; `REPRODUCIBILITY.md` oracle-ceiling withdrawal). CI `repo_hygiene` runs `tests/test_check_published_astex_rates.py` on those surfaces. **Older audits, `CODEBASE_REVIEW_2026-08.md`, `LIVE_QUEUE.md`, and `docs/FLEXAID_FAST_DOCKING_PLAN.md` are not in that scanner.**

There is **no** `result.csv` for HEAD in git and none on this VM.

---

### P1-1. `atoi` vs `env_bool` — Softβ `=true` stays OFF

`EnvFlags.h` exists specifically because `atoi("true")==0` inverts A/B arms. `ProtocolConfig::from_env` still uses `env_truthy_int` → `atoi` for Softβ, Shannon-F, seed elitism, etc. (`LIB/ProtocolConfig.cpp` 30–33, 206, 259–265). `FLEXAIDDS_SOFTBETA_ELECTION=true` silently stays **false**. Cluster ACF uses `flexaids::env_bool` (`LIB/cluster.cpp` 380–381). Two parsers, two answers for the same English word.

`FLEXAIDDS_NO_SAS=0` still disables SAS: presence-only `getenv != nullptr` (`LIB/vcfunction.cpp` 49).

JSON unknown keys are silently ignored (`README.md` 172–173; `LIB/json_value.h` Null fallback). A typo is a default.

### P1-2. Fix B frequency gate is the default elector (and it costs singletons)

Default objective (`FLEXAIDDS_POOLED_ELECTION` unset): drop non-finite / degenerate CF; prefer `Frequency>1`; then min head CF (`LIB/DatasetRunner.cpp` 1719–1749). Softβ ON also sets `include_singletons` (1484–1486), so Softβ on all-singleton ensembles is **CF min with freq=1 admitted**, not an entropy ranking.

`LIB/pooled_election.h` 40–57 cites a Mac measurement (Lane B, `SELECTION_RECOVERY_REPORT.md` — **file not in this repo**) on 84×3 stored arms:

| Elector | Targets with a sub-2 Å top-1 pose / 84 (symmcorr) |
|---------|-----------------------------------------------------|
| Arm’s own elected pose (Fix B default) | **40.0** |
| min-CF over pooled restarts | **45.3** (+5.3) |
| Oracle (best pose in written pool) | **67.3** |
| Random pick | **4.3** |

McNemar on mincf vs default: **p = 0.09–0.39** per arm — not a claim endpoint (`pooled_election.h` 96–99). The +5.3 is a change of **objective**, not of pooling: pooling is already on (3/3 restarts on 254/254 receipts for `guard_arm_20260907_062124`, `seed2_arm_20260908_044750`, `seed3_arm_20260908_123008`). **Those trees are not on this VM.** Treat 40.0 / 45.3 / 67.3 as **comment-cited Mac measurements**, not HEAD docking power.

`FLEXAIDDS_POOLED_ELECTION` is **not** in `flexaidds_flags.cpp` dump.

### P1-3. CF ties: no declared max-tie rule; unstable sort; documented demotion lie

Default elector `std::sort` on score only (`LIB/DatasetRunner.cpp` 1755–1756) — equal CF is unspecified. Engine CF re-sort uses strict `<` so the **earlier cluster index** wins (`LIB/cluster.cpp` 548–552). Consensus `rank0_demoted` uses first min-CF index; “an elected pose tied on CF reads as demoted — rare, and documented” (`LIB/DatasetRunner.cpp` 8388–8389). `mincf` mode ties are enumeration-order (`pooled_election.h` 83–84). BindingMode `EnergyComparator` is strict `<` with no secondary key (`LIB/BindingMode.h` 264–273).

NaN CF into QuickSort is gated **OFF** (`FLEXAIDDS_NAN_RANK_GUARD`, `LIB/gaboom.h` 59–70). A non-finite CF can still participate in `QS_ASC` as `a-b`.

### P1-4. CI gate ≠ campaign DatasetRunner ≠ classic red-pair

`METHODOLOGY.md` §0.1 is still the right *shape*. Line numbers have drifted (it cites `DatasetRunner.cpp:6012` for `mif_enabled`; the write is now ~7244). Timeout row is **stale**: Python `TimeoutExpired` now increments the crash counter (`python/flexaidds/dataset_runner/runner.py` 1808–1824). Methodology still claims it does not (`METHODOLOGY.md` 145–154).

Live divergences that still hold:

| Knob | CI `python -m benchmarks.run` | Campaign `benchmark_datasets` |
|------|-------------------------------|-------------------------------|
| config | compiled-in defaults | writes `dock_config.json` |
| MIF | OFF unless JSON/env | **ON hardcoded** |
| coarse_init | OFF | **ON** |
| retained poses | 10 | **50** |
| permeability / normalize_area / clash ratio | engine defaults | campaign JSON |

Tier-1 docks **4** Astex targets against an **aspirational** `expected_baselines.docking_power_top1: 0.70` (`benchmarks/datasets/astex_diverse.yaml` 210–237) — a 2026-04-19 Claude re-estimate, not a measured campaign rate. JCIM copy in the same file is 0.452. A green tier-1 tick is not Astex-85 docking power.

CI skill job `flexaid_docking_datasetrunner` is **`--dry-run` only** (`.github/workflows/ci.yml` 491–497). Tier-2 push/cron is **`--dry-run`** until Zenodo IDs exist.

### P1-5. Five RMSD implementations, two live on the campaign path

`METHODOLOGY.md` §0. Still true. Claim metric is **symmcorr**; STRICT `success_rmsd` still inherits **serial**. Hungarian (engine REMARK vs DatasetRunner CSV) can disagree (SYBYL type vs element). Permissive `score_offline.py` still in-tree; 1HP0 is the documented false-success.

Do not mix REMARK RMSD, `result.csv` RMSD, and Python CI RMSD in one table.

### P1-6. PoseBust / kekulize: fail-closed pass, fail-open continue

Official `bust_cli` is the default claim backend. Missing `bust` → `pb_ran=false`, `pb_pass=false`, native QC does **not** mint `pb_pass` (`LIB/PoseBust/Engine.cpp` 700–712). Good.

`FLEXAIDDS_KEKULIZE_LIGAND_SDF` default **off** (`LIB/DatasetRunner.cpp` 3027–3029). Geometry-reconstructed ligands keep MDL aromatic order 4; RDKit `bust -l crystal` can `KekulizeException` (receipt comment, `Engine.cpp` 732–739). Measured historical cost: 7/84 Astex wrote 0-byte `bust_raw.csv`. That is a **validator abort**, not a chemistry fail, unless the receipt is read.

Python glob `r"_(\d+)\.pdb$"` (`runner.py` 2048–2050) collides FO `prefix_minPts_N.pdb` with CF `prefix_N.pdb`. C++ enumerator is careful (`DatasetRunner.cpp` 156–231). Mixed leftover FO files in a CF directory **enter** C++ election.

### P1-7. GPF / NATURaL silent zeros on live calls

`TargetServer` registration: `best_center` TODO still ships origin coordinates; `log_Z` falls back to `-predicted_dG / kT` when `ensemble_log_Z==0` (`LIB/DatasetRunner.cpp` 9414–9428). `predicted_dG` itself falls back to CF (`8089–8104`). Grand-canonical occupancy on that path is a **CF-scaled fiction**.

NATURaL calorimetric gaps call `set_inc(..., {0,0})` (`LIB/NATURaL/PoseLocalThermoRewrite.cpp` ~301–349). Experimental surface (`docs/SUPPORT_MATRIX.md`), but zeros are **wired**, not skipped.

GIST is hard-disabled on JSON (`LIB/config_parser.cpp` 137–148). Classic `.inp` still honors `USEGIS` (`LIB/read_input.cpp`). Two config languages.

`THERMO_SCORE` prints `enforced_in_final_election=0` (`LIB/gaboom.cpp` ~1684). Not a physics filter.

`best_score=0.0` when CF parse fails (`DatasetRunner.cpp` 8089) looks like a real proxy score.

### P1-8. Seed elitism header default vs claim-path override

`ProtocolConfig.h` 103: `seed_elitism{true}`. `from_env` default **true** (`ProtocolConfig.cpp` 206). Claim modes including **UNSET** force override 0 (`DatasetRunner.cpp` 8059–8066). A new `BenchmarkMode` that forgets the override will elect `_INI.pdb`. Oracle leak.

### P1-9. Windows / Eigen / FetchContent

Eigen **5.0.1** is hash-pinned (`cmake/FlexAIDEigen.cmake` 68–69, SHA-256 `e9c326dc…`). CI does **not** apt/brew Eigen. This closes the 3.x-vs-5.x eigensolver split that moved tENCoM zero-mode **counts**. Recent `classify_spectrum` magnitude rule is the complementary science fix (`LIB/tENCoM/tencm.cpp`; tests `test_svib_invariants.cpp`). **Do not revert the pin.**

GoogleTest: FetchContent `GIT_TAG v1.18.0`, **no commit hash** (`CMakeLists.txt` 1532–1539). pybind11 in CI is unpinned pip.

Windows: **no** MSVC engine job (`ci.yml` 119–120). Python bindings Windows is `allow_failure: true`. `docs/SUPPORT_MATRIX.md` line 27 still says Windows MSVC is **Supported**. `REPRODUCIBILITY.md` says native MSVC is **not** supported. The matrix is the lie.

AVX-512 Linux: `allow_failure: true`, `BUILD_TESTING=OFF`. Metal hosted smoke: missing compiler → **exit 0**.

`actions/cache@v6` unpinned in benchmark-tier workflows; Claude workflows use unpinned `checkout@v7`.

### P1-10. Stale comparative hub vs HEAD

`docs/implementation/COMPARATIVE_SCIENCE_README.md` last snapshot **2026-07-25**. P2 still HOLD (native CF oracle). Arm pins file `updated: 2026-07-25`. Pilot8 historical: **S1 = S_top10 = BCR = 0/8**, wrong matrix 72d7 (`CAMPAIGN_STATUS_2026-07-25.md` line 35). That is a **science-gate fail**, not a ranking paper.

`docs/FLEXAID_FAST_DOCKING_PLAN.md` §4 still treats Astex-85 success-rate non-regression as an acceptance gate without a receipted baseline.

---

## 5. What closed since the 2026-08-16 audit (do not re-open)

| August finding | Status on `959dd8eb` |
|----------------|----------------------|
| `docs/BENCHMARK.md` / `REPRODUCIBILITY.md` publishing 91.8% / 94.1% as current | **CLOSED** — withdrawn; scanner on those surfaces |
| `scripts/run_dataset.py` advertised as the reproduce path | **CLOSED** — file absent; BENCHMARK says so |
| SHARESCL 0.20 production default | **CLOSED** — `DEFAULT_SHARESCL = 10.0` |
| Softβ S1 default ON | **CLOSED** — default OFF |
| FO dual-suffix enumeration missing | **CLOSED** |
| Inclusive `latm = atm_cnt` (last HETTYP atom dropped) | **CLOSED** — `LIB/read_lig.cpp` 210–215 |
| Eigen 3.x Linux vs 5.x macOS | **CLOSED** — pin 5.0.1 + hash |
| tENCoM zero modes by sign | **CLOSED** on tENCoM path — magnitude classifier; `encom.cpp` still old rule |
| macOS ctest `allow_failure` | **CLOSED** — `macos_cpu_tests` blocking |
| Python timeout indistinguishable from empty dock | **CLOSED in code**, **OPEN in METHODOLOGY.md** (stale §0.1) |
| BindingMode global Z / missing SoftBeta header | **CLOSED** (already on 2026-08-16) |

---

## 6. Benchmark truth table

**Rule:** a number is **REAL** only if this session read a receipt or a literature table. Comment-cited Mac arms are **NAMED, NOT RE-DERIVED**.

| Figure | Where it lives | Status | What it actually is |
|--------|----------------|--------|---------------------|
| **No published FlexAIDdS Astex-85 rate** | `README.md`, `docs/BENCHMARK.md`, `REPRODUCIBILITY.md` | **CURRENT TRUTH** | Unverified / pending receipt |
| 91.8% (78/85) | withdrawn text | **WITHDRAWN** | Former “record”; #437 |
| 94.1% (80/85) | `REPRODUCIBILITY.md` historical block | **WITHDRAWN — oracle ceiling** | `SEED_ELITISM=1` + `NATIVE_SEED_FRAC=0.90` at `8196829f` |
| 88.2% (75/85) as FlexAID 2015 S1 | `docs/BENCHMARK.md` | **WITHDRAWN MISQUOTE** | Not JCIM Table 2 |
| ITC r = 0.93 / CNS 92% | BENCHMARK additional | **WITHDRAWN** | No `analyze_affinity.py`; no receipt |
| JCIM 2015 native FLRP **top-1 45.2% / top-10 66.7%** | `METHODOLOGY.md` §3 | **LITERATURE ANCHOR** | Gaudreault & Najmanovich 2015 Table 2. **Not this engine.** Do not swap labels. Arm A is 2015-era **or** master CF — not the 2015 binary (`3dsig_red_pair_protocol.md` §2) |
| 3Dsig 2017 ~0.66 / ~0.69 | protocol + YAML `success_rate_top10` | **HISTORICAL DECK** | **S_top10** bootstrap medians, 10 sims × 2e6 evals. Skill: not a current receipted FlexAIDdS rate. Benchmarking **not closed** |
| CI gate 0.70 top-1 | `astex_diverse.yaml` expected_baselines | **ASPIRATIONAL** | Claude re-estimate 2026-04-19; 4-target Ubuntu subset |
| 20/79 = 25.3% genuine; BCR 22/79 = 27.8% | `BASELINE_GENUINE_2026-07-24.md` | **OPS SESSION — not publishable** | Denominator 79 scored / 80 finished; predates later PRs; no §0 receipt |
| 17.9% / 31.0% / 48.8% Astex-84 frozen | swarm CSV + `SWARM_ORCHESTRATION.md` | **NAMED REFEREE, other SHA** | min-CF / pool ceiling on frozen poses; engine pin ≠ HEAD; 48.8 is **pool ceiling not S1** |
| 19/84 = 22.6% | `poster_metric_reference.csv` | **STALE ORPHAN** | spyrmsd **pooled min** RMSD; pose tree gone; engine sha `9f47c1bb` |
| `poster_metric_results.csv` 3/19 | truncated permissive CSV | **DO NOT USE** | Incomplete; 1HP0 false-success vs spyrmsd |
| Hungarian ceiling 57.8% | `METHODOLOGY.md` §0 | **INSTRUMENT WARNING** | Same frozen pool, over-permissive assignment |
| Pilot8 B0/B S1=S_top10=BCR=**0/8** | `CAMPAIGN_STATUS_2026-07-25.md` | **SCIENCE GATE FAIL** | Matrix was 72d7; arm A cad-only. Do not Softβ that ensemble |
| Guard / seed2 / seed3 84-target pooling receipts | `LIB/pooled_election.h` comments | **NAMED MAC ARMS, trees absent here** | `guard_arm_20260907_062124`, `seed2_arm_20260908_044750`, `seed3_arm_20260908_123008`; 3/3 pooling; 65.5% elected pose not from r0 |
| Lane B 40.0 / 45.3 / 67.3 | same header | **COMMENT-CITED** | `SELECTION_RECOVERY_REPORT.md` **not in repo**; symmcorr; not HEAD |
| Water ablation `water_ablation_20260906_022153` | drivers under `benchmarks/protocols/drivers/` | **MAC PATH ONLY** | Resume/DONE race documented; 1IGJ contaminated cell (72 targets into one directory) |
| Receptor-strain archive A_rigid **48/84**, B_shrink **32/84** | `PREREG_receptor_strain_2026-09.md` §2 | **ARCHIVE QUOTE, campaign NOT YET RUN** | Autoflex cost −16; PB volume_overlap 0/77 vs 38/77. Pre-registered 2026-09-03. Rigid vs flexed is **not** the DatasetRunner default (autoflex_max=0) |
| VED Stage-0 | — | **NO IN-REPO ARM** | Closest: Wave 0 ACF-strict / `build_wave0` |
| PB “rescoring” as CF re-rank | — | **NOT A LIVE ARM** | In-search `PB_CLASH_WEIGHT` default 0; PB is post-election S2 |
| E2 / E3 tENCoM | strain protocol E2 = elected endpoint; Wave 4 λ ledger | **NOT a tENCoM election arm** | λ inert on election |
| Softβ S1 campaign | policy + code | **FLAG OFF** | Cannot raise S1 if BCR=0 |
| Any HEAD `result.csv` % | — | **DOES NOT EXIST HERE** | Mac/iCloud only |
| Tier-2 nightly | `benchmark-tier2.yml` | **DRY-RUN** | Not a docking result |

### Selection vs sampling (what survives)

On **comment-cited** 84×3 Mac arms: elected ~40% vs pool oracle ~67% ⇒ **sampling ceiling ≫ election**. Fix B vs plain min-CF is +5.3 pp, not statistically clean per arm. On the **July 24 OPS session** (25.3% / BCR 27.8%): election gap ~2 targets, sampling is the wall — **unverified, other SHA**. On **pilot8**: BCR=0, so **no ranking experiment is allowed**. Those three stories are **not the same campaign**. Do not average them.

### FlexAID 2015 / literature comparators (caveats only)

- JCIM Table 2 is **top-1 45.2% / top-10 66.7%** native FLRP. Protocol ≠ current DatasetRunner (MIF, coarse_init, T=300 ACF emission, 50 poses, Fix B, PB, matrix 9dc9 vs whatever 2015 used).
- 3Dsig red bars are **S_top10**, 10×2e6, bootstrap median — not S1.
- Glide / Vina / GNINA rows in `docs/BENCHMARK.md` remain **approximate literature**, with the FlexAIDdS cell correctly blank/withdrawn.
- Repo git history does not contain the 2015 FlexAID tree (`docs/swarm/2026-08-13/SWARM_COMMON_PREAMBLE.md` notes history starts 2026-05-21). Arm A is a **pin or reconstruction**, not a bit-identical 2015 binary unless `arm_pins.json` Mach-O is rebuilt and receipted.
- **Do not build a league table from this audit.**

---

## 7. Integrity issues (harness, not scoring)

Documented, still the resume contract:

1. **DONE / receipt race.** Cell with `RUN_RECEIPT.json` but no `DONE` → resume docked the **entire dataset into one cell** (72 targets, 9,034 poses, 12.4 h). Gate: `run_water_ablation_w4.sh` pre-scrub + post-check rc=90 (`benchmarks/protocols/drivers/README.md` 19–35).
2. **STOP sentinels.** `touch $B/STOP` exists on water-ablation / determinism drivers; original `run_water_ablation.sh` lacked it.
3. **Contaminated pools.** 1IGJ redo (`run_igj_redo.sh`); leftover FO files in CF dirs enter C++ election; Python glob last-integer collision; `posehash()` globbing every `*.pdb` including orphans.
4. **Reconstructed DONE with `rc=unknown` must not count as `ok`** (`docs/run-uniformity/CONVENTION.md`).
5. **Early exit / `FLEXAIDDS_NO_SEC`.** Historical Astex-84 relaunch left security-channel early-exit **on** (`handoff_swe/AUDIT_astex84_vs_3dsig2017_20260809.md`): most restarts never spent 2e6 evals. Claim contract: `FLEXAIDDS_NO_SEC=1`.
6. **DoF budget.** Claim runs scale **population**, not generations (`AGENTS.md`). Freezing 1000×6000 because an old log printed it is the wrong experiment. LIVE_QUEUE still says 6000.

---

## 8. Claim-ready vs experimental (defaults)

| Quantity | Live default? | Enters rank-0? | Physical units? |
|----------|---------------|----------------|-----------------|
| CF / Voronoi | yes | GA search always | CF a.u. |
| Engine ACF (`T>0`) | **yes** on JSON/DR | emission order | CF a.u., T is score-T |
| DatasetRunner Softβ S1 | **no** | no | CF a.u. if ON |
| Fix B freq>1 | **yes** | yes | — |
| FO density clustering | **no** (CF default); arm B yes | defines modes | Å / MinPts |
| Spread guard | **no** | no | — |
| Seed elitism | header true; UNSET **forced off** | `_INI.pdb` if forgotten | oracle |
| Shannon H of CF histogram | yes (GA stop) | no | nats |
| ShannonThermoStack `deltaG` | print if enabled | no | mixed / unearned kcal |
| Cartesian ligand ANM / HVIB | weight 0 in CF; HVIB diagnostic ON | no | model-scale |
| Torsional ΔS_vib | opt-in basis | **no** (fail-closed 0) | model-scale |
| PoseBusters `pb_pass` | `bust_cli` if installed | post-election S2 | boolean |
| Flexed-receptor GA | **no** (`autoflex_max=0`) | n/a | n/a |
| GIST | JSON hard-off; classic inp can enable | n/a | broken evaluator |
| `PB_CLASH_WEIGHT` | 0 | no | — |
| `THERMO_SCORE` | off | **not wired** | diagnostic |

---

## 9. Recommended next actions (ordered, OVAT)

Respect: **do not re-optimize the CF energy matrix.** Conformational entropy = **pose-space density** (FO / occupancy of geometric modes), not a 20-bin CF histogram. Vibrational = **torsional-consistent tENCoM**, not Cartesian ENCoM mixed into CF.

0. **Documentation quarantine (this week, no dock).** One pin (9dc9) in skills + `METHODOLOGY.md` §3 + `LIVE_QUEUE.md` banner. Fix Softβ policy ghost `ACF_STRICT` / ZH-OFF lie. Split README What’s New ranking row. Delete or qualify kcal/mol in `python/flexaidds/__init__.py` docstring. Align `SUPPORT_MATRIX.md` Windows with CI. Refresh `METHODOLOGY.md` §0.1 timeout row and DatasetRunner line numbers. **One change per PR.**

1. **Receipt the Mac trees, do not re-dock first.** If `~/flexaidds_results` still holds `guard_arm_*` / `seed2_*` / `seed3_*` / water ablation, run `scripts/aggregate_claim_metrics.py` + `scripts/check_run_receipt.py` against 9dc9 and the HEAD binary SHA. Publish a **named, SHA-pinned, denominator-85** table or formally withdraw the comment numbers in `pooled_election.h`. OVAT: **no new GA**.

2. **Native CF oracle on any ranking claim.** `scripts/native_cf_oracle_gate.py` fail-closed before Softβ / FO election papers. Pilot8 0/8 is the template for “do not rank”.

3. **Sampling OVAT (only after oracle PASS), generations fixed.** One lever per arm: SHARESCL already 10 — do not touch unless receipted. Next: cleft / coarse_init / MIF **ablation to match classic FlexAID** (campaign ON vs CI OFF is a confound). Then FO@TEMPER21 vs CF@T=0 on the **same** binary and matrix (classic A/B0/B), headline **S_top10**. Do **not** enable DatasetRunner Softβ S1 on the same arm.

4. **Election OVAT, offline.** `FLEXAIDDS_POOLED_ELECTION=mincf` on frozen poses (already designed). Report McNemar. Keep default Fix B until a pre-registered arm wins. Do not also flip Softβ in that arm.

5. **PoseBusters S2, official `bust` only.** Kekulize CCD path as a **prep** A/B (`FLEXAIDDS_KEKULIZE_LIGAND_SDF=1`) with atom-type invariance tests already described in `DatasetRunner.cpp` 3015–3022. Do not treat KekulizeException as docking failure without the receipt argv.

6. **tENCoM validator OVAT.** `FLEXAIDDS_TENCOM_BASIS=torsional` + `_ligtopo.json` on elected poses; never `TENCOM_WEIGHT>0` in CF until a pre-registered ledger study. Do not mix Cartesian HVIB rows with torsional ΔS_vib.

7. **Flexed receptor.** DatasetRunner default is rigid search. The 48/84 vs 32/84 archive is a **different experiment**. Run the pre-registered strain campaign as written (`PREREG_receptor_strain_2026-09.md`) or leave it. Do not “just turn on autoflex_max=5” on Astex redock.

8. **Parser unification.** Migrate `env_truthy_int` to `env_bool`. Fail-closed unknown JSON keys. Register `FLEXAIDDS_POOLED_ELECTION` in the flag dump. Presence-only `NO_SAS` is a footgun — OVAT last.

9. **GPF / NATURaL.** Do not claim occupancy until `best_center` and real `ensemble_log_Z` are filled or the path is compile-gated off.

10. **CI.** Pin googletest by commit. Either add a real MSVC job or stop saying Windows is supported. Do not treat tier-1 0.70 as science.

---

## 10. What looks solid — do not touch

- **Claim firewall on README / BENCHMARK / REPRODUCIBILITY** + `tests/test_check_published_astex_rates.py` on those surfaces
- **Softβ S1 default OFF** + `[SOFTBETA-ELECT]` logs either way
- **`elected_pose.pdb` + SHA256** as the audit identity (not root `_0.pdb`)
- **Empty/missing energy matrix `Terminate`s** (`LIB/read_emat.cpp`)
- **`success_pb` requires `pb_ran`**; NativePoseQC cannot mint official pass
- **FO literature single MinPts + dual-suffix packaging** and C++ enumerator
- **Inclusive `latm`** last-atom emission
- **SHARESCL 10** production default
- **Eigen 5.0.1 hash pin** + tENCoM magnitude zero-mode rule + OpenMP buffer clamp
- **macOS blocking ctest** job
- **Wave 4 tENCoM λ ledger inert on election** (`test_wave4_tencom_ledger_firewall.py`)
- **Astex search receptor rigid** unless `FLEXAIDDS_AUTOFLEX_MAX` — reproducible, if you do not call it induced-fit
- **Seed-echo path-based `_INI.pdb`** (not CF tolerance)
- **Pooled election `mode_snapshot()`** so `best_score` and RMSD cannot diverge mid-run
- **CF vs thermodynamics language** in `AGENTS.md`, `docs/SCORING.md`, skill deception-proof table
- **Do not re-fit `MC_st0r5.2_6.dat`**

---

## 11. CI snapshot (this session)

GitHub `main` after #495 (2026-09-15 ~15:41Z): FlexAID CI **success** (44m55s), Sanitizers success, tsan success, coverage success, CodeQL success, license success, Benchmark Tier 2 success (**dry-run**). That is **build health**, not Astex docking power.

---

## 12. Auditor’s one-line for LP

You have a **defensible CF engine**, a **layered entropy story that is not yet one product**, and **no receipted Astex-85 rate on this SHA**. The next week of work is **documentation pin unification + Mac-tree aggregation**, not another GA. The next science week is **sampling (FO@21, oracle-gated) with generations fixed**, not Softβ, not a new energy matrix, not Cartesian ΔS_vib in the score.
