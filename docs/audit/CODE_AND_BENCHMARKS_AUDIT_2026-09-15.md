# FlexAIDdS code and benchmarks audit — 2026-09-15

**Tip audited:** `959dd8eb60b0ed035e2076426f6eafc32c80d62c` (`main`, merge of PR #495)  
**Remote:** `LeBonhommePharma/FlexAIDdS` (`origin/main` matched this SHA after `git fetch`)  
**Auditor:** Cursor Grok 4.6 (cloud agent), first-hand file / git / test-contract inspection  
**Mode:** Diagnostic. No engine defaults, ranking, or claim-ready Astex numbers were changed. Two trivial firewall clarifications shipped with this PR (`docs/EXPERIMENTAL_CAPABILITIES.md`, `docs/KNOWN_LIMITATIONS.md`).

This document re-audits architecture, the science-claim firewall, NATURaL / DualAssembly / PoseHelix / PoseLocal, benchmarks, reproducibility, correctness risks, and CI vs the support matrix. Every factual statement is anchored to a file read in this session. Live docking-power rates are **not** restated as current FlexAIDdS results.

Companion prior audits (not re-run as campaigns): `docs/audit/2026-08-16_science_and_code_audit.md`, `docs/audit/26h-swarm/SCIENCE_AUDIT.md`, `docs/audit/AUDIT_2026-07-23.md`.

---

## 1. Executive summary (one page)

FlexAIDdS on this tip is a **mature claim-culture repository with an immature claim package**. The engine, StatMech ledger, DatasetRunner harness, and CI firewalls are real. A receipted FlexAIDdS Astex docking-power rate is not.

**What is honest today**

- README, `docs/BENCHMARK.md`, and `REPRODUCIBILITY.md` refuse a live Astex-85 success rate. The former 80/85 (94.1%) figure is labelled a withdrawn oracle ceiling (`SEED_ELITISM` / native-pose seeding), not docking power.
- DatasetRunner Softβ S1 (`FLEXAIDDS_SOFTBETA_ELECTION`) is **OFF** by default. PoseHelixThermoRewrite (PR #471) and PoseLocalThermoRewrite (PR #472) are **OFF** by default and do **not** write `G_natural`.
- Source-audit tests (`tests/test_thermo_claim_firewall.py`, `python/tests/test_cf_naming_clarity.py`, `tests/test_check_published_astex_rates.py`, `tests/p0_claim_contract/`) actually enforce CF ≠ ΔG and block several withdrawn percentages from public vehicles.
- Literature FlexAID numbers in YAML are labelled published-external: Astex native top-1 **0.452**, Astex non-native **0.385**, HAP2 **0.220** (Gaudreault & Najmanovich 2015 JCIM Table 2); 3Dsig 2017 top-10 bootstrap medians **0.66 / 0.69** (FlexAID / FlexAID_dS).

**What is not claim-ready**

- There is **no** in-tree campaign `result.csv` / `RUN_RECEIPT` for a blind Astex-85 / FlexADS STRICT package. `benchmarks/astex_diverse/expected/` holds a README placeholder only.
- CI Tier-1 docks a 4-target subset against an **aspirational** `expected_baselines.docking_power_top1: 0.70` (YAML itself says this is not a measured FlexAIDdS rate). Tier-2 push/schedule is **`--dry-run` until Zenodo IDs exist**.
- STRICT / `claim_ready` still requires RMSD ≤ 2.0 Å **and** official PoseBusters (`bust_cli`) **and** tENCoM/Eigen on the elected pose SHA, plus protocol eligibility. That conjunction is specified, not deposited.

**The most important science-path surprises (defaults, not flags)**

1. **NATURaL DualAssembly growth is default-ON** for any receptor with residues on the FlexAID `GA()` epilogue unless `--folded` / `advanced.assume_folded`. It writes `FA->natural_deltaG`, which `BindingMode::compute_energy()` adds. Intra-target mode **order** is unchanged (pose-independent constant). Absolute ledger values, `has_natural`, and wall-clock are not. Docs still present NATURaL as an experimental *workflow* you opt into.
2. CSV column `predicted_dG` can be Helmholtz F **or literally CF** when StatMech is off (`LIB/DatasetRunner.cpp`). The header looks like experimental ΔG. Comments and tests document this; downstream plots will not read the comments.
3. **Five RMSD implementations** exist (`METHODOLOGY.md` §0). CI gates on positional Python RMSD; claim aggregator can use symmetry-corrected spyrmsd. Unlabelled “RMSD” is not reportable.
4. `docs/SUPPORT_MATRIX.md` still says Linux **GCC ≥ 10 / Clang ≥ 10** and **Windows MSVC 2022 Supported**. CMake fatals below **GCC 14 / Clang 18**; `ci.yml` **excludes** the FlexAID engine on MSVC. The matrix overclaims the release surface.
5. `benchmarks/datasets/astex_nonnative.yaml` is **deprecated in-file** (2026-09-11): the old roster is largely cross-*protein*, not cross-dock. Using it as “Astex non-native” is a scientific error. Replacement is `astex_nonnative_v2.yaml`.

**Claim-ready vs experimental (this SHA)**

| Surface | Status on `959dd8eb` |
|---------|----------------------|
| FlexAID / FlexAIDdS CLI docking (CF search) | Core 1.0 execution surface (`docs/VALIDATED_CAPABILITIES.md`) |
| Python `flexaidds` analysis APIs (labelled) | Validated for reporting when provenance is `proxy_only` |
| StatMech canonical F/H/S/Cv ledger | Validated for **analysis**, not experimental ΔG_bind |
| tENCoM CLI / vib diagnostics | Validated as a tool; STRICT still needs Eigen receipts |
| FlexAIDdS Astex S1 / S2 / STRICT rates | **Not claim-ready** — pending receipt |
| FlexAID 2015 / 3Dsig 2017 published numbers | **External published** — cite the paper/YAML, not this binary |
| DatasetRunner Softβ S1, THERMO_SCORE, USE_SHANNON | Experimental / diagnostic; default OFF |
| NATURaL DualAssembly growth | **Experimental code on the default GA epilogue** |
| PoseHelix / PoseLocal | Experimental, default OFF, no GA atom-coord hook |
| Swift / Fleet / iCloud / CUDA / Metal runtime | Experimental |
| Grand-canonical Ξ fixtures | Math harness only — “do not claim real ΔG” |
| DualAssembly Aβ42 scientific track | Spec + synthetic CLI; `--real-ga` refuses |

**Recommended next five actions** are in §10. None of them rewrite Astex numbers.

---

## 2. Provenance of this audit

Inspected in this session (non-exhaustive): `AGENTS.md`, `METHODOLOGY.md`, `docs/VALIDATED_CAPABILITIES.md`, `docs/EXPERIMENTAL_CAPABILITIES.md`, `docs/SUPPORT_MATRIX.md`, `docs/CI_RELEASE_GATES.md`, `docs/KNOWN_LIMITATIONS.md`, `docs/thermodynamics.md`, `docs/BENCHMARK.md`, `docs/BENCHMARKS.md`, `docs/REPRODUCIBILITY.md`, `REPRODUCIBILITY.md`, `docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md`, `docs/POSE_HELIX_THERMO_REWRITE.md`, `docs/POSE_LOCAL_THERMO_REWRITE.md`, `docs/FLEET_BENCHMARKS.md`, `docs/ICLOUD_BENCHMARK_STORAGE.md`, `workorders/CLAIM_LANGUAGE_FREEZE.md`, `benchmarks/README.md`, `benchmarks/protocols/{admission_metrics_contract,science_exclusions,three_engine_entropy_comparison}.md`, `benchmarks/datasets/{astex_diverse,astex_nonnative,hap2}.yaml`, `.github/workflows/{ci,benchmark-tier1,benchmark-tier2}.yml`, `LIB/{gaboom,BindingMode,DatasetRunner,NATURaL/*,tENCoM/tencm_cpu_fallback}.cpp`, `tests/p0_claim_contract/test_fixed_denominator.py`.

`PRODUCT.md` is **referenced** (`docs/EXPERIMENTAL_CAPABILITIES.md` promotion rule 4, `docs/KNOWN_LIMITATIONS.md`) and **absent** from the tree at this SHA.

PRs #471 / #472 are merged: `a4096559` (2026-09-08, PoseHelix + Zhao JCP 2011 lineage), `7bcb33e0` (2026-09-08, PoseLocal cited tables).

---

## 3. Architecture map — Core 1.0 vs experimental

### 3.1 Layer model (do not conflate)

| Layer | What it does | Units | Elects poses? |
|-------|----------------|-------|----------------|
| **L1 Search** | GA fitness on Voronoi **CF** (Vcontacts / `vcfunction`). Default fitness model **SMFREE** blends niche-share with soft-β weights on CF (`LIB/gaboom.cpp`, `ProtocolConfig` default `"SMFREE"`). | CF a.u. | Explores; SMFREE changes *selection*, not the CF operand |
| **L2 Cluster** | CF / FO / DP modes (`cluster.cpp`, FastOPTICS, DensityPeak) | geometry + CF | Groups |
| **L3 Rank** | Classic Softβ ACF \(\tilde G=\tilde H-T\tilde S\) on **CF** when `TEMPER>0` and `force_cf_rank_emission` is false (`config_defaults.h`: `classic_entropy_ranking: true`). DatasetRunner Softβ S1 is a **separate** post-cluster election and is OFF. | CF a.u.; T is score temperature, not \(k_B T\) kcal | Yes, engine path when T>0 |
| **L4 Ledger** | `StatMechEngine` / `BindingMode::get_thermodynamics()` F, H, S, Cv; optional `G_vib`, `G_natural` | kcal when \(k_B\) is used; still `proxy_only` without full validated path | No, unless an explicit ranking flag is on |
| **L5 Experimental** | PoseHelix / PoseLocal / Swift / Fleet / GPU backends / DualAssembly *rewrites* | mixed | No on defaults for rewrites |

`AGENTS.md` scientific guardrails remain the claim law: never sell L1/L3 as experimental ΔG_bind.

### 3.2 Binaries and `LIB/`

Both `FlexAID` and `FlexAIDdS` link `flexaid_core` (`LIB/CMakeLists.txt` `FLEXAID_CORE_SOURCES`) and share `LIB/top.cpp`. `FlexAIDdS` additionally compiles DatasetRunner / PoseBust redock sources when `BUILD_FLEXAIDDS_FAST` is on (default ON).

Notable `LIB/` modules in core:

- Scoring: `Vcontacts.cpp`, `vcfunction.cpp`, `cffunction.cpp`, `VoronoiCFBatch.h`, HBond/GIST evaluators
- GA: `gaboom.cpp`, `GAContext`, `ga_constants.h`, optional `cmaes_search` (`FLEXAIDDS_SEARCH=cmaes`)
- Thermo: `statmech.cpp`, `BindingMode.cpp`, `encom.cpp`, `ThermodynamicEngine` (score gate via `FLEXAIDDS_THERMO_SCORE`, default OFF)
- Shannon: `ShannonThermoStack/`, `shannon_ga.cpp`
- tENCoM: `LIB/tENCoM/`
- NATURaL: compiled into `flexaid_core` (not a separate optional object for FlexAID)
- PoseBust, CavityDetect, LigandRingFlex, ChiralCenter, DiFT, GrandPartitionFunction (`FLEXAIDS_GRAND_CANONICAL` default ON)

`ENABLE_DUAL_ASSEMBLY_TOOL` defaults **ON** (standalone `dual_assembly` CLI). `ENABLE_NATURAL_HAMMERHEAD` defaults OFF. `BUILD_SWIFT_BRIDGE` defaults OFF.

### 3.3 Python

`python/flexaidds/` is a Core 1.0 surface. `_core` is optional (`BUILD_PYTHON_BINDINGS` default OFF). Dual binding sources exist: `python/bindings/core_bindings.cpp` (CMake) vs `python/flexaidds/_core.cpp` (setuptools). Pure-Python fallbacks cover StatMech/ENCoM when `_core` is missing. `grand_canonical` still documents C++ wiring as future (`HAS_GRAND_BINDINGS = False` in package init at last inspection of the architecture map).

High-level `docking.py` can `rank_by_free_energy()`; examples now say “F-like” / `claim_validity`, which is better than older ΔG copy, but the method name remains easy to misread.

### 3.4 Swift / fleet / GPU

`swift/` is experimental (`docs/SUPPORT_MATRIX.md`). CUDA / Metal / ROCm / AVX-512 / WebGPU are experimental backends. Metal **CMake default is ON on Apple+ObjCXX** while the support matrix still marks Metal experimental — hosted CI only compiles shaders (`macos_metal_compile_smoke`); full link is `metal-self-hosted.yml` (`workflow_dispatch`, label `self-hosted-m3`).

### 3.5 Docs vs code: `G_natural` “validated”

`docs/VALIDATED_CAPABILITIES.md` lists additive `G_vib` / `G_natural` with presence flags as **validated for reporting/analysis**. `docs/EXPERIMENTAL_CAPABILITIES.md` lists NATURaL workflows as experimental. Both can be true only if `G_natural` is inert unless a validated path filled it. On this tip DualAssembly **does** fill `FA->natural_deltaG` on ordinary protein docks (§5). That is a documentation/firewall split, not a PoseHelix leak.

---

## 4. Science-claim firewall

### 4.1 What works

| Control | Evidence |
|---------|----------|
| Public README badge | `Astex-85-unverified \| pending receipt` (`README.md`) |
| CF naming contract | `docs/thermodynamics.md` DatasetRunner CSV table; `docs/KNOWN_LIMITATIONS.md` |
| Softβ S1 default OFF | `LIB/ProtocolConfig.h` `election_shannon_free_energy{false}`; `tests/test_protocol_config.cpp` |
| THERMO_SCORE does not elect | `enforced_in_final_election=0` in `gaboom.cpp`; README protocol-flags table |
| Shannon GA monitor default OFF | `FLEXAIDDS_USE_SHANNON`; DatasetRunner comments at `DatasetRunner.cpp` ~7371 |
| Source firewall tests | `tests/test_thermo_claim_firewall.py` (forbidden phrases; `proxy_only` runtime strings); `python/tests/test_cf_naming_clarity.py`; `tests/test_check_published_astex_rates.py` (blocks 91.8% / 94.1% etc. without withdrawal language) |
| Claim aggregator | `scripts/aggregate_claim_metrics.py` + `tests/p0_claim_contract/test_fixed_denominator.py` (N=**85**, missing targets = failures) |
| DOI cache | `scripts/check_dois.py` + `tests/fixtures/doi_resolution/doi_resolution_cache.json` (CI `repo_hygiene`) |

### 4.2 Residual misread surfaces (not silent ranking bugs)

**P0 — `predicted_dG` name vs CF fallback**

`LIB/DatasetRunner.h` documents the historical header. Assignment:

```text
result.predicted_dG = have_free_energy ? free_energy_F
                    : (best_dG != 0) ? best_dG : best_cf;
```

When Post-GA StatMech is off, a ΔG-shaped column can hold CF. Affinity Pearson code then does `-predicted_dG / 1.3636` and comments “proxy score, monotone rescale only” (`DatasetRunner.cpp` ~9688–9706). That is disciplined in-source and still the sharpest **plot-level** claim hole.

**P0/P1 — `best_score`, `G_bind`, `thermo_G_bind`**

`best_score` is elected CF (`DatasetRunner.h`). Logs print `G_bind=` with `claim_validity=proxy_only energy_domain=cf_arbitrary_units` (`gaboom.cpp`). Identifiers remain claim-hostile if tags are stripped.

**P1 — Engine Softβ when TEMPER>0**

Product default: `classic_entropy_ranking: true`, `force_cf_rank_emission: false` (`LIB/config_defaults.h`). Arm B (FO + TEMPER 21) elects by ACF on CF a.u. That is **documented classic FlexAID**, not DatasetRunner Softβ S1, and still not kcal ΔG_bind. Campaign logs that omit “CF a.u. soft-β” will be misquoted.

**P1 — Induced-fit language**

CCBM receptor entropy in the **bound** ensemble is firewalled (`LIB/BindingMode.h` S_receptor high → induced-fit / population shift; `tests/test_ccbm.cpp`). Marketing/figure strings can still say “Fully Flexible Induced-Fit Docking” (`python/flexaidds/figures.py`). Keep CCBM estimands and poster language apart.

**P2 — getenv-only scoring knobs**

`METHODOLOGY.md` §0.2: `permeability`, `pb_pocket_weight`, `pb_clash_weight` are LOST unless a sidecar records them. Completing a run without that sidecar is not reproducible scoring.

### 4.3 Flag cheat sheet (runtime)

| Knob | Default | Ranking? |
|------|---------|----------|
| `FLEXAIDDS_SOFTBETA_ELECTION` | OFF | DatasetRunner S1 reorder only if ON |
| `FLEXAIDDS_THERMO` | OFF | Ledger / CSV |
| `FLEXAIDDS_THERMO_SCORE` | OFF | Print only; `enforced_in_final_election=0` |
| `FLEXAIDDS_USE_SHANNON` | OFF | Diagnostic; must not enter CF/fitness/election (contract) |
| `FLEXAIDDS_HVIB` | ON unless `=0` | Vib channel / validator, not Astex election by itself |
| `FLEXAIDDS_SEED_ELITISM` | OFF | Oracle ceiling if 1 |
| `FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE` | OFF | DualAssemblyRunner diagnostic mixture only |
| `FLEXAIDDS_POSE_HELIX_THERMO_REWRITE` | **does not exist** | Flag is `NATURaLConfig::enable_pose_helix_rewrite` (false) |
| `--folded` / `advanced.assume_folded` | **false** | **Opt-out** of DualAssembly growth |

---

## 5. NATURaL / DualAssembly / PoseHelix / PoseLocal

### 5.1 Completeness

| Piece | Path | Status |
|-------|------|--------|
| DualAssemblyEngine | `LIB/NATURaL/NATURaLDualAssembly.{h,cpp}` | Growth loop, partial CF, Shannon, bursts, nucleation, TM; omittable only via compile macro `FLEXAIDS_OMIT_DUAL_ASSEMBLY_ENGINE` on the **standalone CLI / some tests**, not FlexAID |
| DualAssemblyRunner | `DualAssemblyRunner.{h,cpp}` | Injected GA callbacks; PoseLocal diagnostic hook |
| Ribosome / RNAP ODE | `RibosomeElongation.{h,cpp}` | Implemented |
| Translocon | `TransloconInsertion.{h,cpp}` | Hessa 2007; gated by `has_membrane_topology` |
| PoseHelix | `PoseHelixThermoRewrite.{h,cpp}` | Pure function complete; sidecar only |
| PoseLocal | `PoseLocalThermoRewrite.{h,cpp}`, `NnMotifTables.h`, `data/*.tsv` | Pure function complete; `{0,0}` labelled gaps (`TODO(burgundy)`) |
| CLI | `dual_assembly_main.cpp` | **Synthetic backend only**; `--real-ga` prints refuse and exits 2 |
| DISCO | `disco_natural.hpp` | Header comments; not called from `run()` |

`DualAssemblyEngine::run()` explicitly does **not** apply PoseHelix (`NATURaLDualAssembly.cpp` ~879–882).

### 5.2 Default-OFF vs default-ON (the split)

**Rewrites — OFF (as advertised)**

- `NATURaLConfig::enable_pose_helix_rewrite = false`
- `DualAssemblyConfig::enable_pose_local_thermo_rewrite = false`
- Env `FLEXAIDDS_POSE_LOCAL_THERMO_REWRITE` opt-in
- Docs correctly say rewrites do not feed `G_natural` / Astex / FlexADS contracts

**Growth — ON (not obvious from “experimental workflow” language)**

`auto_configure` (`NATURaLDualAssembly.cpp` ~172–178):

```text
// Enable for any receptor with residues.
if (n_residues > 0) {
    cfg.enabled                 = true;
    cfg.co_translational_growth = true;
```

Caller (`gaboom.cpp` ~1920–1942): if `!FA->assume_folded` and ligand atom range is valid, construct `DualAssemblyEngine`, `engine.run()`, then `FA->natural_deltaG = engine.final_deltaG()`.

`assume_folded` defaults **0 / false** (`top.cpp`, `config_defaults.h`, `config_parser.cpp`). There is no `FLEXAIDDS_NATURAL=0`.

`BindingMode::compute_energy()` adds `nat_dg` on both classic Softβ and StatMech branches (`BindingMode.cpp` ~541–564). Comments at `get_thermodynamics()` (~587–591) correctly note the term is pose-independent and **cancels in intra-run ΔΔG**. Consequences:

- Elected **rank order within one dock** should not flip solely from this constant.
- Absolute `compute_energy()` / breakdown `G_natural_kcal_mol` / `has_natural` **do** change.
- `compute_partial_cf()` temporarily mutates `FA_->num_optres`, calls `vcfunction`, then restores (`NATURaLDualAssembly.cpp` ~905–924). Residual side-effect risk if `vcfunction` caches more than CF.
- Wall-clock: growth steps scale with residue count on **every** unfolder protein dock.

This is **P0 for claim contamination / experimental bleed onto the default path**, not P0 for “PoseHelix silently reranks Astex.”

### 5.3 Missing GA atom-coord hooks

Documented in `docs/POSE_HELIX_THERMO_REWRITE.md` and `docs/POSE_LOCAL_THERMO_REWRITE.md`:

- No GA → `ReceptorNtCoord` / `LigandPose` snapshots
- `DualAssemblyRunner` does not call PoseHelix (RMSDs, not atom coords)
- PoseLocal needs pre-labelled `pose_rewrite_elements` + `pose_rewrite_poses`
- No Arrhenius `k_fold` coupling of rewritten ΔG
- Aβ42 track (`docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md`): canonical PDB **5OQV** or **2NAO**; sequence given; scientific nucleation match is a **plan**, not a deposited rate
- `docs/DUAL_ASSEMBLY_COTRANSLATIONAL.md` §8 still cites `tests/test_dual_assembly_runner.cpp` — **file absent**; coverage lives in `tests/test_nascent_chain_scheduler.cpp`

Tests that do exist: `test_pose_helix_thermo_rewrite.cpp`, `test_pose_local_thermo_rewrite.cpp`, `test_natural.cpp`, `test_nascent_chain_scheduler.cpp`. No GoogleTest that runs `DualAssemblyEngine` against real docked atoms.

### 5.4 Citation hygiene

| Citation | In-repo locator | This-session check |
|----------|-----------------|-------------------|
| Zhao, Zhang, Chen, *J. Chem. Phys.* **135**, 245101 (2011) | doi:10.1063/1.3671644 | Cache: registered; title “Cotranscriptional folding kinetics of ribonucleic acid secondary structures”; **OK** |
| Zhao et al., *J. Phys. Chem. B* **115**, 3987 (2011) | Repeated in docs/headers; **DOI `10.1021/jp109255g` in `RibosomeElongation.cpp:5`** | Cache (`tests/fixtures/doi_resolution/doi_resolution_cache.json`, checked_utc 2026-09-11): that DOI is **Werner et al., *J. Phys. Chem. C* 115, 5063–5072 (2010)** — gold nanoparticles. **Wrong paper.** `check_dois.py` only asserts the identifier *resolves*, so CI stays green. |
| Xia 1998 | doi:10.1021/bi9809425 | Cache present; used for RNA NN |
| SantaLucia 1998 PNAS Table 2 | doi:10.1073/pnas.95.4.1460 | DNA table; tests assert not aliased to RNA |
| Scholtz 1991 | doi:10.1073/pnas.88.7.2854 | Helix–coil |
| Meier & Seelig 2008 | doi:10.1021/ja077231r | Labelled experimental sheet midpoint |

Docs correctly **warn against mixing JCP vs JPCB roles**. They do not catch the fabricated/wrong DOI on the protein-side Zhao citation. This audit does **not** invent a replacement DOI.

“FlexADS” in PoseLocal/experimental docs means **Astex Diverse Set claim contracts**, not a sibling `FlexADS_Benchmarks` repository (no matches for `FlexADS_Benchmarks` / `FlexAID_Benchmarks` at this SHA).

---

## 6. Benchmark inventory

### 6.1 Astex Diverse (ADS) / “FlexADS”

| Asset | Path | Maturity |
|-------|------|----------|
| Structures | `benchmarks/astex_diverse/astex_diverse/<PDB>/` | In-tree |
| Canonical YAML | `benchmarks/datasets/astex_diverse.yaml` (`docking_mode: self_docking`) | Identity-gated |
| SHA CSV | `benchmarks/datasets/astex_diverse_sha256.csv` | Frozen checksums |
| JSON pair list | `benchmarks/datasets/benchmark_astex_native_85.json` | Campaign entry |
| Frozen claim roster | `benchmarks/protocols/astex85_target_manifest.json` **N=85**, includes **2HR7** | Aggregator denominator |
| 2HR7 science exclusion | `benchmarks/protocols/science_exclusions.md` wants **N=84**; roster file `state/astex85_codes_84.txt` | **File missing** |
| Expected poses | `benchmarks/astex_diverse/expected/README.md` only | No metric snapshot |
| Campaign `result.csv` | (none under workspace) | Not deposited |

**Metrics (normative, not rates):**

- **S1** — elected top-1 in-place RMSD ≤ 2.0 Å (diagnostic)
- **S2 / `success_pb`** — S1 **and** PoseBusters on the **same** elected pose (`pb_backend=bust_cli`)
- **STRICT / `claim_ready`** — S2 + tENCoM/Eigen OK + pose SHA pins + `protocol_claim_eligible=1` + `seed_echo=0` + `native_pose_seeded=0` + matrix pin (`admission_metrics_contract.md`)
- **S3 / BCR** — pool ceiling — never headline success
- Classic 3Dsig red-bar is **S_top10** (any of ranks 0..9), 10 sims × 2e6 evals, bootstrap median — **do not mix** with JCIM 2015 top-1 (`astex_diverse.yaml` comments)

**Numbers allowed to cite from this checkout (external / withdrawn / aspirational — never as live FlexAIDdS STRICT):**

| Number | Meaning | Source |
|--------|---------|--------|
| 0.452 | FlexAID 2015 JCIM Table 2 Astex native FLRP **top-1** | `astex_diverse.yaml` `published_baselines` |
| 0.66 / 0.69 / 0.78 / 0.82 / 0.88 | 3Dsig 2017 **top-10** bootstrap medians (flexaid / flexaid_ds / flexx / vina / rdock) | `deck_2017_baselines` |
| 0.385 | FlexAID 2015 Astex non-native FLRP top-1 | `astex_nonnative.yaml` `published_baselines` (file **deprecated** — see §6.2) |
| 0.220 | HAP2 native FLRP top-1 | `hap2.yaml` |
| 0.70 | Aspirational CI / YAML `expected_baselines.docking_power_top1` | Same YAML; labelled **not measured** |
| 80/85 = 94.1% | **Withdrawn** oracle ceiling | `REPRODUCIBILITY.md` |

`tests/p0_claim_contract/test_fixed_denominator.py` asserts aggregator **N=85**. `METHODOLOGY.md` §0.2 audit-repair contract: an 84-target execution list does not change the claim denominator. `science_exclusions.md` disagrees (2HR7 out → 84). Until that is resolved, **do not publish N/84 or N/85 interchangeably**.

### 6.2 Astex non-native

`benchmarks/datasets/astex_nonnative.yaml` header (2026-09-11): **DO NOT USE FOR NEW RUNS**. Measured defects: roster is largely **CROSS_PROTEIN_INVALID**; C++ `astex_nonnative_targets()` is a second independent defect (only 24/63 native codes are Astex Diverse). Successor: `astex_nonnative_v2.yaml`. Published 0.385 remains a **paper** number, not a receipt for v2.

### 6.3 Other tracks

| Track | YAML / dir | Claim posture |
|-------|------------|---------------|
| CASF-2016 | `casf2016.yaml` | Preliminary / pending receipt (`docs/BENCHMARKS.md`) |
| ITC-187 | `itc187.yaml`; Pearson **0.93**, RMSE **1.4**, ranking **78%** labelled **target / preliminary** | Not a live claim |
| DUD-E / psychopharm / HAP2 live FlexAIDdS | various | Preliminary |
| Grand synthetic Ξ | `benchmarks/grand_synthetic/` | Exact math; T=298 K, c°=1 M; “Do not claim real ΔG” |
| Competition | `benchmarks/datasets/competition_example.yaml` | Harness; ΔΔG_ref **1.364** kcal/mol at 298 K is a **literature-inspired** fixture, not a dock |
| DualAssembly Aβ42 | docs + `scripts/run_dual_assembly_cotranslational.sh` | Experimental; CLI `--real-ga` unimplemented |
| Smoke | `benchmarks/smoke/` | Build/ctest/pytest replay — not docking power |
| Poster CSVs | `benchmarks/astex_repro/poster_metric_*.csv` | Historical/local; `score_offline.py` documented over-permissive (`METHODOLOGY.md` §0) |

### 6.4 CI benches vs campaign

| Harness | Workflow | What it is |
|---------|----------|------------|
| PR sanity dock | `benchmark-tier1.yml` | ~4 Astex targets; seed `FLEXAIDDS_TIER1_SEED=20260816`; `SEED_ELITISM=0`; gate vs aspirational 0.70; PoseBust-only PRs can skip via `tier1_pr_needs_dock.py` |
| Nightly/full | `benchmark-tier2.yml` | Push + weekly cron **force `--dry-run`** until Zenodo records exist; MPI job `if: false` |
| Skill dry-run | `ci.yml` `flexaid_docking_datasetrunner` | `--dataset astex_diverse --tier 1 --dry-run --resume --package` |
| Claim firewall | `ci.yml` `repo_hygiene` | `validate_thermo_claims.py`, published-rate tests, aggregator, skill validate |

`METHODOLOGY.md` §0.1: **CI gate and campaign are not the same experiment** (permeability, normalize_area, clash ratio, coarse_init, MIF, retained poses, timeout accounting). Do not compare an RMSD from `python -m benchmarks.run` to DatasetRunner `result.csv` without stating those divergences.

**MIF drift to re-check before quoting §0.1 literally:** `config_defaults.h` now sets `seeding.mif_enabled: true` (authorised comment, 2026-08-02). `config_parser.cpp` fallback if the key is **absent** is still `false`. DatasetRunner **hardcodes** `"mif_enabled": true` in `dock_config.json` (~7244). Campaign path is MIF ON. Whether hosted Tier-1 still matches the “compiled-in MIF OFF” paragraph depends on whether that job writes JSON defaults — treat §0.1 as **possibly stale** on MIF.

### 6.5 Hardware microbench tables

`docs/BENCHMARKS.md` reports Shannon stack speedups including CUDA A100 **3,575×** and Metal M2 Ultra **412×**, with “re-run `benchmark_dispatch` before quoting.” `docs/USERGUIDE.md` repeats **3575×/412×** as a tip without that caveat. No receipt package for those factors was found in this checkout. Treat as **unreceipted microbench language**, not a Core 1.0 performance claim.

---

## 7. Reproducibility

| Piece | Status |
|-------|--------|
| `docs/REPRODUCIBILITY.md` maturity levels | Replayable / preliminary / published-external — sound |
| `REPRODUCIBILITY.md` blind script | `scripts/reproduce_astex85.sh` defaults `SEED_ELITISM=0`; `--oracle-ceiling` labelled |
| `RUN_RECEIPT` contract | `docs/run-uniformity/RUN_RECEIPT_CONTRACT.md`; fixtures only under `tests/fixtures/receipts/` |
| Seed | Canonical `FLEXAID_SEED`; METHODOLOGY determinism seed **12345**; `FLEXAIDDS_SEED_BASE` is **not** that seed |
| Matrix | `MC_st0r5.2_6.dat`; skill pin MD5 `72d7c7396702331d96ff12d18f831796` (`AGENTS.md` / flexaidds skill) |
| DoF budget | METHODOLOGY §0: **2000 gen × pop 1000** product 2e6; AGENTS.md: claim runs **fix generations, scale population** via `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1`. Skills must cite METHODOLOGY, not fork |
| Restarts | Published Astex protocol 10; `ops/reference_config.env` pins `PARALLEL_RESTARTS=0` but **not** `FLEXAIDDS_RESTARTS` (compiled default 5) — METHODOLOGY §0.2 trap |
| iCloud / fleet | Live OUT **local-first** (`docs/ICLOUD_BENCHMARK_STORAGE.md`); iCloud thin mirror; Fleet experimental (`docs/FLEET_BENCHMARKS.md`, `docs/EXPERIMENTAL_CAPABILITIES.md`) |
| Windows replay | `REPRODUCIBILITY_WINDOWS.md` → WSL2; not native MSVC engine CI |

A claim is repository-reproducible only with bundle layout under `benchmarks/` (README, manifest, run, expected, environment). Astex **expected/** is not populated. Therefore Astex docking power cannot be “replayable from repository artifacts” today.

---

## 8. Correctness risks

### 8.1 Experimental / optional bleed into ranking

| Mechanism | Default dock? | Elects? |
|-----------|---------------|---------|
| DualAssembly → `natural_deltaG` → `compute_energy` | **Yes** unless `--folded` | Constant: **order unchanged**; absolute scores change |
| PoseHelix / PoseLocal | No | No |
| THERMO_SCORE | No | No |
| DatasetRunner Softβ S1 | No | Only if env ON |
| CMA-ES | No (`FLEXAIDDS_SEARCH=cmaes`) | Search change if ON |
| Classic Softβ ACF | **Yes if T>0** | Yes (CF a.u.) |
| SMFREE fitness | **Yes** (ProtocolConfig / DatasetRunner default) | Changes GA *selection* on CF; stagnation tracks CF (`gaboom.cpp` ~898–899). Historical comment: overflow to fitness cap 1000 made GA a random sampler (`DatasetRunner.cpp` ~7380–7384) — treat as **known sharp edge**, verify on current CF scale before claiming SMFREE = “entropy-aware optimization” |
| JSON MIF | Defaults file **true**; parser missing-key **false**; DatasetRunner **true** | Seeding / poses (measured 1MQ6 comment in `config_defaults.h`) |
| `reference_ligand.seed_fraction` parser default **0.25**, `pose_seed_enabled` default **true** | If that JSON block is applied without zeros | Claim-hostile vs blind protocol. DatasetRunner comments say campaign emits `seed_fraction: 0.0` — do not assume every CLI path does |

### 8.2 UB / buffers (legacy `.inp`)

Still live on untrusted / long config lines (not closed by `tests/test_buffer_safety.cpp`, which covers GIST helpers):

- `LIB/read_input.cpp` `strcpy` into `clf_file` (LOCCLF); `sscanf %s` into `field[7]`
- `LIB/gaboom.cpp` `sscanf %s` into `field[9]`; unbounded `reflig_file`
- `docs/SECURITY_HARDENING_ROADMAP.md` still marks H-1…H-7 PENDING (April 2026 inventory). `KNOWN_LIMITATIONS.md` already warns not to treat that roadmap as closure.

**P0 for untrusted configs; P2 for signed in-repo benchmark bundles.**

### 8.3 OpenMP

- PR #492 clamped OpenMP team size in `LIB/tENCoM/tencm.cpp`. **`tencm_cpu_fallback.cpp` still does `#pragma omp parallel` without `num_threads(n_threads)`** against vectors sized from `omp_get_max_threads()` (~33–54, Hessian ~148+). P0 **if** fallback+OpenMP is used with a larger team than the pre-sized buffer (nested parallel / `omp_set_num_threads` after capture).
- `gaboom.cpp` / `Vcontacts.cpp` use `thread_local` / `default(none)` — residual complexity, not a new finding.

### 8.4 NaN

- StatMech empty/non-finite: throws (`statmech.cpp`) — good.
- Softβ empty / no finite members: `G = +∞` (`SoftBetaFreeEnergy.h`) — fail-soft; can make a mode unrankable.
- CF (`vcfunction` / `cffunction`): no `isfinite` guards found — bad coords → NaN CF → Softβ drop.
- `FLEXAIDDS_FIXED_ORDER_LSE` default **false** (`log_sum_exp.h`) — OpenMP reduction order can change ln Z (reproducibility P1).
- DualAssemblyRunner empty engine diagnostic can emit NaN (`DualAssemblyRunner.cpp`).

---

## 9. CI / release gates vs `SUPPORT_MATRIX.md`

`docs/CI_RELEASE_GATES.md` matches `.github/workflows/ci.yml` jobs reasonably well:

| Job | Blocking? |
|-----|-----------|
| `pure_python_results` | Yes |
| `cxx_core_build` linux-gcc / linux-clang / mpi / asan + ctest | Yes (avx512 `allow_failure`) |
| `macos_cpu_tests` + ctest, Metal OFF | Yes, no `continue-on-error` |
| `macos_metal_compile_smoke` | Yes; exit 2 (no metalc) → soft-skip success |
| `python_bindings_smoke` | Yes (Linux); Windows bindings historically `allow_failure` per CI_RELEASE_GATES |
| `typescript_claim_firewall` | Yes |
| `repo_hygiene` | Yes |
| `flexaid_docking_datasetrunner` | Yes (dry-run) |

**Splits vs support matrix / CLAUDE.md:**

| Doc | Reality at this SHA |
|-----|---------------------|
| SUPPORT_MATRIX Linux GCC ≥ 10 / Clang ≥ 10 | `CMakeLists.txt` **FATAL** GCC < 14, Clang < 18, AppleClang < 16 |
| SUPPORT_MATRIX Windows MSVC 2022 Supported (CLI) | `ci.yml`: “windows-msvc excluded: deep MSVC C++20/C++26 incompatibilities”; CMake caps MSVC to C++20 and warns only `_core` is expected to build |
| SUPPORT_MATRIX Python 3.9–3.11 | CI uses 3.11; 3.9 not gated here |
| Metal experimental | Consistent with CI_RELEASE_GATES; do not claim Metal production without green `metal-self-hosted.yml` |

`PRODUCT.md` missing means Core 1.0 “product boundary” is only the three capability docs — and those disagree on `G_natural` (§3.5, §5.2).

---

## 10. Prioritized findings

### P0 — claim contamination / silent science-path mismatch

1. **DualAssembly growth default-ON** on FlexAID `GA()` epilogue (`auto_configure` + `gaboom.cpp` + `BindingMode::compute_energy`). Experimental NATURaL writes `G_natural` on ordinary protein docks. Intra-target order likely unchanged; ledger and “experimental = off the default path” story are wrong. Opt-out is `--folded` only.  
2. **`predicted_dG` CF fallback** under a ΔG header (`DatasetRunner.cpp`). Highest practical misclaim vector for plots/ITC Pearson.  
3. **No FlexAIDdS Astex docking-power receipt** while aspirational 0.70 and USERGUIDE hardware × factors still sit in docs. Public vehicles correctly say unverified; **do not** fill the gap from memory.  
4. **Astex N=85 vs N=84** (`astex85_target_manifest.json` + aggregator tests vs `science_exclusions.md` + missing `state/astex85_codes_84.txt`). Publishing either rate without resolving 2HR7 is a denominator cheat waiting to happen.  
5. **`astex_nonnative.yaml` is scientifically invalid for new runs** (in-file deprecation). Easy to miss if a skill still points at the old slug.

### P1 — correctness / contract bugs

6. **`tencm_cpu_fallback.cpp` OpenMP team OOB** — PR #492 did not land here.  
7. **Legacy `%s` / `strcpy` overflows** in `read_input.cpp` / `gaboom.cpp` GA `.inp` parsers.  
8. **Zhao JPCB DOI `10.1021/jp109255g`** resolves to Werner 2010 JPC C in the committed DOI cache. CI does not fail title mismatch.  
9. **SUPPORT_MATRIX compiler/Windows overclaim** vs CMake + `ci.yml`.  
10. **Five RMSD instruments**; CI vs claim path still differ (`METHODOLOGY.md` §0.0).  
11. **JSON `seed_fraction` / `pose_seed_enabled` parser defaults** are claim-hostile if a path applies `reference_ligand` without zeros.  
12. **SMFREE overflow / “random sampler” historical behaviour** still documented next to the default fitness model — needs a current fail-closed test that SMFREE remains a function of CF, not a flat 1000.  
13. **METHODOLOGY §0.1 MIF OFF on CI** may be stale vs `config_defaults.h` MIF true.  
14. **`PRODUCT.md` missing**; DualAssembly doc cites a missing `test_dual_assembly_runner.cpp`.  
15. **`FLEXAIDDS_POSE_HELIX_THERMO_REWRITE` env does not exist** — operators following experimental-docs naming will think they toggled helix rewrite.

### P2 — gaps

16. PoseHelix/PoseLocal: no real-GA atom snapshots; `{0,0}` gaps; `--real-ga` refuse.  
17. Tier-2 dry-run until Zenodo; smoke `expected/` empty.  
18. Hardware 3575×/412× unreceipted; USERGUIDE states them as fact.  
19. Dual Python binding sources; grand-canonical Python `HAS_GRAND_BINDINGS = False`.  
20. Buffer-safety tests vs roadmap PENDING list.  
21. `FLEXAIDDS_FIXED_ORDER_LSE` off by default.  
22. Figure/marketing induced-fit copy outside CCBM tests.

### P3 — nice-to-have

23. Rename CSV `predicted_dG` → `ensemble_F_or_cf_fallback` in a **new** schema version (keep old header for live campaigns).  
24. Rename log field `G_bind` → `proxy_G_bind_cf_units`.  
25. Index `docs/audit/` so this file is discoverable next to 2026-08-16.  
26. Align `docs/BENCHMARKS.md` ITC/CASF “target” tables with the same “do not quote” banner as README Astex.

---

## 11. Recommended next five actions (ranked)

1. **Fail-closed DualAssembly on the docking default path** (follow-up PR, not this one): require `--folded` inverted, or `FLEXAIDDS_NATURAL=0` default, or skip `auto_configure` unless an explicit NATURaL flag is on. Keep PoseHelix/Local OFF. Add a ctest that `FA->natural_deltaG == 0` on a tiny protein dock without the flag. Update `VALIDATED_CAPABILITIES.md` so `G_natural` is not “validated additive correction” until that path is labelled.  
2. **Stop shipping ΔG-shaped CF:** emit an explicit `cf_fallback` / `has_free_energy` column next to `predicted_dG`, and fail `test_cf_naming_clarity` if assignment can copy `best_cf` without a sibling flag.  
3. **Pick one Astex denominator** (85 including 2HR7-as-failure vs 84 excluding 2HR7), commit the missing roster file or delete the 84 language, and make `science_exclusions.md` and `test_fixed_denominator.py` agree. Do **not** rewrite any success *numerator*.  
4. **Fix the Zhao DOI comment** in `RibosomeElongation.cpp` (remove `10.1021/jp109255g` or replace only after a human finds the real resolving DOI). Extend `check_dois.py` to flag DOI/title mismatch against nearby citation strings.  
5. **Clamp OpenMP in `tencm_cpu_fallback.cpp`** the same way as `tencm.cpp` (#492), and add a test that a nested/larger team cannot index past `thread_contacts.size()`.

Then, separately: a **blind** Astex campaign with `RUN_RECEIPT` + `result.csv` under the admission contract — the only way a FlexAIDdS rate becomes claim-ready. Do not use `astex_nonnative.yaml` for that work.

---

## 12. What this audit did not do

- Did not run Astex, PoseBusters, or DualAssembly Aβ42 campaigns.  
- Did not rebuild FlexAID or run full `ctest` (documentation PR; firewall tests run on the doc delta).  
- Did not independently resolve a correct Zhao JPCB DOI.  
- Did not change SMFREE, MIF, TEMPER, Softβ S1, or Astex YAML numbers.  
- Did not treat `benchmarks/astex_repro/poster_metric_*.csv` as STRICT evidence.

---

## 13. Claim-ready TODAY vs experimental (copy table)

**Claim-ready / allowed from this repository today**

- Core CLI / Python execution surfaces as *software*, not as a docking-power percentage.  
- StatMech ledger mathematics (log-sum-exp, empty-ensemble throw) for **labelled analysis**.  
- CF ≠ ΔG language in README / KNOWN_LIMITATIONS / thermodynamics CSV contract.  
- Published **external** FlexAID 2015 / 3Dsig 2017 figures cited from YAML (with statistic type).  
- Grand-synthetic Ξ arithmetic fixtures.  
- “Astex-85 unverified / pending receipt” as the live FlexAIDdS docking-power statement.

**Experimental or not claim-ready**

- Any FlexAIDdS Astex S1/S2/STRICT percent, including 0.70, ~70%, 81.2%, 25.9%, 94.1%.  
- NATURaL DualAssembly (including default-ON growth), PoseHelix, PoseLocal, Aβ42 nucleation match.  
- Softβ S1 election, THERMO_SCORE as a science result, Shannon-in-fitness.  
- ITC Pearson 0.93 / CASF r≈0.88-class tables.  
- Swift, Fleet, iCloud live GA, CUDA/Metal **runtime** docking, AVX-512.  
- Metal “production GPU” without self-hosted full gate.  
- `G_bind` / `predicted_dG` as experimental ΔG_bind.  
- RMSD-only “docking success” without PoseBusters on modern/PB tables.

**One-line verdict:** The firewall around *quoting* Astex rates is in better shape than in the 2026-08-16 audit. The firewall around *running experimental thermo on the default GA epilogue* (DualAssembly growth → `G_natural`) and around *ΔG-named CF columns* is not. Do not claim FlexAIDdS docking power from this SHA.
