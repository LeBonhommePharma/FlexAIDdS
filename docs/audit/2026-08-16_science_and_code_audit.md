# FlexAIDdS — Full Science and Code Audit

**Date:** 2026-08-16
**Tip audited:** `5c891e6eb033455bdf0ba03356b292a49a4835d9` (`main`)
**Auditor:** Cursor Grok 4.6 (cloud agent), first-hand source inspection
**Companion prior audits:** `docs/audit/26h-swarm/SCIENCE_AUDIT.md` (2026-07-15),
`docs/audit/AUDIT_2026-07-23.md`, `docs/audit/2026-08-13_multithreading_reproducibility_science_robustness.md`,
`CODEBASE_REVIEW_2026-08.md`
**Mode:** Diagnostic. No engine, ranking, or thermodynamic behaviour was changed.

This document re-verifies the July/August science and code findings against the
current tip and records what is still open. Every factual statement below is
anchored to a file read in this session. Numbers that would require a live Astex
campaign are **not** restated as current results.

---

## 1. Executive scientific verdict

| Question | Answer on this tip |
|----------|--------------------|
| Does the **math** of soft-β \(\tilde G\) on CF make sense? | **YES** — one shared identity in `LIB/SoftBetaFreeEnergy.h` |
| Is that math **one objective** across cluster / BindingMode / DatasetRunner? | **MOSTLY** — formula and local \(Z\) now agree; **CF operand still splits** (`evalue` vs `app_evalue`) and FO emission still refuses the re-scored pose |
| Can you claim “3Dsig Shannon free energy election” on **defaults**? | **NO** — DatasetRunner Softβ S1 is correctly **OFF** by default |
| Can you claim physical binding free energy \(\Delta G\)? | **NO** (correct if language is disciplined; several user-facing docs still are not) |
| Can you quote a current Astex-85 success rate from this repository? | **NO** — README is honest (“unverified”); `docs/BENCHMARK.md` and `REPRODUCIBILITY.md` still publish 91.8% / 94.1% without a live receipt on this tip |
| Is FO production policy (single literature MinPts + dual-suffix packaging) on main? | **YES** — dual-suffix enumeration is on `main` with tests |
| Is DoF budget (scale pop, fix gen) the right claim contract? | **YES** |
| Is the CI gate the same experiment as the campaign path? | **NO** — `METHODOLOGY.md` §0.1 still holds |

**One-line verdict:**
The scientific *idea* (search on the Voronoi CF proxy; rank modes by local soft-β
\(\tilde G=H-TS\) on CF; density-based FO modes; fixed-gen pop×DoF) is coherent and
much closer to a single implementation than it was on 2026-07-15. The **claim
surface is not**. User-facing benchmark pages still advertise rates this tip
cannot defend; the 94.1% reproduction script is still a 90% native-seeded
oracle; CI and campaign still dock different physics; and several getenv-only
scoring knobs still vanish when the shell exits.

Treat entropy / ΔG / Astex-headline language from this tip as **unproven**
unless a current `result.csv` + `RUN_RECEIPT.json` for that exact binary is in
hand (AGENTS.md deception-proof contract; METHODOLOGY.md §0–§3).

---

## 2. Scientific layer model (must not be conflated)

| Layer | What it optimizes / reports | Units | Elects poses? |
|-------|----------------------------|-------|----------------|
| **L1 Search** | GA fitness = Voronoi **CF** (contact-function scoring proxy) | a.u. | No (explores) |
| **L2 Cluster** | Binding modes (CF / FO / DP) | geometry + CF | Groups poses |
| **L3 Rank (classic / 3Dsig)** | Soft-β \(\tilde G=\tilde H-T\tilde S\) on **CF** | CF a.u., \(T\) is a score temperature (\(\beta=1/T\)), **not** \(k_B T\) | Yes, when TEMPER>0 (engine) or Softβ S1 ON (DatasetRunner) |
| **L4 Thermo ledger** | StatMech / BindingMode physical F, H, S, Cv (+ vib if a real eigenvalue channel exists) | kcal when \(k_B\) is used | No, unless an explicit ranking flag is on |

**AGENTS.md contract:** never sell L1/L3 as experimental \(\Delta G\). L3 is a
ranking objective on a scoring proxy.

Default DatasetRunner S1 is **L1 min-CF** (`FLEXAIDDS_SOFTBETA_ELECTION=0`).
Arm B classic FlexAID is engine FO + TEMPER 21 ACF — **not** DatasetRunner Softβ.

---

## 3. What closed since the 2026-07-15 science audit

Re-read against `5c891e6e`. These P0 items from `docs/audit/26h-swarm/SCIENCE_AUDIT.md`
are **closed on main**.

| July finding | Status now | Anchor |
|--------------|------------|--------|
| S1 — DatasetRunner `soft_T` hardcoded 298, never reads TEMPER | **CLOSED** | `LIB/DatasetRunner.cpp` now takes `dock_temperature_K`; env → dock → 298. Shared `soft_beta::resolve_soft_T` |
| S2 — FO BindingMode never writes `.mcf` | **CLOSED** | `LIB/BindingMode.cpp` writes member-CF sidecar after `write_pdb` |
| S3 — BindingMode used **global** \(Z\) | **CLOSED** | `compute_enthalpy/entropy/energy` use local SoftBeta over mode members |
| S4 — `LIB/SoftBetaFreeEnergy.h` missing | **CLOSED** | Header on main; included by cluster, BindingMode, DatasetRunner |
| S5 — Shannon S1 default ON | **CLOSED** | `ProtocolConfig::election_shannon_free_energy` default **false**; both env aliases default OFF |
| FO dual-suffix packaging incomplete | **CLOSED** | `enumerate_emitted_cluster_heads()` on main; gtests in `tests/test_dataset_runner.cpp` |
| Python RMSD prefix truncation | **CLOSED** | `python/flexaidds/dataset_runner/runner.py` returns `-1.0` on shape mismatch |
| `docking_power` dropped failed targets | **CLOSED** | Denominator is attempted targets; sentinels are never successes (`python/tests/test_metrics_docking_power.py`) |
| CI timeout indistinguishable from empty dock | **CLOSED on Python path** | `TimeoutExpired` now increments crash count (`runner.py`); campaign already recorded timeouts |
| SHARESCL 0.20 production default | **CLOSED** | `scripts/generate_flexaid_inp.py` `DEFAULT_SHARESCL = 10.0`; `tests/test_ga_sharescl_default.py` |
| Python RMSD not symmetry-corrected (#365) | **CLOSED as element-blocked Hungarian** | `python/flexaidds/benchmark.py::_symmetry_permutation` (PR #365, `2959d587`). **Not** spyrmsd graph-isomorphism |

Soft-β identity tests exist and are real gates: `tests/test_classic_entropy_ranking.cpp`
(`SoftBetaIdentity.*`, `GEqualsHminusTS`, `GEqualsEminMinusTlnZ`, duplicate-invariant
`free_energy_strict`).

---

## 4. Remaining science findings

### C1 — CRITICAL: user-facing benchmark pages still publish rates this tip cannot defend

`README.md` is correctly fail-closed: Astex-85 badge is `unverified | pending receipt`
and the prose says this repository publishes no receipted success rate.

That honesty is contradicted by two still-live documents:

1. **`docs/BENCHMARK.md`** still titles a section “The 91.8% Record” and tables
   **“FlexAID∆S (current) 91.8% (78/85)”** against Glide/Vina/GNINA. It also
   claims ITC-187 Pearson **r = 0.93** for `G_bind` vs experimental \(\Delta G\),
   and gives reproduce commands that invoke **`scripts/run_dataset.py` and
   `scripts/analyze_affinity.py` — neither file exists** in this tree.
2. **`REPRODUCIBILITY.md`** still publishes **80/85 (94.1%)** RMSD_hungarian < 2.0 Å
   at commit `8196829f…`, with `FLEXAIDDS_NATIVE_SEED_FRAC=0.90` and
   `FLEXAIDDS_SEED_ELITISM=1` in `scripts/reproduce_astex85.sh`. Multiple
   targets in that protocol elect the crystal ligand (0.00 Å). METHODOLOGY.md §0
   forbids reporting seed-elitism / `_INI.pdb` RMSD as the result.
   `reproduce_astex85.sh` then “verifies” against a ±3% band around 94.1 rather
   than measuring a blind rate.

This is the same CRITICAL finding as `CODEBASE_REVIEW_2026-08.md` §1. It has
**not** been retired. A reader who opens `docs/BENCHMARK.md` instead of
`README.md` is sold a current SOTA number.

**Science rule:** until a blind, unseeded, PoseBusters-gated Astex-85
`result.csv` + `RUN_RECEIPT` for *this* tip exists, the only honest headline is
the README’s: unverified.

---

### C2 — CRITICAL: CI gate and campaign path remain different experiments

`METHODOLOGY.md` §0.1 is still an accurate description of the code.

| Knob | CI / `python -m benchmarks.run` | Campaign `DatasetRunner` |
|------|----------------------------------|---------------------------|
| config | compiled-in defaults (`top.cpp`) | writes `dock_config.json` |
| `permeability` | 1.0 (`top.cpp`) | 0.9 |
| `normalize_area` | 0 | true |
| `intermolecular_clash_ratio` | 0.0 | 0.75 |
| `coarse_init.enabled` | OFF | ON (hardcoded) |
| `mif_enabled` | OFF unless JSON/env | ON (hardcoded, `DatasetRunner.cpp`) |
| retained poses | 10 | 50 / restart |

`top.cpp` still returns immediately from MIF init when
`!(FA->mif_enabled || FA->grid_prio_percent < 100.0f)`. The gate therefore still
docks **without building a pocket field**. A green tier-1 tick is not evidence
about campaign physics.

**Do not compare an RMSD from one path to an RMSD from the other** without
stating these divergences. Changing either default is a methodology change and
must land in `METHODOLOGY.md` first.

---

### H1 — HIGH: CF operand is still not one quantity

Two aggregators exist in `LIB/ic2cf.cpp`:

```text
get_apparent_cf_evalue = com+wal+sas+elec+hbond+gist+metal+entropy+pb_clash
get_cf_evalue          = apparent + con + (optional tencom_weight * h_rep)
```

| Writer | `REMARK CF=` | `.mcf` / ranking members |
|--------|--------------|--------------------------|
| `cluster.cpp` (CF clustering) | **re-scored** `emitted_cf` = `get_cf_evalue` | `app_evalue` unless pose-inconsistent, then recomputed apparent |
| `BindingMode.cpp` (FO) | stored `chrom->evalue` | stored `chrom->app_evalue` |

FO **does** call `ic2cf` to rebuild coordinates, then **discards** the re-score
for the REMARK line, with an explicit comment that Vcontacts can return an
uncapped clash at emission (`BindingMode.cpp`). CF clustering did the opposite
fix (1HP0 phantom CF): it **substitutes** the re-scored pose when search CF
disagrees.

DatasetRunner CF-rank-0 (the **default** elector) parses `REMARK CF=` (so
`evalue` / `get_cf_evalue` on FO; re-scored `get_cf_evalue` on CF clustering).
Softβ, when ON, reads `.mcf` (`app_evalue`). BindingMode classic ranking uses
`pose.CF`, constructed from `app_evalue`.

**Scientific conclusion:** even with SoftBetaFreeEnergy.h shared, the three
callers can rank different numbers. Constraint term `cf.con` and optional
`h_rep` are in L1 search CF and in `REMARK CF=`, but not in apparent CF / `.mcf`.
A comment at `DatasetRunner.cpp` still claims `REMARK CF=` is `app_evalue`. That
comment is false.

---

### H2 — HIGH: Grand-canonical `log_Z` still has a zero-sentinel and a CF→kcal fallback

`LIB/DatasetRunner.cpp` (TargetServer registration):

```text
if (result.ensemble_log_Z != 0.0)
    sess.log_Z = result.ensemble_log_Z;
else
    sess.log_Z = -predicted_dG / (kB_kcal * T);
```

`ensemble_log_Z` defaults to `0.0`. A real \(\log Z = 0\) (Z = 1) is therefore
indistinguishable from “not emitted” and is replaced by the proxy.

`predicted_dG` itself is Helmholtz F **when present**, else `best_dG`, else
**best CF** (`DatasetRunner.cpp`). Dividing a CF arbitrary-unit score by
\(k_B T\) in kcal/mol produces a number with no thermodynamic meaning. AGENTS.md
requires real ensemble `log_Z` from `BindingPopulation.get_global_ensemble()`
for GPF paths; this fallback violates that when the ledger line is absent.

Single-ligand canonical docks that never touch TargetServer are unaffected
(pre-GPF bit-identity). Multi-ligand / `--conc` / competition YAML is where
this becomes a science bug.

---

### H3 — HIGH: getenv-only scoring knobs are still LOST provenance

`METHODOLOGY.md` §0.2: anything that goes through `getenv` and nowhere else is
gone when the shell exits. Examples still in that class:

| Knob | Read at | In `dock_config.json`? | In `RUN_RECEIPT` / `ProtocolConfig`? |
|------|---------|------------------------|--------------------------------------|
| `FLEXAIDDS_PB_CLASH_WEIGHT` | `top.cpp` / `config_parser.cpp` (env overrides JSON) | only if DatasetRunner wrote it — **it does not** | **no** |
| `FLEXAIDDS_PB_POCKET_WEIGHT` | same | **no** | **no** |
| `FLEXAIDDS_PB_CLASH_ELECT_WEIGHT` | `BindingMode.cpp` only | **no** | **no** |
| `FLEXAIDDS_SAS_WEIGHT` | DatasetRunner → JSON `sas_weight` | yes, when DR writes config | not as env snapshot |
| `FLEXAIDDS_COM_FLOOR` | engine getenv | **no** | **no** |
| `FLEXAIDDS_VCT_NORM` | presence getenv | DR can emit `vct_normalize_contacts` | ProtocolConfig yes |

`LIB/RunReceipt.cpp` embeds `protocol_config` and hashes, but **no
`scoring_env` object**. `scripts/check_run_receipt.py --require-scoring-env`
looks for keys the writer never emits. Default receipt check does **not**
require `scoring_env`, so campaigns pass the gate without recording the
weights that actually scored the poses.

`ProtocolConfig` still does not carry `pb_clash_weight` / `pb_pocket_weight`.

---

### H4 — HIGH (science-claim): BindingMode vibrational correction is a production no-op

`BindingMode::compute_vibrational_correction()` (`LIB/BindingMode.cpp`) now
documents that the previous eigenvalue channel was a misread of eigenvector
storage, and **returns 0.0 on every real path**, caching that zero. Ranking,
clustering, and output order are unchanged by tENCoM/ENCoM on this function.

DatasetRunner HVIB is a **different** post-GA ligand-ANM diagnostic over emitted
PDBs; comments say eigenvalues never enter CF or selection. `ProtocolConfig`
defaults `hvib_enabled{true}` (`FLEXAIDDS_HVIB=0` to disable), while
`DatasetRunner.cpp` still comments “Gated by FLEXAIDDS_HVIB=1 (default OFF)” —
the comment is wrong; the typed default is ON.

**Do not claim “tENCoM validated the pose ranking” from this tip.** tENCoM as a
STRICT *validator binary* (`ligand_tencom_pose` / `tencom_entropy_diff`) is a
separate artifact path; the in-engine \(-T S_\mathrm{vib}\) term does not move
S1.

`tests/test_binding_mode_vibrational.cpp` still describes itself as validating
“ENCoM-based -T*S_vib integration into BindingMode free energy”. If those tests
pass by constructing a fabricated atom-0 eigenvalue layout, they are the
false-confidence pattern called out in `docs/audit/2026-08-13_test_coverage_quality.md`.

---

### H5 — HIGH: five RMSD instruments remain; METHODOLOGY.md is stale on #365

METHODOLOGY.md §0 still says the in-repo Python metric is **not**
symmetry-corrected and that #365 is open. On this tip #365 **has merged**
(`2959d587`): `compute_rmsd(..., elements=)` runs an **element-blocked
Hungarian** assignment. It is still **not** spyrmsd graph-isomorphism.

The five instruments are still live:

1. **CI / Python gate** — `python/flexaidds/benchmark.py::compute_rmsd` (in-place,
   element-Hungarian when elements are passed).
2. **Offline reference** — `benchmarks/astex_repro/score_reference.py` (spyrmsd).
   Still not wired into CI.
3. **Offline permissive** — `benchmarks/astex_repro/score_offline.py`. Still in
   tree; still documented as over-permissive (1HP0).
4. **Engine Hungarian** — `LIB/calc_rmsd.cpp::calc_Hungarian_RMSD` (groups by
   **SYBYL type**; writes `REMARK RMSD`).
5. **DatasetRunner Hungarian** — `dataset::hungarian_rmsd` (groups by
   **element**; writes `result.csv`).

A cross-check **now exists** (`tests/test_dataset_runner.cpp` `RmsdCrossCheck.*`)
and pins that the two C++ solvers **must disagree** when type and element
partitions differ (split carbon types → dataset 0, engine 10; dummy type
spanning N/O → engine 0, dataset 10). Mixing a REMARK RMSD and a CSV RMSD in
one table is still invalid.

Python `docking_power` uses **`rmsd < 2.0`**. DatasetRunner success uses
**`rmsd <= 2.0`**. AGENTS.md / benchmarking skill say **≤ 2.0**. A pose at
exactly 2.000 Å is a Python miss and a C++ hit.

---

### H6 — MEDIUM–HIGH: protocol defaults vs published Astex protocol

| Axis | Published / METHODOLOGY | Code default |
|------|-------------------------|--------------|
| Restarts | 10 (Astex protocol) | `ProtocolConfig::restarts{5}` |
| Softβ S1 | OFF for claim (correct) | OFF (correct) |
| Arm B TEMPER | 21 (ranking hyperparameter, not Kelvin) | engine CONFIG; DatasetRunner S1 only if Softβ ON |
| Seed elitism | OFF for claim | `ProtocolConfig::seed_elitism{true}`; DatasetRunner **forces false** for AUTONOMOUS / DEFINED_CLEFT / UNSET |
| Eval scale | pop×DoF, gen fixed | `eval_scale_dihedral{1}` (correct) |
| RNG stream fix | documented defect | `FLEXAIDDS_RNG_STREAM_FIX` default **OFF** (parity with frozen numbers) |
| NaN rank guard | F5 still present | `FLEXAIDDS_NAN_RANK_GUARD` default **OFF** |

`ops/reference_config.env` pinning `FLEXAIDDS_PARALLEL_RESTARTS` but not
`FLEXAIDDS_RESTARTS` remains a METHODOLOGY §0.2 footgun: a campaign against it
silently takes 5 restarts, not 10.

---

### M1 — MEDIUM: `is_hydrogen_atom` still misclassifies mercury (and helium)

`LIB/hbond_potential.h`:

```cpp
return atom.element[0] == 'H' ||
       (atom.element[0] == ' ' && atom.element[1] == 'H') ||
       atom.name[0] == 'H';
```

`"HG"` and `"HE"` match clause 1. PDB hydrogens named `1HB` miss clause 3.
`tests/test_hbond_potential.cpp` does **not** mention mercury. The 2026-07-23
audit finding is **still open**. Psychopharmacology ligands are not the usual
Hg hosts, but heavy-atom derivative structures and the element table are.

Virtual-H still collects at most two heavy neighbors; protonated tertiary
amines (the common CNS basic center) still cannot get a correct third-neighbor
direction. Hydroxyl azimuth is still `perp3` frame-dependent. These were
MAJOR/MINOR in July and were not closed.

---

### M2 — MEDIUM: METHODOLOGY.md has drifted from the tree

| METHODOLOGY claim | Tree on this tip |
|-------------------|------------------|
| §0 Python RMSD “NOT symmetry-corrected”; #365 open | #365 merged; element-Hungarian is on |
| §4 ctest expect **11/11** | This configure registered far more than 11 binaries (see §8) |
| cmake path `<host-specific-cmake-path>` (Darwin Homebrew in this session) | Host-specific; Linux CI uses distro/kitware cmake |

The file remains the right *procedural* source of truth for parity / Astex-85
A/B / merge gates. The **instrument table in §0** needs a methods edit (there,
not in a skill fork) so agents stop repeating a closed #365.

---

### M3 — MEDIUM: hardcoded Darwin paths in tests / scripts

`tests/test_campaign_methodology_gates.py::test_g4_2_niche_distance_drives_shipped_cpp_binary`
defaults `SCRATCH` to a Darwin `<temporary-directory>` under the `/var/folders/`
prefix (host-specific UUID path omitted). On this Linux host that raises
`PermissionError` on that Darwin prefix rather than skipping. Same prefix
appears in `scripts/patch_bcr_from_poses.py`. This is an AGENTS.md hygiene miss
(machine-specific absolute path) and a false-red / false-environment test.

---

### M4 — MEDIUM: three C++ test files still never run

On disk but **not mentioned in `CMakeLists.txt`:**

| File | What it would test |
|------|--------------------|
| `tests/test_ga_population.cpp` | GA operators / fitness_stats |
| `tests/test_binding_mode_io.cpp` | BindingMode I/O |
| `tests/test_production_blockers.cpp` | Named “production blockers”; not in any target |

Unchanged since the 2026-08-13 coverage audit.

`tests/test_aggregate_oracle_ceiling.py` still `return`s (not `pytest.skip`)
when `~/flexaidds_results/v43_…` is absent. If collected, that is a silent
**PASS**. It lives under `tests/` next to C++ sources, so default `pytest
python/tests/` may not collect it — dead theater rather than green CI, unless
someone runs `pytest tests/`.

---

### M5 — LOW–MEDIUM: license scan still swallows scancode failure

`.github/workflows/license-scan.yml` pins `scancode-toolkit>=32,<33` (schema
fix vs the August review) and `check_licenses.py` now reads
`license_detections[].matches` with a score floor. The scan step still ends
with `|| true`. If scancode dies before writing JSON, the next step fails
closed (`FileNotFoundError`). If it writes a truncated-but-valid empty-ish
report with files, behaviour depends on contents. The `|| true` is still a
smell; the checker is no longer a total no-op.

---

### M6 — LOW: CMake 3.28.3 cannot configure C++26 + OpenMP on GNU

`cmake_minimum_required(VERSION 3.28)` plus `CMAKE_CXX_STANDARD 26` makes
CMake 3.28.3’s `FindOpenMP` `try_compile` fail:
“requires the language dialect CXX26 … GNU does not support this, or CMake
does not know the flags.” This session configured successfully with **CMake
4.4.2**. Those are the only versions tested here (3.28.3 fail, 4.4.2 pass).
This audit does **not** claim CMake ≥3.31 is sufficient — that was not
build-tested. A host stuck on 3.28.3 cannot configure this tree. Worth
documenting the tested pair in INSTALLATION / CI (out of scope for this
docs-only change).

---

## 5. Code quality (non-ranking)

The August 13 multithreading audit’s F1 (`lazy_thread_rng` stream collapse)
is **still the default path**, now behind `FLEXAIDDS_RNG_STREAM_FIX` default
OFF, with a comment that frozen reference numbers were produced under the
defect. That is the correct parity posture. The defect is not gone; it is
quarantined.

F5 (NaN CF elected rank-0) is likewise gated OFF
(`FLEXAIDDS_NAN_RANK_GUARD`). Production ranking can still elect NaN if a
non-finite CF appears.

FO MinPts single-pass literature composite (`fo_choose_minpts`) remains the
right production rule vs a MinPts ladder.

`-ffast-math` + SIMD-width dispatch still means **not bit-reproducible across
binaries/architectures** (August F4). Claim comparisons must pin binary
SHA256.

`predicted_dG` remains a historical CSV column name. Docs in
`DatasetRunner.h` and `docs/thermodynamics.md` now label it correctly as an
ensemble F estimate / CF fallback. The column name itself will keep generating
overclaims in downstream notebooks.

---

## 6. What you *can* defend scientifically today

If the methods section is precise and receipts exist:

1. **Search** optimizes Voronoi CF (proxy), matrix `MC_st0r5.2_6.dat`.
2. **Modes** from single-pass FO with literature-inspired MinPts, or CF
   clustering, depending on CLUSTA.
3. **Engine emission** at TEMPER>0 uses local ACF / SoftBeta (duplicate-invariant
   unless `FLEXAIDDS_ELECT_LEGACY_ACF=1`).
4. **DatasetRunner S1 default** is min finite head CF (Softβ OFF).
5. **DoF budget** scales population, generations fixed (claim path).
6. **S1/S2/S3** definitions are encoded in the claim aggregator; S2 still
   requires PoseBusters on the elected pose.
7. **Apo strip / ligand integrity / native CF oracle** scripts exist as
   fail-closed prep gates.
8. **Soft-β algebra** \(G = E_{\min}-T\ln Z = H-TS\) is unit-tested.

You **cannot** yet defend, from this tip alone:

- Any numerical Astex-85 success rate as “current”
- “DatasetRunner S1 = 3Dsig Shannon free energy on FO modes” (flag OFF; and
  even ON, CF operand / FO re-score policy still differ)
- “Entropy arm and CF arm differ only by TEMPER”
- “tENCoM changed which pose was elected”
- “Computed thermodynamic binding free energy / ITC r = 0.93”
- Equating a tier-1 CI RMSD with a DatasetRunner campaign RMSD

---

## 7. Prioritized remediation

### P0 — claim hygiene (no engine change)

| # | Action |
|---|--------|
| 1 | Relabel `docs/BENCHMARK.md` 91.8% / ITC r=0.93 as **historical, unreproduced on this tip**; delete or rewrite the reproduce commands that point at missing scripts; match README’s unverified badge |
| 2 | Relabel `REPRODUCIBILITY.md` / `reproduce_astex85.sh` 94.1% as **oracle-ceiling (90% native seed + elitism)**, not docking power; stop verifying against a ±3% band around 94.1 |
| 3 | Update `METHODOLOGY.md` §0 for #365 merged (element-Hungarian, still not spyrmsd) and §4 ctest count |

### P1 — ranking identity (opt-in, parity-gated)

| # | Action |
|---|--------|
| 4 | One CF operand for REMARK / `.mcf` / BindingMode `pose.CF` / DatasetRunner election; decide `get_cf_evalue` vs apparent and apply on **both** CF-cluster and FO emission |
| 5 | Align FO emission with cluster.cpp pose-score consistency (or document FO as search-CF-only and stop calling it the same CF) |
| 6 | `ensemble_log_Z` optional/NaN instead of `!= 0.0`; never form \(\log Z\) from CF/\(k_B T\) |
| 7 | `scoring_env` (or ProtocolConfig fields) for every getenv scoring knob; make `--require-scoring-env` the default receipt check |
| 8 | Pin `FLEXAIDDS_RESTARTS=10` in published Astex receipts; default 5 must be labeled |

### P2 — hardening

| # | Action |
|---|--------|
| 9 | Exact element match for hydrogen; tertiary-amine virtual H; rotor-donor closed-form angle |
| 10 | Wire `score_reference.py` or retire `score_offline.py`; unify `<` vs `≤` 2.0 Å |
| 11 | Register or delete the three unbuilt `tests/test_*.cpp` files; `pytest.skip` oracle-ceiling |
| 12 | Replace Darwin `/var/folders/` scratch defaults with `tempfile` |
| 13 | Drop `\|\| true` on scancode; document a **tested** CMake minimum for CXX26+OpenMP (this session: 3.28.3 fail, 4.4.2 pass; do not claim 3.31 without a configure test) |
| 14 | Opt-in `FLEXAIDDS_RNG_STREAM_FIX` / `NAN_RANK_GUARD` only after a new baseline campaign |

---

## 8. Execution evidence (this session)

Host: Linux x86-64 cloud agent, GCC 14.2, CMake 4.4.2, Eigen 3.4.0, OpenMP 4.5.
Configure: `-DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_CUDA=OFF -DFLEXAIDS_USE_METAL=OFF`.
System CMake 3.28.3 **cannot** configure this tree (CXX26 + FindOpenMP `try_compile`); see M6.

| Suite | Result |
|-------|--------|
| `ctest --test-dir build --output-on-failure` | **94/94 passed** (4.47 s). METHODOLOGY.md §4 still says expect 11/11. |
| `test_classic_entropy_ranking --gtest_filter='SoftBeta*'` | **15/15 passed** (identity, strict duplicates, gated election, T resolve) |
| `test_dataset_runner --gtest_filter='RmsdCrossCheck*'` | **3/3 passed** (type vs element partitions) |
| `test_protocol_config` election/Softβ filter | **4/4 passed** including `ElectionShannonDefaultOffOptInAndLegacyZh` |
| `FlexAID --help` | exit 0, usage text present |
| Python science-gate subset (excluding Darwin-scratch test) | **79 passed, 1 skipped** |
| Same subset **including** `test_g4_2_niche_distance_drives_shipped_cpp_binary` | **1 failed** — `PermissionError` on Darwin `<temporary-directory>` prefix (finding M3) |
| `python3 scripts/check_repo_hygiene.py` | **OK** (does not flag the Darwin scratch default in `tests/`) |

No Astex dock was run. No success rate is reported.

Logs: `ctest_2026-08-16_science_audit.log`, `gtest_softbeta_identity.log`, `gtest_rmsd_crosscheck.log`, `pytest_science_gates.log`.

---

## 9. Evidence anchors (code, this tip)

| Claim | Anchor |
|-------|--------|
| Softβ shared header + local Z | `LIB/SoftBetaFreeEnergy.h`; `LIB/BindingMode.cpp` `compute_energy`; `LIB/cluster.cpp` ACF block |
| Softβ S1 default OFF | `LIB/ProtocolConfig.cpp` `from_env` |
| Dock T into DatasetRunner election | `LIB/DatasetRunner.cpp` `select_pose_freq_gated_pooled` |
| FO `.mcf` | `LIB/BindingMode.cpp` `write_mcf_sidecar` |
| FO discards re-score | `LIB/BindingMode.cpp` comment at ic2cf / `REMARK CF=% evalue` |
| CF vs apparent | `LIB/ic2cf.cpp` `get_cf_evalue` / `get_apparent_cf_evalue` |
| Dual-suffix enum | `LIB/DatasetRunner.cpp` `enumerate_emitted_cluster_heads` |
| RMSD refuse-on-mismatch | `python/flexaidds/dataset_runner/runner.py` `_pose_rmsd_vs_reference` |
| Element-Hungarian Python | `python/flexaidds/benchmark.py` `_symmetry_permutation` |
| Type vs element Hungarian | `tests/test_dataset_runner.cpp` `RmsdCrossCheck` |
| Vib correction zero | `LIB/BindingMode.cpp` `compute_vibrational_correction` |
| GPF log_Z sentinel | `LIB/DatasetRunner.cpp` TargetServer block |
| Receipt has no scoring_env | `LIB/RunReceipt.cpp` `build_run_receipt_json` |
| 94.1% native seed | `scripts/reproduce_astex85.sh`; `REPRODUCIBILITY.md` |
| 91.8% / r=0.93 | `docs/BENCHMARK.md` |
| README unverified | `README.md` Astex-85 badge |
| Hydrogen predicate | `LIB/hbond_potential.h` `is_hydrogen_atom` |
| RNG defect default | `LIB/RngSeed.h` `lazy_thread_rng` legacy branch |
| NaN sort default | `LIB/gaboom.h` `flexaids_nan_rank_guard_flag` |
| CI vs campaign | `LIB/top.cpp` defaults; `LIB/DatasetRunner.cpp` JSON emission |
| Unbuilt tests | `tests/test_ga_population.cpp`, `test_binding_mode_io.cpp`, `test_production_blockers.cpp` |

---

## 10. Relationship to prior audits

This audit does **not** replace METHODOLOGY.md. It does not restate parity /
Astex-85 A/B / ctest *procedures* — those stay in METHODOLOGY.md §1–§4.

It **does** record that several July P0 implementation holes are closed, that
the August RNG/NaN findings are gated rather than fixed, and that the
**unresolved load-bearing problem is no longer “Softβ is three formulas”** —
it is **“the claim surface and the two harnesses still do not describe one
experiment, and CF is still two aggregators.”**
