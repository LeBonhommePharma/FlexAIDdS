# FlexAIDdS — Past-week science and code audit (hard METHODOLOGY.md §1 bar)

**Date:** 2026-08-18  
**Window:** 2026-08-11 01:42Z → 2026-08-18 01:42Z  
**Week base:** `0971bd1a` (`actions/setup-node` bump, 2026-08-10)  
**Tip audited:** `bde7908cc600b396a957d22dcedc72f0c11251fd` (`main`, `#455`)  
**Companion (do not rewrite):** `docs/audit/2026-08-18_last24h_science_and_code_audit.md`  
**Prior full-tree audit on this week’s mid-point:** `docs/audit/2026-08-16_science_and_code_audit.md` (tip then `5c891e6e`)  
**Auditor:** Cursor Grok 4.6 (cloud agent), first-hand source + tests this session  
**Mode:** Diagnostic. This PR does not change engine, ranking, or thermodynamics.

Every factual statement below was re-read from files, `git log`/`git diff`, or test output in this session. No Astex dock was run. **No success rate is reported.**

The 24-hour write-up was too willing to treat “correctness bugfixes” as compatible. This document uses the letter of **METHODOLOGY.md §1**: an intended behaviour change must be opt-in behind an env flag that **defaults OFF**, and parity (1G9V, seed 12345, elected CF equal, 10 poses byte-identical) must hold with the flag off. Ungated default-path CF, ranking, election, or sampling changes are **epoch breaks**, even when the commit message says “fix” or “Science-Impact: none.”

---

## 1. One-line verdict

Across this week the tree gained real claim hygiene (`#437` withdrew 91.8%/94.1%) and a large **gated-OFF** acceleration / RNG / LUT stack. It also landed **two default-path epoch breaks** that a `benchmark_datasets` binary from `0971bd1a` does not share: **`#454` (ungated campaign `CF.hbond` topology)** and **`#416` (DatasetRunner parallel-restart fan-out now auto-throttled unless `FLEXAIDDS_MAX_CONCURRENT_RESTARTS=0`)**. Search still optimizes the Voronoi **CF scoring proxy**. Nothing here computes experimental **ΔG**. This repository still publishes **no receipted Astex-85 docking-power rate** on `bde7908c`. Do not mix pre-week and post-week CF tables.

---

## 2. The bar (not optional)

METHODOLOGY.md §1 (reproducibility / parity gate, run before **any** merge):

1. Build candidate and baseline separately; record both md5s.  
2. Dock 1G9V, `FLEXAID_SEED=12345`, `OMP_NUM_THREADS=1`, 2000 gen / pop 1000, crystal-pose seed off.  
3. PASS = elected CF equal **and** all 10 elected poses byte-identical.  
4. Any intended behaviour change must be opt-in behind an env flag that defaults OFF; parity must hold with the flag OFF.

This week’s merge messages often asserted “default-flag behaviour unchanged” without a 1G9V §1 receipt on disk. That assertion is not evidence. Below, “epoch” means: a scientist cannot treat `result.csv` / elected CF from before the merge as comparable to after, unless a named restore knob is proven off and §1 (or an equivalent completed-restart identity argument) is actually run.

---

## 3. Scope

`git log` on `origin/main` in the window: **60** non-merge commits, **~30** merge commits, **95** commits `0971bd1a..origin/main`. Engine/CI/python/tests/scripts/workflows: **101 files, +8954 / −541**.

### 3.1 Default-path / scoring-adjacent (must not be waved through)

| PR / commit | Subject | Hard §1 class |
|-------------|---------|----------------|
| `#415` `4183560c` | Contacts-epoch counter moved into stamp buffer | **Gated OFF.** Default memset path. Parser landmine (`getenv != nullptr`). Not a CF epoch if the var is **unset**. |
| `#416` `5fa020a5` | Restart concurrency auto-cap (`-1` default) | **Default-on scheduling epoch** for DatasetRunner `restarts > 1`. Restore: `FLEXAIDDS_MAX_CONCURRENT_RESTARTS=0`. |
| `#411` | Launcher fail-closed if `*.def` missing | No scoring. Prevents null 0.0% campaigns. |
| `#420` `5e294994` / `1824ba0e` / `99a43c7e` | RNG stream map, NaN CF rank, MI sign, StatMech NaN throw / Cv clamp | RNG/NaN **gated OFF**. StatMech API **ungated** (ledger / abort, not DatasetRunner election). |
| Accel `#429`/`#430` / chunks 1–6 / `#438` | LUT, rigid fastpath, niche hash, two-stage screen, keyed jitter, overlay | **Advertised gated OFF.** Overlay cannot enable LUT/epoch in `vcfunction.cpp` (TU-static snapshot). Fastpath is live `getenv`. |
| `#445` | Live `register_result` | Pose ranking no. **GPF occupancy now live.** Aug 16 **H2** is load-bearing. |
| `#453` | `native_score` sums all CF channels | Oracle `cf_native` / `probe_cf` only (`FLEXAIDDS_SCORE_NATIVE`). |
| `#454` | H-bond roles from topology, not charge | **Ungated campaign `CF.hbond` epoch.** See 24h H24-1. |
| `#455` | `FLEXAIDDS_FITNESS_MODEL` | Unset DatasetRunner = SMFREE (old DR hardcode). `FlexAID -c` missing JSON key still **PSHARE**. |
| `#448` | Ledger tENCoM λ | Default OFF; `inert_on_election`. Vib correction still 0.0. |
| `#446` | tENCoM mode↔structure pairing after sort | Writers follow `structure_index`. Docking election no. |
| `#443`/`#444`/`#449` | STRICT ∩ 85; Python `<= 2.0`; five-way RMSD | Claim / metric surface, not search CF. |
| `#452` | Exclude `1of6`/`2bys` from Tier-1; pin seed `20260816` | **CI subset only.** Campaign N still 85. Deepens C2. |
| `#450` | Blind Astex-85 receipt protocol | No dock. **Matrix pin fork 72d7 vs 9dc9.** |
| `#437` | Withdraw unverified 91.8%/94.1% | Closed Aug 16 **C1** on user-facing pages. |

### 3.2 Hygiene / docs / site (not CF)

`#414` compile_commands, `#417` setup-node, `#419` MT audit docs, `#421` A/B artifact guard, swarm pack `52d2566d`, `#423`/`#432` gitignore, `#431` site, `#433` site-stats off main, `#434` coverage ratchet, `#436`/`#451` audit + Tier-2 rebuild, site-stats auto commits, `5e67e6b7` output-prefix init (defined a previously uninitialized path; campaign `.grand.txt` prefixes already matched PDB ids).

### 3.3 Open, not merged

| PR | Subject | Ship? |
|----|---------|-------|
| `#439` | Always emit vib REMARK 0.0; proxy-only GPF/CF labels | **Hold** — rebase; ranking already 0-vib on main |
| `#441` | ParallelDock last-atom copy + `thread_local` ParEvalWS | **Hold** — real UB; tests do not exercise `create_workspace` |

---

## 4. Does a campaign binary from this tip match week-base `0971bd1a`?

| Change | Compiled CI (`FlexAID -c`, hbond JSON default **false**) | DatasetRunner campaign (hbond JSON **hardcoded true**; `restarts` default **5**; parallel restarts on when `restarts > 1`) | Election formula |
|--------|----------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------|------------------|
| `#454` H-bond topology | No (hbond off) | **Yes — `CF.hbond` in GA + rank** | Soft-β still \(\tilde G+\)vib(0); **inputs to \(\tilde G\) change** |
| `#416` restart throttle | n/a (single process) | **Yes — default auto-cap** (`-1`). Completed restart CF *can* match if every child finishes with the same seed/OMP. Pool composition / timeouts need not. | Unchanged formula |
| `#415` epoch rewind | No if env **unset** | No if env **unset**. `=0`/`=off` **enables** the epoch path (presence parser). | Unchanged if unset |
| Accel / LUT / RNG / NaN rank | No (default OFF) | No (default OFF). Overlay does **not** turn LUT/epoch on in `vcfunction.cpp`. | Unchanged |
| `#453` native_score | Diagnostic | `cf_native` only | Pose CF already had all channels |
| `#455` fitness_model | Env unused by `FlexAID -c` | Unset = SMFREE (old DR) | Unchanged unless `PSHARE` |
| `#445` GPF `register_result` | n/a | Post-election occupancy now live | Pose ranking no |
| StatMech NaN throw / MI sign / Cv clamp | If `StatMechEngine::compute*` is called | DatasetRunner election does **not** call it | Ledger / API |

**METHODOLOGY.md §1 against week-base is expected to FAIL on DatasetRunner defaults** because of `#454`. A 1-restart 1G9V with hbond off can still fail §1 if the test uses DatasetRunner’s hardcoded hbond JSON. A `FlexAID -c` JSON that leaves `hbond_enabled: false` should not see `#454`. Nobody ran that dock in this session.

Do not stitch Astex `result.csv` from before `#454` to after `#454`. Do not stitch 10-restart parallel campaigns from before `#416` to after `#416` without proving every restart completed.

---

## 5. Finding W1 (High) — `#454` is an ungated campaign scoring epoch

Fully specified in the 24h audit (H24-1). Restated here so the week document does not go soft:

- 3-arg `classify_hbond_donor` still returns **false** for every `N_sp3` / `O_sp3`. 3-arg acceptor for `N_sp3` is `partial_charge < 0.3`. Charge-0 SDF → ligand sp3 N was **acceptor-only**.
- `#454` passes topology when `bond[0] > 0`. Proteins get bonds via `residue_conect`. DatasetRunner writes `"hbond_enabled": true` and `"hbond_search_enabled": true` (`LIB/DatasetRunner.cpp`). CI / `config_defaults.h` keep hbond **false**.
- **`CF.com` unchanged; `CF.hbond` + virtual-H change** on the campaign path. No `FLEXAIDDS_HBOND_TOPOLOGY` flag.
- Ammonium rescue is **incomplete** on charge-0 SDF: `N.3` implicit H valence is 3 only if `charge < 0.3` (`LIB/top.cpp`). Tests call `encode_from_sybyl` with explicit `n_hydrogens` and **never** `top.cpp`. `HbondTopology.formal_charge` is never filled from `atoms[i].charge`.
- Do **not** claim “fixed CNS ammonium on Astex.” What moved: 1°/2° ligand amines, hydroxyls, and protein Ser/Thr/Tyr/Lys when bonds exist.

This is not a compatible bugfix. It is a new scoring epoch that happened to be scientifically motivated.

---

## 6. Finding W2 (High) — `#416` moved the default restart scheduler

`ProtocolConfig::max_concurrent_restarts` default is **`-1`** (auto). `restart_concurrency_cap()` then returns `max(1, cpu_budget / (omp_per_worker × workers))`. Legacy unlimited fan-out is **`0`**, not the unset default.

The merge claimed “Science-Impact: none. Scheduling only.” Under §1 that sentence is not available. The intended behaviour **did** change, and the restore knob is not the default.

What is true, and what is not:

- **True:** `OMP_NUM_THREADS` per child is unchanged. Parent-side `dock_config.json` / ligand / receptor are written before the fork, in launch order. Each child is an independent process with its own seed. A **completed** restart with the same argv/seed/OMP can be bit-identical to the unthrottled child.
- **True:** remaining timeout is `per_job_timeout_s` minus time since **that** child’s fork, not a shared parent budget stolen by earlier siblings.
- **Not proven:** campaign `result.csv` identity. Default `FLEXAIDDS_RESTARTS=5` (ProtocolConfig) and published Astex protocol **10**. Default parallel restarts = `(restarts > 1)`. Auto-cap **serializes** forks on a typical 11-core / 2-worker split (cap often 1–2). Later restarts start later. Outer campaign wall clocks, host load, and “how many restarts produced cluster heads for Fix B pooling” can move. Incomplete restarts change the election pool.
- **Not a CF formula change.** It is a **default-on sampling-completeness / schedule** change.
- `FLEXAIDDS_PARALLEL_RESTARTS=true` is still parsed with **`std::atoi` ≠ 0** (`LIB/ProtocolConfig.cpp`). `true`/`on`/`yes` → 0 → parallel restarts **off**. That parser predates this week (`0971bd1a` already had it). `#416` stacks a new default cap on top of a parser `EnvFlags.h` was written to kill.

Restore for mixable 10-restart tables: `FLEXAIDDS_MAX_CONCURRENT_RESTARTS=0` **and** every restart completed. §1 1G9V with one restart does not exercise this code.

---

## 7. Finding W3 (Medium) — `#415` is gated, but the gate is the wrong parser

`FLEXAIDDS_CONTACTS_EPOCH` remains **unset ⇒ OFF**. The rewind bug (counter in `FA_Global`, buffer resident across `tl_fa[t] = *FA`) only existed on the flag-ON path. Default `memset(FA->contacts, 0, MAX_ATOM_NUMBER * sizeof(int))` is still there. Allocations grew to `CONTACTS_BUFFER_SIZE` (`MAX_ATOM_NUMBER + 1`); the extra slot is unused on the memset path. `FA_Global` lost the `contacts_epoch` field — layout-only on default flags.

That is **not** a default-path CF epoch **if the variable is unset**.

It **is** a gate-integrity failure:

```185:186:LIB/vcfunction.cpp
static const bool contacts_epoch_mode =
    (std::getenv("FLEXAIDDS_CONTACTS_EPOCH") != nullptr);
```

`EnvFlags.h` exists specifically because `getenv != nullptr` treats `=0` / `=off` / `=false` as **ON**. `flags::active("FLEXAIDDS_CONTACTS_EPOCH")` uses `env_bool` (those values are OFF). The hot path does not.

TU-static snapshot runs before `top.cpp` calls `apply_to_environ()`. Overlay `FLEXAIDDS_FLAGS=epoch` **cannot** enable `contacts_epoch_mode` in `vcfunction.cpp`. Conservative for default-off; broken for the advertised overlay.

No in-tree campaign script exports this flag (only tests set `"1"`). Frozen Astex-84 (`docs/swarm/2026-08-13/`) did not rely on it. Still: anyone who exported `FLEXAIDDS_CONTACTS_EPOCH=0` “to be sure it is off” was already on the (then-corrupt, now-fixed) epoch path. Those runs are a separate epoch and were silently wrong before `#415`.

---

## 8. Finding W4 (Medium) — overlay vs hot-loop snapshots

`c5d5aa40` pushes resolved flags into `getenv` so legacy call sites honour `FLEXAIDDS_FLAGS=…`. That only works if the call site **re-reads** getenv after `apply_to_environ()`.

| Site | Reader | Overlay can enable? | Default |
|------|--------|---------------------|---------|
| `FLEXAIDDS_RIGID_FASTPATH` in `vcfunction` / `Vcontacts` | `rigid_fastpath_requested()` live `env_bool` | **Yes** | OFF |
| `FLEXAIDDS_GET_YVAL_LUT` in `vcfunction.cpp` | TU-static `get_yval_lut_enabled_cached()` | **No** (snapshotted before `top()`) | OFF |
| `FLEXAIDDS_CONTACTS_EPOCH` in `vcfunction.cpp` | TU-static `getenv != nullptr` | **No** | OFF |
| `FLEXAIDDS_NAN_RANK_GUARD` | process-static `env_bool` in `gaboom.h` | Only if set before first use of the magic static | OFF |
| `FLEXAIDDS_RNG_STREAM_FIX` | `env_bool`, re-read on seed-epoch | Yes if re-seeded after overlay | OFF |

Accel chunks 1–6, `#438` LUT hoist, keyed Voronoi jitter, niche hash, two-stage cube screen (`--two-stage` / `--coarse-prefilter`): **default OFF / CLI opt-in**. Tests this session: `test_get_yval_lut` 7/7, `test_rigid_fastpath` 6/6, `test_flexaidds_flags` 12/12, `test_contacts_epoch` 6/6. That is not a 1G9V §1 dock. It is enough to say the advertised defaults are OFF **when env is unset**.

RNG stream fix stays default OFF **on purpose**: the broken single-generator interleave (GA 0x9A800D / Vcontacts 0x0C0A11 / FOPTICS 0xF0701C5) is what frozen numbers were produced under (`LIB/RngSeed.h`). Enabling it is its own epoch. NaN rank guard stays OFF: finite-vs-finite unchanged; non-finite CF can still be elected rank-0.

---

## 9. Finding W5 (Low–Medium) — ungated StatMech API, not DatasetRunner election

`1824ba0e` / `99a43c7e`:

- `StatMechEngine::compute` / `compute_at_temperature` **throw** on non-finite sample energy and on all-zero multiplicity.
- `Cv` clamped `max(0, var) / (kB T²)`.
- `compute_joint_ensemble` MI sign flipped from \(S_J - S_R - S_L\) to \(S_R + S_L - S_J\), then clamped ≥ 0.

DatasetRunner does not call `StatMechEngine::compute` or `compute_joint_ensemble`. BindingMode’s `ligand_receptor_mutual_information()` already used \(I = S_L + S_R - S_J\) before this week. Default election remains min finite head **CF**; Softβ S1 remains OFF.

These are still **ungated** numeric/control-flow changes for any caller of the C++ StatMech API (Python bindings, CCBM tests, future thermo ranking). A previously completing NaN ensemble now aborts. That is fail-closed, not bit-identical. It did not follow “flag, default old behaviour.”

---

## 10. Finding W6 — `#445` made Aug 16 H2 live

`register_result` was on the same line as `// TODO` (commented out). It now runs. Pose ranking unchanged. Grand-canonical occupancy / `p_bind` / `.grand.txt` now see real `log_Z` **or** the `log_Z==0` → CF/`kT` fallback (H2). Single-ligand canonical paths with no TargetServer must remain pre-GPF; this session did not re-dock that identity.

---

## 11. Dual defaults and pin forks (week still split)

| Split | A | B | Why it matters |
|-------|---|---|----------------|
| H-bond | `config_defaults.h` / CI JSON **false** | DatasetRunner JSON **true** | `#454` is a campaign-only CF epoch |
| Fitness model | `config_parser.cpp` missing key → **PSHARE** | `config_defaults.h` + ProtocolConfig + DR → **SMFREE** | `FlexAID -c` ≠ DatasetRunner |
| Energy matrix | Claim / `generate_flexaid_inp.py` **`9dc93717…` (9dc9)** | `#450` receipt protocol **`72d7c739…`** | Mixing pins is mixing physics |
| Restarts | ProtocolConfig / many launchers **5** | METHODOLOGY.md published Astex **10** | Pool ceiling vs as-run |
| CI vs campaign | Tier-1 excludes `1of6`/`2bys` (`#452`) | Campaign N=85 | C2 still open; CI got *more* unlike campaign |
| RMSD cutoff | Python `docking_power` **`<= 2.0`** (`#444`) | Several summarize / offline scorers still **`< 2.0`** | Same pose, two verdicts |
| `native_score.h` | Comment: RMSDST **ignored** | `native_score.cpp` **loads** RMSDST | Header lies |

`#437` closed the user-facing 91.8%/94.1% C1. README remains unverified. That honesty does not license quoting the frozen Astex-84 17.9%/31.0%/48.8% as a current tip rate — those numbers are a **named referee CSV** on engine pin `dfc065ac…` / repo pin `aa15464e`, not `bde7908c`.

---

## 12. vs 2026-08-16 audit (C/H/M) at week tip `bde7908c`

| Finding | Status |
|---------|--------|
| **C1** unverified 91.8%/94.1% | **Closed in tree by `#437`** |
| **C2** CI ≠ campaign | **Open.** `#452` made CI’s pool smaller. Physics table in METHODOLOGY.md §0.1 unchanged |
| **H1** CF operand `evalue` vs `app_evalue` | **Open.** `#453` only fixed the native oracle sum |
| **H2** `log_Z==0` sentinel | **Open and live** (`#445`) |
| **H3** getenv-only scoring provenance | **Open.** Overlay does not reach LUT/epoch hot loops. New `FLEXAIDDS_FITNESS_MODEL` at least lands in ProtocolConfig JSON |
| **H4** vib correction production no-op | **Still 0.0** |
| **H5** five RMSDs; `<` vs `<=` | **Partial** (`#444`, `#449`) |
| **H6** restarts 5 vs 10 | **Open**, now plus `#416` default throttle |
| **M1** hydrogen predicate / tertiary virtual-H | **Open**, and **moved** by `#454` |
| **New W1 / H24-1** | Ungated `#454` `CF.hbond` |
| **New W2** | Ungated `#416` default restart cap |
| **New W3** | `CONTACTS_EPOCH` presence parser ≠ `env_bool` |
| **New** | 9dc9 vs 72d7 in `#450` |

---

## 13. Code review (non-science)

**Well done**

- `#415` makes the epoch/buffer split unrepresentable; `test_contacts_epoch` has a `WILL_FAIL` prefix-layout twin.
- Accel stack kept LUT/fastpath/RNG **default OFF** in `EnvFlags.h` / `get_yval.h` / `gaboom.h` / `RngSeed.h`.
- `#411` fail-closed launcher is the correct shape for “missing `AMINO.def` looks like a 0% science result.”
- `#437` is the only user-facing-rate fix that actually closed C1.
- `#455` rejects unknown / `pshare`; empty string does not invent a third arm.

**Gaps**

- `#454` shipped without a dock/CF test with `hbond_enabled: true` showing Lys NZ / ligand ammonium `cf.hbond` before vs after, and without a gate.
- `#416` shipped “Science-Impact: none” with ctest **not** run at commit time; no §1; no proof that 10/10 restarts still complete under the same outer wall clock.
- `contacts_epoch_mode` and `get_yval_use_lut` are TU-static; overlay is theatre for those two.
- `ProtocolConfig` still uses `atoi` for `FLEXAIDDS_PARALLEL_RESTARTS` after `EnvFlags.h` documented why that is poison.
- `native_score.h` still says RMSDST is ignored; `.cpp` loads it.
- Open `#441`: last-atom copy (`assign(atoms_, atoms_ + atm_cnt)` vs gaboom `natm + 1`) remains real UB.

---

## 14. Tests run this session

Build: CMake **4.4.2**, GCC **14.2**, `BUILD_TESTING=ON`, OpenMP ON, CUDA/Metal OFF, out-of-tree build directory. System CMake 3.28.3 still cannot enable CXX26+OpenMP (same as 2026-08-16).

| Suite | Result |
|-------|--------|
| `test_contacts_epoch` | **6/6 PASSED** |
| `test_flexaidds_flags` | **12/12 PASSED** |
| `test_get_yval_lut` | **7/7 PASSED** |
| `test_hbond_amine_roles` | **13/13 PASSED** |
| `test_protocol_config` | **15/15 PASSED** (includes fitness-model env + `max_concurrent_restarts` −1/0/3) |
| `test_statmech` filter `*NaN*:*NonFinite*:*Mutual*:*Joint*:*HeatCapacity*:*Cv*` | **11/11 PASSED** |
| `test_rigid_fastpath` | **6/6 PASSED** |
| pytest: published-rate check, blind receipt, campaign methodology gates, claim aggregator, run-receipt, Python `docking_power` | **56 passed**, 2 deselected (`test_g4_2_niche_distance_drives_shipped_cpp_binary` hits a Darwin `/var/folders/…` path and is not runnable here) |
| `python3 scripts/check_repo_hygiene.py` | **OK** |

No full `ctest`. No 1G9V §1 parity dock. No Astex-85. No live `result.csv`.

---

## 15. Claim-language bans on this tip

Do not say: current Astex-85 rate, genuine docking power, oracle ceiling as claim, S_top10, ITC r=0.93, “tENCoM elected,” “computed ΔG,” “`#454` fixed ammonium on Astex,” “`#416` cannot change science,” “`#415` is a no-op even with `CONTACTS_EPOCH=0`,” or that CI Tier-1 equals campaign. Seed-elitism / `_INI.pdb` is still forbidden as a result (METHODOLOGY.md §0). Default DatasetRunner election is still min finite head **CF**; Softβ S1 remains OFF.

Refuse S1 / S_top10 / STRICT % until all of the following exist **for this SHA’s binary**: `resolve_build.py --check` pin; **one** matrix MD5 (do not mix 9dc9 and 72d7); METHODOLOGY.md §1 1G9V vs `0971bd1a` (expect FAIL on DatasetRunner defaults because of `#454`); blind N=85, restarts=10, seed-echo=0, rank-0 in-place RMSD `<= 2.0 Å` with a named instrument; PoseBusters + tENCoM/Eigen on the same pose SHA-256; `result.csv` + `RUN_RECEIPT.json` with `fitness_model=SMFREE`, hbond flags, and `max_concurrent_restarts` recorded; `aggregate_claim_metrics.py --headline strict`.

Until those artifacts exist, the only honest numeric statement is: **this repository publishes no receipted Astex-85 docking-power rate on `bde7908c`.**

---

## 16. Recommended next (not done here)

1. **Gate `#454`** behind `FLEXAIDDS_HBOND_TOPOLOGY` default OFF, or accept a named scoring epoch and run a receipted A/B. Do not call it compatible.  
2. **Default `#416` restore to `0` (unlimited)** until a §1 + 10-restart completeness receipt exists, **or** keep auto-cap but document it as an epoch and require `FLEXAIDDS_MAX_CONCURRENT_RESTARTS` in every `RUN_RECEIPT`.  
3. Parse `FLEXAIDDS_CONTACTS_EPOCH` with `env_bool` (default false); re-read after overlay like fastpath. Same for `get_yval_use_lut` in `vcfunction.cpp`.  
4. Drive implicit H / `topo.formal_charge` from the same protonation evidence the `#454` classifier claims to use.  
5. Rebase `#439` (labels) and `#441` (ParallelDock copy + real workspace test).  
6. Align `#450` matrix pin with 9dc9, or name 72d7 as a packing fork in the receipt and never let `claim` accept either silently.  
7. Replace `atoi` on `FLEXAIDDS_PARALLEL_RESTARTS` with `env_bool`.  
8. Delete the stale “RMSDST ignored” paragraph in `LIB/native_score.h`.  
9. Do not merge pre-week and post-`#454` CF tables. Do not merge pre/post `#416` 10-restart pools without a completed-restart proof.
