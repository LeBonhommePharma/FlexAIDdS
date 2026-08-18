# FlexAIDdS — Last-24h science and code audit

**Date:** 2026-08-18  
**Window:** 2026-08-17 00:38Z → 2026-08-18 00:38Z (plus the tightly coupled Monday-night wave ~25 min older)  
**Tip audited:** `bde7908cc600b396a957d22dcedc72f0c11251fd` (`main`, `#455`)  
**Parent of window:** `f9e1d046` (`#444` merge)  
**Auditor:** Cursor Grok 4.6 (cloud agent), first-hand source + tests this session  
**Companions:** `docs/audit/2026-08-16_science_and_code_audit.md`, `docs/audit/2026-08-18_past_week_science_and_code_audit.md` (hard §1 bar, full week), `METHODOLOGY.md`, `AGENTS.md`  
**Mode:** Diagnostic. This PR does not change engine, ranking, or thermodynamics.

Every factual statement below was re-read from files, `git log`/`git diff`, or test output in this session. No Astex dock was run. **No success rate is reported.**

---

## 1. One-line verdict

The last 24 hours improved claim hygiene and CI honesty. They also landed **one ungated, default-on campaign scoring change** (`#454` H-bond topology). A `benchmark_datasets` binary from this tip is **not** METHODOLOGY.md §1 bit-identical to a binary from 24 h earlier. Search still optimizes the Voronoi **CF scoring proxy**. Nothing here computes experimental **ΔG**. This repository still publishes **no receipted Astex-85 docking-power rate** on `bde7908c`.

---

## 2. Scope

### 2.1 Strict 24 h (12 commits on `origin/main`)

| PR | Subject | Science class |
|----|---------|----------------|
| `#451` | Tier-2: no stale `build/` cache; push + Sunday cron `--dry-run` | CI hygiene |
| `#452` | Tier-1: exclude `1of6`/`2bys` by receptor size; pin seed `20260816` | CI subset only |
| `#453` | `native_score` sums all ten `get_cf_evalue` channels | Oracle diagnostic |
| `#454` | Ligand/protein sp3 N/O H-bond roles from topology, not charge | **Default-on CF.hbond** |
| `#436` | Land 2026-08-16 full science/code audit + host-path sanitize | Docs |
| `#455` | `FLEXAIDDS_FITNESS_MODEL` (`SMFREE` default, `PSHARE` opt-in) | Gated; unset = old DR |

Diff vs `f9e1d046`: 20 files, +1257 / −128.

### 2.2 Adjacent Monday-night wave (same session, ~25 min before the cutoff)

`#443` STRICT ∩ frozen 85, `#444` Python RMSD `<= 2.0`, `#445` live `register_result`, `#446` tENCoM mode–structure pairing, `#448` ledger-only tENCoM λ (OFF), `#449` five-way RMSD harness, `#450` receipt-gated blind Astex-85 protocol. `#437` (withdraw 91.8%/94.1%) is slightly older still and is what actually closed the Aug 16 **C1** user-facing-rate finding.

### 2.3 Open, not merged

| PR | Subject | Ship? |
|----|---------|-------|
| `#439` | Always emit vib REMARK 0.0; proxy-only GPF/CF labels | **Hold** — rebase; ranking already 0-vib on main |
| `#441` | ParallelDock last-atom copy + `thread_local` ParEvalWS | **Hold** — real UB, tests do not exercise `create_workspace` |

---

## 3. Does a campaign binary score differently?

| Change | Compiled CI (no JSON, hbond off) | DatasetRunner campaign | Election formula |
|--------|----------------------------------|------------------------|------------------|
| `#454` H-bond topology | No | **Yes — CF.hbond in GA + rank** | Soft-β still \(G̃+\)vib(0); **inputs to \(G̃\) change** |
| `#453` native_score channels | Diagnostic | `cf_native` / `probe_cf` only | Pose CF already had all channels |
| `#455` fitness_model | Env unused by `FlexAID -c` | Unset = SMFREE (old DR hardcode) | Unchanged unless `PSHARE` |
| `#451` / `#452` | CI pool / cron | No | No |
| `#436` | Docs | Docs | No |
| `#445` (adjacent) | n/a | Post-election GPF occupancy now live | Pose ranking no |

**METHODOLOGY.md §1 (1G9V, elected CF equal, 10 poses byte-identical) is not in this window and is expected to FAIL on DatasetRunner defaults across `#454`.** `#455` is the only change that followed “flag, default old behavior, fail-closed typo.”

---

## 4. Finding H24-1 (High) — `#454` is default-on campaign scoring

### 4.1 What the old classifier did

The 3-arg `classify_hbond_donor` still returns **false for every `N_sp3` and `O_sp3`**. The 3-arg acceptor for `N_sp3` is `partial_charge < 0.3f`. PDB-derived SDF charges are 0, so **every ligand sp3 N was acceptor-only**.

```145:154:LIB/atom_typing_256.h
        case N_sp3:
            // Covers 1°/2°/3° amine. Tertiary amine (no N-H) dominates drug
            // scaffolds and n_hydrogens cannot distinguish it → conservative
            // acceptor-only default (no donor role).
            return false;
        case O_sp3:
            // Covers ether (no O-H) and hydroxyl. Ether dominates drug
            // scaffolds and n_hydrogens cannot distinguish it → conservative
            // acceptor-only default (no donor role).
            return false;
```

### 4.2 What production now does

`top.cpp` always passes topology when `bond[0] > 0`. Proteins get bonds from `residue_conect` before typing. Ligands from SDF/MOL2 carry bond lists.

```1238:1240:LIB/top.cpp
		if (FA->is_protein) {
			residue_conect(FA, atoms, residue, deftyp);
		}
```

```1684:1692:LIB/top.cpp
						atom256::HbondTopology topo;
						topo.n_heavy_neighbors = heavy_neighbor_count(i);
						topo.known             = (atoms[i].bond[0] > 0);
						atoms[i].type256 = atom256::encode_from_sybyl(
							atoms[i].type,   // SYBYL type 1–40
							atoms[i].charge, // partial charge (MOL2 or AMBER ff14SB)
							n_hydrogens,     // explicit + conservative implicit H
							topo             // heavy-atom substitution evidence
						);
```

With `topo.known`, 1°/2° amines become **amphoteric**, hydroxyls (`O_sp3`, ≤1 heavy) become **donors**, 3° amines stay acceptor-only, quaternary/ammonium (`coordination ≥ 4` or `N.4`) become donor-only.

This is live on the **claim path**. DatasetRunner hardcodes H-bond search on:

```6053:6057:LIB/DatasetRunner.cpp
                   << "    \"hbond_enabled\": true,\n"
                   << "    \"hbond_search_enabled\": true,\n"
                   << "    \"hbond_rank_enabled\": "
                   << (hbond_rank ? "true" : "false") << ",\n"
                   << "    \"metal_coord_enabled\": true,\n"
```

CI / no-JSON engine defaults still have `hbond_enabled: false` (`LIB/config_defaults.h`). A green Tier-1 tick is **not** this change’s campaign effect (METHODOLOGY.md §0.1).

**CF.com is unchanged** (still indexes SYBYL `atom.type`). **CF.hbond and virtual-H geometry change.** Soft-β election consumes those CF values.

### 4.3 The advertised ammonium rescue is incomplete on charge=0 SDF

Production implicit H for SYBYL `N.3` is still charge-gated:

```1602:1605:LIB/top.cpp
					case 8: { // N.3
						const int valence = (a.charge >= 0.3f) ? 4 : 3;
						const int h = valence - heavy_bonds;
						return h > 0 ? h : 0;
```

On the PDB-derived SDF path (`charge == 0`), a heavy-atom-only protonated 3° amine still looks like R3N with 0 H → **still acceptor-only**. Tests that pass `n_hydrogens` into `encode_from_sybyl` never see this (`tests/test_hbond_amine_roles.cpp` does not call `top.cpp`).

`HbondTopology.formal_charge` is also never filled from `atoms[i].charge` in `top.cpp`.

Do **not** claim “fixed CNS ammonium pharmacophore on Astex.” What actually moved on charge=0 SDF: **1°/2° ligand amines and hydroxyls**, plus **protein Ser/Thr/Tyr OH and Lys NZ** whenever `residue_conect` populated `bond[]`.

### 4.4 METHODOLOGY.md §1

Intended scoring changes must be env-gated **OFF** with parity holding when off. `#454` has **no flag**. Suggested name: `FLEXAIDDS_HBOND_TOPOLOGY` (unset = 3-arg classifier). Then §1 on 1G9V with the flag off, then a documented A/B with it on, before any campaign number is compared to pre-`#454`.

---

## 5. Other 24 h changes (ranked)

### H24-2 (Medium) — `#453` `cf_native` baseline moved

GA `get_cf_evalue` already summed ten channels + optional `tencom_weight * h_rep` (`LIB/ic2cf.cpp`). `score_native_pose` previously omitted `elec` / `gist_desolv` / `metal_coord` / `entropy`. It now copies those four. Campaign JSON has `tencom_weight: 0.0` and GIST/elec structurally ~0, so the practical oracle delta is **`metal_coord` + VCT `entropy`**. Pose ranking is unchanged. Post-`#453` native-vs-elected CF tables are a **new baseline**; do not mix with pre-`#453` oracle gaps.

`h_rep` is still not copied. The unit “test” is a source-string grep (`tests/test_dump_pop_refstructure_wiring.py`), not `native_total == get_cf_evalue`.

### H24-3 (Medium) — `native_score.h` contradicts `native_score.cpp`

The header still says `FLEXAIDDS_RMSDST` is **intentionally ignored** because raw SDF coordinates are in the wrong frame. The `.cpp` **loads that SDF and injects crystal XYZ** before `vcfunction`. Easy to “fix” the diagnostic by deleting RMSDST and silently score the blinded pose.

### H24-4 (Medium) — `#455` is fail-closed in ProtocolConfig, unused on the engine CLI

Unset/empty → `SMFREE`. Only exact `SMFREE` or `PSHARE`. Unknown throws. DatasetRunner writes `ga.fitness_model` from `protocol_cfg_`. Unset DatasetRunner path is bit-identical to historical SMFREE.

The engine does **not** `getenv("FLEXAIDDS_FITNESS_MODEL")`. `LIB/config_parser.cpp` still defaults a **missing** JSON key to **`PSHARE`**, which disagrees with `config_defaults.h` (`SMFREE`). Exporting the env var and running `FlexAID -c` without DatasetRunner is a silent no-op. SMFREE remains a **CF Boltzmann fitness**, not ΔG, and not Softβ S1 (still OFF).

### H24-5 (Medium) — `#451` weekly cron is fail-closed; `workflow_dispatch` is not

Push and Sunday cron force `--dry-run` until Zenodo IDs exist. Manual dispatch still defaults `dry_run: "false"`. A click against empty `benchmark_data/` is the remaining live path. `FILTER_ARGS` is also always `--all` even when `DATASET_FILTER` is set (pre-existing, edited in this file).

### L24-1 (Low) — `#452` does not change campaign N=85

`1of6` / `2bys` leave the **Tier-1 eligible pool** only (apo ATOM+HETATM > 15098). `DatasetRunner::astex_diverse_codes()`, YAML `targets:`, and the frozen STRICT manifest still have 85 including those two. Comment still calls `2c3i` “newly pinned”; seed `20260816` actually draws `1gpk 1mq6 1xm6 2cet`. `1n1m` (12843 atoms) stays eligible. Do not quote a 4-target CI draw as Astex-85.

### L24-2 (Low) — `#436` audit text is already stale on C1

The merged 2026-08-16 audit still says `docs/BENCHMARK.md` / `REPRODUCIBILITY.md` publish 91.8% / 94.1%. **`#437` withdrew those rates** a few hours earlier. Treat that C1 paragraph as a snapshot of `5c891e6e`, not of `bde7908c`.

---

## 6. Adjacent Monday wave (needed to read the 24 h)

| PR | Moves ranking? | Note |
|----|----------------|------|
| `#443` STRICT ∩ frozen 85 | No | Fail-closed aggregator. Off-manifest extras cannot inflate n. Missing targets stay failures. |
| `#444` Python `<= 2.0` | No | Aligns Python `docking_power` with METHODOLOGY.md §0 and C++ `success_rmsd`. Exclusive `< 2.0` still lives in `validate_benchmark_results.py`, `score_offline.py`, `score_reference.py`, several summarize scripts, and leftover docs. |
| `#445` live `register_result` | Pose no | Call was on the same line as `// TODO` (commented out). Now every DatasetRunner receptor registers GPF occupancy. **Aug 16 H2 (`log_Z==0` → CF/`kT` fallback) is now load-bearing, not dormant.** |
| `#446` tENCoM pairing | Docking no | Writers follow `structure_index` after `sort_by_free_energy()`. |
| `#448` ledger λ | No | `FLEXAIDDS_LEDGER_TENCOM_LAMBDA` default OFF; `inert_on_election=1`. `compute_vibrational_correction()` still returns 0.0. |
| `#449` five-way RMSD | No | Tests/docs. METHODOLOGY.md §0 still says “#365 … still open” even though `python/flexaidds/benchmark.py::_symmetry_permutation` exists. |
| `#450` blind receipt protocol | No dock | `claim` refuses a % without `RUN_RECEIPT.json`. **Pin fork:** this script and `docs/implementation/BLIND_ASTEX85_RECEIPT_PROTOCOL.md` cite matrix MD5 `72d7c739…`; claim aggregator / `arm_pins.json` / `generate_flexaid_inp.py` cite **`9dc93717…` (9dc9)**. |

---

## 7. vs 2026-08-16 audit (C/H/M)

| Finding | On `bde7908c` |
|---------|----------------|
| **C1** unverified 91.8% / 94.1% on user-facing pages | **Closed in tree by `#437`**, not by `#436` |
| **C2** CI ≠ campaign | **Open.** `#452` makes CI’s target pool smaller; §0.1 physics table unchanged |
| **H1** CF operand split (`evalue` vs `app_evalue`) | **Open.** `#453` only fixed the native oracle sum |
| **H2** `log_Z==0` sentinel | **Open and now live** (`#445`) |
| **H3** getenv-only scoring provenance | **Open.** New `FLEXAIDDS_FITNESS_MODEL` at least lands in ProtocolConfig JSON |
| **H4** vib correction production no-op | **Still 0.0.** `#448` tags inert; `#439` would make the zero visible |
| **H5** five RMSDs; `<` vs `<=` | **Partially closed** (`#444`, `#449`). Exclusive scripts remain |
| **H6** restarts 5 vs published 10 | **Open** |
| **M1** hydrogen predicate / tertiary virtual-H | **Open** |
| **M2** METHODOLOGY drift (#365 “open”, §4 ctest 11/11) | **Open** |
| **New H24-1** | Ungated `#454` CF.hbond role change |
| **New** | 9dc9 vs 72d7 matrix pin in `#450` |

---

## 8. Code review (non-science)

**Well done**

- `topo.known == false` reproduces the 3-arg verdict (`LIB/atom_typing_256.h`); PDB atoms with empty `bond[]` stay byte-identical on H-bond flags.
- Second `#454` commit (`fa824812`) uses coordination number = heavy + H so R–NH3+ is not amphoteric *when nH is actually 3*.
- `#455` rejects unknown / `pshare`; empty string does not invent a third arm; `ClearProtocolEnv` unsets the var.
- `#451` kills the proven-stale `build/` cache (content-addressed ccache only).
- `#452` pins the draw in YAML + workflow + a test that both `astex_diverse.yaml` copies agree.

**Gaps**

- No dock/CF test with `hbond_enabled: true` showing Lys NZ / ligand ammonium `cf.hbond` before vs after.
- H-bond tests never call `conservative_implicit_h_count`.
- `#455` does not assert DatasetRunner JSON `ga.fitness_model` equals `protocol_cfg_.fitness_model`, or that the engine `getenv`s the var.
- `#441` (open): last-atom copy is real (`assign(atoms_, atoms_ + atm_cnt)` vs gaboom `natm + 1`); `thread_local` under OpenMP is unproven.

---

## 9. Tests run this session

Build: CMake **4.4.2**, GCC **14.2**, `BUILD_TESTING=ON`, OpenMP ON, CUDA/Metal OFF, out-of-tree `build_audit24`. System CMake 3.28.3 still cannot enable CXX26+OpenMP (same as 2026-08-16).

| Suite | Result |
|-------|--------|
| `python/tests/test_tier1_random_subset.py` + native_score wiring + blind receipt + aggregate claim metrics | **55 passed** |
| `./test_hbond_amine_roles` | **13/13 PASSED** |
| `./test_protocol_config` | **15/15 PASSED** (includes `FitnessModelFromEnv`) |

No full `ctest`. No 1G9V §1 parity dock. No Astex-85.

---

## 10. Claim-language bans on this tip

Do not say: current Astex-85 rate, genuine docking power, oracle ceiling as claim, S_top10, ITC r=0.93, “tENCoM elected,” “computed ΔG,” “#454 fixed ammonium on Astex,” or that CI Tier-1 equals campaign. Seed-elitism / `_INI.pdb` is still forbidden as a result (METHODOLOGY.md §0). Default DatasetRunner election is still min finite head **CF**; Softβ S1 remains OFF.

Refuse S1 / S_top10 / STRICT % until all of the following exist **for this SHA’s binary**: `resolve_build.py --check` pin; one matrix MD5 (do not mix 9dc9 and 72d7); METHODOLOGY.md §1 1G9V vs pre-`#454` (expect FAIL on defaults); blind N=85, restarts=10, seed-echo=0, rank-0 in-place RMSD `<= 2.0 Å` with a named instrument; PoseBusters + tENCoM/Eigen on the same pose SHA-256; `result.csv` + `RUN_RECEIPT.json` with `fitness_model=SMFREE` and hbond flags; `aggregate_claim_metrics.py --headline strict`.

Until those artifacts exist, the only honest numeric statement is: **this repository publishes no receipted Astex-85 docking-power rate on `bde7908c`.**

---

## 11. Recommended next (not done here)

1. Gate `#454` behind `FLEXAIDDS_HBOND_TOPOLOGY` default OFF, or accept a new scoring epoch and run a receipted A/B.  
2. Drive implicit H / `topo.formal_charge` from the same protonation evidence the classifier claims to use.  
3. Rebase `#439` (labels) and `#441` (ParallelDock copy + real workspace test).  
4. Align `#450` matrix pin with 9dc9, or document 72d7 as a named packing fork in the receipt.  
5. Default Tier-2 `workflow_dispatch.dry_run` to `true`.  
6. Delete the stale “RMSDST ignored” paragraph in `LIB/native_score.h`.
