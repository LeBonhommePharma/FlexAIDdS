# WORKORDER — Clustering P2: CF-independent within-cluster representative (medoid / density-peak center)

**Priority:** P2 (cheapest principled selection fix)
**Domain:** pose clustering + within-cluster rank-0 representative election
**Author:** OPS/Benchmarker (Bonhomme) · paste-in spec for the committing engineer (Claude Code / Grok / Codex)
**Repo:** `/Users/lp.more/Projects/FlexAIDdS`
**HEAD when written:** `3e674479c66b460605f2267c5dc89a749fc0391e` on branch `fix/fo-cluster-emit-low-budget` (⚠ NOT on `origin/main`)
**Constraints honored:** open-source, dependency-minimal, energy-matrix-based (MC_st0r5.2_6.dat / Vcontacts), env-gated DEFAULT-OFF (defaults bit-identical to current), reproducible (`FLEXAID_SEED=12345`), rank-0 symmetry-corrected Hungarian RMSD @ 2.0 Å, benchmark norm pop=1000 / 2000 gen / 10 restarts.

> **You are the committer.** OPS does not commit or build (host is RAM-constrained ~5 GB free; a single 2000-gen dock OOMs during emission). OPS owns all measurement/canary runs. Implement exactly what is below, prove flag-OFF bit-identity, and stop — do NOT run the benchmark.

---

## 0. TL;DR — read this first, the fix already half-exists and is broken

The HEAD commit `3e674479c` (title: *"Fix: DP cluster head election + pose_score_consistent REMARK; add CF medoid refinement"*) **already attempted this P2**, but ships three defects that violate the hard constraints and/or invert the intended change:

1. **`cluster.cpp` medoid is DEFAULT-ON** (`FLEXAIDDS_MEDOID_REFINE` unset ⇒ `true`, `cluster.cpp:213-216`). This makes **defaults NOT reproduce prior behavior** — a hard-constraint violation. Must flip to default-OFF.
2. **The medoid is Boltzmann/CF-weighted** (`w_n = exp(−β·(E_n−E_min))`, `cluster.cpp:235-254`), re-injecting the CF signal into what must be a **CF-independent geometric** selector. The whole premise (verified: within-target Spearman(CF,RMSD) ≈ −0.002) is that CF is orthogonal to pose correctness, so weighting the medoid by CF partially defeats the purpose. Must offer a **pure unweighted geometric medoid**.
3. **`DensityPeak_Cluster.cpp` went the WRONG way**: HEAD flipped `#define OUTPUT_CLUSTER_CENTER true → false` (`DensityPeak_Cluster.cpp:27`), i.e. it **removed** density-peak-center election and reverted DP to lowest-CF — the *opposite* of this P2. Its stated motive ("score_pose_consistent=0") optimizes a **CF-consistency diagnostic that is ⊥ RMSD** — the exact anti-goal. And it is a hardcoded `#define`, not an env gate.

This work order **supersedes** that commit's clustering changes: it (a) makes the CF/leader medoid pure + default-OFF, (b) re-exposes the DP center as an env-gated option (default-OFF = current lowest-CF), (c) unifies all backends under one gate, and (d) specifies the **full-population dump prerequisite** — which, I discovered, *partially already exists* as the `.rrd` files, so this is cheaper than the domain report assumed.

**Net verdict:** worth_it = **yes**, but eyes-open. On the *current emitted representative pool* the medoid swap is provably **neutral** (measured Δ = −0.005 Å; see §8): emitted clusters are 2 Å-tight so all members have near-equal native-RMSD. Its value is (i) correcting two live constraint violations already committed, (ii) delivering the validation instrument (dump), and (iii) real RMSD payoff **only once near-native poses survive to the pool** (that is the SEARCH/representation work order's job). Do not expect this alone to move the 0/19 success rate.

---

## 1. Current state audit — every anchor verified at HEAD `3e674479c`

| backend | dispatch | within-cluster representative NOW | file:line |
|---|---|---|---|
| **CF / leader** (default, `clustering_algorithm="CF"`) | `top.cpp:2956-2957` → `cluster()` | greedy head = lowest-CF seed, THEN Boltzmann-CF-weighted medoid **(default-ON)** overwrites head | `cluster.cpp:190-304` (medoid), emit at `:391-563` |
| **FO** (FastOPTICS, `"FO"`) | `top.cpp:2944-2948` → `FastOPTICS_cluster()` → `BindingMode::output_BindingMode()` | **OPTICS density-center** (min reachDist), fallback lowest-CF; + gated `FLEXAIDDS_PB_AWARE_PROMOTION` hook | `BindingMode.cpp:561-571`, gate `:574-618`, elector `:907-914` |
| **DP** (DensityPeak, `"DP"`) | `top.cpp:2949-2953` → `DensityPeak_cluster()` | **lowest-CF** (`Representative`/`BestCF`); density `Center` computed then discarded | `DensityPeak_Cluster.cpp:27` (`#define`), `:463-473`, `:567-568` |
| **InStream** (online GA-time) | `InStreamClustering.cpp` | tracks `best_score` = lowest CF (class is misnamed "ClusterMedoid" — it is NOT a geometric medoid) | `InStreamClustering.cpp:130-148` |

Key structural facts the implementation must respect:
- **CF/leader** stores clusters as `Clus_GAPOP[k]=head_idx`, `Clus_TOP[cl]=head_idx`, `Clus_FRE[cl]`, `Clus_TCF[cl]`, `Clus_ACF[cl]`. Between-cluster rank-0 is by **ACF** (soft-β free energy) when `T>0 && !force_cf_rank_emission` (`cluster.cpp:306-321`), NOT by the representative's CF. The medoid change must therefore **not touch ACF** (it is computed over all members and is representative-independent) — only which member is emitted. The HEAD commit gets this right (`cluster.cpp:285` leaves ACF alone); preserve that.
- Coordinates for all `num_chrom` chromosomes are pre-cached in `coord_cache` (`cluster.cpp:100-120`), stride `nAtoms_clus*3`, distances via `flexaids::sum_sq_distances_f`. A pure medoid reuses this cache — **zero extra geometry rebuilds**.
- `.mcf` sidecar (member CFs, `cluster.cpp:546-562`) enumerates members by `Clus_GAPOP[k]==Clus_TOP[j]`; the HEAD medoid remaps `Clus_GAPOP` (`:277-283`) so the sidecar stays consistent. Preserve.
- **DP** already has both `pCluster->Center` (isCenter chromosome) and `pCluster->Representative` (lowest-CF) populated (`DensityPeak_Cluster.cpp:459-473`); the `.cad` and inter-cluster RMSD already reference the center. So exposing the center as the emitted pose is a **one-branch flip** already scaffolded at `:567-568`.

---

## 2. Integration points (exact, file:line at HEAD `3e674479c`)

**IP-1 — CF/leader medoid: make pure + default-OFF.** `LIB/cluster.cpp:212-304`.
Replace the default-ON Boltzmann-CF-weighted block with a gate that (a) defaults OFF, (b) supports a pure unweighted geometric medoid, (c) keeps the Boltzmann variant only as an explicit opt-in mode.

**IP-2 — DensityPeak center: env-gate, default-OFF = current lowest-CF.** `LIB/DensityPeak_Cluster.cpp:27` (the `#define OUTPUT_CLUSTER_CENTER false`) and the election branch `:567-568`. Replace the compile-time `#define` with a runtime read of the same env gate; default (unset) keeps `false` (lowest-CF, current HEAD behavior).

**IP-3 — FastOPTICS/BindingMode: add medoid as an election mode (optional, low-value here).** `LIB/BindingMode.cpp:561-571` + elector `:907-914`. This backend already elects a consensus center (min reachDist), so it already satisfies the P2 intent. Add a true-medoid branch under the same gate only for cross-backend consistency; it is NOT required to close the finding.

**IP-4 — InStream (secondary).** `LIB/InStreamClustering.cpp:130-148`. The online path tracks lowest-CF under a "medoid" name. Out of scope for the benchmark (InStream is not the benchmark clustering path), but flag in a code comment that it is CF-tracking, not geometric, to avoid future confusion.

**IP-5 — Full-population dump prerequisite (the binding item).** `LIB/write_rrd.cpp` (`write_rrd` and `write_DensityPeak_rrd`) + call sites `cluster.cpp:566`, `DensityPeak_Cluster.cpp:676`, both gated on `FA->refstructure==1` (set by DatasetRunner via `FLEXAIDDS_RMSDST`, `DatasetRunner.cpp:6234` / `native_score.cpp:115`). See §5 — the `.rrd` already IS a full-population pre-emission dump; augment it, don't rebuild.

---

## 3. Env gate design

Single primary gate, default-OFF, with a mode selector so all backends move together and defaults stay bit-identical:

```
FLEXAIDDS_CLUSTER_REP   (unset|"lowcf"  = current behavior, DEFAULT)
                        ("medoid"        = pure unweighted geometric medoid)
                        ("bmedoid"       = Boltzmann-CF-weighted medoid, the HEAD-commit variant)
                        ("center"        = density-peak / OPTICS center where the backend has one)
```

Rules:
- **Unset or `lowcf` ⇒ every backend reproduces current HEAD-minus-defects behavior bit-for-bit.** Specifically: CF/leader emits the lowest-CF head (NO medoid), DP emits `Representative`/BestCF, FO emits its OPTICS center (its existing default).
- `FLEXAIDDS_MEDOID_REFINE` (the HEAD env var) is **retired/aliased**: if present and non-zero, treat as `FLEXAIDDS_CLUSTER_REP=bmedoid` (backward-compat for any script already using it), but log a deprecation line. The default-ON behavior is removed.
- Parse once into an enum on `FA_Global` (e.g. `FA->cluster_rep_mode`) in `config_parser.cpp` alongside the existing `FA->beta` read (`config_parser.cpp:158`), so all three backends read one resolved value and the gate is testable. If touching `config_parser` is off-limits per house rule, read `getenv` locally in each backend with an identical helper — but centralize the string→enum map in one header to avoid drift.

---

## 4. Algorithm — pure geometric medoid (mode `medoid`)

For each cluster with ≥3 members (`Clus_FRE[cl] ≥ 3`; singletons and pairs keep the head — a 2-member medoid is degenerate):

```
members = { k : Clus_GAPOP[k] == old_head }              # includes head
medoid  = argmin_{m ∈ members}  Σ_{n ∈ members, n≠m}  sqRMSD(x_m, x_n)
```
- `sqRMSD` uses the cached `coord_cache` and `flexaids::sum_sq_distances_f` already in scope (`cluster.cpp:147`); work in squared-distance units (skip the sqrt — monotone, so argmin is identical). This is the **unweighted** version of the cost the HEAD commit already computes at `cluster.cpp:256-273`; delete the weight vector.
- **No CF, no Boltzmann weight, no temperature dependence** — that is the point (CF ⊥ RMSD). The medoid is a function of geometry alone.
- Complexity **O(k²·A)** per cluster (k = members, A = ligand heavy atoms), k small (median emitted cluster size here ≈ 10–25, cap by `num_of_results`). Negligible vs the 2M-eval GA.
- On tie (equal summed sqRMSD within 1e-9), keep the lower array index (deterministic, seed-independent).
- Remap identically to HEAD (`cluster.cpp:277-287`): repoint `Clus_GAPOP[k]` from `old_head`→`medoid` for all members, set `Clus_TOP[cl]=medoid`, `Clus_TCF[cl]=chrom[medoid].app_evalue`; **leave `Clus_ACF[cl]` untouched** (ranking invariant). Emit the `[MEDOID_REFINE]` stdout line so OPS can audit which clusters moved.

Mode `bmedoid` = the HEAD algorithm verbatim (`w_n = exp(−β·(E_n−E_min))/Z`) — kept as an opt-in for ablation only. Mode `center` (DP/FO) = emit `pCluster->Center` / min-reachDist pose (already computed).

REMARK additions on every emitted pose (all backends): `REMARK cluster_rep_mode=<lowcf|medoid|bmedoid|center>` and, when the pick moved, `REMARK cluster_rep_shifted=1 head_cf=<..> rep_cf=<..> wRMSD2=<..>` so the elected-pose provenance is auditable from the PDB alone.

---

## 5. THE PREREQUISITE — full-population pre-clustering dump (mostly already exists)

The domain report states no full-population dump exists and that no selection change is validatable without one. **Correction from on-disk inspection:** the `.rrd` files already emitted by `write_rrd`/`write_DensityPeak_rrd` (gated `refstructure==1`, active in benchmark mode) ARE a per-restart full-population dump. Verified format (`write_rrd.cpp`):

```
idx  cluster_head(GAPOP)  rmsd_to_head  rmsd_raw  rmsd_sym  evalue(CF)  [ gene_0 … gene_npar ]
```
one row per chromosome (1000 rows/restart), with **symmetry-corrected RMSD-to-native (`rmsd_sym`), CF, and cluster label** — exactly the join the report asked for. Confirmed present for all 19 ablation targets under `…/three_engine/A/3dsig_full85_scratch_3b2fa57cc/<PDB>/<PDB>_r0.rrd`.

**What `.rrd` is missing (augment `write_rrd` to add these — small, additive, and independently useful):**
1. **CF components** `com` / `wal` per chromosome (currently only total `evalue`). Needed to reconcile with `cf_pose_components_19targets.csv` and to test component-based selectors. Source: the per-chromosome `cfstr` already populated during scoring — write `cf.com`, `cf.wal` alongside `evalue`.
2. **A stable `pose_id`** that joins each population row to its emitted representative PDB when it is one (else `-1`). Lets OPS confirm "the near-native population pose was/was not the elected representative".
3. **The elected flag** — 1 if this chromosome is the emitted rank-j representative, else 0.

**Format proposal (backward-compatible — append columns after the existing 6, before the `[genes]` block, OR add a sibling `.pop.tsv`):** prefer a sibling `<prefix>_rN.pop.tsv` with a header row so downstream parsers don't break on the fixed-width `.rrd`:
```
idx	cluster	rmsd_to_head	rmsd_raw	rmsd_sym	cf_total	cf_com	cf_wal	pose_id	is_elected
```
**Where to emit:** the natural site is inside `cluster()` / `DensityPeak_cluster()` right after clustering and before/around the existing `write_rrd` call (`cluster.cpp:566`), since `Clus_GAPOP`, `Clus_TOP`, member CFs and (with `refstructure==1`) per-chrom RMSD are all in scope there. Gate the extra columns on `refstructure==1` (already the `.rrd` gate) plus an opt-in `FLEXAIDDS_DUMP_POP=1` so default benchmark output is unchanged. DatasetRunner already sets `FLEXAIDDS_RMSDST`; it should additionally set `FLEXAIDDS_DUMP_POP=1` when an audit dump is requested (one line near `DatasetRunner.cpp:6234`).

This prerequisite is the highest-value part of the work order: it is what makes the medoid (and every future selection experiment) measurable, and it is ~cheap because the `.rrd` scaffolding already exists.

---

## 6. New dependencies / licensing

**None.** No Eigen, no external potential, no new data file. Pure C++ over the existing `coord_cache` + `flexaids::sum_sq_distances_f` (already linked). The energy matrix (`MC_st0r5.2_6.dat`) and Vcontacts base are untouched — the medoid is geometry-only and does not read the scoring matrix at all. Clean-room: nothing vendored; all code is original to this repo. No license impact.

---

## 7. Acceptance gates (what OPS runs to verify — you must pass 1–3 before handing back)

1. **Flag-OFF bit-identity (MANDATORY).** With `FLEXAIDDS_CLUSTER_REP` unset (and `FLEXAIDDS_MEDOID_REFINE` unset), a dock of the canary set **1G9V, 1SJ0, 1OPK, 1M2Z, 2HB1** (`FLEXAID_SEED=12345`, pop=1000/2000 gen/10 restarts, matrix pinned via `FLEXAIDDS_DATA_DIR`) must produce **byte-identical emitted `_j.pdb` poses and `.cad`/`.mcf`** vs a build of HEAD with the HEAD medoid **forced OFF** (`FLEXAIDDS_MEDOID_REFINE=0`). This proves the default path reproduces prior behavior. (NB: it will NOT match HEAD's *default* output, because HEAD's default is ON — that is the defect being fixed; compare against HEAD+`=0`.)
2. **`ctest` green** (existing suite) with the gate unset.
3. **Mode plumbing smoke:** on 1G9V, `FLEXAIDDS_CLUSTER_REP=medoid` emits `REMARK cluster_rep_mode=medoid` and a `[MEDOID_REFINE]` line for ≥1 cluster; `=lowcf` emits neither and equals gate #1 output; `=center` under `clustering_algorithm=DP` emits the density center (differs from lowcf on ≥1 cluster). No crash, no NaN RMSD, `score_pose_consistent` REMARK still written.
4. **Elected-pose RMSD delta (OPS-owned, informational).** OPS re-elects on the canary and on the 19-target ablation: report median rank-0 symmetry-corrected RMSD for {lowcf, medoid, bmedoid, center}. **Expected on the current emitted pools: Δ ≈ 0** (see §8) — this gate documents neutrality, it is not a pass/fail on RMSD.
5. **Prerequisite dump audit.** With `FLEXAIDDS_DUMP_POP=1`, the `.pop.tsv` for 1OF6 contains ≥1 sub-3 Å (`rmsd_sym<3`) population row (verified to exist: 6 such rows), and OPS can confirm from `is_elected`/`pose_id` whether any sub-3 Å pose was elected. This is the instrument that makes future selection work measurable.

---

## 8. Offline validation — MEASURED, no engine build/run (this is the evidence)

I reconstructed the engine's leader-clustering + ACF election **offline** from on-disk data and ran the exact operator swap. Data: the frozen `3dsig_full85_r10_cf_fix` pool (510 representative PDBs/target **with Cartesian coordinates** + full CF-component REMARKs) and the `3dsig_full85_scratch_3b2fa57cc` `.rrd` full-population dumps, all 19 ablation targets. Method faithful to `cluster.cpp` (2 Å leader capture seeded by ascending CF; ACF cluster ranking; within-cluster selector swapped). RMSD = symmetry-corrected, matching the report.

**Baseline reproduction is sound:** offline lowest-CF election median = **7.08 Å**, matching the report's engine `current_lowCF` median 7.29 Å and elected-mean 7.99 Å. This validates the harness.

| finding | value | meaning |
|---|---|---|
| medoid moves the pick | **334 / 822** multi-member clusters (41%) | the operator is active, not a no-op |
| ΔRMSD medoid − lowCF | mean **−0.005 Å**, median 0.00 Å | **neutral on the emitted pool** |
| ΔRMSD bmedoid − lowCF (HEAD variant) | mean −0.002 Å | even more conservative (CF weight pulls toward lowCF) |
| why neutral | emitted cluster diameter median **2.36 Å** (≈ 2 Å leader cutoff) | tight clusters ⇒ all members ≈ equal native-RMSD |
| emitted-pool ceiling | median **4.72 Å**, 0/19 sub-2 Å | no hit to recover by re-selection |
| **full-population `.rrd` ceiling** | median **3.79 Å**; 1OF6 down to **2.08 Å**; 3/19 have sub-3 Å | **near-native poses exist upstream, lost before emission** |
| within-cluster oracle headroom (full pop) | mean **1.35 Å**, p90 4.56 Å | a real selection prize exists — but only in the full population |

**Interpretation (the honest read):** the medoid swap is **correct and safe** (strictly geometric, reversible, changes the pick sensibly in 41% of clusters) but its measurable RMSD benefit on *today's emitted representative set is ~0 Å*, because those clusters are already 2 Å-tight. The lever only pays off once the pipeline stops discarding the near-native poses that DO exist in the raw GA population (full-pop ceiling 3.79 Å vs emitted 4.72 Å) — i.e. the payoff is gated on the representation/sampling fixes, not on this operator alone. The strongest independent value here is the **dump prerequisite** (§5), which turns "is selection the problem?" from unmeasurable into measurable.

Figure: `clustering_medoid_offline_validation.png` (panel a: ΔRMSD distribution, 41% picks change yet Δ̄≈0; panel b: per-target elected-vs-ceilings showing near-native lives in the full population, not the emitted pool). Tables: `offline_per_cluster_medoid_vs_lowcf.csv` (822 clusters), `offline_fullpop_rrd_headroom_19targets.csv`, `offline_election_emitted_pool_19targets.csv`, `offline_validation_summary.json`.

---

## 9. Risks & caveats

- **Neutral-on-current-pool is a feature-flag risk, not a correctness risk.** If a reviewer expects this to move 0/19, that expectation is wrong (see §8) and must be reset — the fix is a prerequisite/correctness change, not a success-rate change in isolation. Ship it for the constraint fixes + the instrument, and measure the RMSD payoff only after near-native survives clustering.
- **HEAD-commit regression surface.** HEAD's `cluster.cpp` medoid is default-ON, so anyone who built HEAD and ran a benchmark since `3e674479c` has NON-default-reproducing output already. Flipping to default-OFF will change HEAD's *default* output back — that is intended, but note it in the commit so a prior HEAD benchmark is re-run, not compared stale.
- **DP `#define` → runtime gate** slightly changes a hot constant into a branch; negligible (one `getenv` cached once). Keep the `Center`/`Representative` NULL-guards (`DensityPeak_Cluster.cpp:566`).
- **Medoid on non-superposed coordinates.** All poses share the receptor frame (redocking), so direct (non-superposed) heavy-atom RMSD is valid — matching how `cluster.cpp` already computes intra-cluster RMSD. Do not add Kabsch; it would be inconsistent with the clustering metric.
- **Symmetry.** The medoid uses direct heavy-atom RMSD (element order stable within a target), NOT Hungarian — consistent with `cluster.cpp`'s own clustering distance. Final rank-0 RMSD reporting stays Hungarian/symmetry-corrected as today. Do not conflate the two.
- **`.rrd` augmentation must stay behind `FLEXAIDDS_DUMP_POP`** so default benchmark artifacts are unchanged (bit-identity gate #1 covers PDB/.cad/.mcf; keep `.rrd` default-format identical unless the opt-in is set).

---

## 10. Effort estimate

- **IP-1 (CF/leader pure medoid + default-OFF gate):** ~40 LOC net (mostly deleting the weight vector from the existing block + inverting the default + adding the enum read). **Low complexity** — the geometry loop already exists and is correct.
- **IP-2 (DP center env-gate):** ~15 LOC (`#define`→runtime read + reuse existing `Center` branch). **Low.**
- **IP-3 (FO/BindingMode medoid mode):** ~30 LOC, **low-medium**, OPTIONAL (skip for v1; FO already elects a center).
- **IP-5 (dump augmentation):** ~50–80 LOC in `write_rrd.cpp` + 1 line in DatasetRunner. **Medium** (touches file format; keep opt-in). Highest value.
- **Gate centralization (`config_parser` enum):** ~20 LOC. **Low.**

**Total for v1 (IP-1, IP-2, IP-5, gate):** ~130–160 LOC, **low-medium** complexity, single-pass implementable. IP-3/IP-4 deferrable. No new deps, no build-system change.

---

## 11. Commit hygiene

Branch `feat/cluster-rep-medoid` off HEAD `3e674479c`, **no merge to main** until OPS clears gates 1–3. The commit must (a) remove the default-ON `FLEXAIDDS_MEDOID_REFINE` behavior, (b) NOT alter `Clus_ACF` / between-cluster ranking, (c) leave `.rrd` default format byte-identical unless `FLEXAIDDS_DUMP_POP=1`. State in the message that this supersedes the clustering portion of `3e674479c` and that flag-OFF was proven bit-identical to `3e674479c` + `FLEXAIDDS_MEDOID_REFINE=0`.
