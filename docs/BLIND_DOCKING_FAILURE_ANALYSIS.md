# Why autonomous blind docking collapses while the oracle ceiling stays high

**FlexAID∆S — root-cause analysis, HW/SW optimization plan, and pre-emptive bug list**
Author aid: investigation on `master`, 2026-07-08. All claims below carry `file:line` or committed-artifact evidence. Where a number is from a run I could not re-execute here, it is tagged **[snapshot]** (committed repo file) or **[live]** (result CSV read from disk this session). Nothing is estimated.

---

## 0. First-principles framing (the DOF ledger)

A rigid-plus-torsion docking search over 6 rigid DOF (3 translation, 3 orientation) + N torsions. The scoring/entropy machinery only matters *after* the search has placed the ligand in the right cubic ~10³ Å³ voxel of a ~10⁶ Å³ protein. So decompose "success" multiplicatively:

```
P(sub-2Å) = P(right pocket) · P(right translation | pocket) · P(right orientation+torsion | translation) · P(rank the good pose #1 | it was sampled)
```

The oracle ceiling and the blind result are **not measuring the same factors**. That is the whole story, and the code proves it.

---

## 1. The measured decomposition (what actually fails)

| Regime | Site handling | N | success <2Å | mean RMSD | Evidence |
|---|---|---:|---:|---:|---|
| Historical oracle ceiling (v31) | site given | 85 | **65 (76.5%)** | — | `CODEX_HANDOFF_v39_investigation.md` **[snapshot]** |
| Current oracle single-GA (v112) | site given | 85 | **38 (44.7%)** | 4.21 | `MULTI_CLEFT_RESTORATION_COMPARISON.md` **[snapshot]** |
| **Blind multi-cleft, score-selected** | autonomous | 85 | **4 (4.7%)** | 19.76 | same **[snapshot]** |
| **Blind multi-cleft, oracle-best cleft** | pick best cleft by RMSD | 85 | **5** | — | same **[snapshot]** |

Two independent gaps, and the second is the surprising one:

**Gap A — oracle regression: 76.5% → 44.7%.** Even *with the site handed in*, the engine lost ~20 targets. This is the documented v34→v39 cliff (cluster election + native-seed round-trip), see §3.5/§3.6.

**Gap B — blind collapse: 44.7% → 4.7%.** This is the user's question.

The decisive number is **oracle-best-cleft = 5/85** vs **score-selected = 4/85** (`summarize_multicleft_astex.py:98-101` computes both). If the failure were *pocket ranking*, giving the system a free oracle over which cleft to pick would recover most of the oracle 38 — it recovers **one** target (4→5). Therefore:

> **Blind docking is not failing at pocket *selection*. It is failing at pose *generation and registration inside pockets it already has* — including the correct one.** Even when the true cleft is among the tried clefts, a sub-2Å pose is essentially never produced.

That inverts the usual "blind docking = bad pocket picker" assumption. The pocket picker is fine; the machinery that turns a cleft into a search space is broken.

---

## 2. Why the oracle ceiling is "high" — it barely tests translation

The oracle path does **not** blind the translational DOF. Two mechanisms pin the ligand centroid to the true site:

- `write_blinded_ligand()` applies a uniformly-random SO(3) rotation **about the heavy-atom centroid and keeps that centroid exactly where the crystal put it** — `LIB/DatasetRunner.cpp:1480-1483` ("keep the heavy-atom centroid exactly where it is (site preserved) and apply a … rigid ROTATION about that centroid"). The comment is explicit that translating the ligand "strands the search in the wrong pocket."
- `SITE-CONFINE` then overrides the search centroid with the oracle/ligand centroid and sets the translational cutoff `rmax2` from the ligand extent — `LIB/top.cpp:1756-1759` (`using_oracle → cx,cy,cz = oracle_cx,…`) and the expanding-radius confinement immediately below (~`top.cpp:1770+`).

So the oracle "ceiling" is really an **orientation + torsion + scoring ceiling with translation ≈ given**. It is a legitimate number for a *re-docking* benchmark, but it is **not** an upper bound on blind performance and should never be quoted as one. The 76.5%→4.7% "gap" is, to first order, *the cost of actually having to solve the 3 translational DOF and localize the pocket* — which the blind path does through a completely different, and broken, code route.

**Exact-science caveat to bake into every writeup:** report the oracle number as "oracle (translation-pinned) ceiling," not "oracle ceiling." Otherwise the 76.5% is an apples-to-oranges denominator.

---

## 3. Root causes, ranked by confidence

### 3.1 [HIGH] The blind path skips site-confinement *and* keeps the IC frame anchored at the whole-receptor centroid

In `top.cpp`, the entire SITE-CONFINE block is gated `if (!using_explicit_cleft)` — `LIB/top.cpp:1720` and the reflig/oracle overrides at `:1735-1759`. In blind multi-cleft mode (`using_explicit_cleft == true`, set at `:1655`), none of that runs. The cleft becomes the grid verbatim from `generate_grid(FA, spheres, …)` (`top.cpp:1710`), then `calc_cleftic()`.

But the internal coordinate origin is **`FA->ori` = geometric centroid of *all* protein atoms**, computed in `calc_center()` (`LIB/calc_center.cpp:22-38`; it averages every ATOM over every residue). The cleftgrid ICs (`calc_cleftic`) and `buildcc`'s GPA1/GPA2 grandparent reference are all encoded relative to that whole-protein centroid — the code says so at `top.cpp:1527-1532`: *"the cleftgrid IC (calc_cleftic) are encoded relative to this receptor-center ori … gene[0] translates GPA0 … Overwriting … breaks that reference frame when gene[0] translates GPA0 far from the ligand starting position."*

Consequence, from first principles: for an off-center cleft (the common case — surface pockets sit far from the protein's center of geometry), the translational gene `gene[0]` must express a large displacement from `FA->ori` to reach the cleft. In oracle mode, `rmax2`/`rcut` re-anchor and bound that displacement around the true centroid. In blind mode that bounding is **absent**, so:
- the translational search sphere is effectively centered on the protein centroid, not the cleft;
- sampling density (poses per Å³) in the actual cleft voxel drops by the ratio of (search volume)/(cleft volume);
- the GA spends its budget roaming, and the correct translation is rarely sampled — matching the observed mean RMSD of 19.76 Å (i.e., poses land ~halfway across the protein).

**This is the primary Gap-B driver.** Fix in §5.1.

### 3.2 [HIGH] Native-seed IC↔Cartesian round-trip is broken → search space is misregistered

`CODEX_HANDOFF_v39_investigation.md` smoking gun #3: `NATIVE-SEED-RMSD` = 4.62 Å (1HNN), 5.46 Å (1HP0), 2.87 Å (1GPK) in v35+, but **0.00 Å in v34_ctrl**. That value is the RMSD after encoding the crystal pose into gene space and decoding it back. A non-zero round-trip means the gene↔Cartesian map is not an identity on the true pose — the coordinate frame the GA searches in is *rotated/reflected/translated* relative to the crystal frame. Emitted from `LIB/native_score.cpp` (the file flagged dirty — read only).

Combine 3.1 and 3.2: even if a cleft perfectly overlaps the true site, a pose that is sub-2Å in Cartesian space does **not** correspond to a reachable, low-energy point in the (mis-registered) gene space. That is exactly why **oracle-best-cleft is only 5/85** — the machinery cannot represent the right answer, independent of pocket choice. The handoff points at the latent reflection bug in the `ictogene`/`genetoic` crossover descendants and the GPA0/grid-centroid construction as the suspects. This must be fixed before any blind number is meaningful.

### 3.3 [HIGH] Thermodynamic / entropy ranking channel is inert (zeroed)

Live confirmation this session, `results/astex_jcim2015_fair_20260708_0002/1HNN/result.csv` **[live]**:
`predicted_TdS = 0.0000`, `shannon_entropy = 0.0000`, and `predicted_dH == predicted_dG == best_score == -189.86` (dH is just being aliased to the raw CF, not an enthalpy). So the ∆S engine that is supposed to break ties between clusters is contributing **zero signal**. This reproduces smoking gun #1 from the handoff on a fresh 2026-07-08 build. When ∆S = 0, `G_bind = H − TΔS` collapses to `H`, and cluster election reduces to raw complementarity function — the entropy-driven part of "FlexAID∆S" is switched off in practice.

### 3.4 [HIGH] Cluster election reports a worse pose than it found

Same 1HNN row **[live]**: `best_cluster_rmsd = 7.35 Å` but the **elected** `rmsd_to_crystal = 14.54 Å`. The system generated a 7.35 Å cluster and then handed back a 14.54 Å pose as rank-0. This is handoff smoking gun #2, still live. `cluster.cpp:230-241` now re-sorts clusters ascending by representative `chrom[Clus_TOP[b]].evalue` (the lowest-CF fix, commit `cd9004d`), and that sort *is* present — so the surviving election error is upstream: the entropy-augmented `ACF` (`cluster.cpp:154-155`) still feeds `QuickSort_Clusters` at `:216` before the CF re-sort, and with §3.3 zeroing the entropy term, the ACF ranking is degenerate/noisy. Net: on this target the rank-0 pose is neither the lowest-CF nor the lowest-RMSD cluster. (Note: 7.35 Å is itself a failure, so §3.1/§3.2 dominate here; election error is an independent, additive loss visible on targets where a good cluster *does* exist.)

### 3.5 [MED] Pocket detection returns the *largest* cleft, not the *ligandable* one

`CleftDetector.cpp:23` — "the largest cluster is returned as the binding cleft." SURFNET/GetCleft largest-void ≠ cognate site for a large fraction of Astex (well-known: the biggest cleft is often an inter-domain groove or crystallographic artifact). The multi-cleft campaign mitigates this by running top-3 clefts, which is why *ranking* costs only 1 target (§1) — but coverage still caps the achievable set. This is a secondary contributor and only becomes rate-limiting **after** 3.1–3.2 are fixed.

### 3.6 [MED] `read_spheres` column parsing is fragile — silent grid corruption

`LIB/read_spheres.cpp:44-48` reads the sphere radius from fixed columns `buffer[61..65]` (the B-factor field) with **no length check** on `buffer` and no validation that the field is numeric. GetCleft writes radius there, so the happy path works, but: (a) a sphere file with radius in the occupancy column (55-60), or short lines, yields `radius = 0` → degenerate/collapsed spheres → `generate_grid` produces a near-empty or single-point grid → the GA searches nothing. There is no post-parse sanity gate (`num_grd`/radius distribution) before docking. This is both a robustness bug and a reproducibility landmine across GetCleft versions.

---

## 4. Cross-cutting: a memory/perf bug that is already half-patched (keep the guard)

`DatasetRunner.cpp:2065-2069`: mmCIF read through the fixed-column PDB parser misplaces coordinates into the 0–999 Å range → Vcontacts box `dim³ = 333³ = 37M` → ~2 GB per worker. The current code prefers `.pdb` and only falls back to CIF, which dodges it. **Do not remove that preference**, and add the explicit guard in §5.6 so a future CIF-only target cannot silently OOM the box or, worse, run with garbage coordinates and report confident nonsense.

---

## 5. Fixes (software) — ordered by expected Gap-B recovery per unit effort

**5.1 Re-anchor and bound the translational search in blind mode (biggest win).**
In the `using_explicit_cleft` branch of `top.cpp`, after `generate_grid`/`calc_cleftic`, run the *same* confinement used in oracle mode but with the **cleft centroid** as the anchor: compute `(cx,cy,cz)` = mean of the cleft sphere centers, set `rmax2` from the cleft sphere extent (max sphere-center distance + max radius), and apply the expanding-radius grid trim. Concretely: lift the SITE-CONFINE block out of `if(!using_explicit_cleft)` (`top.cpp:1720`) and feed it the cleft centroid instead of the reflig/oracle centroid. This makes the translational gene range O(cleft size), not O(protein size), restoring sampling density in the correct voxel. Expected to move blind from ~5% toward the pocket-coverage ceiling.

**5.2 Fix the native-seed round-trip before trusting any blind number (§3.2).**
Bisect the v31→8c0c840 range restricted to: GPA0 construction, grid centroid, `ic2cf`/`ictogene`/`genetoic`, native-seed placement, cleft-relative transforms (the handoff's exact list). Acceptance gate: `NATIVE-SEED-RMSD < 0.1 Å` on 1HNN/1HP0/1GPK, restored to the v34_ctrl 0.00 Å. Add this as a **hard CI assertion** (see §7) so it can never silently regress again — this is the bug that already cost 20 oracle targets.

**5.3 Un-zero the ∆S channel (§3.3).**
Trace where `predicted_TdS`/`shannon_entropy` are written in `DatasetRunner.cpp` back through the Boltzmann weighting; the entropy stack is producing values (`TdS_shannon = 4.35`, `TdS_vib = -0.0`, `H_pop`, `H_rep_*` columns are populated in the 1HNN row) but they are not propagating into `predicted_TdS`. This is a wiring bug (populated intermediate, zeroed output field), not a physics bug — cheap and high-value because it re-enables the entropy tie-break that is FlexAID∆S's entire thesis.

**5.4 Make cluster election deterministic on the intended objective (§3.4).**
Decide the election objective explicitly: if the product is ∆G-ranked, election must sort on `G_bind = H − TΔS` *after* 5.3 restores TΔS; if it is CF-ranked, drop the ACF `QuickSort` (`cluster.cpp:216`) entirely and keep only the evalue re-sort (`:230-241`). Right now you sort twice on two different objectives with a degenerate second key. Add a unit test asserting `rank0.rmsd ≤ min(cluster.rmsd) + ε` on a synthetic ensemble.

**5.5 Pocket coverage: replace largest-cleft with druggability-ranked top-K (§3.5).**
Rank clefts by a volume × enclosure × hydrophobic-contact score (cheap, already have the atoms) rather than raw volume, and keep top-3–5. This lifts the coverage ceiling once 5.1/5.2 make in-pocket docking work again.

**5.6 Harden `read_spheres` (§3.6).**
Add: length check (`strlen(buffer) >= 66`), `strtod` with error check for the radius, reject `radius <= 0.3 Å`, and a post-parse gate `if (n_spheres < MIN || radius_p95 < 0.5) Terminate(...)` with a real message. Fail loud, never dock an empty grid.

---

## 6. Hardware / performance optimizations (tie each to the physics, not to hype)

The blind problem is *sampling-limited*, so the right HW play is **more independent search at fixed wall-clock**, plus removing the memory cliff. In rough order of leverage:

1. **Grid-decomposed multi-start across clefts (embarrassingly parallel).** Blind mode already runs one process per cleft. Make that a first-class parallel-dock path (`top.cpp` advertises `--parallel-dock` at `:367`) with a shared read-only receptor grid via the Strategy-A on-disk cache (`FLEXAIDDS_GRID_CACHE_DIR`, `top.cpp:107-223`). One receptor grid built once, K clefts × M restarts fanned across cores. This directly raises `P(right translation|pocket)` at fixed time — the factor that is currently starving.

2. **SIMD the Vcontacts/CF inner loop (AVX2 today, AVX-512 opt-in).** `VoronoiCFBatch.h` (std::span batch interface) + `flexaids_configure_simd` already exist. Vectorize the per-pose contact accumulation over the batch axis; the CF evaluation is the GA hot loop, and every doubled eval throughput converts 1:1 into more generations / more restarts. Verify with `benchmark_vcfbatch`.

3. **GPU histogram for the Shannon configurational-entropy stack.** `ShannonThermoStack/` already has CUDA (`shannon_cuda.cu`) and Metal (`shannon_metal.metal`) histogram kernels behind a dispatch layer. Once §5.3 re-enables ∆S, the per-generation entropy binning becomes a real cost; keep it off the CPU critical path. This is throughput-neutral for correctness but lets you afford entropy every generation instead of only gen500/gen1000 (note the `TdS_shannon_gen500/1000 = NA` columns — entropy is currently sampled sparsely).

4. **Kill the box-dim memory cliff (§4).** Add the coordinate-range guard so a mis-parsed structure can never allocate a 2 GB `333³` box. This is a stability/throughput fix: one OOM worker stalls a whole benchmark shard.

5. **tENCoM vibrational entropy on GPU only when it changes ranking.** `tENCoM` has CUDA/Metal bridges; gate the expensive Hessian diagonalization behind a "does ΔS_vib actually reorder the top-K clusters?" check so you pay for it only on the ~top few poses, not the whole population.

**Do not** reach for an ML pocket-predictor to paper over §3.1/§3.2 — that would hide a coordinate-frame bug behind a black box and make the eventual failure unreproducible. Fix the mechanics first; the mechanics are cheap.

---

## 7. Pre-emptive bug / exploit / footgun list (with the fix)

| # | Class | Where | Risk | Fix |
|---|---|---|---|---|
| B1 | Silent wrong-frame search | `native_score` IC round-trip §3.2 | Every blind result invalid; already cost 20 oracle targets | CI gate: `NATIVE-SEED-RMSD < 0.1 Å` blocks merge |
| B2 | Zeroed output field | `predicted_TdS`/`shannon` §3.3 | ∆S product ships with ∆S=0 | Wire intermediate→output; assert non-zero on ≥1 known target |
| B3 | Unbounded translational search | `top.cpp:1720` blind branch §3.1 | Blind ≈ random over protein | Cleft-centroid confinement (§5.1) |
| B4 | Parser OOB / silent 0-radius | `read_spheres.cpp:44` §3.6 | Empty grid, docks nothing, no error | Length+numeric+range checks (§5.6) |
| B5 | Coordinate misparse → 2 GB box | `DatasetRunner.cpp:2065` §4 | Worker OOM or confident garbage | Range guard + keep PDB-preferred |
| B6 | Double-objective election | `cluster.cpp:216` vs `:230` §3.4 | Reports worse-than-found pose | Single objective + unit test (§5.4) |
| B7 | Benchmark-integrity | oracle quoted as blind ceiling §2 | Overstated headline number | Label "translation-pinned oracle" everywhere |
| B8 | Reproducibility drift | GetCleft sphere column convention §3.6 | Same input, different grid across tool versions | Pin sphere-file schema; assert on read |
| B9 | 0-oracle silent run | `DatasetRunner.cpp:3814` | (already fixed — keep) fails loud if `FLEXAIDDS_ORACLE_SITE_DIR` unset | Leave the `[FATAL]` guard in place |
| B10 | Float32 centroid accumulation | `calc_center.cpp` sums into `float` | Large multimers lose precision in `ori`, shifting the IC frame origin | Accumulate in `double`, divide once |

B9 is already correctly implemented (`DatasetRunner.cpp:3809-3816`) — an explicit abort when an Astex run has 0 oracle sites, which is exactly what silently caused the v35–v39 cliff. Preserve it and add the B1/B3 analogues for blind mode (abort if cleft grid < MIN points or cleft centroid confinement was skipped).

---

## 8. Bottom line

The oracle ceiling is high because oracle mode **pins the 3 translational DOF** (`DatasetRunner.cpp:1480`, `top.cpp:1756`) — it measures orientation + torsion + scoring, not localization. Blind mode has to solve translation and pocket ID through a separate route that is broken in three compounding ways: the translational search is **unbounded and anchored at the whole-protein centroid** in the cleft branch (§3.1, `top.cpp:1720`/`calc_center.cpp`), the **native-seed IC↔Cartesian frame is mis-registered** so the correct pose isn't even representable (§3.2, `NATIVE-SEED-RMSD` 2.9–5.5 Å), and the **∆S ranking channel is zeroed** (§3.3, live 1HNN `TdS=0`). The proof that it's these and not pocket-picking: giving the system a free oracle over *which* cleft to use recovers exactly **one** target (4→5 of 85). Fix §5.1 and §5.2 first — they are cheap, mechanistic, and gate everything downstream — then re-enable ∆S (§5.3) and clean up election (§5.4). Only then does a pocket-coverage or ML play (§5.5) buy anything real.
