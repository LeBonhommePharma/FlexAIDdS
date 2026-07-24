# Fair Comparative Benchmarking Methodology — FlexAIDdS vs FlexAID-2015 vs FlexAID-first-entropy

**Goal.** Compare, fairly and objectively, three points on the FlexAID lineage on the Astex Diverse 85 redocking benchmark, isolating the effect of the entropy concept and of the FlexAIDdS extensions.

**Status.** Draft protocol (2026-07-24). Version pins finalized below; execution pending.

---

## 1. The three engines (pinned to the FlexAID lineage)

Reference repo: `LeBonhommePharma/FlexAID` (mirror of `NRGlab/FlexAID`), local clone at `../flexaid` / `../FlexAID`. 377 commits, 2012–2020, branches `master` and `entropy` (originally `entropize`).

| Engine | Definition | Ref | Notes |
|---|---|---|---|
| **V1 — FlexAID 2015 (JCIM baseline)** | Published CF + GA + FastOPTICS methodology, **no entropy scoring** | `master` branch, paper-era commit (mid-2015) | The Gaudreault & Najmanovich 2015 methodology (doi:10.1021/acs.jcim.5b00078). C code. |
| **V2 — FlexAID first-entropy** | First working `BindingPopulation::Entropize()` + `SolvatedPartitionFunction` | `entropy` branch, first stable-entropy commit (~2016) | The direct conceptual ancestor of FlexAIDdS's Shannon/partition-function pose entropy. C/C++. |
| **V3 — FlexAIDdS** | Current entropy-driven fork | `9dc9-production` (PR #300, clean 9dc9 matrix + validated fixes) | C++26. Adds ΔG_eff, ΔS_vib (tENCoM), Shannon pose entropy, PoseBusters QC. |

> **Pin discipline.** Each engine is built once from an immutable commit SHA recorded in the run receipt. V1 must be the *published* state — before entropy was wired into scoring — so exact commit selection is verified by confirming `Entropize()` is NOT in the CF/selection path (it may exist as dead scaffolding on `master`, which is acceptable only if unused).

---

## 2. Fairness controls (the crux of objectivity)

A comparison is only meaningful if every difference in the *result* is attributable to the *engine*, not to prep, protocol, scoring, or compute noise. Controls:

1. **Identical dataset + preparation.** ONE canonical Astex-85 input set (same apo receptors, same ligand SDFs, same defined binding-site definition) fed to all three engines. This is non-negotiable: the session already showed prep confounds (un-stripped HEM on 1G9V/1Q4G) can dominate. All engines dock the *same* atoms.
2. **Identical search protocol.** Blind redock, defined-cleft (no native-pose seeding, `seed=OFF`), identical search budget (population, generations, restarts) mapped to each engine's equivalent parameters. Document any parameter that cannot be matched 1:1.
3. **N random seeds (≥3, ideally 5).** Docking GAs are stochastic; a single seed is not a measurement. Report mean ± dispersion across seeds; this is what makes the comparison a *measurement* rather than an anecdote.
4. **Uniform EXTERNAL scoring.** Success is computed by ONE external RMSD scorer (symmetry-corrected, e.g. `spyrmsd`/Hungarian) applied identically to every engine's emitted rank-0 pose vs the crystal ligand. **Never** trust an engine's self-reported success column (the session showed FlexAIDdS's `success`/`claim_ready` are QC-gated and not comparable). Denominator frozen at 85 (a target that crashes or emits nothing = failure, never removed).
5. **Matrix control — TWO arms.** The interaction matrix differs across versions and is a confound:
   - **Arm A (as-published):** each engine uses its OWN native/validated matrix. Answers "which released engine docks best."
   - **Arm B (matrix-held-constant):** all engines use the SAME matrix (9dc9). Isolates the *algorithm/entropy* changes from the *matrix* changes.
   Report both; their difference quantifies how much of any gap is matrix vs algorithm.
6. **Serialized compute (no contention).** Runs execute ONE engine at a time on the same hardware, pinned threads, determinism flags where available. Box contention was the single largest confound this session — it must be eliminated, not managed.
7. **Genuine-blind enforcement.** All engines verified to run without native-pose seeding / seed-echo (rmsd≈0 with best_score≈native = disqualified). The metric is genuine top-1, not seed-inflated.

---

## 3. Primary metric & statistics

- **Primary:** genuine top-1 sub-2Å success rate = (# targets whose rank-0 elected pose has symmetry-corrected RMSD < 2.0 Å) / 85, with **binomial 95% CI** (Wilson).
- **Paired analysis:** per-target win/loss matrix across engines; **McNemar's test** for each pairwise comparison (V1↔V2, V2↔V3, V1↔V3) — the correct test for paired binary outcomes on the same targets.
- **Averaging seeds:** per target, success = majority (or mean sub-2Å fraction) across the N seeds; report both the pooled rate and the per-seed spread.
- **Secondary/diagnostic:** full RMSD distribution (not just the 2 Å threshold — report 1/2/2.5/3 Å bins), near-native sampling frequency, wall-clock per target, and pose validity (PoseBusters) as a *separate* axis (not folded into the success metric).

---

## 3b. Shared-search design for V1↔V2 (paired; do NOT search twice)

**Finding (verified against the code, 2026-07-24).** Between `master` (V1) and `entropy` (V2):
- `Entropize()` is invoked from `FOPTICS.cpp` (clustering) — the entropy is a **post-GA re-ranking of the converged population**, never a search-time fitness.
- `gaboom.c` differs by +418 added lines, but those are **new standalone functions** (`generate_genetic_variants`, `generate_single_gene_mutants`, `generate_multiple_genes_mutants`, `generate_true_positive_cluster`, `generate_true_negatives_clusters`) — clustering *validation scaffolding*, not the search loop. Only **1** added line mentions entropy/BindingMode. Core operators are untouched (`mutate`/`boom`/`populate` diffs are `genlim*` → `const genlim*` signature propagation).
- V1 elects via `cluster.c` (CF-ordered); V2 adds `entropy_cluster.c` (+134 lines) and reorders clusters "after considering entropy" (`505d764`).

**⇒ Consequence: run the GA ONCE, score it TWICE.** Because the search is common and the entropy acts only at election, running two independent GAs would make part of any V1↔V2 difference pure stochastic noise. Instead:

1. Run the search **once** per (target, seed) using the V2 binary.
2. Emit the **full final population** with per-pose coordinates and per-pose CF (FlexAID already writes ranked results; ensure the whole population, not just rank-0, is dumped).
3. Apply **both election rules to that same population, offline in the analysis layer**:
   - **E_CF** (V1-style): elect lowest-CF cluster head.
   - **E_S** (V2-style): entropy-weighted cluster election (`Entropize()` semantics).
4. Score both elected poses with the same external RMSD scorer.

**Why this is better, not just cheaper:** it is a **paired** design — identical search, identical population, identical prep, identical scorer. The *only* difference is the election rule, so the measured Δ is attributable to the entropy concept alone. It also halves compute and removes GA-seed variance from the most important comparison. Applying the election rules offline (rather than in-engine) keeps both rules transparent, auditable, and identically implemented.

**Scope limits (stated honestly):**
- This arm answers **"does the entropy concept improve pose election?"** — the core scientific question.
- It does **not** reproduce "V1 exactly as published" (V2's binary still carries the const-correctness/scaffolding deltas). The as-published engine comparison remains a separate, lower-priority arm requiring V1's own binary.
- **V3 (FlexAIDdS) cannot share this search** — it is an independent C++26 engine with its own GA (coarse-init, elitism, restarts) and its own entropy that *does* enter fitness (SMFREE/ACF sharing). V3 must run its own search and be compared unpaired (hence §2.3 N seeds + §3 CIs).

**Risk de-scoped:** because the entropy question no longer depends on a working V1 binary, the hardest feasibility risk (building 2015-era C on macOS/arm64) is removed from the critical path.

## 3c. VERIFIED build + output facts (2026-07-24) — engines are executable on darwin/arm64

**⚠ LANDMINE — never build from the clones' working trees.** Both `../flexaid` and `../FlexAID` have **33 uncommitted modified files** (a prior session's WIP) that **comment out the entropy dispatch in `top.c` and strip ~455 lines from `gaboom.c`**. A build from the dirty tree silently produces an entropy-less binary (283 KB, zero `Entropize` symbols). **Always build from the pinned SHA via `git worktree add`.**

| | V1 (2015 baseline) | V2 (first-entropy) |
|---|---|---|
| Pinned SHA | **`b555e0e`** (Makefile portability fix; master HEAD `9aa7995` is byte-identical apart from a README) | **`1a6ae0b`** (entropy HEAD) — methodology first functional at **`312f0c9`** (2017-02-23, "reranked according to CFdS instead of CF in entropy_cluster()"), built at HEAD only because the 2017 commits predate the Makefile portability fixes |
| Entropy present? | **No** — `BindingMode.cpp/h` exist in-tree but are *not compiled or linked*; `top.c` calls classic `cluster()` only. Binary: **0** `Entropize` symbols | **Yes** — links `BindingMode/FOPTICS/ColonyEnergy/entropy_cluster`. Binary: **1** `Entropize` symbol |
| Built / docks | yes / yes (283 KB) | yes / yes (467 KB) |

**Build (both, ~10 s, warnings only):**
```bash
git worktree add /tmp/flexaid_v2_pristine 1a6ae0b     # or b555e0e for V1
cd /tmp/flexaid_v2_pristine/BIN
make -f Makefile.Linux64 CXX=g++ BOOST_INCLUDES=/opt/homebrew/include -j4
```
**Only macOS adaptation** (CLI override, no source edits): `BOOST_INCLUDES` → `/opt/homebrew/include` (Apple-Silicon homebrew; Boost is header-only here, so no link libs). Apple clang aliases `g++` and targets arm64 natively.
Staged binaries: `~/flexaid_ref_binaries/FlexAID_V{1_master,2_entropy}`.

**Full-population output (enables the offline paired election).** Run V2 with `TEMPER 0 / CLUSTA CF` to get the raw un-elected population from a single search:
- **`<out>.rrd`** — every GA survivor (`NUMCHROM` lines): `idx  clusterID  clusterRMSD  RMSD_ref  CF  [genes]`. Requires `RMSDST` set.
- **`<out>.cad`** — per cluster: `TOP=<idx> TCF=<lowest CF> ACF=<avg CF> freq=<size>` (TCF order **is** the E_CF election order).
- **`<out>_N.pdb`** — Cartesian representative (lowest-CF member) per cluster + CF breakdown (`CF.app/com/sas/wal/con`).
- Config: **⚠ `MAXRES 0` writes ZERO results, not "all"** — empirically verified on identical inputs (`MAXRES 0` → `num_of_results=0`, only `_INI.pdb`; `MAXRES 500` → 500 results/PDBs). **Use `MAXRES = NUMCHROM`.** `NUMCHROM` = population size; `TEMPER 0`→CF election, `TEMPER>0`→entropy; `CLUSTA CF|FO`. Also set `OUTGENER 50` (at `OUTGENER 1` logs are ~35% of the disk footprint). Invocation: `FlexAID CONFIG.inp ga.inp <out_prefix>`.
- *Caveat:* `_N.pdb` holds only the per-cluster representative; all-pose Cartesians require `.rrd` (internal coords) or enabling `output_dynamic_BindingMode()` (`BindingMode.cpp:50`).

### The two election rules (reference semantics, to be reimplemented offline)
From `BindingMode.cpp`, per pose *i* with `CF_i = chrom->app_evalue` and `T = TEMPER` (integer):

```
w_i = exp(−(1/T)·CF_i)          # β = 1/T  — NOT 1/(k_B·T)
Z   = Σ_{i ∈ whole population} w_i
per binding mode (cluster) m:
  H_m = Σ_{i∈m} (w_i/Z)·CF_i
  S_m = −Σ_{i∈m} (w_i/Z)·ln(w_i/Z)        # Shannon, no reference state
  G_m = H_m − T·S_m                        # compute_energy()
```
- **E_S (V2-style):** `Entropize()` sorts modes **ascending by `G_m`**; `mode[0]` wins. Written structure = the mode's lowest-CF member (`elect_Representative(false)`).
- **E_CF (V1-style):** rank modes by raw CF, lowest wins (= `.cad` TCF order).

> **Note the normalization:** `w_i/Z` is normalized over the **whole population**, so within-mode probabilities do not sum to 1 — this must be reproduced exactly. **β = 1/T (no k_B)** matches FlexAIDdS's convention (adding k_B collapses the weights → S≡0), and the reference's `TEMPER 21` is the same calibration as FlexAIDdS's "ISMB 2017 T=21".

## 3d. Integration with the existing harness (`scripts/generate_flexaid_inp.py`)

A mature generator already exists and defines the classic-FlexAID arms. **Reuse it**; do not rebuild.

| Arm | Binary | TEMPER | CLUSTA | Meaning |
|---|---|---|---|---|
| **A** | 2015 pin (`b555e0e`) | 0 | CF | FlexAID 2015 as published |
| **B0** | master/entropy build | 0 | CF | entropy build with entropy OFF |
| **B** | entropy (`1a6ae0b`) | 21 | FO | entropy ON (as shipped) |
| **C** | entropy | 298 | FO | high-T variant (gated) |
| **C0** | FlexAIDdS | — | — | separate runner (not this generator) |

It already encodes good fairness practice: `DEFAULT_MAXRES = 50` deliberately matches the FlexAIDdS cluster-emit ceiling, one canonical prep per target, and clean-apo/ligand-integrity gates.

**Two issues resolved / flagged:**

1. **The "A == B twin" hazard is now resolved.** The generator warns that B0 is *"twin of A when bin A==B SHA — not an independent control"*. The build audit settles it: `b555e0e` (master) links **no** entropy — `BindingMode.cpp/h` are present but never compiled, `top.c` calls classic `cluster()` only, binary has **0** `Entropize` symbols (283 KB) — while `1a6ae0b` (entropy) has **1** (467 KB). Pinning A→`b555e0e`, B→`1a6ae0b` gives a genuine, verifiable binary difference.

2. **⚠ Arm B changes TWO variables at once.** B differs from A/B0 in *both* `TEMPER` (0→21, entropy ranking) **and** `CLUSTA` (CF→FO, Fast-OPTICS density clustering). A B-vs-A difference therefore cannot be attributed to entropy alone. Hence the two complementary readouts:
   - **Q1 — "the entropy concept" (isolated):** the §3b paired offline election — ONE search, fixed CF clustering, `E_CF` vs `E_S` on the same population (`scripts/elect_paired.py`). Only the election rule varies.
   - **Q2 — "the entropy pipeline as shipped":** arm **B** vs arm **A/B0** — entropy ranking *plus* FO clustering, i.e. the engine as released. Legitimate, but reported as a pipeline effect, not an entropy effect.

   Both are reported; their difference measures how much of B's effect is FO clustering rather than entropy.

3. **Seeding caveat (affects replicate design).** The generator notes the staged Mach-O binaries **seed from `time(0)` and may ignore `STRTSEED`**. Consequence: runs are *not* bit-reproducible, but independent invocations are genuinely independent replicates — so N-seed replication works by simply re-invoking, while exact reproduction requires archiving outputs (record this in the receipt).

**Validation status of the election layer:** `scripts/elect_paired.py` was checked against real engine output (1GPK smoke, 100 poses → 21 modes): its `E_CF` winner reproduces the engine's own `.cad` rank-0 (`Cluster 0 TOP=0 TCF=-30.09`) exactly, and log-sum-exp weighting absorbs the clash-sentinel poses (CF ≈ 10⁴, visible as `ACF=19608`) without underflow.

## 3e. PILOT RESULTS + blockers (2026-07-24, 8 targets × 3 seeds, V2 @ 500 chrom × 500 gen)

**Sizing (measured):** median **30 s/run** single-threaded (min 18, max 42; n=32). Per arm: 85×3 seeds ≈ 2.1 CPU-h (~42 min wall @3 concurrent); 85×5 seeds ≈ 3.5 CPU-h (~70 min). Full matrix (3 arms × 85 × 5) ≈ **5 CPU-h, <2 h wall, ~8 GB disk**. **Disk is the binding constraint, not CPU** (~19 MB/run; 12 GB free — free ≥8 GB before launching).

**Early signal (pooled 24 obs):** `E_CF` 3/24, `E_S` 3/24. The elected pose **differs on 12/24 (50%)** — the entropy rule is genuinely active, *not* a no-op — but discordant outcomes are 1 vs 1, **McNemar exact p = 1.00**. ⇒ *At this n, entropy election changes the pose but not the outcome.*

**Two facts that constrain any conclusion:**
- **Seed noise dominates.** Per-target E_CF RMSD spread across 3 seeds reaches **7.84 Å** (1P62), 6.66 (1SJ0), 3.93 (1GPK). A single seed is uninterpretable → **N=5 seeds is mandatory**, and the V1-vs-V2 engine table at n=1 carries **no information** (do not quote it).
- **Sampling ceiling:** a sub-2 Å pose exists in the population for only **3–4 of 8** targets. Studying election rules on populations where the answer is absent ~60% of the time has little headroom — **raise `NUMCHROM` before drawing election conclusions** (cheap at 30 s/run).

**Blindness re-verified:** `POPINIMT RANDOM`, "generated 500 randomized individuals", and **zero poses < 0.05 Å across all 32 runs** — no seed echo. Binary pins re-verified by `nm -a`: **5** `entrop*` symbols in V2, **0** in V1.

### ⛔ P0 blockers — the entropy claim is NOT yet defensible
1. **Population mismatch (fidelity gap).** The `.rrd` exposes only the final `NUMCHROM` chromosomes, but `Entropize()` elects over the **whole GA snapshot**: `.cad` cluster frequencies are 139–25,172 vs 500 `.rrd` rows (1GPK cluster 0: **209 members in-engine, 4 in the `.rrd`**). Since `S_m` depends on mode membership, the offline `E_S` is reconstructed from **~2%** of the poses the reference actually uses ⇒ **the §3b paired election is currently a PROXY, not a reproduction.** *Fix:* enable the reference's own `output_dynamic_BindingMode()` (`BindingMode.cpp:50`, one-line uncomment — instrumentation only, writes each Pose as a MODEL with its CF), rebuild V2, and re-derive `E_S` from the full snapshot.
2. **Internal RMSD is not symmetry-corrected.** V1's Hungarian column quantifies the cost: 1MEH 3.99→**2.01**, 1L7F 4.78→**2.94**, 1OWE 9.83→**5.90**. Internal RMSD systematically overstates error and would flip success calls. **§2.4's external symmetry-corrected scorer (spyrmsd) on the emitted `_N.pdb` is mandatory — never use the `.rrd` column.**

### Fixed in this round
- `elect_paired.py` is now **column-layout aware** (V1 writes 6 numeric cols with CF **last** and a Hungarian RMSD at col 4; V2 writes 5 with CF at col 4 — the previous hardcoded `cf=parts[4]` silently read V1's *RMSD* as its CF), and it **drops `cluster < 0`** poses (MAXRES-overflow poses were fabricating one spurious mega-mode of up to 376/500 poses).
- `MAXRES 0` corrected above; `OUTGENER 50` recommended.
- `SOFTWA 0.40` in the canonical prep is **inert for both classic binaries** (token absent from both) — harmless for V1↔V2 fairness, but a documented non-1:1 vs FlexAIDdS.

## 4. Execution plan

1. **Build** V1 (`master`, C, Makefile adapted for this platform), V2 (`entropy`, C/C++), V3 (FlexAIDdS, CMake C++26). Record each build's compiler, flags, and commit SHA.
2. **Freeze inputs**: one canonical Astex-85 prep (receptors/ligands/sites), checksummed, shared read-only by all engines.
3. **Run**: each engine × 85 targets × N seeds, blind, serialized. Isolated output dirs, local storage.
4. **Score**: uniform external RMSD scorer over every rank-0 pose. Build the 3 × 85 × N outcome table.
5. **Analyze**: success ± CI per engine (Arms A & B), McNemar pairwise, RMSD distributions, per-target win/loss heatmap.
6. **Report**: a single results table + significance + the "how much is matrix vs algorithm vs entropy" decomposition.

---

## 5. Threats to validity (and their controls)

| Threat | Control |
|---|---|
| Prep differences (HEM, protonation, site) | §2.1 one canonical prep |
| Matrix differences | §2.5 two arms (native vs constant) |
| Seed-echo inflation | §2.7 blind enforcement + §2.4 external scoring |
| Stochastic single-seed noise | §2.3 N seeds + §3 CI |
| Box contention / CPU starvation | §2.6 serialized runs |
| Self-reported-metric incomparability | §2.4 uniform external scorer |
| Old-C vs C++26 numerical drift | document compiler/flags; the metric is geometric (RMSD), not score-value, so scale differences don't bias it |
| Parameter non-equivalence across engines | §2.2 document every non-1:1 mapping |

---

## 6. Open decisions (require sign-off before execution)

1. **V1 exact commit** — the paper-era `master` SHA (verify entropy not in scoring path).
2. **V2 exact commit** — the first *stable* `Entropize()` state to pin.
3. **N seeds** — 3 (faster) vs 5 (tighter CI).
4. **Which arm first** — Arm A (as-published) is the headline; Arm B (matrix-constant) is the mechanistic decomposition. Recommend both, Arm A first.
5. **RMSD scorer** — spyrmsd (symmetry-correct) recommended over the engines' internal RMSD.
