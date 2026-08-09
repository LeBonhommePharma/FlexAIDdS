# Science audit — Astex-84 dG campaign vs Morency 3Dsig 2017

**Auditor:** Claude Code (Fable 5) · **Date:** 2026-08-09
**Baseline document:** `~/Downloads/Morency_LP_3Dsig_2017.pdf` (23 slides, ISMB/ECCB 2017)
**Campaign under audit:** `/Users/lp.more/flexaidds_results/astex84_dG_20260809_141245` (live, run by Claude Science)
**Method:** 6 independent dimension auditors -> 9 adversarial verifiers -> synthesis.
9 findings confirmed, 0 refuted, 27 lower-severity findings passed through unverified.

## Two corrections the verifiers made to first-pass readings (recorded so they are not repeated)

1. `ensemble_log_Z = 0.0000` and `report_T = 21.0000` on every row are **struct defaults that were
   never written** (the thermo engine is off in this build path). They are NOT measurements.
   Do not read `log Z = 0` as `Z = 1` / "the entropy term is inert" - that inference is wrong.
   `report_T = 21.0` is a documented reporting-only constant (`LIB/flexaid.h:692`, "kT_ISMB,
   ISMB 2017 calibration") that by design never feeds `G_bind` or CF scoring.
2. `cf_top1_* == entropy_top1_*` on 20/20 targets is NOT evidence that the entropy election is
   inert. Soft-beta fired on every target and moved the elected pose off the min-CF head on
   13-14 of 19. The columns are equal because the CF-arm branch is dead code
   (`LIB/DatasetRunner.cpp:7415`). It is a silent-zero artifact.

---

## Addendum — root cause of the `rmsd = -1.0000` rows (found after the audit workflow ran)

Section 3 / C4 below records the `-1.0000` sentinel as a false-positive trap but not its cause.
Traced afterwards, and independently confirmed by Claude Science:

`LIB/DatasetRunner.cpp:6222` writes the pose-blinded ligand to a **target-level shared path**
`run/<PDB>/<PDB>_dockin.sdf`, **inside** the restart loop (brace-matched body 5845–6446, with the
fork at 6422). Unlike the sibling artifacts `_pruned.pdb` (:6280) and `_prepped.pdb` (:6296),
which have exists+mtime write guards (:1674-1676, :6298-6302) and are written once on `ri=0`,
`write_blinded_ligand()` is called **unconditionally every iteration** — 10 rewrites per target.
Concurrently forked siblings therefore read a file a later iteration is rewriting:
`truncated at line 2` → `Failed to read ligand file` → `exit code 2`.

`ret` then aggregates any-failure-wins (:6440 sequential, :6467 parallel), so
`docking_completed = (ret == 0 && n_poses > 0 && !result.stuck)` (:6830) goes false and the whole
crystal-RMSD + elected-pose-persist block at :6866 is skipped — even though the surviving restarts
docked fine and the pooled election had already chosen and logged a winner.

**Why it stayed invisible:** both log lines that should surface it report success unconditionally.
`:6472-6474` prints `n_restarts` (the *configured* count, never `ri_ret`), and `:6444` pushes every
prefix into the pool outside any success check, so `:6496-6499` reports `10/10`. For 1HQ2 claim.log
says `10 restarts completed` then `pooling 10/10` while two restarts were dead. The only trace is in
per-restart `stderr.log`.

**Rate:** 2 of 22 completed targets (9.1%) — 1HQ2 (r6,r7) and 1N2V (r1 alone) — projecting ~7–8 of 84.
It is a race, so it is nondeterministic across campaigns.

**Recovery:** `scripts/recover_voided_targets.py` (branch `tools/recover-voided-targets`), read-only,
reconstructs RMSD from the `[3DSIG-RANK] rank=0` pose in claim.log. Recovers 1HQ2 → 8.4112 Å and
1N2V → 2.3241 Å. Note 1N2V is the third-best result in the run to date, so the race is voiding
competitive targets. Recovered top-1 is 0/22; a naive `rmsd < 2.0` filter reports 2/22 because
`-1.0 < 2.0`.

**Fix scope** (three items, not two): per-restart `_dockin` path; gate `docking_completed` on
poses-existing rather than all-restarts-succeeding; and emit a `restarts_ok` column — there is
currently no restart-count field in result.csv (`elected_restart` is the winner's index, and
`restarts_supporting` exists only as a comment at :1371), so dropping any-failure-wins without it
would trade a loud failure for silent per-target sampling-budget heterogeneity.

---

## Science Audit Verdict — `astex84_dG_20260809_141245` vs Morency 2017 3Dsig, slides 15–16

Run root: `/Users/lp.more/flexaidds_results/astex84_dG_20260809_141245`
Build: `repo_sha = 54882666 (main)`, `engine_md5 = 0aaee5b9267210ed2be42540672f490d` (`provenance.txt`)
State at audit: 22 target dirs created, 20 `result.csv` written, 84 targets contracted.

---

## 1. VERDICT

**No — not as designed; only-under-restatement, and even then only for the absolute FlexAIDdS number, never for the delta.** The campaign measures a top-1 point estimate at N=84 with a ~0.31× evaluation budget; slides 15–16 report the bootstrapped median of a top-10 success rate at N=85 with 10 × 2,000,000 evals per case. Worse for the headline claim of slide 23: the campaign runs **one arm only** — `cf_top1_*` and `entropy_top1_*` are byte-identical on 20/20 targets because the CF-control branch is dead code — so the +0.03 FlexAID→FlexAIDdS delta that the 2017 conclusion rests on is **structurally unmeasurable** from this run's CSVs. The good news: nothing is lost. Every restart's 50 cluster heads, their `.mcf` member sidecars, and the crystal ligands are on disk, so the 2017 estimand and both arms are recoverable post hoc without re-docking.

---

## 2. BLOCKERS
*(ordered by how badly each breaks the comparison to slides 15–16)*

### B1. The campaign has no control arm. `cf_top1_*` ≡ `entropy_top1_*` on 20/20 targets.

Slide 23's claim is a **paired delta under one protocol**. The CSV columns designed to carry it are hard-wired identical.

- `result.pose_source` has exactly one write site: `LIB/DatasetRunner.cpp:7122` → `result.pose_source = is_ini ? "ini_elitism" : "ga_cluster";`. Grep over `LIB/` returns no other assignment.
- The dual-estimand guard `LIB/DatasetRunner.cpp:7415` is `if (!softbeta_on || result.pose_source != "softbeta")` → **always true** → lines 7416–7423 copy `elected_*` into *both* `cf_top1_*` and `entropy_top1_*`. The CF-arm computation at 7425–7467 never executes.
- Data confirms: `cf_top1_pose_sha256 == entropy_top1_pose_sha256` on 20/20; `election_mode` empty on 20/20; `consensus_count = -1`; `rank0_demoted = 0`.
- This is **not** evidence that entropy is inert. Softβ fired on every target (39× `Softβ S1 ON`, 0× OFF; 117 `[3DSIG-RANK]` lines in `claim.log`) and moved the elected pose off the min-head-CF pose on **13–14 of 19** targets, giving away up to 55.2 CF units (1IGJ −99.18 elected vs −154.34 min-head; 1K3U 42.9 units; 1KZK 30.8; 1G9V 28.2).

**Consequence:** Δ = mean(`entropy_top1_rmsd`<2) − mean(`cf_top1_rmsd`<2) = **0.000 by construction**. Any analyst reading the CSV would conclude entropy does nothing. That is a silent-zero artifact, not a result.

### B2. The metric is top-1; 2017 scored top-10. Measured gap on the partial sample: 0.000 → 0.121.

- Stop condition (`handoff_swe/HANDOFF_CLAUDE_SCIENCE_benchmark_relaunch_20260809.md:16-18`): "the top-1 <2 Å rate computed by you from `rmsd_to_crystal`". Slide 14: success if RMSD < 2.0 Å **in the 10 best results predicted by the method**. The repo's own frozen protocol agrees — `docs/implementation/3dsig_red_pair_protocol.md:15`: "S_top10: any of top 10 ranked modes has RMSD ≤ 2.0 Å".
- Current top-1 on completed targets: **0/19 valid** (1HQ2 is the `-1.0000` sentinel). Minimum `rmsd_to_crystal` observed is 1M2Z at **2.1067**. Nothing is under 2 Å.
- 2017-faithful per-simulation S_top10 on the same data (each of 10 restarts scored on its own 10 lowest-CF emitted heads): **24/198 = 0.1212**, stable at 0.121–0.126 under every ranking key tried (CF, soft_beta_G, free_energy, emitted index, cluster frequency).
- Pooled-across-restarts top-10 = 3/20 = 0.150. Best-of-all-~500-poses (crystal-informed ceiling) = 7/20 = **0.350** (`best_cluster_rmsd` ≤ 2.0: 1G9V 1.5614, 1GPK 1.8793, 1HNN 1.0180, 1IA1 0.8877, 1JD0 1.1312, 1L2S 1.9919, 1MQ6 1.0061).
- `result.csv` cannot express S_top10: header at `LIB/DatasetRunner.cpp:7752-7776` is 83 columns (verified against `run/1G9V/result.csv`: 83 header, 83 data), and the only crystal-RMSD fields are `rmsd_to_crystal`, `rmsd_hungarian`, `cf_top1_rmsd`, `entropy_top1_rmsd`, `best_cluster_rmsd`. `rmsd_to_crystal` is recomputed from `result.elected_pose_path` alone (`DatasetRunner.cpp:7386-7394`); `pose_ledger/*.json` holds a single `role="elected"` record.
- Fail-closed confirmed empirically: `python3 scripts/bootstrap_3dsig_s_top10.py --arm-dir <OUT>/run` → `error: S_top10 fail-closed: 20 result.csv lack mode_rmsd_*; Refusing BCR/top1 fallback.`

**Consequence:** the emitted number is strictly ≤ the 2017-comparable number by a target-dependent, unbounded amount. On this partial sample the gap is 12.1 points, not a rounding error.

### B3. Effective budget is ~0.32× of the 2017 per-simulation budget.

- Configured: 1000 pop × 1000 gen = **1,000,000 evals/restart** = 0.50× of slide 14's 2,000,000. Repo protocol requires 2M: `METHODOLOGY.md:29-30` ("2000 generations, population 1000 … Do not change for accuracy runs"), `docs/implementation/COMPARATIVE_BENCHMARK_METHODOLOGY.md:35`, `3dsig_red_pair_protocol.md:17`, `arm_pins.json:19 "evals_per_sim": 2000000`.
- Eval accounting is exact, not inferred: `rep_model="BOOM"`, `boom_fraction=1.0` (`LIB/config_parser.cpp:222,225`; neither key present in `run/1G9V/dock_config.json`) → `nnew = 1000` offspring/gen (`LIB/gaboom.cpp:2209`), offspring get `status=' '` (2330, 2367), eval loop skips `'n'` (3252) → exactly 1000 CF evals/gen. Initial pop 1000 (`gaboom.cpp:522`). `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=fixed` → `eval_scale_dihedral=-1` (`LIB/ProtocolConfig.cpp:257-260`) → fixed pop/gen (`DatasetRunner.cpp:5798-5802`); `FLEXAIDDS_BUDGET_SCALE=0` → no multiplier (5810). `claim.log` line 11: `GA: pop=1000 gen=1000 (1000k evals/complex)`.
- **Early exit was left on** in violation of `METHODOLOGY.md:32` / `COMPARATIVE_BENCHMARK_METHODOLOGY.md:41` (both mandate `FLEXAIDDS_NO_SEC=1`): `run/provenance.json` has `no_sec=false`, `benchmark_mode=false`. Termination census over the 198 restarts that had printed `Done.`: **103 entropy-convergence, 79 fitness-stagnation, 16 reached gen 1000** — the cap binds on ~8%.
- Actual spend from the last `Generation:` line of 195 GA logs: min 258, p25 558, **median 618**, p75 868, max 999 → median ≈ 618,000 + 1,000 init + ≤6,400 coarse-init ≈ **625k evals ≈ 0.31×**.

**Do not** call this a "strict lower bound." The sign of the budget effect on top-1 success is unmeasured here, and it interacts with the election rule. Fixing it requires **both** `--ga-generations 2000` and `FLEXAIDDS_NO_SEC=1`; the generation count alone changes nothing for 92% of restarts.

### B4. Two mutually inconsistent temperatures in one row; the reported dG̃ is not the objective that elected the pose. *(not adversarially verified — lower severity tier)*

The election forms G̃ at `T_soft = 300` in **dimensionless CF units** (`FLEXAIDDS_ELECTION_SOFT_T` unset → `election_soft_T=0.0` in `RUN_RECEIPT.json` → dock branch at `DatasetRunner.cpp:1016-1021`; `claim.log:76` `T=300.0000 source=dock (soft-β CF a.u., not k_B·T)`), while `predicted_dG` / `predicted_dH` / `predicted_TdS` — the columns the campaign is *named* after — use `kB·300 = 0.596 kcal/mol`, a 503× different temperature. In the election, entropy dominates enthalpy ~10:1; in the reported ledger it is ~1%. Also flagged in this tier: `ensemble_log_Z = 0.0000` and `report_T = 21.0000` on all 20 rows are **struct defaults, never written** (thermo engine off) — do not read log Z = 0 as Z = 1.

---

## 3. CONFOUNDS
*(contaminate the number/delta; do not void the run)*

**C1 — The election is a pose-count vote, not H − TS with a working energy filter.** At T = 300 CF a.u. the formula is implemented correctly (max |G − (H − T·S)| = 0.0149 over 117 `[3DSIG-RANK]` rows) but sits deep in the entropy-dominated regime: S̃ is within 2.3e-3 nats (median) of `ln(n_members)` (1G9V rank0 S=5.8827 vs ln(359)=5.8833; 1GM8 rank0 S=5.7289 vs ln(308)=5.7301). Within-target sd: `T·ln n` = 343.4 CF units vs Emin 36.1 and spread 14.1 — a ~10:1 variance ratio. `T·ln2 = 207.9` CF units per doubling of member count, against a **full** within-target head-CF range of only 102–250 units (1J3J 102.2, 1G9V 138.6, 1KE5 190.0, 1MQ6 249.8). The winner is `argmax(n_members)` in 16/20 (the 4 exceptions are ≤2% count ties broken on enthalpy) and differs from the min-CF head in 18/20. Concretely, 1G9V's elected mode has H̃ = −67.36 vs the loser's −116.49: 49.1 units of enthalpy surrendered because 359 > 145 members (300·ln(359/145) = 271.9). Slide 11 says the winning wide well contains multiple poses **with favorable energy**; at this T the favorability filter is gone.

**C2 — GA clones cast the votes.** The elected `.mcf` member counts are heavily duplicated: `1G9V r8/1G9V_1.mcf` n=359 / 46 unique, one CF value ×203 (that clone sits 32 CF units *worse* than its own head); `1KE5 r7/1KE5_0.mcf` n=1572 / 112 unique, one CF ×956 (92.9%); clone fraction 87.4% (1MQ6), 87.2% (1G9V), 82.6% (1K3U), 72.6% (1JD0). Note the member count 1572 **exceeds the population of 1000** — Z is summed over a GA snapshot with generation/restart structure, not an independent conformational census. `DatasetRunner.cpp:1122` calls multiplicity-sensitive `soft_beta::free_energy()`; `cluster.cpp:206-208` and `BindingMode.cpp:390-391` call `free_energy_strict(UniqueGeometry)`. Strict flips the elected head on 8–9/19 targets. **Fidelity note, important:** slide 12 sums over *all* N_poses and slide 11 names multiplicity as the signal, so the non-strict path is arguably the slide-12 formulation and `UniqueGeometry` is a 2026 repo addition with no 2017 counterpart. What is unambiguously wrong is the **provenance inconsistency**: the elected pose's own REMARK reads `Cluster 1: … Average CF:-1206.50789` (strict, written by `cluster.cpp`) while DatasetRunner elected it on G = −1832.18 (classic), and nothing records which variant ranked the modes. Success is 0/9 under both on the flipped targets, so today's headline number is unchanged either way.

**C3 — 2HR7 silently dropped; the denominator lives in `/tmp`.** Live process: `--only-codes /tmp/codes84.list` (84 lines, verified). `astex_diverse_codes()` at `LIB/DatasetRunner.cpp:4042-4055` returns 85. Set diff: *in canonical, not in file*: `['2HR7']`; *in file, not in canonical*: `[]`. The list is byte-identical (sorted) to `astex84_v2_20260808_223645/codes84.list` — inherited, never justified. 2HR7's cache is complete and prepares cleanly (`cache_v2/astex_diverse/2HR7/{2HR7.cif,.pdb,_apo.pdb,_ligand.sdf,_ligand_centered_site.pdb}`; `astex84_v2 prep.log:137-139`). The bias is **directionally favourable but small**: 2HR7 is a top-1 failure in the last full-85 run (`rmsd_to_crystal=11.5524`), so 30/85 = 0.3529 → 30/84 = 0.3571, i.e. **+0.0042**, ~1/13 of one bootstrap SE. The real defect is attributability: `RUN_RECEIPT.json` records no case list and `OUT/` holds no copy. (Separately: the repo's *other* 84-list, `benchmarks/astex_repro/astex84_no1hq2.txt`, excludes 1HQ2 and **keeps** 2HR7 — two different 84-sets are floating around.)

**C4 — `-1.0000` sentinel is a false-positive trap.** 1HQ2 has `rmsd_to_crystal = -1.0000`, `rmsd_hungarian = -1.0000` (no elected pose) yet carries `best_score = -255.7534` and `predicted_dG = -213.7031`. A literal reading of the stop condition ("<2 Å from `rmsd_to_crystal`") scores it as a **success**. N silently drifts below 84.

**C5 — Scoring-function deltas vs slide 10.** *(this cluster was not adversarially verified — treat as leads, each individually checkable)*
- Steric wall: `soft_wall_cutoff=0.4` routes CF.wal through `soft_wall_fitness_energy()`, a Hermite-cubic/quadratic in overlap `(cr−d)` hard-clamped at `WAL_CONTACT_CAP=50`/contact — a different functional form from slide 10's unbounded `K_wall·[(1/d)^12 − (1/(Pe(ri+rj)))^12]`.
- `sas_weight = 0.40` (engine default 1.0) attenuates the slide-10 ligand-solvent term 2.5×, and the source comment states it was chosen because 1.0 collapsed Astex genuine top-1 from 52.2% → 6% — i.e. **fitted on this test set**.
- `vct_dist_weight_r0 = 7` multiplies every complementarity contact by `exp(−d/7)`; slide 10 modulates ε by contact **surface area** only. The value was selected by Astex success rate.
- Matrix `OUT/bin/MC_st0r5.2_6.dat` md5 `9dc93717dfed0698006d88dd6a9627bc` is hand-patched in 7 cells off the trained `204b75ef…`; three edits (13-40, 14-40, 15-40) are the SOLVENT column, i.e. slide 10's implicit-water ε(i′,w).
- CF total = com + wal + sas + elec + **hbond** + gist_desolv + **metal_coord** + entropy + pb_clash — six channels beyond slide 10's three; the H-bond Gaussian is active on 87.6% of emitted poses and flips top-1 on 2/21 targets.
- Water-typed-as-carbon is **checked and cleared** for this build (fix `739bc3d3` is an ancestor of `54882666`) — but 84/85 cached apo receptors still carry 271–508 explicit HOH records, scored as O.3 contacts *on top of* the implicit-water channel.

**C6 — Search-protocol additions with no 2017 counterpart.** *(not adversarially verified)* Oracle SITE-CONFINE restricts translation to (ligand extent + 0 Å margin) around a centroid a median 1.88 Å from the true ligand centroid — far tighter than 2017's GetCleft cleft. `coarse_init` spends up to 6,400 **uncounted** CF evals/restart and injects the 25 best into gen 0; MIF seeding adds an interaction map over the confined grid. "CF" clustering is greedy fixed-radius leader clustering at a hard 2.0 Å ordered RMSD — a tiling, not slide 13's density-based clustering — and `cluster_rmsd = 2.0` makes a "binding mode" exactly one success ball (23% of head pairs sit within two success radii). `max_results = 50` truncates clustering, so unassigned chromosomes never enter any Z or S̃.

**C7 — Blinding holds; PoseBusters is not diagnostic here.** *(not adversarially verified)* All completed rows: `native_pose_seeded=0`, `seed_fraction=0.0000`, `seed_echo=0`, and the docked input is genuinely orientation-blinded — **no seed-echo inflation**, unlike earlier campaigns. Caveat: `native_pose_seeded` is a mode constant, not a measurement. PoseBusters pass is 3/20, but `native_qc_pass` is 9/20 — the **crystal** pose fails the same suite on 11/20 (modal failure `minimum_distance_to_waters`, 13/20), so PB here reads receptor prep, not docking quality. Also: the `success` column is now an alias for `success_rmsd` (RMSD ≤ 2.0 **AND** seed-echo **AND** pose-hash gates), not "docking ran" as the handoff claims — it is stricter than a hand-computed RMSD count, and neither is the 2017 statistic.

---

## 4. WHAT THE CAMPAIGN CAN LEGITIMATELY CLAIM

As designed, exactly this and nothing more:

> "On 84 of the 85 Astex Diverse targets (2HR7 excluded; case list not preserved in the campaign record), FlexAIDdS at `repo_sha 54882666`, run with an oracle-confined binding site, 10 restarts of ≤1,000,000 CF evaluations each (median actual spend ≈625,000 owing to entropy-convergence and stagnation early exits), and a soft-β free-energy election at T = 300 CF arbitrary units, placed its **single elected pose** within 2.0 Å of the crystal ligand in X of 84 cases."

Mandatory riders in the same paragraph:
1. **"This is a top-1 election-quality metric, not the 2017 top-10 success rate, and is not comparable to the slide 15/16 boxplots."**
2. **"No FlexAID control arm was run; no FlexAID→FlexAIDdS delta is reported."**
3. Budget is ~0.31× of the 2017 per-simulation budget.
4. Point estimate only; no bootstrap interval.
5. If the number is used at all, state whether `-1.0000` rows were excluded and the resulting N.

What LP **cannot** say from this run: any sentence containing "0.66", "0.69", "+0.03", "improves on FlexAID", "entropy helps", "reproduces 3Dsig", or any placement of this number on the 2017 boxplot axis.

For calibration on how weak a single point estimate is here: binomial SE at N=84 is 0.0517 at p=0.66, 0.0505 at p=0.69, 0.0500 at p=0.30, 0.0327 at p=0.10. The 2017 delta of +0.03 is **0.58 SE** paired-optimistic and **0.42 SE** unpaired (SE of an unpaired 0.66-vs-0.69 difference = 0.0722). One 84-case arm has roughly one-third the resolution needed. The 2017 FlexAID and FlexAIDdS boxes already visibly overlap each other in the deck.

---

## 5. WHAT IT WOULD TAKE TO ACTUALLY REPRODUCE SLIDE 15/16

Six deltas. All six are required; any five gives a number that still cannot be placed on the boxplot.

1. **Top-10 criterion, with the pool defined in writing before scoring.** Score success as: RMSD ≤ 2.0 Å for **any** of a simulation's 10 best-ranked emitted results. The current campaign retains 10 restarts × 50 heads ≈ 500 candidates per target, so "top 10" is ambiguous and the two readings give different numbers (measured: per-simulation top-10 = 0.121; pooled-across-restarts top-10 = 0.150). Declare one — (a) top-10 heads of a single simulation, matching 2017's "the 10 best results predicted by the method" for one of the 10 runs, or (b) global top-10 by ranking key across restarts — and never mix. Ranking must be **crystal-blind**; `best_cluster_rmsd` / `conditional_scanned_pool_ceiling` is min-RMSD-over-all-~500-poses chosen **using the answer** (0.350 on the partial set) and is an oracle upper bound, not a top-10 proxy.
2. **10 × 2,000,000 evaluations per case, actually spent.** `--ga-population 1000 --ga-generations 2000` **and** `FLEXAIDDS_NO_SEC=1` / `benchmark_mode=true`. Generations alone do not fix it: 187/198 restarts currently exit early. Decide explicitly whether `coarse_init`'s ≤6,400 evals and MIF seeding count against the budget or are declared as protocol additions absent from 2017.
3. **10,000-iteration bootstrap, reported as median + IQR + whiskers.** Per `3dsig_red_pair_protocol.md:18`. `scripts/bootstrap_3dsig_s_top10.py --bootstraps 10000` already implements it. Report a boxplot, not a scalar — the deck's own presentation is a boxplot precisely because the point estimate is not informative at N≈85.
4. **N = 85, including 2HR7.** Its cache is complete and prepares cleanly, so there is no technical reason to exclude it. If it must be excluded, publish a pre-registered crystal-blind criterion in `provenance.txt` and report both N=84 and N=85. Move the code list out of `/tmp` into the campaign directory — the denominator of a scientific claim must not live in a file that reboots away.
5. **THE PAIRED CONTROL ARM — the single most important item.** Slide 23 asserts a delta, not a level. That requires two arms scored under one identical protocol: **arm A = FlexAID**, election = min CF among the same emitted heads; **arm B = FlexAIDdS**, election = soft-β G̃. Same seed, same site, same budget, same clustering, same emission ordering, same pose pool. Then bootstrap the **paired** difference: resample case indices once, score *both* arms on that same resample, take the difference per resample, report the median and interval of the difference distribution. Unpaired differencing of two independent boxplot medians has SE 0.0722 at N=84 and cannot resolve +0.03; only the paired design removes the between-case variance. Note the emission-order contamination that must be disclosed: with `classic_entropy_ranking=true` → `force_cf_rank_emission=false`, the 50 retained heads are truncated in **entropy** rank order, so an offline CF arm approximates rather than exactly reproduces a `force_cf_rank_emission=true` FlexAID baseline. A clean paired run needs both emission orders, or a retention limit high enough that truncation is not binding.
6. **One temperature, declared with units.** Either the election and the reported ledger share a T, or the report states plainly that the pose was elected at T = 300 CF a.u. and the ΔG̃/ΔH̃/TΔS̃ ledger is a separate kB·T = 0.596 kcal/mol quantity that did not participate in selection. Additionally: fix and declare `free_energy` vs `free_energy_strict`, since the elected pose's own REMARK currently reports the *other* variant's G than the one that elected it.

Secondary, needed for a defensible claim of "same scoring function": the slide-10 deviations in C5 (capped soft-core wall, `sas_weight=0.40` tuned on Astex, `exp(−d/7)` distance decay, 7 hand-patched matrix cells including three solvent-column entries, four extra CF channels). These do not block a 2017-shaped *protocol* comparison but they mean the engine being benchmarked is not the 2017 CF, which must be stated.

---

## 6. CHEAP MEASUREMENTS AVAILABLE RIGHT NOW

All read-only, no re-docking, on data already on disk. Ordered by value per unit effort.

**M1 — Recover S_top10 and its bootstrap from the existing poses (recovers B2 entirely).**
Poses: `<OUT>/run/<code>/<code>_{0..49}.pdb` and `<OUT>/run/<code>/r{1..9}/<code>_{0..49}.pdb` — 504–510 pose PDBs per target on 1G9V/1GM8/1KE5. Each carries `REMARK CF=`, `REMARK Cluster <i>: Rank (top):<r> Average CF:<x> Frequency:<f>`, `REMARK binding_mode = <i>`, `REMARK pose_rank = <n>`. Crystal reference: `/Users/lp.more/flexaidds_results/cache_v2/astex_diverse/<code>/<code>_ligand.sdf`. `scripts/parse_flexaid_arm_results.py` already emits `mode_rmsd_0..9` in crystal-blind emitted-rank order; pipe into `scripts/bootstrap_3dsig_s_top10.py --bootstraps 10000`. Declare the pool convention first (item 5.1 above). Expected: ~0.12 per-simulation, ~0.15 pooled, vs 0.00 top-1.

**M2 — Reconstruct the CF control arm offline (recovers B1 without re-running anything).**
Per target take `argmin` over the head CF (first line of every `<OUT>/run/<code>/**/*.mcf`, verified to be the head CF), RMSD that pose against `<code>_ligand.sdf` with the same ordered element-matched routine, report top-1 and S_top10 for both arms as a paired table plus the paired bootstrap of the difference. Already partially done in this audit: elected_cf differs from global min head CF on 13/18 targets with the deltas listed in B1. Disclose the entropy-truncation caveat from item 5.5.

**M3 — Classic vs strict re-election A/B on the same ensemble.**
`scripts/acf_strict_offline_reelect.py` implements `free_energy_strict(UniqueGeometry)` over the written `.mcf` sidecars. Report both rates and, per target, `n_members` and `n_unique` of the elected cluster so clone inflation (1G9V 359/46, 1KE5 1572/112) is auditable in the record.

**M4 — T_soft sweep over the `.mcf` sidecars.**
Recompute the election at T ∈ {0.5, 1, 2, 5, 10, 30, 100, 300} CF a.u. and plot top-1 and S_top10 vs T. Choose T from observed intra-mode CF spread (the elected 1G9V mode spans 64.1 units) so weights actually discriminate. Report the S̃ − ln(n) residual per target as the diagnostic that the entropy term carries information beyond a headcount; today the median residual is 2.3e-3 nats.

**M5 — Bootstrap the top-1 vector you already have.**
Zero-cost sanity floor: run the 10,000-resample bootstrap on the final `rmsd_to_crystal` vector, excluding `-1.0000` sentinels, and publish median + p05/p95 + the explicit N. This does not make the number comparable to 2017, but it stops a bare scalar from being read as precision it does not have.

**M6 — Two one-line integrity checks before any analysis.**
(a) Count `-1.0000` rows and state N explicitly (currently 1/20: 1HQ2). (b) Copy `/tmp/codes84.list` into `<OUT>/` next to `provenance.txt` — the file that defines the denominator of the claim currently exists only in a temp directory.

**M7 — Budget-slope bound on a subset (only non-free item, ~20 targets).**
Paired A/B, same seed: arm 1 = current settings, arm 2 = `--ga-generations 2000` + `FLEXAIDDS_NO_SEC=1`. This is the only way to bound how much of the deficit in B3 is real. Do not do this until the current campaign finishes — the audit constraint is read-only.