# PRE-REGISTRATION — receptor-strain / wall-cap / solvation-reference campaign

**Status: PRE-REGISTERED, NOT YET RUN. Frozen 2026-09-03.**

Nothing in this document was measured by the campaign it describes. Every number
in §2 is quoted from the ARCHIVE. Every number in §5 is a THRESHOLD chosen before
any treated cell exists. If a threshold here is edited after a treated cell has
been scored, the campaign is no longer a test and must be reported as exploratory.

Companion driver: `scripts/run_prereg_receptor_strain.sh`.
Engine/harness gates under test: stages 1–3 of `/tmp/flexaidds_patches/` (see §4).

---

## 0. One-paragraph summary

FlexAID scores with IMPLICIT solvation: the absence of an atom means that space
is water. That is coherent for a RIGID receptor. For a FLEXIBLE receptor it is
not — when the GA rotates a side chain away, the vacated volume reads as bulk
solvent, and bulk solvent is cheap to displace. The search can therefore
manufacture a pocket by evicting a side chain and is then paid for filling it.
Two further terms fail to object: `Pose::receptor_strain` is declared, wired into
the ensemble, and identically zero because it is never computed; and the
per-contact steric ceiling `WAL_CONTACT_CAP = 50.0` bounds what burial into the
vacated hole can cost. This campaign tests whether pricing those three things
removes the observed cost of flexibility.

---

## 1. HYPOTHESIS (falsifiable)

> **H1.** The accuracy penalty that receptor side-chain flexibility imposes on
> FlexAIDdS is caused by **eviction-then-burial that the objective does not
> price**: volume vacated by a moved side chain is scored as bulk solvent
> (`eps(i,w)*S(i,w)` against a receptor envelope the search itself invented), the
> conformational cost of the eviction is charged at exactly 0.0
> (`receptor_strain` never computed), and the steric objection to burying into the
> vacated volume is capped at 50 per contact.

> **H1 predicts a coupled signature**, not two independent ones: the
> interpenetration failures and the top-1 loss are the SAME event seen twice —
> once by the validator (in the crystal frame the ligand overlaps protein that the
> search moved) and once by the ranking (a physically wrong pose outscores the
> near-native one). Pricing the three unpriced quantities must therefore move
> **both**.

**H0 (what we are trying to reject).** The flexibility penalty is not caused by
unpriced eviction; it is caused by something else (dimensionality of the enlarged
chromosome, GA convergence, election, site definition), and the three gates will
move at most one of the two signatures, or neither.

---

## 2. CONTROL — the archive. Quoted, not re-derived.

Archived arms (`NO_SEC` campaign): seed base **12345**, **3000** generations,
`FLEXAIDDS_SCORED_ONLY=1`, `--clustering CF`, `--mode defined-cleft-redock`,
oracle site directory, CF election.

| arm | receptor | top-1 (S1) | notes |
|---|---|---|---|
| `A_rigid` | rigid | **48 / 84** | no autoflex |
| `B_shrink` | flexible | **32 / 84** | `FLEXAIDDS_AUTOFLEX_MAX=5`, `FLEXAIDDS_AUTOFLEX_METAL_SHRINK=1` |
| net | | **−16** | McNemar exact **p = 1.45e-04** |

PoseBusters, **77 evaluated cells per arm**, crystal-frame receptor:

| check | rigid | flexible |
|---|---|---|
| `volume_overlap_with_protein` | **0 / 77** fail | **38 / 77** fail |
| `minimum_distance_to_protein` | **19** fail | **55** fail |
| any internal-geometry check | **0** fail | **0** fail |

**The internal-geometry row is the load-bearing one.** The flexible arm's poses
are not distorted molecules. They are well-formed ligands sitting inside the
protein. That is what makes "the receptor moved and nothing charged for it" the
leading explanation rather than "the GA broke the ligand".

### 2.1 Reconstructed discordance — MUST BE CONFIRMED BEFORE THE CAMPAIGN STARTS

`p = 1.45e-04` two-sided exact-binomial with a net of −16 is uniquely consistent
with **b = 17 recoveries, c = 1 loss, D = b + c = 18**
(`2·P(X ≤ 1 | n = 18, ½) = 1.4496e-04`). This is a **reconstruction from the
reported p**, not a reading of the archive.

> **PRECONDITION.** Before any treated cell runs, `(b, c)` must be read directly
> from the archive's per-target table. If it is not `(17, 1)`, §6 (power) is
> wrong and this document must be corrected and re-frozen first.

### 2.2 What the control is NOT

`A_rigid` is not "the right answer". It is the reference the flexible arm lost 16
targets against. A treated flexible arm that merely matches `A_rigid` has removed
a defect; it has not demonstrated that flexibility helps. **This campaign does not
test whether receptor flexibility improves docking.** It tests whether the
flexibility penalty is the unpriced-eviction defect. Any claim of the former from
these data is out of scope.

---

## 3. POPULATION, ARMS, AND EVERYTHING HELD FIXED

* **Targets.** The **84-target list and the 77-cell PoseBusters subset are taken
  VERBATIM from the archive**, as files. They are not re-derived from
  `astex85_target_manifest.json` and not recomputed from a success mask.
  (Astex-85 in this repository is not the Hartshorn set; re-deriving the list
  would silently change the population.)
* **Held fixed in every arm** (any deviation voids the comparison):
  seed base 12345 · restarts 3 · 3000 generations · population 1000 ·
  `--mode defined-cleft-redock` · `--clustering CF` · the same oracle site
  directory · the same cache root · `FLEXAIDDS_SCORED_ONLY=1` · the same
  `--omp-threads` value · the same engine SHA-256 in every cell.
* **`FLEXAIDDS_SCORED_ONLY=1` is kept** precisely because the archive used it.
  It changes which atoms appear in the pose PDB and therefore what the validator
  and the RMSD both consume. Changing it here would confound the treatment.
* **Engine provenance precondition.** The build must contain the cleft-grid
  determinism fix (`8dc88b4e`). Before it, the SURFNET probe-merge order permuted
  the cleft-grid indices that become GA gene 0, so multi-threaded arms differ for
  a reason that has nothing to do with any gate here. The driver asserts this
  with `git merge-base --is-ancestor`.

### 3.1 Arms

| arm | autoflex | stage-1 | stage-2 | stage-3 | role |
|---|---|---|---|---|---|
| `A_rigid`  | off | off | legacy | dynamic | archive replication, rigid |
| `B_shrink` | 5 + metal-shrink | off | legacy | dynamic | archive replication, flexible = **the comparator** |
| `T1_strain`  | 5 + metal-shrink | **on** | legacy | dynamic | exploratory decomposition |
| `T2_walcap`  | 5 + metal-shrink | off | **flex** | dynamic | exploratory decomposition |
| `T3_solvref` | 5 + metal-shrink | off | legacy | **crystal** | exploratory decomposition |
| `T4_all`     | 5 + metal-shrink | **on** | **flex** | **crystal** | **CONFIRMATORY** |

**Only `T4_all` is confirmatory.** `T1`–`T3` are a mechanistic decomposition and
are reported with **unadjusted p values explicitly labelled EXPLORATORY**. This
is the multiplicity control: one confirmatory contrast (`T4_all` vs `B_shrink`),
declared here, before any cell runs.

### 3.2 Instruments, on in EVERY arm including the rigid one

* `FLEXAIDDS_CONTACT_PROFILE=1` — pure diagnostic sidecar, adds nothing to CF.
* `FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1` — pure output; writes the receptor as
  scored. **Subject to the disk pilot in §9.**

An instrument that is only on in the treated arm is not an instrument, it is a
confound. Both are on everywhere or off everywhere.

### 3.3 Explicitly NOT in this campaign

`FLEXAIDDS_PB_RECEPTOR=flexed` is **not an arm here**. It changes the
**measurement frame**, not the model. A flexed-frame overlap check cannot tell
whether the receptor STATE is physical: if the search evicted a side chain, the
overlap disappears *because the protein moved*, and the number improves for the
wrong reason. Crystal-frame and flexed-frame are **different estimands**. They may
be reported side by side; they may **never** be pooled, differenced, or used
interchangeably in one column. Every number in §5 is **crystal-frame**.

---

## 4. TREATMENTS — the exact gates

| stage | variable | default | treated value |
|---|---|---|---|
| 1 | `FLEXAIDDS_RECEPTOR_STRAIN` | off | `1` |
| 1 | `FLEXAIDDS_RECEPTOR_STRAIN_T` | `298.15` | `298.15` (unchanged) |
| 2 | `FLEXAIDDS_WAL_CAP_MODE` | `legacy` | `flex` |
| 2 | `FLEXAIDDS_WAL_CAP_FLEX` | unset (`0` = no ceiling on the soft-core wall) | unset |
| 3 | `FLEXAIDDS_SOLVATION_REF` | `dynamic` | `crystal` |
| 4 | `FLEXAIDDS_CONTACT_PROFILE` | off | `1` in **all** arms |
| 5 | `FLEXAIDDS_WRITE_FLEXED_RECEPTOR` | off | `1` in **all** arms |
| 5 | `FLEXAIDDS_PB_RECEPTOR` | `crystal` | `crystal` in **all** arms |

Each is echoed into every per-case `dock_config.json`:
`scoring.receptor_strain`, `scoring.receptor_strain_temperature_K`,
`scoring.wal_cap_mode`, `scoring.wal_cap_flex`, `scoring.solvation_ref`,
`scoring.autoflex_max`, `output.contact_profile`, `output.write_flexed_receptor`,
`output.pb_receptor`.

### 4.1 The precondition that voids everything if unmet

Stages 1–3 are **structurally inert without `FLEXAIDDS_AUTOFLEX_MAX > 0`**: with
a rigid receptor no residue is ever off rotamer 0, no contact is attributable to a
flexed side chain, and no volume is vacated. At least eight past features in this
repository were structurally absent because the harness never supplied the
precondition **while every log looked clean**. The driver therefore refuses to
score any arm until it has observed, on disk:

* `REMARK CF.strain` present and **non-zero on at least one pose** in `T1`/`T4`;
* at least one `[WAL_CAP]` attribution event in `T2`/`T4`;
* at least one `[SOLV_REF]` occlusion event in `T3`/`T4`;
* `REMARK n_residues_off_input_rotamer` **> 0** on at least one companion in
  every flexible arm.

A gate that is set but never fires is a **null measurement reported as a result**.
That is the failure this section exists to make impossible.

---

## 5. PRIMARY PREDICTION — declared before any cell runs

**H1 predicts BOTH of the following. They are conjunctive. Either one alone is
NOT a confirmation.** All numbers are crystal-frame, endpoint E1 (§7),
seed operator §8, tier V0/V1 as noted.

### P1 — interpenetration (validity signature)

> `T4_all`, `volume_overlap_with_protein` failures ≤ **8 of 77** evaluated cells
> (from 38/77; a ≥ 79 % reduction),
> **and** the paired McNemar `T4_all` vs `B_shrink` on that indicator shows
> **b ≥ 30** fail→pass with **c ≤ 2** pass→fail (two-sided exact p ≤ 2.5e-07 at
> b=30, c=2).

### P2 — net top-1 (accuracy signature)

> `T4_all` S1 ≥ **45 / 84** — i.e. net vs `A_rigid` ≥ **−3**, recovering ≥ 13 of
> the 16 targets flexibility cost —
> **and** the paired McNemar `T4_all` vs `B_shrink` on S1 shows **b ≥ 10**
> recoveries with **c ≤ 2** new losses (two-sided exact p ≤ 0.0386 at b=10, c=2;
> see the §6 table for the general (b, c) requirement).

**CONFIRMATION = P1 ∧ P2.** Anything else is §6 or §7 below.

### Ambiguous zone, declared in advance so it cannot be re-drawn later

If `volume_overlap_with_protein` lands in **(8, 30]** of 77, or S1 lands in
**(38, 45)** of 84, the result is **INCONCLUSIVE**. It is neither a confirmation
nor a refutation, no threshold may be moved to make it one, and the correct
report is "the effect is real but partial; the mechanism is not fully priced."

---

## 6. FALSIFICATION — explicit

| # | outcome | verdict |
|---|---|---|
| **F1** | **P1 holds, P2 fails** (interpenetration ≤ 8/77 but S1 ≤ 38/84, i.e. net ≤ −10) | **H1 REFUTED as a single cause.** The unpriced-eviction defect is real and now fixed, but it is **not** what cost the 16 targets. These are **two separate problems**. Report both, do not merge the gates into the default, and open a separate investigation into the accuracy loss (dimensionality of the enlarged chromosome and election are the next suspects). Publishing the validity improvement as an accuracy result is forbidden. |
| **F2** | **P2 holds, P1 fails** (S1 ≥ 45/84 but ≥ 30/77 still interpenetrate) | **H1 REFUTED.** Accuracy improved by some route other than the named mechanism — most plausibly the gates shrank the reachable search space. The poses still interpenetrate, so the model still does not know the receptor moved. This result may **not** be reported as a solvation/strain finding. |
| **F3** | neither P1 nor P2 | **H1 REFUTED outright.** |
| **F4** | the OFF-arm bit-identity canary fails (§10, G1) | **CAMPAIGN VOID.** The patches are not default-preserving; no number from any arm may be reported. |
| **F5** | any liveness check in §4.1 fails | **CAMPAIGN VOID for that arm.** The treatment never reached the model. Escalate; do not patch around it and do not report the arm as a null. |
| **F6** | `A_rigid` or `B_shrink` replication differs from the archive by **> 3 targets** | **STOP before scoring any treated arm.** The environment is not the archive's; no comparison against 48/84 or 32/84 is licensed. |
| **F7** | any arm has < 84 completed cells, or any cell has `rmsd_to_crystal < 0` | that arm is **INCOMPLETE**, not partially scored. `-1.0` is the shared-`_dockin.sdf` race sentinel and passes a naive `rmsd < 2` filter; every filter in the analysis must guard `rmsd >= 0`. |

**F1 is the outcome most likely to be mis-reported**, because both halves look
like good news. It is written first for that reason.

---

## 7. ENDPOINT — declared in advance

The engine emits a pose elected by **soft-beta free energy**. The archived
campaign reported the **minimum-CF** pose. **These are different poses.**
Reporting one against a control measured on the other is a category error.

| id | definition | role |
|---|---|---|
| **E1** | **argmin apparent CF** over the emitted pose set of the cell, read from `REMARK CF.app` in the pose PDBs (excluding `*_INI.pdb`) | **PRIMARY** — the only endpoint comparable to the 48/84 and 32/84 control |
| **E2** | the pose the engine actually emits as elected (`elected_pose_path` / `result.csv → rmsd_to_crystal`) | **CO-PRIMARY, reported separately, never substituted for E1** |

Both are computed for every arm and both appear in the results table. **Neither may
be swapped for the other after the data are seen.** A disagreement between E1 and
E2 is itself a reportable finding about election, and must be reported rather than
resolved by choosing the better number.

Success on either endpoint: `rmsd_to_crystal >= 0.0 AND <= 2.0` Å, symmetry-aware
(Hungarian) where the harness provides it, using the **same** RMSD column for
every arm.

---

## 8. SEED OPERATOR — frozen in advance, applied to EVERY endpoint

The previous campaign used a **majority** operator for top-1 and a **union**
operator for the ceiling in the same table, and had to retract. That cannot recur.

1. **PRIMARY:** a single pre-registered seed base, **`FLEXAIDDS_SEED_BASE=12345`,
   `FLEXAIDDS_RESTARTS=3`**. The 3 restarts are internal to the engine's own
   election and are **not** a seed operator.
2. **The primary operator is "single frozen seed base", and it is applied
   identically to S1, to E1, to E2, to the pool ceiling, to every PoseBusters
   indicator, and to every validity tier.** No exceptions.
3. **SENSITIVITY (optional, non-confirmatory):** seed bases 22345 and 32345.
   If run, they are combined by **MAJORITY-OF-3, applied identically to every
   endpoint in the table**. They may never upgrade an INCONCLUSIVE to a
   confirmation.
4. **UNION IS BANNED** at every endpoint, including the pool ceiling. A ceiling
   computed as a union over seeds is not comparable to a top-1 computed as a
   majority, and a table mixing the two is void.
5. Any table in which two endpoints used different seed operators is **void** and
   must be recomputed, not annotated.

---

## 9. VALIDITY TIERS — all three reported, always together

| tier | definition |
|---|---|
| **V0 (S1)** | RMSD-only success on the declared endpoint. No validity filter. **This is the archive's 48/84 and 32/84** and the only tier comparable to it. |
| **V1 (PoseBusters-as-run)** | `success_pb` = V0 ∧ `pb_pass`, with the harness's **current** check set exactly as it runs today — including both water checks (`LIB/PoseBust/BustCli.cpp` deliberately selects them, matching upstream `redock.yml`). |
| **V2 (protocol-scoped)** | V0 ∧ every PoseBusters check passes **except** `minimum_distance_to_waters` and `volume_overlap_with_waters`, recomputed from `pb_failed_keys` by removing exactly those two keys. |

### 9.1 Why V2 exists — the scientific justification, not a convenience

**FlexAID models solvent IMPLICITLY. The model's definition of "solvent" is
"space with no atom in it".** It never places a water molecule, never scores one,
and has no term whose value depends on where a water is. The engine strips waters
on **504 of 504 cells** of this campaign, so the receptor the pose was generated
against contains **no waters at all**.

A water check therefore evaluates one of two things, and neither is attributable
to the docking:

1. against a water-free receptor, it is **vacuous** — the harness itself marks it
   skipped (`ChecksProtein.cpp` emits `vacuous(...)` for both water keys when the
   crop has no waters); or
2. against waters re-imported from the deposited structure, it charges the pose
   for overlapping molecules **the model was never shown and cannot represent**.

In case (2) a failure is a property of the *protocol's receptor preparation*, not
of the *scoring function under test*. It is out of protocol scope.

### 9.2 The three rules that keep V2 honest

1. **V2 is reported ALONGSIDE V1, never instead of it.** A table containing V2
   without V1 is a violation of this pre-registration.
2. **V2 is not the headline.** The confirmatory endpoint (§5) is V0 for P2 and a
   single named check for P1. V2 is context.
3. **The previously observed water-dominated validity figure must not be
   published as this model's validity**, in either direction — neither as a
   damning number nor as one V2 conveniently removes.

---

## 10. POWER CHECK

Two-sided exact-binomial (McNemar) at p = ½. Minimum recoveries `b` needed for
two-sided p < 0.05, as a function of new losses `c`:

| `c` (new losses) | minimum `b` | `D = b + c` | p at that point |
|---|---|---|---|
| 0 | 6 | 6 | 0.0313 |
| 1 | 8 | 9 | 0.0391 |
| 2 | 10 | 12 | 0.0386 |
| 3 | 12 | 15 | 0.0352 |
| 4 | 13 | 17 | 0.0490 |
| 5 | 15 | 20 | 0.0414 |
| 6 | 17 | 23 | 0.0347 |
| 8 | 20 | 28 | 0.0357 |

**Minimum discordance for significance is D = 6**, and only at perfect imbalance
(`b = 6, c = 0`). The archived contrast had **D = 18** (`b = 17, c = 1`).

**Can the planned cell count reach it?** Yes.

* **P2.** 84 paired targets; the treated arm can in principle flip any of them.
  Flexibility cost 16 targets, so `b ≤ 16 + (targets rigid also missed)`, giving a
  realistic ceiling well above the required `b ≥ 10` at `c ≤ 2`. The archive
  itself observed D = 18 on this exact population. **Adequately powered.**
* **P1.** 77 paired cells with 38 discordant failures already observed
  (`b = 38, c = 0`, p = 7.3e-12). The required `b ≥ 30, c ≤ 2` is far above D = 6.
  **Adequately powered.**

**Therefore the campaign MAY run.**

### 10.1 The honest limit of that power

If `T4_all` recovers **fewer than 6 targets with zero new losses**, no exact test
on 84 paired targets can reach p < 0.05, and the correct report is a **null with
a stated detection floor** — "an effect smaller than ~6/84 is not detectable by
this design" — not a trend, not a "direction of improvement", and not a larger
campaign justified post hoc on this data.

### 10.2 Not powered, and declared as such

The **per-target** claims (`T1` vs `T2` vs `T3` decomposition, per-target
attributions from the contact-profile sidecars) are **not powered** and are
labelled exploratory wherever they appear. No p value is computed for them.

---

## 11. ANALYSIS PLAN (fixed here; the driver implements exactly this)

1. Build **once**. Every arm uses the **same** binary; the engine SHA-256 is
   recorded in every cell receipt and asserted equal across arms.
2. Gate G1 (inertness) → G2 (wiring) → G3 (liveness) → G4 (control replication)
   → scoring. **Each gate writes a PASS/FAIL verdict file; the next stage refuses
   to start without PASS.** An exit code of 0 is not a gate.
3. Per cell, from `result.csv`: `rmsd_to_crystal`, `rmsd_hungarian`,
   `success_rmsd`, `pb_pass`, `success_pb`, `pb_failed_keys`, `pb_volume_overlap`,
   `pb_min_lig_prot_dist`, `pb_ran`, `elected_pose_path`, `cf_native`,
   `pose_sha256`. Guard `rmsd >= 0` on every filter.
4. E1 recomputed offline from `REMARK CF.app` over the cell's pose PDBs; E2 read
   from the harness. Both in the table.
5. Paired 2×2 tables `T4_all` vs `B_shrink` for P1 and P2; two-sided exact
   binomial. `T1`–`T3` the same, labelled EXPLORATORY.
6. `A_rigid` reported for context; the confirmatory contrast is against
   `B_shrink`, because H1 is a statement about what flexibility broke.
7. Report V0, V1, V2 as three columns for every arm.
8. No arm is dropped, no target is dropped, no threshold is moved.

---

## 12. RISKS AND KNOWN CONFOUNDS, stated before the fact

* **Dimensionality.** The flexible arms carry one extra gene per flexible residue.
  `B_shrink` carries it too, so the confirmatory contrast is balanced on it — but
  `A_rigid` is not, which is exactly why `A_rigid` is context and not the
  comparator.
* **Determinism.** OMP thread count changes results by up to ~2 Å per target.
  Every arm uses the same `--omp-threads`; it is recorded per cell and asserted
  equal.
* **`receptor_rotamer_prep`.** If on, stage 3's frozen solvent reference is the
  **pre-relaxed** conformer, not the deposited crystal one. The harness prints
  this; the value is recorded per cell and must be identical across arms.
* **Disk.** `FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1` writes one receptor-sized PDB per
  emitted pose. **A single-target pilot must size this before the full grid**; if
  the projected total exceeds the free-space budget, the instrument is turned off
  in **all** arms, not just some.
* **Cost.** 6 arms × 84 targets × 3 restarts × 3000 generations. This is a large
  campaign. It should not start until G1–G3 have passed on one target.

---

## 13. AMENDMENTS

Any change to §5, §6, §7, §8 or §9 after the first treated cell completes makes
this an exploratory study. Amendments are appended below with a UTC timestamp and
the reason; nothing above is edited in place.

*(no amendments)*
