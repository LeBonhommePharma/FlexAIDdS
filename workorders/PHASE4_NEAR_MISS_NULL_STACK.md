# Phase-4 near-miss sampling stack — frozen null results (publication table)

**Status:** FREEZE (offline gate; no dock)  
**Date:** 2026-07-28  
**Defers to:** `BENCHMARK_SELF_EVAL_CONTRACT.md` (status enum, residual path),  
`a_posteriori_gate_ledger.md` (canonical SoT), per-gate posteriori workorders.  
**Purpose:** Single SI/Methods table so Phase-4 near-miss negatives cannot be re-labeled as ACCEPT.

Panel class for all rows: **NEAR_MISS** (`1N1M`, `1L7F`) unless noted.  
Matrix pin: **9dc9** (`md5 9dc93717dfed0698006d88dd6a9627bc`).  
Phase-4 docks: **NO_SEC=1**. Full-85: **BLOCKED** until sampling ACCEPT.

---

**SUPERSEDED 2026-07-31.** The next-step order recorded in this file
(`new_search_arch` first, `scoring_locked_decoy_work` behind it) was
**reversed** by decision in the Buzz channel *Benchmarking FlexAIDdS*
(`716e79b1-6a4f-4cde-8851-22836ba8738c`), authorized by Bonhomme.

**Current order:** instrument fix → frozen-gate replication → scoring-locked
decoy work → `new_search_arch` *(deferred, not cancelled)*.

**Rationale.** The sampling hypothesis is a closed null: six one-variable
gates, best magnitude −0.019 Å against a −0.5 Å floor. Nothing came close.
The burden of proof therefore sits on `new_search_arch` and it has not met
it. Scoring-locked decoy work probes the other candidate wall
(discrimination) at a fraction of the cost.

`new_search_arch` gets budget only if the decoy work fails to move the
needle. Deferred is not cancelled — do not delete it from this stack.

Provenance: see `PLANS/FLEXAIDDS_PHASE5_DECISION_RECORD.md` for the cited
event IDs.

---

## ⚠ `pb_pass` is dominated by water compatibility — reporting rules

Measured over all 36 real PoseBusters 0.6.5 CSVs in
`/Users/lp.more/flexaidds_benchmark_results/` (34 scored poses). Computed
independently by Bumble and by Opus, with separate reimplementations of the
evaluator; every cell agreed.

| | |
|---|---|
| `pb_pass` TRUE as the gate stands | **10/34 — 29%** |
| `pb_pass` TRUE with the two water columns excluded | **24/34 — 71%** |
| Poses failing on water **and nothing else** | **14 — 41% of the panel** |

Failing columns, most frequent first: `minimum_distance_to_waters` 22/34,
`volume_overlap_with_waters` 12/34, `minimum_distance_to_protein` 6/34,
`internal_steric_clash` 3/34, `minimum_distance_to_inorganic_cofactors` 1/34,
`internal_energy` 1/34. Zero NaN/blank cells — this is not uncomputed-value
poisoning, it is the gate doing what it was built to do.

**Every non-excluded column is ANDed into `pb_pass`.** So 41% of scored poses
are chemically clean against protein, cofactors, internal geometry and energy,
and fail the claim gate solely because they displace crystallographic water.

### ⚠ SCOPE — this is a property of the receptor tier, not of the docking

All 36 bust receipts record a runtime-cache receptor. Exact audit over the
31 unique targets:

```
receipts using ~/.flexaidds runtime cache       36/36
cache receptor == repository canonical           0/31
cache receptor contains HOH                     30/31
repository-canonical receptor contains HOH       0/31
cache receptor == deprecated tracked prep       24/31 unique targets
cache receptor is a THIRD, unidentified variant  7/31 unique targets
```

**Under the canonical preparation these columns cannot fail** — there are no
crystallographic waters to displace, so both water checks are vacuous by
construction.

### State the interaction precisely — an earlier draft of this block did not

It is **wrong** to say "an artifact of the receptor tier, not a property of
the docking." Receptor preparation makes a water failure *representable*; the
elected pose determines whether it *actually occurs*. The correct statement is
a **protocol-tier × pose interaction**: the legacy tier defines a
water-sensitive endpoint, and the docking output supplies the outcome.

### Frozen boundaries

- 29% → 71% is valid **for the legacy / noncanonical runtime-cache protocol
  only.**
- It is **not** an independent finding from receptor identity.
- It does **not** estimate the canonical-prep strict pass rate.
- **Do not claim all runtime-cache receptors are the deprecated tracked prep.**
  Seven targets — `1GPK, 1HQ2, 1L7F, 1MEH, 1N1M, 1N2J, 1OPK` — are a third
  variant whose provenance is **unresolved**. (`1GPK` appears in both arms, so
  28/36 *receipts* match the deprecated tree exactly, not 36/36.)

### The vacuity is not limited to water — generalize the caveat

The seven third-variant targets keep **every water** but have their non-water
heteroatoms stripped:

```
       non-water HETATM   stripped relative to deprecated prep
1GPK   cache 0 / dep 56    NAG
1HQ2   cache 2 / dep 17    CL, MG, PH2
1L7F   cache 2 / dep 127   CA, BMA, MAN, NAG
1MEH   cache 1 / dep 8     K, CSO
1N1M   cache 0 / dep 403   HG, BMA, FUC, MAN, NAG, NDG
1N2J   cache 0 / dep 12    BAL
1OPK   cache 0 / dep 16    MYR
```

**The vacuity here is asymmetric, and saying "all four cofactor columns on
all seven" would overstate it.** Checked against PoseBusters 0.6.5's own
element classification (`posebusters/tools/protein.py`), which counts Mg, Ca
and K as inorganic-cofactor elements:

- both **organic**-cofactor checks are vacuous on **all seven**;
- both **inorganic**-cofactor checks are vacuous on **four of seven**
  (`1GPK, 1N1M, 1N2J, 1OPK`);
- `1HQ2` (2 × MG), `1L7F` (2 × CA) and `1MEH` (1 × K) retain inorganic
  cofactors, so those checks remain representable there.

The mechanism is the same as for water: a receptor stripped of a species
makes a clash with that species unrepresentable.

**Scale, so this is not over-read:** 8 of 34 scored poses come from these
seven targets, and the one cofactor failure in the whole set (`1JJE`) is not
among them. **This changes no number today.**

**The general rule it establishes, which is the durable part:** the strict
endpoint is **vacuous on different column subsets for different target
subsets**, depending on what each receptor preparation retained.

`pb_pass` itself remains **formally well-defined** — it is the same AND over
27 booleans regardless of target. What varies is the *opportunity profile*:
which components were capable of failing at all. So: **aggregate `pb_pass` is
well-defined only relative to the exact fixed receptor-preparation manifest,
and it is not transportable or comparable across receptor tiers/manifests
without reporting which checks were non-vacuous per target.** The evaluator is
not undefined; the aggregate is not portable.

The evidence set contains **at least three receptor tiers**, not two.

Every strict `pb_pass` figure the campaign produced is a figure against a
noncanonical receptor tier that the repository does not sanction for new work.
Codex's reporting semantics above are unaffected and still correct.

### Why this outranks the pre/post-fix boundary

The boundary in Block 2 separates two eras. **This one applies to both of
them.** Pinning never changed what is scored (the pinned list is a presence
assertion on the header; the evaluator walks every non-excluded column
regardless). So the water dominance predates #310 entirely and is unaffected
by it. Any `claim_ready` figure from any era is a water-displacement-filtered
figure.

### Water is a treatment-sensitive mediator, NOT a confound

The two arms in this tree are `armA_mincf` and `armA_smartwater_rawcom`. One
of them varies **water handling**, so water compatibility sits *downstream* of
the treatment.

```
armA_mincf              n=29  pass=9   pass_if_water_excluded=21
armA_smartwater_rawcom  n=5   pass=1   pass_if_water_excluded=3
```

**Terminology matters here and an earlier draft of this note got it wrong.**
A confound is a *pre-treatment* variable correlated with both treatment and
outcome. Water compatibility is not that: the smart-water arm changes water
handling, so water compatibility is a **treatment-sensitive outcome/mediator.**
It becomes **selection bias** only if we filter on full `pb_pass` and then
compare pose accuracy among the survivors. That is a narrower and more
accurate hazard than "confound," and it points at a different remedy —
do not subset, rather than do not measure.

**No arm effect is claimed and none can be** — n=5 on one side supports
nothing.

### The scientific call — DECIDED, with upstream basis

**Should displacing a crystallographic water fail a redock claim? Yes — in
the strict endpoint.** Not by default and not by accident: PoseBusters 0.6.5
`redock.yml` explicitly selects both `minimum_distance_to_waters` and
`volume_overlap_with_waters` (`redock.yml:199-215, 268-283`). Full upstream
`pb_pass` includes retained-water compatibility **by definition.** Dropping
those fields would not repair the metric; it would create a custom metric that
is no longer "PoseBusters pass."

The 29% → 71% counterfactual does not show the upstream check is wrong. It
shows that **water compatibility dominates this dataset and protocol** — and
that the campaign had been calling a composite environment-compatibility
endpoint "chemistry."

### Locked reporting semantics

1. **Strict headline.** Fixed denominator, identical elected pose, direct
   whole-ligand RMSD ≤2 Å **and every upstream 0.6.5 redock boolean including
   water.** The names `pb_pass` / "strict" are reserved for exactly this
   endpoint.
2. **Diagnostic decomposition**, on the same fixed denominator: report
   `pb_pass_nonwater`, `water_distance_pass`, `water_overlap_pass`, and
   `water_only_failure` alongside it. **Never promote the 71% non-water rate
   to the strict headline.**
3. **Arm comparison.** Report each endpoint over the same predeclared target
   denominator. **Do not subset to `pb_pass==1` before comparing arms** — that
   subsetting is the selection bias.
4. **Protocol changes.** Stripping waters, marking waters displaceable, or
   using a custom PoseBusters config is a **new protocol tier and its own
   comparability boundary** — not a silent repair of historical numbers.

Decided by Codex, 2026-07-31, on upstream `redock.yml` evidence.

### Scope of this caveat — audited, not assumed

It applies to anything gated on `pb_pass` / `claim_ready`. It does **not**
reach the Phase-4 sampling conclusions.

**This was checked, not inferred.** At `main` `11ce273c`, these eight
workorder documents contain zero occurrences of `pb_pass` or `claim_ready`:

```
workorders/PHASE4_GATES_ACTUALIZED.md    pb_pass=0  claim_ready=0
workorders/CAMPAIGN_GATE_SUMMARY.md      pb_pass=0  claim_ready=0
workorders/a_posteriori_gate_ledger.md   pb_pass=0  claim_ready=0
workorders/G4_2_NICHE_CART.md            pb_pass=0  claim_ready=0
workorders/G4_4_EARLY_STOP.md            pb_pass=0  claim_ready=0
workorders/NEW_SEARCH_ARCH_APRIORI.md    pb_pass=0  claim_ready=0
workorders/INVERSION_MAP.md              pb_pass=0  claim_ready=0
workorders/COARSE_ORIENT_MATCHED_AB.md   pb_pass=0  claim_ready=0
```

Every Phase-4 acceptance criterion in `PHASE4_GATES_ACTUALIZED.md` is stated
in RMSD/BCR terms — `mean dRMSD <= 0 on >=4/5 SEARCH-MISS`,
`mean dBCR <= -0.5 A or >=1 target crosses BCR<2 A`, `n_niches occupied must
increase`, generations-reached distribution. **No Phase-4 gate reads
`pb_pass`.** The sampling null therefore stands independent of everything in
this block.

One artifact mentions `pb_pass` — `S4_PHENOTYPE_UNIQUE_STATUS_REPORT.md` —
as a reported side column carrying its own warning not to promote it to a
STRICT claim without a full receipt and bust_cli audit, and it lists any
`claim_ready` implication as out of scope. The exception confirms the rule.

**Scope limit on this audit itself:** it covers the eight documents named
above at commit `11ce273c`. It is not a claim about every artifact in the
repository. Audited by Opus, 2026-07-31, static read.

Provenance: Opus and Bumble, 2026-07-31, static analysis over existing
campaign outputs. No docking or rebuild was run to produce these numbers.

---

## ⚠ Comparability boundary — 2026-07-31

**Numbers produced before this date are not comparable to numbers produced
after it.** The discontinuity has **two independent axes**, and a reader who
knows about only one of them will still draw a wrong conclusion.

| Axis | Before | After |
|---|---|---|
| **Backend** | `native_pose_qc_fallback` — a silently degraded fallback | real PoseBusters 0.6.5 |
| **Check definition** | pinned `no_protein_clashes` — a column upstream 0.6.5 **never emitted** | the **canonical 27** PoseBusters 0.6.5 redock booleans, water included |
| **Gate authority** | scoring iterated *the CSV's* headers; the pin only asserted presence | scoring iterates **the pinned list**; set equality enforced in both directions |

The third axis is new as of PR #310 (`0af9fe82`, refined at `3734ee7d`) and
was **verified a no-op on today's data** before the semantics changed — 34
CSVs, one header layout, exactly the 27, no extras, identical verdicts. So it
introduces no numeric discontinuity on 2026-07-31. It changes what the metric
*is* going forward: a schema change in either direction now fails closed and
demands a deliberate pin bump, where previously an added column would have
been silently ANDed in and a removed one silently dropped.

Gate membership is decided by a **literal four-name metadata exclusion** —
`file`, `molecule`, `position`, `rmsd_≤_2å` — not by substring match. 31
columns minus those four is exactly the 27.

**VERIFIED against upstream — no sharp edge.** The RMSD column name is a
literal string in the 0.6.5 config, not computed from the threshold:

```yaml
# posebusters/config/redock.yml:286-292
function: rmsd
rmsd_threshold: 2.0
rmsd_within_threshold: "RMSD ≤ 2Å"
```

Under stock upstream 0.6.5 the name is fixed and cannot drift. It changes
only if someone edits that label in a custom config — which is already a new
protocol tier requiring its own comparability boundary. Failing closed there
is the guard working, not a false alarm.

**Incidental finding, recorded because it is a trap.** `rmsd_threshold` and
the display label are **independent fields**. Setting `rmsd_threshold: 5.0`
without editing the string emits a column still labelled `RMSD ≤ 2Å` carrying
a 5 Å verdict. Nothing consumes that column today, so it cannot reach a
number — but **the label is not evidence of the threshold**, and anything
that starts consuming it must read `rmsd_threshold`, not the header.



Either axis alone changes what `claim_ready` means. Together they mean a
pre-fix `pb_pass` and a post-fix `pb_pass` are **different measurements that
share a name.**

**Rule:** no table, plot, or claim may place a pre-fix and a post-fix number
in the same column without carrying this note. If the two must appear
together, label the axis explicitly per row.

**Additional caution — the fallback failed silently.** It degraded and went
unnoticed for an entire campaign. That is the specific failure mode this
boundary exists to prevent recurring: not a wrong number, but a wrong number
that looked fine.

**Why the old pin was worse than a wrong name.** `no_protein_clashes` was not
a check that changed — it was a column upstream never emitted. The parser
fails closed on a missing mandatory header, so every real bust run returned
`pb_pass=0` for a *schema* reason and the campaign silently fell back to
`native_pose_qc_fallback`. Pre-fix `pb_pass` therefore does not encode a
chemistry verdict at all. Do not read pre-fix values as weak evidence; read
them as no evidence.

Provenance: PR #310, commit `5a5b7a55`. Column evidence re-derived from a real
0.6.5 bust CSV (`astex85_threearm_20260722_224149/.../1IA1_..._bust.csv`).

### RESOLVED — the pre-fix era is five bands, not one regime

This note marks where the pre-fix era *stops*. It does not follow that the
pre-fix era was internally uniform, and it was not: `main` accumulated the
four stacked CI failures over **17 days**, not all at once.

A failure becomes reachable at `max(defect introduced, checker introduced)` —
a defect with no checker is inert; a checker with no defect is silent.

| Failure | Defect | Checker | **First reachable** |
|---|---|---|---|
| bare `g.log_Xi` | `0da7c25a` 2026-07-08 | `6e3d17bd` 2026-07-08 | **2026-07-08** |
| stale `test_palette_cycles` | `2fb4eb76` 2026-07-09 | predates | **2026-07-09** |
| orphan `metal_microbench_enhanced.cpp` | `69aa0fab` 2026-07-14 | `17ae1d15` 2026-05-25 | **2026-07-14** |
| `competition_example.yaml` no `docking_mode` | `b27889c0` 2026-07-08 | `2fffe6fe` 2026-07-25 | **2026-07-25** |

Two asymmetries defeat a naive "when was the bug written" reading, in
opposite directions. The YAML file was written 2026-07-08 and sat harmless
for 17 days — it became a failure only when the validator turned fail-closed
on 2026-07-25, making the oldest file the newest failure. Conversely the
orphan guard existed from 2026-05-25, seven weeks before anything tripped it.

**The five bands:**

```
before 2026-07-08   0 of 4 present
2026-07-08          1 of 4   (log_Xi)
2026-07-09          2 of 4   (+ palette)
2026-07-14          3 of 4   (+ orphan guard)
2026-07-25          4 of 4   (+ docking_mode)   <- the state diagnosed in #311
```

**Rule: two pre-fix numbers are comparable to each other only if they fall in
the same band.** Check the date before comparing.

**Read this as narrowing, not widening.** Before this was dated, "unknown
internal structure" meant *any* two pre-fix numbers might be incomparable.
Now it means four specific dates, and most pairs are fine. The era is
tractable.

**Two limits on the dating itself:**

1. These are dates of **reachability derived from code history** — when each
   failure *could* fire. They are **not** a reconstruction of observed CI
   history. Nobody has established which days CI actually ran red, and per
   the stacked-failure lesson a red matrix reports only the first failure, so
   the observed record would understate the count regardless.
2. Scoped to **these four failures** — the ones enumerated in #311 — at
   `11ce273c`. No search was made for a fifth.
3. The orphan-guard date is scoped to **the 11 orphans that were silenced**,
   not to the guard's whole life. No one has shown that no other orphan
   existed between 2026-05-25 and 2026-07-14 and was later resolved.

### Two different durations wear the same units — say which one you have

Birth order was **log_Xi (07-08) → palette (07-09) → orphan guard (07-14) →
docking_mode (07-25)**. Discovery order was almost exactly inverted: the
orphan guard and the dataset-schema failure were found first, the two
~3-week-old Python defects last.

The reason is mechanical. The C++ configure failure and the schema failure
**abort their jobs before the Python defects are reached at all.** A red
matrix does not surface its oldest defect first; it surfaces its
**earliest-aborting** one.

**Rule: a duration read off CI is a lower bound on defect age, never an
estimate of it.** The durations in the table above are **ages**, derived from
each defect's introducing commit — they are citable as ages precisely because
they were not read off CI.

The two quantities wear the same units and are not the same thing:

| Quantity | Source | What it is |
|---|---|---|
| "`main` has been red for N days" | observed CI | a **floor** — masking hides everything behind the earliest-aborting job |
| "this defect has been in `main` for N days" | introducing commit | an **age** — knowable only after the archaeology |

They coincide here only because the birth commits were found. Before that work
existed, all anyone had was the floor.

This is the same family as the `ctest` count and the column names, but a
distinct member of it. A name can be wrong (`no_protein_clashes`) or
overstated (a test name, the `RMSD ≤ 2Å` label). A duration is neither: it is
**true, and misleading relative to a question it cannot answer.** So the
discipline is not "distrust the number" — it is **state which question the
number answers.** Both durations above are correct; only one of them answers
"how old is this bug."

All commits verified as ancestors of `main` at `11ce273c`
(`git merge-base --is-ancestor`, all yes). Archaeology by Opus, 2026-07-31;
full detail including all 11 orphan dates in
`RESEARCH/MAIN_CI_FOUR_STACKED_FAILURES_ARCHAEOLOGY.md`.

---

## Baseline provenance requirement

A frozen number is only frozen if it is regenerable. Every pinned baseline
must cite `{commit SHA, binary SHA, config path}` **and every one of those
must be tracked in git.** A baseline bound to an untracked config is not
reproducible from the state it names.

This is not hypothetical: the deletion of a pilot leaf is what made the
Phase-4 SCORING-LOCKED `cf_native` values unreproducible, and the
replacement configs were themselves untracked when the re-pin began.

### Input identity matters as much as config identity

The SCORING-LOCKED baseline exists in **two protocol tiers that are not
interchangeable**, and the difference is chemically material.

| Tier | Receptor prep | 1OQ5 | 1SQ5 | 1YGC |
|---|---|---|---|---|
| **Repository-canonical Astex** | `benchmarks/astex_diverse/astex_diverse/{PDB}/{PDB}_apo.pdb` | `-34.229803` | `-73.790961` | `+7.303774` |
| **Legacy production-runtime** | `benchmarks/astex_diverse/data/astex_diverse/{PDB}/` — a **deprecated historical second prep** | `-34.229803` | `-73.413683` | `-0.871396` |

The legacy tier is what the runtime cache (`~/.flexaidds/benchmarks/...`)
holds, byte-identical, and it is what reproduces the historical workorder
numbers. **It is not canonical.** `benchmarks/datasets/CANONICAL.md` and
`benchmarks/astex_diverse/README.md` define the canonical tree and classify
`data/astex_diverse/` as deprecated — "do not use for new work."

Immutable snapshots of the deprecated prep are retained at
`ops/gates/configs/legacy_runtime_receptors/`. They are **snapshots of a
deprecated prep, not canonical inputs**, and the directory name says so.

**Why the numbers move:** the deprecated prep **retains waters** —
233 / 885 / 260 across the three targets — plus Zn/Ca where present. The
canonical apo copies do not. That is the mechanism behind the 1YGC sign flip
from `-0.871396` to `+7.303774`.

Ligand hashes matched throughout. Only the receptors differed, which is why
this hid until hash-level inspection.

**Rule:** any baseline citing SCORING-LOCKED numbers must name its tier.
A figure reproducing the historical workorder is a *legacy production-runtime*
figure, not a canonical one, and the two must never share a column.

### RESOLVED — the receptor tier *causes* the water dominance

This was raised as "the two may be entangled." They are not merely entangled.
**They are cause and effect**, and the causal direction runs from this block
to Block 4.

Every bust run records its own argv. All 36 receipts name a runtime-cache
receptor:

```
36/36  -p /Users/lp.more/.flexaidds/benchmarks/astex_diverse/<PDB>/<PDB>_apo.pdb
```

**The cache is not a single tier.** Across 31 unique targets: 24 are
byte-identical to the deprecated tracked prep; **7 are a third variant**
(`1GPK, 1HQ2, 1L7F, 1MEH, 1N1M, 1N2J, 1OPK`) matching neither canonical nor
deprecated, **provenance unresolved and unowned.** None match canonical.

Water content across the full evidence set:

```
runtime cache        : 30 of 31 receptors carry waters (28-931 HOH; 1IGJ has 0)
repository-canonical : 31 of 31 carry ZERO waters
```

**The water columns cannot fail against a water-free receptor.** Under the
repository-canonical preparation there are no crystallographic waters to
displace, so `minimum_distance_to_waters` and `volume_overlap_with_waters` are
vacuous by construction.

**What this does to the 29% → 71% figure is stated once, in Block 4, and is
not restated here.** See "SCOPE — protocol-tier × pose interaction" and the
frozen boundaries beneath it. Do not paraphrase that framing into this block:
an earlier draft carried its own copy, Block 4's copy was corrected and this
one was not, and the document contradicted itself for one commit. **One
claim, one home.**

**KNOWN GAP — the binary is not commit-stamped.** The receipt records the
binary SHA256 and its build date (2026-07-28) but **not** the source commit,
because the build did not capture one (`source_commit: null`). Do not assume
it corresponds to `11ce273c`. So the baseline is replayable only while a
binary matching the recorded SHA is retained: **its inputs are tracked; it is
not source-regenerable from tracked state.** Closing this needs a rebuild
with commit stamping.

---

## Stack table

| Gate | One variable | R | Status | Magnitude | L4 | OUT | Workorder |
|------|--------------|--:|--------|-----------|----|-----|-----------|
| G4.1 BOOM | `BOOM_FRAC` ∈ {0.05,0.1,0.2} vs unset | 2 | **FAIL (null mag)** | best mean_dBCR **−0.019** (frac010); floor ≤−0.5 or BCR&lt;2 | control 0 / tx 236 [BOOM] **PASS** | `~/flexaidds_results/g4_1_boom_near_miss_20260726_200953` | `G4_1_NEAR_MISS_POSTERIORI.md` |
| ELECTION_V135 | `ELECTION_V135=1` (τ=25) vs unset | 5 | **FAIL (null mag)** | elect identical 6.3999 / 3.9907; gap shrink 0 | protocol markers live both arms | `~/flexaidds_results/election_v135_near_miss_20260726_225823` | `ELECTION_V135_POSTERIORI.md` |
| G4.3 MUTATION | `MUTATION_GRANULAR=1` vs unset | 2 | **PASS_LIVENESS** | mean_dBCR **+0.118**; 1L7F elect 3.99→6.25 | control 0 / tx 8 [MUT-GRAN] **PASS** | `~/flexaidds_results/g4_3_mut_gran_near_miss_20260727_122215` | `G4_3_MUTATION_POSTERIORI.md` |
| S4 A PHENO_UNIQUE | `PHENOTYPE_UNIQUE=1` vs unset | 2 | **PASS_LIVENESS** | mean_dBCR **−0.057**; 1N1M elect still 6.40 | control 0 / tx 4 [NEW-SEARCH-ARCH] **PASS** | `~/flexaidds_results/s4_pheno_unique_near_miss_20260727_211213` | `S4_PHENOTYPE_UNIQUE_POSTERIORI.md` |

### G4.1 best arm detail (from OUT flip_order_decision)

| arm | 1L7F BCR | 1N1M BCR | mean_dBCR |
|-----|---------:|---------:|----------:|
| control | 3.9907 | 4.5515 | — |
| frac005 | 4.0834 | 4.5515 | +0.0464 |
| frac010 | 3.9523 | 4.5515 | **−0.0192** |
| frac020 | 4.0834 | 4.5515 | +0.0464 |

### G4.3 detail (from evidence/g4_3_posteriori.txt)

| arm | 1L7F elect/BCR | 1N1M elect/BCR |
|-----|----------------|----------------|
| control | 3.9907 / 3.9907 | 6.3999 / 4.1954 |
| mut_gran | 6.2458 / 4.1128 | 6.0053 / 4.3088 |

Treatment binary SHA256 (post dup-clear relaunch):  
`19f300d9798d4985423fb501697ee3b397cc57596040ccfe7c84a8f8165225f6`  
Git tip (fix): `b8b19468` (dup-clear) / posteriori commit lineage on `fix/dump-pop-refstructure-autonomous`.  
**Match caveat:** control used pre-relaunch binary path; treatment used post-`b8b19468` binary—document in SI; magnitude null so no false PASS.

---

## Flip residual (after stack)

```
rule: G4.3_null_phase4_sampling_exhausted
priority_order: [new_search_arch, scoring_locked_decoy_work, full85_still_blocked]
```

Machine source: `g4_3_mut_gran_near_miss_20260727_122215/evidence/flip_g4_3.json`.

---

## Publication residual path (contract)

1. ~~Phase-4 sampling ACCEPT on near-miss~~ → **not met** (this freeze).  
2. Optional non-burial scoring residual for **SCORING-LOCKED** (class-matched only).  
3. Full-85 — **blocked**.  
4. Claim language: CF/contact-function scoring proxy; no true ΔG without STRICT (PB + tENCoM).

---

## Pins required on every future closed gate

matrix 9dc9 · per-arm binary SHA256 · git tip · FLEXAIDDS_* env · R · pop/gen · NO_SEC · Sol #9 · L4 stderr+r* · status enum · `evidence/accept.txt`.

Enforcers: `scripts/benchmark_self_eval.py`, `scripts/campaign_flip_order.py`, `scripts/benchmark_coord.py`.

---

## Explicit non-claims

- Not genuine top-1 / PoseBusters rates.  
- Not permission to enable BOOM_FRAC, ELECTION_V135, or MUTATION_GRANULAR in claim recipe.  
- Not evidence for memetic unlock or burial re-panel.
