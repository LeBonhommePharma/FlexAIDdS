# Admission + Metrics Contract (Normative)

**Comparative science hub:** `docs/implementation/COMPARATIVE_SCIENCE_README.md`  
**Status:** Normative for claim-table aggregation and abstract / headline rates.  
**Aligned with:** `benchmarks/protocols/three_engine_entropy_comparison.md` §1.4–§5, `AGENTS.md`, audit 2026-07-17.  
**Enforcement:** `scripts/aggregate_claim_metrics.py` (fail-closed).

---

## 1. Success metrics (report all; headline is claim_ready)

| ID | Definition | Role |
|----|------------|------|
| **S1** | Elected pose **ordered direct** `rmsd_to_crystal` ≤ 2.0 Å and `seed_echo==0` | RMSD-only diagnostic of election |
| **S2** | S1 ∧ `pb_pass` on the **same** elected pose (`success_pb`) | Intermediate (RMSD∧PB) |
| **STRICT** / **claim_ready** | Engine `claim_ready==1`: S2 ∧ official PoseBusters ∧ tENCoM/Eigen complete on exact pose SHA-256 ∧ protocol eligibility ∧ score–pose consistency | **Primary / claim KPI** |
| **S3** | `conditional_scanned_pool_ceiling` / `best_cluster_rmsd` ≤ 2.0 Å | **Diagnostic only** (scanned heads/members, **not** any-pose) |

### Field mapping (DatasetRunner `result.csv`)

| Metric | Required / preferred fields |
|--------|------------------------------|
| Elected RMSD (S1) | **`rmsd_to_crystal` only** (ordered direct, whole-ligand). **Never** `rmsd_hungarian` for S1/admission |
| Hungarian | `rmsd_hungarian` — diagnostic symmetry metric only |
| PoseBusters | `pb_pass`, `success_pb`, `pb_backend==bust_cli` for claims |
| Strict claim | **`claim_ready==1`** (and receipts: pose hashes, `tencom_status==ok`, `eigen_status==ok`) |
| Pool ceiling (S3) | `conditional_scanned_pool_ceiling` → else `best_cluster_rmsd` (scanned emission pool only) |

### Hard rules

1. **Headline claim rate = claim_ready only** (not S1, not S3).
2. **S1 uses ordered direct RMSD only** — never Hungarian for claim S1.
3. **Never report S3 as abstract / headline success.** Label as diagnostic conditional scanned-pool ceiling.
4. **Always report S1, S2, STRICT (claim_ready), and S3 separately.**
5. **Election gap** = S3=1 and S1=0 (scanned pool had ≤2 Å; elector missed). Diagnostic only.
6. Generator **CF top-1** and **entropy/consensus reranked top-1** are separate estimands (paths, scores, hashes, RMSDs) when both are present.

---

## 2. Claim admission (row must pass all)

A row is **claim-eligible** only if:

| Check | Required value |
|-------|----------------|
| `seed_echo` | **0** (explicit; missing fails closed) |
| `native_pose_seeded` | **0** (explicit; missing fails closed) |
| `matrix_md5` | equals campaign matrix pin when present |
| `protocol_claim_eligible` | true when column present |
| **`claim_ready`** | **1** when column present (strict admission) |
| Hash receipts (when present) | `rmsd_pose_sha256` / `posebusters_pose_sha256` / `tencom_pose_sha256` match `pose_sha256` |
| Validator status (when present) | `tencom_status==ok`, `eigen_status==ok`, `pb_backend==bust_cli` |

Rows missing `claim_ready` but satisfying seed gates are admitted only for **legacy diagnostic** tables; they **must not** contribute to STRICT headline. The aggregator reports them under `N_legacy_no_claim_ready` and excludes them from STRICT rates.

### Default matrix pin

| Field | Value |
|-------|--------|
| Canonical matrix | campaign-dependent (see RUN_RECEIPT) |

---

## 2b. Solvent condition (NORMATIVE — added 2026-09-04)

**Every PoseBusters-derived rate in this contract is conditional on which waters the
validator saw, and until now that was not the set the engine scored against.**

### The engine's shipped policy

`config_parser.cpp:362-369` defaults to `remove_water=true`,
`keep_structural_waters=true`, `structural_water_bfactor_max=20.0`.
`modify_pdb.cpp:158-166` implements that combination as **B-factor-thresholded
conserved-water retention**: an `HOH` is dropped only when its B-factor exceeds
the threshold. It is *not* full stripping. The banner at `modify_pdb.cpp:59`
prints whenever `exclude_het || remove_water` and is **not** a count.

Measured with `FLEXAIDDS_WRITE_FLEXED_RECEPTOR=1` on 1JD0: the receptor as scored
holds 4325 atoms including 156 `HOH`, every one at B ≤ 20 and none above; 504
waters in, 156 retained. The loader's 4338 reconciles exactly as
4169 non-water + 156 conserved + 13 ligand. Across the 84-target roster,
**2822 of 28608 crystallographic waters (9.9%) reach the search**; 18 targets have
waters and retain none; 1 target (1IGJ) has none at all.

### Hard rules

7. **Every PB-derived rate MUST name its receptor condition.** A bare
   `pb_pass`/`success_pb`/`claim_ready` rate is incomplete without it.
8. **The claim-bearing condition is ENGINE-MATCHED**: the conditioning receptor
   passed to `bust -p` must contain exactly the waters the engine retained
   (`<T>_apo.pdb` minus `HOH` with B > `structural_water_bfactor_max`).
9. **Full-receptor PB rates are SUPERSEDED, not merely noisier.** They judge the
   pose against waters the search never saw — 348 per target on 1JD0, ~25,786
   across the set. Measured on 77 paired rigid cells (seed 12345), STRICT rule
   (§2b rule 14): all-waters 41.6%, engine-matched **68.8%**, no-waters 72.7%.
   The full-receptor tier is biased downward by ~27 points and must not be
   reported as a headline.
10. **Water removal can only remove failures**, so an engine-matched rate that
    exceeds a full-receptor rate is expected, not suspicious. Zero losses across
    the correction is the internal check; a loss indicates a defect.

### The PB pass rule (NORMATIVE)

14. **A pose passes only if EVERY check column in `bust_raw.csv` is literally
    `True`**, excluding the three non-check columns (`file`, `molecule`,
    `position`) and the `rmsd_≤_2å` column (`BustCli.h:42` records but excludes
    it). **A blank is a FAILURE, not a pass** — a blank means the check did not
    run. This is 27 check columns.

15. **Do NOT derive the check set by testing which columns are boolean-valued
    across the corpus.** A predicate that admits `''` into the boolean set
    silently drops every sometimes-blank column. Measured consequence on the
    154-cell NO_SEC corpus: such a predicate kept 16 columns and dropped 11 —
    `internal_steric_clash`, `bond_lengths`, `bond_angles`,
    `aromatic_ring_flatness`, `non-aromatic_ring_non-flatness`,
    `double_bond_flatness`, `double_bond_stereochemistry`,
    `tetrahedral_chirality`, `molecular_formula`, `molecular_bonds`,
    `internal_energy`. **`internal_steric_clash` is `False` on 3 of 154 cells**,
    so dropping it makes the rule LOOSER: poses with intramolecular clashes pass.
    Ten of the eleven are vacuous in practice (always `True`, blank on the same
    2 cells, both of which fail other checks); the strict rule changes the rigid
    engine-matched count by exactly 1 cell (54 → 53). The looseness is small but
    it is in the direction that flatters the engine.

### PB denominator: 8 targets are NOT ASSESSABLE under the current invocation

16. **A PB rate MUST state its own denominator, and a target PoseBusters could
    not score is NOT ASSESSED — never a failure.** Measured on the 84-target
    roster: 7 targets wrote a **0-byte** `bust_raw.csv` with `pb_pass=false`, and
    an 8th (`1TW6`) produced a row whose chemistry columns are blank. So a rate
    quoted `x/84` counts 7–8 unmeasurable cells as failures, and the honest
    denominator is **76–77**.

    | Target | Symptom |
    |--------|---------|
    | `1K3U` `1N2V` `1U1C` `1U4D` `1XOZ` `1Y6R` `2BSM` | 0-byte `bust_raw.csv`, `pb_pass=false` |
    | `1TW6` | row present, chemistry columns blank |

    **Cause, established by re-invocation.** The ligand SDF for these targets
    carries an aromatic bond block RDKit cannot kekulize
    (`KekulizeException: Can't kekulize mol`), so `MolFromMolFile` returns
    `None` and `SDMolSupplier` yields 0 non-`None` molecules. This is **in the
    cache input, not written by the engine** — the native `<T>_ligand.sdf` and
    the written pose SDF have the *identical* bond-order histogram
    (1K3U: 9×single, 3×double, 10×aromatic) and both fail identically, while
    control targets parse in both. It is therefore **not** the pose→SDF writer
    defect reported elsewhere.

    **Why it is fatal only in the harness.** The harness passes
    `-l <crystal_sdf>` (`BustCli.cpp:243-244`), which puts PoseBusters on the
    RMSD path and that path requires a sanitized reference. Measured:
    `bust <pose> -p <apo>` returns **1 row** for all 8 targets, while
    `bust <pose> -l <crystal> -p <apo>` returns **0 rows** on 1K3U and 1 row on
    1JD0. Scored without `-l`, three of the eight (`1K3U`, `1XOZ`, `1Y6R`) pass
    **every physical check**; all eight report `sanitization=False` and
    `inchi_convertible=False`. So the exclusion is **not neutral** — it removes
    targets that would have passed the chemistry.

17. **Do not "fix" this by dropping `-l`.** Without a reference, PoseBusters
    cannot compute RMSD and the result cannot satisfy S2 (§1). The two admissible
    repairs are: (a) repair the cache SDF aromaticity so the reference sanitizes,
    or (b) take RMSD from the engine's own `rmsd_to_crystal` and run PoseBusters
    without `-l`, declaring the split provenance. Either way the 8 targets must
    be reported as a named, counted exclusion until repaired.

### Scope limit on the numbers above

These figures are an **S2-grade PB tier, NOT `claim_ready`** — `claim_ready`
additionally requires tENCoM/Eigen completion, protocol eligibility and
score–pose consistency per §1. `claim_ready` *inherits* the solvent condition
through its PB component and must be recomputed under rules 8 and 14 before it
is quoted. The 16-column rule reproduces the harness's stored `pb_pass` on
154/154 cells, so the harness's own `pb_pass` carries the same looseness.

### Ablation status

The retention policy is **the shipped default and has never been ablated.** It
was never chosen: the `protein` config block was not emitted until commit
`70dd2ed4`, so every key took its parser fallback on every cell this project has
ever run. All three conditions are now reachable
(`FLEXAIDDS_REMOVE_WATER`, `FLEXAIDDS_KEEP_STRUCTURAL_WATERS`,
`FLEXAIDDS_STRUCTURAL_WATER_BFACTOR_MAX`) and default-off is proven byte-identical
at production settings. **Methods text must state the policy** rather than imply
a dry receptor.

---

## 2c. Scored endpoint (NORMATIVE — added 2026-09-04)

**The two candidate endpoints disagree on a majority of cells, so a rate without a
declared endpoint is unreadable.**

| Endpoint | Definition | Status |
|----------|------------|--------|
| `argmin(ACF)` | Cluster free energy over members; what the engine **ships** | Primary — declare explicitly |
| `argmin(CF)` | Lowest contact-function pose in the emitted pool | Secondary / sensitivity |

Measured on 83 rigid NO_SEC cells (seed 12345, ordered-direct RMSD): ACF 40/83 =
48.2%, CF 44/83 = 53.0%. Per-cell endpoint agreement is 34.9% rigid and 39.8%
flexible.

### Hard rules

11. **Every claim table MUST carry an endpoint column**, and every reported rate
    MUST name its endpoint. The two are separate estimands (§1 rule 6).
12. **Banner-to-pose joins key on the restart directory path**, never on the CF
    value — two restarts share a rank-0 CF to 3 d.p. on 22% of cells, making a
    value-join non-deterministic.
13. **Multi-seed aggregation is majority-of-seeds, never union.** A union
    systematically flatters the noisier arm.

    **WHY, so this is not "simplified" back into a union.** A p-value is a
    function of the number of INDEPENDENT observations. Three seeds on one target
    are three looks at the same target, not three targets: they share the
    receptor, the pocket, the ligand and the reference pose. Majority-of-seeds
    collapses them to ONE per-target verdict, so the headline denominator is the
    number of TARGETS (84), never the number of cells (84 x 3 = 252). Reporting
    252 would be pseudo-replication — the test would compute significance for a
    sample size the experiment does not have, and it fails in the flattering
    direction, so nothing in the output flags it.

    Union is worse than merely non-independent: taking the best of three seeds is
    a maximum over noise, which is biased upward by construction.

    **The same rule binds any per-item replication, not just seeds:** restarts
    within a cell, multiple poses of one target, multiple structures of one
    receptor. MEASURED elsewhere in this project: an n=11 two-state result was 11
    pairs drawn from only 5 receptors (3+3+2+2+1), so its reported Wilcoxon p was
    computed over non-independent observations; the honest n was 5. Report the
    count of independent units and name what the unit is.

    **Residual, to be stated as a limitation rather than fixed:** the 84 Astex
    targets are not 84 fully independent draws either — the set contains related
    proteins (several kinases, several proteases). That affects interval width,
    not the point rate, and it must be declared, not corrected away.
| **Fallback MD5** | `9dc93717dfed0698006d88dd6a9627bc` (aggregator default; receipt wins) |

---

## 3. Aggregator CLI

```bash
python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]
python3 scripts/aggregate_claim_metrics.py --c0-full85
python3 scripts/aggregate_claim_metrics.py <dir> --headline strict   # default
python3 scripts/aggregate_claim_metrics.py <dir> --headline s1       # RMSD-only diagnostic
python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only
```

`--headline s3` requires `--diagnostic-only` or exit 2.

---

## 4. PoseBusters receipts (engine)

Upstream `bust` results must preserve **raw CSV before schema failure returns**, and receipt
(`<pdb>_bust_receipt.json` under the per-target PoseBusters sidecar):

- resolved binary path, file SHA-256, version string (if available)
- full argv, exit status, raw CSV SHA-256 / path
- version-pinned mandatory check columns — **every** listed check header required;
  duplicate headers and header/value column-count mismatches fail closed

Shared strict finite PDB coordinate decoder: `LIB/PoseBust/PdbCoords.h`
(used by DatasetRunner RMSD and PoseBust loaders). Topology transfer is
graph-identity only (element sequence + full atom counts including H/Du);
positional remapping of repeated elements is rejected.

## 5. Dual election estimands (DatasetRunner)

When Softβ / entropy election is active, persist **separately**:

| Estimand | Fields |
|----------|--------|
| Generator CF top-1 | `cf_top1_pose_path`, `cf_top1_score`, `cf_top1_rmsd`, `cf_top1_pose_sha256` |
| Entropy/consensus top-1 | `entropy_top1_pose_path`, `entropy_top1_score`, `entropy_top1_rmsd`, `entropy_top1_pose_sha256` |

Elected headline still uses the Softβ/consensus path; CF top-1 never silently
overwrites `seed_echo` / elected provenance.

## 6. Pool ceiling naming

`conditional_scanned_pool_ceiling` (= legacy `best_cluster_rmsd`) is the min
**ordered** RMSD over **actually enumerated** emitted cluster heads and
existing `_member*` files. It is **not** any-pose / full emission census.
Pool ceiling must never clear or mutate `seed_echo` or `pose_source`.
