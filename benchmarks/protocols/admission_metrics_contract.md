# Admission + Metrics Contract (Normative)

**Status:** Normative for claim-table aggregation and abstract / headline rates.  
**Aligned with:** `benchmarks/protocols/three_engine_entropy_comparison.md` §1.4–§5, `AGENTS.md`, audit 2026-07-17.  
**Enforcement:** `scripts/aggregate_claim_metrics.py` (fail-closed).

---

## 1. Success metrics (report all; headline is claim_ready)

| ID | Definition | Role |
|----|------------|------|
| **S1** | Elected pose **ordered direct** `rmsd_to_crystal` ≤ 2.0 Å and `seed_echo==0` | RMSD-only diagnostic of election |
| **S2** | S1 ∧ `pb_pass` on the **same** elected pose (`success_pb`) | Intermediate (RMSD∧PB) |
| **STRICT** / **claim_ready** | Engine `claim_ready==1`: S2 ∧ official PoseBusters with version/config receipt ∧ tENCoM/Eigen complete on exact pose SHA-256 ∧ no seed/oracle contamination ∧ protocol eligibility ∧ score–pose consistency | **Primary / claim KPI** |
| **S3** | `conditional_scanned_pool_ceiling` / `best_cluster_rmsd` ≤ 2.0 Å | **Diagnostic only** (scanned heads/members, **not** any-pose) |

### Field mapping (DatasetRunner `result.csv`)

| Metric | Required / preferred fields |
|--------|------------------------------|
| Elected RMSD (S1) | **`rmsd_to_crystal` only** (ordered direct, whole-ligand). **Never** `rmsd_hungarian` for S1/admission |
| Hungarian | `rmsd_hungarian` — diagnostic symmetry metric only |
| PoseBusters | `pb_pass`, `success_pb`, `pb_backend==bust_cli` for claims |
| Strict claim | **`claim_ready==1`** plus complete same-pose hashes, `pb_backend==bust_cli`, `pb_pass==1`, upstream PB version/config receipt, `tencom_status==ok`, and `eigen_status==ok` |
| Pool ceiling (S3) | `conditional_scanned_pool_ceiling` → else `best_cluster_rmsd` (scanned emission pool only) |

### Hard rules

1. **Headline claim rate = claim_ready only** (not S1, not S3).
2. **S1 uses ordered direct RMSD only** — never Hungarian for claim S1.
3. **Never report S3 as abstract / headline success.** Label as diagnostic conditional scanned-pool ceiling.
4. **Always report S1, S2, STRICT (claim_ready), and S3 separately.**
5. **Election gap** = S3=1 and S1=0 (scanned pool had ≤2 Å; elector missed). Diagnostic only.
6. Generator **CF top-1** and **entropy/consensus reranked top-1** are separate estimands (paths, scores, hashes, RMSDs) when both are present.

---

## 2. Claim denominator and strict numerator

The claim denominator is fixed before outcomes are observed. Every row with an
explicit preregistered `protocol_claim_eligible==1` contributes to `N_claim`,
including docking failures, RMSD failures, PoseBusters failures, tENCoM/Eigen
failures, missing receipts, and `claim_ready==0`.

Outcome and receipt fields never remove a preregistered target from that
denominator. They decide only whether the target enters the STRICT numerator:

| STRICT numerator check | Required value |
|------------------------|----------------|
| `protocol_claim_eligible` | **1** (explicit preregistration) |
| `seed_echo` / `native_pose_seeded` | **0** (explicit; missing fails closed) |
| `matrix_md5` | equals campaign matrix pin when present |
| `claim_ready` | **1** |
| Ordered RMSD / PB | S1 true, `pb_pass==1`, `pb_backend==bust_cli` |
| Same-pose hashes | `pose_sha256`, `rmsd_pose_sha256`, `posebusters_pose_sha256`, and `tencom_pose_sha256` are all nonempty and identical |
| Upstream PB receipt | receipt exists and has nonempty binary path/hash, version, full argv, raw CSV path/hash, `backend==bust_cli`, `pb_pass==1`, and exit 0 |
| Validator status | `tencom_status==ok`, `eigen_status==ok` |

Rows missing explicit `protocol_claim_eligible` are legacy/incomplete rows. An
observed legacy row remains in `N_claim` (so mixed schema can never turn one
success plus one incomplete row into 1/1), is reported under
`legacy_diagnostics`, and cannot enter the STRICT numerator. The strict
headline is marked suppressed until the schema is complete.

The preregistered target-ID list is loaded from `--expected-targets` (JSON,
CSV, or newline/comma text) or `expected_target_ids` / `target_ids` / `targets`
in `RUN_RECEIPT.json`. Expected IDs with no result row are synthesized as
strict failures in the denominator. Missing IDs, duplicate result IDs,
duplicate manifest IDs, and unexpected IDs are exposed under `completeness`
and suppress the strict headline; a partial row set is never labelled a full
campaign result.

### Default matrix pin

| Field | Value |
|-------|--------|
| Canonical matrix | campaign-dependent (see RUN_RECEIPT) |
| **Fallback MD5** | `9dc93717dfed0698006d88dd6a9627bc` (aggregator default; receipt wins) |

---

## 3. Aggregator CLI

```bash
python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]
python3 scripts/aggregate_claim_metrics.py <campaign_dir> --expected-targets expected_ids.txt
python3 scripts/aggregate_claim_metrics.py --c0-full85
python3 scripts/aggregate_claim_metrics.py <dir> --headline strict   # default
python3 scripts/aggregate_claim_metrics.py <dir> --headline s1       # RMSD-only diagnostic
python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only
```

`--headline s3` requires `--diagnostic-only` or exit 2.

---

## 4. PoseBusters receipts (engine)

Upstream `bust` results must preserve **raw CSV before schema failure returns**.
The nested `<pdb>_bust_receipt.json` is pinned to:

- schema `posebusters-0.6.5-redock-csv-v1`, exactly 27 required checks
- package `posebusters==0.6.5` and launcher output `bust 0.6.5`
- config `redock` with SHA-256
  `4d551d898ff29a404f16e02ad5a7a2d4235e6b7b14e9a3e27f7c66b4d16b2da9`
- nonempty path plus 64-hex SHA-256 for package RECORD, launcher, config,
  three inputs, raw CSV, and validated CSV
- nonempty argv, exit 0, `backend==bust_cli`, `ran==true`, `pb_pass==true`

When receipt output files are still accessible, the aggregator re-hashes them
and fails STRICT on a mismatch. Duplicate headers and header/value
column-count mismatches fail closed in the engine parser.

Shared strict finite PDB coordinate decoder: `LIB/PoseBust/PdbCoords.h`
(used by DatasetRunner RMSD and PoseBust loaders). Topology transfer is
graph-identity only (element sequence + full atom counts including H/Du);
positional remapping of repeated elements is rejected.

## 5. Election estimands (DatasetRunner)

Persist each election boundary **separately**, without reconstructing it from
the final elected pose. Generator CF and SoftBeta/entropy are frozen from the
same eligible GA-emitted cluster-head census; `_INI`, FREQSEL, seed elitism,
and later consensus reranking cannot overwrite either snapshot:

| Estimand | Fields |
|----------|--------|
| Generator CF top-1 | `cf_top1_pose_path`, `cf_top1_score`, `cf_top1_rmsd`, `cf_top1_pose_sha256` |
| SoftBeta/entropy top-1 | `entropy_top1_pose_path`, `entropy_top1_score`, `entropy_top1_rmsd`, `entropy_top1_pose_sha256` |
| Final consensus top-1 | `consensus_top1_pose_path`, `consensus_top1_score`, `consensus_top1_rmsd`, `consensus_top1_pose_sha256` |

The elected headline uses the final post-reranker/seed/consensus path. CF and
entropy top-1 never silently inherit that state or overwrite `seed_echo` /
elected provenance.

## 6. Pool ceiling naming

`conditional_scanned_pool_ceiling` (= legacy `best_cluster_rmsd`) is the min
**ordered** RMSD over **actually enumerated** emitted cluster heads and
existing `_member*` files. It is **not** any-pose / full emission census.
Pool ceiling must never clear or mutate `seed_echo` or `pose_source`.
