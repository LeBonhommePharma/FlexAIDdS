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
| **Fallback MD5** | `72d7c7396702331d96ff12d18f831796` (aggregator default; receipt wins) |

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
