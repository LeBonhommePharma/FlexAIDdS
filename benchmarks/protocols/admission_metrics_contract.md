# Admission + Metrics Contract (Normative)

**Status:** Normative for claim-table aggregation and abstract / headline rates.  
**Aligned with:** `benchmarks/protocols/three_engine_entropy_comparison.md` §1.4–§5, `AGENTS.md` scientific guardrails, `benchmarks/BENCHMARK_STANDARD.md` (no-seed, seed_echo).  
**Enforcement:** `scripts/aggregate_claim_metrics.py` (fail-closed claim filters; S3 never primary).

---

## 1. Success metrics (report all; headline uses S1)

| ID | Definition | Role |
|----|------------|------|
| **S1** | Elected (top-1) pose RMSD ≤ 2.0 Å | **Primary / claim KPI** |
| **S2** | S1 ∧ PoseBusters pass on elected pose | Modern secondary |
| **S3** | Any emitted cluster pose RMSD ≤ 2.0 Å (BCR / sampling ceiling) | **Diagnostic only** |

### Field mapping (DatasetRunner `result.csv`)

| Metric | Preferred fields (first hit) |
|--------|------------------------------|
| Elected RMSD | `rmsd_hungarian` → else `rmsd_to_crystal` |
| S1 flag (optional) | `success_s1` → else recompute from elected RMSD ≤ 2.0 and `seed_echo==0` → else `success_rmsd` |
| PoseBusters | `pb_pass` / `success_pb` (`success_pb` should equal S1 ∧ `pb_pass`) |
| S3 / BCR | `best_cluster_rmsd` ≤ 2.0 Å |

**Hungarian preferred** for S1 when present (symmetry-corrected elected pose vs crystal).  
`rmsd_to_crystal` is serial-order and is a fallback / diagnostic, not the preferred claim RMSD.

### Hard rules

1. **Never report S3 as abstract / headline success.** Always label S3 “diagnostic (BCR / any-pose ceiling).”
2. **Always report S1, S2, and S3 separately** when aggregating claim rows.
3. **Election gap** = targets with S3=1 and S1=0 (sampling found a near-native pose the elector missed). Report counts; do not fold into S1.

---

## 2. Claim admission (row must pass all)

A row is **claim-eligible** only if:

| Check | Required value |
|-------|----------------|
| `seed_echo` | **0** |
| `native_pose_seeded` | **0** |
| `matrix_md5` | equals campaign **matrix pin** |

Optional column `protocol_claim_eligible` (engine) is treated as an additional gate when present: claim rows require it true **or** missing (missing → recompute from seed flags). Rows with `protocol_claim_eligible=0` are excluded from claim aggregates even if seed flags are clean, to honour runner metadata.

### Default matrix pin

| Field | Value |
|-------|--------|
| Canonical matrix | `MC_st0r5.2_6.dat` |
| **Default MD5 pin** | `72d7c7396702331d96ff12d18f831796` |

Override pin sources (first available):

1. CLI `--matrix-md5`
2. Campaign `RUN_RECEIPT.json` → `matrix_md5`
3. Campaign `provenance.json` → `matrix_md5`
4. Default pin above

Per-row `matrix_md5` (if present) must match the pin. Campaign-level pin is applied when rows omit the column.

**Seeded / oracle-ceiling campaigns** (native inheritance) are a **separate science track**. Do not mix their rates into three-engine or TIER-1 claim tables. Use `scripts/aggregate_oracle_ceiling.py` for that track only.

---

## 3. Aggregator CLI (enforceable)

```bash
# Default claim aggregation (S1 primary)
python3 scripts/aggregate_claim_metrics.py <campaign_dir> [--json out.json]

# C0 full85 via env (after source ~/.flexaidds_env)
python3 scripts/aggregate_claim_metrics.py --c0-full85

# Explicit pin / receipt
python3 scripts/aggregate_claim_metrics.py <campaign_dir> --matrix-md5 72d7c7396702331d96ff12d18f831796

# FAIL: S3 as headline without diagnostic-only
python3 scripts/aggregate_claim_metrics.py <dir> --headline s3          # exit 2
python3 scripts/aggregate_claim_metrics.py <dir> --headline s3 --diagnostic-only  # OK, still labels diagnostic
```

Exit codes:

| Code | Meaning |
|------|---------|
| 0 | Claim rows present; aggregation OK |
| 1 | No claim-eligible rows (or empty campaign) |
| 2 | Usage / contract violation (`--headline s3` without `--diagnostic-only`, bad path) |

---

## 4. C0 full85 path

After `source ~/.flexaidds_env`:

```text
$FLEXAIDDS_RESULTS/campaigns/C0_full85_defined_cleft_nativeseed_forbidden
```

Expect per-target `<PDB>/result.csv`, plus `RUN_RECEIPT.json` / `provenance.json` with `matrix_md5`.

---

## 5. What this contract does **not** change

- Pose ranking, clustering, election, or docking engine code.
- GA search or CF/contact-function scoring.
- Seeded oracle-ceiling science (separate aggregator).

---

## 6. Cross-references

| Doc | Role |
|-----|------|
| `benchmarks/protocols/three_engine_entropy_comparison.md` | Multi-arm protocol, schema, hypotheses |
| `benchmarks/BENCHMARK_STANDARD.md` | Tiers, seed_echo, BCR computation notes |
| `AGENTS.md` | Scientific guardrails (proxy vs thermo; no silent ranking change) |
| `CLAUDE_BENCHMARK_HANDOFF.md` | Live C0 full85 ops handoff |
| `scripts/aggregate_claim_metrics.py` | Enforcement implementation |
| `scripts/aggregate_oracle_ceiling.py` | Seeded ceiling track only |
