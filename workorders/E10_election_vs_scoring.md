# E10 — election vs scoring / sampling (STEP 1)

**Date:** 2026-07-25  
**Phase:** W0.1 / E10 offline (no re-dock)  
**Script:** `scripts/e10_election_vs_scoring.py`  
**Campaign:** `~/flexaidds_results/v_autonomous_20260724_160919`  
**Materialize:** named rsync + thin result.csv copy → **85** `result.csv`; rank-0 PDBs still partial (~10) so size-bias uses REMARK when present.

## Results (full result.csv set)

| Metric | Value |
|--------|------:|
| Targets analyzed | **85** |
| seed_echo=0 | **85** |
| Elected RMSD < 2.0 (proxy genuine) | **21/85 = 24.7%** |
| BCR < 2.0 | **24/85 = 28.2%** |
| Election-gap (BCR≤2.5, elected>2.0, seed_echo=0) | **16/85 = 18.8%** |
| Size-bias suspects (REMARK soft_β when PDBs present) | **10** (only targets with local pose PDBs) |

Reference baseline (OPS doc): genuine **20/79=25.3%**, BCR **22/79=27.8%**. E10 elected-rmsd rate on full 85 may differ if seed_echo / sentinel columns differ — use for **gap structure**, not to re-derive 25.3% without `aggregate_claim_metrics`.

## Interpretation

- Election-gap fraction **18.8%** — material but **not** "election is the whole wall"; many failures have BCR≫2 (sampling).
- Size-bias on materialized heads confirms **legacy ACF multiplicity** on this frozen run (soft_β_G ≪ CF with high freq). **Do not** cite this campaign as post-`free_energy_strict` proof.
- **Independent scorer prefers near-native head over elected:** uncommon in notes; election gap often = BCR near-native pool not electing rank-0, not necessarily a better-CF head.

## ACCEPT vs STEP 1

| Criterion | Result |
|-----------|--------|
| E10 on local leaf | PASS (N=85) |
| Large election-dominated failure? | **NO** (~19% gap; sampling still primary for BCR≪2 failures) |
| STOP before sampling? | **NO** — continue; wall STEP 2 already FAIL on production panel |

## Cadence

- Phase: STEP 1 E10  
- One variable: n/a  
- Genuine proxy: 21/85; BCR: 24/85; election gap: 16  
- **PASS** (continue to wall redesign / sampling diagnosis)

Machine: implementer `e10/e10.{json,csv,md}`
