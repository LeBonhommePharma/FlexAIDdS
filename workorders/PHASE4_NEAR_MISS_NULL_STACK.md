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
