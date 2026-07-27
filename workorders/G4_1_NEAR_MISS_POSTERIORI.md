# G4.1 BOOM_FRAC near-miss — a posteriori (FINAL)

**OUT:** `/Users/lp.more/flexaidds_results/g4_1_boom_near_miss_20260726_200953`  
**Completed:** 2026-07-27T02:56:32Z (`ALL_ARMS_DONE`)  
**Contract:** `workorders/G4_1_NEAR_MISS_APRIORI.json` + `BENCHMARK_SELF_EVAL_CONTRACT.md`

## L4 BOOM injection
| arm | 1L7F | 1N1M |
|-----|------|------|
| control | 0 | 0 |
| frac005/010/020 | 158 | 78 |

Control zero + treatment live → L4 **PASS**.

## Magnitude (best_cluster_rmsd / elect)
| arm | 1L7F BCR | 1N1M BCR | mean ΔBCR |
|-----|----------|----------|-----------|
| control | 3.9907 | 4.5515 | — |
| frac005 | 4.0834 | 4.5515 | **+0.0464** |
| frac010 | 3.9523 | 4.5515 | **−0.0192** (best) |
| frac020 | 4.0834 | 4.5515 | **+0.0464** |

Floor mean_ΔBCR ≤ −0.5 or BCR&lt;2 → **FAIL (null magnitude)**.

## Decision
- `accept_g4_1 = False`
- flip rule: `G4.1_null_or_l4_fail`
- merged with 1N1M election_offline → **`election_fix_P0`**

## Next (launched)
`election_v135_near_miss_20260726_225823` — `FLEXAIDDS_ELECTION_V135=1` vs control, R=5, codes 1N1M+1L7F, matrix 9dc9, NO_SEC=1.
