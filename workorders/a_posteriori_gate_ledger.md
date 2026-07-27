# A posteriori gate ledger (canonical)

Status enum: **PASS / FAIL / VOID / INVALID / MISSING_OUT / IN_FLIGHT**.  
Source of floors: `METHODOLOGY.md`, `workorders/PHASE4_GATES_ACTUALIZED.md`, `workorders/BENCHMARK_SELF_EVAL_CONTRACT.md`.  
Do not average wrong-class targets; no silent rewrites.

| Gate | Status | One variable / check | Evidence |
|------|--------|----------------------|----------|
| E10 + M2 triple | **PASS** | offline election vs scoring split | `workorders/E10_election_vs_scoring.md` |
| Native–Elected CF inversion | **PASS** | pose role fixed LOCCLF; SEARCH-MISS=5 / SCORING-LOCKED=3 | `workorders/INVERSION_MAP.md` |
| WAL_COERCIVE | **VOID** | structural no-op (B3) | `workorders/WALL_ORACLE.md` |
| pb_clash SEARCH-MISS (legacy) | **VOID** | wrong panel / ROADMAP_v2 | `workorders/PB_CLASH_ORACLE.md` |
| pb_clash SCORING-LOCKED 2b′ | **FAIL** | weight 1/5/10; 0 sign flips | `workorders/PB_CLASH_SCORING_LOCKED.md` |
| Burial lever | **VOID/retired** | empty weight window | CAMPAIGN_GATE_SUMMARY |
| COARSE 64 vs 256 matched | **FAIL** | COARSE_ORIENTATIONS only; mean ΔBCR +3.44 | `MATCHED_AB_GATE.md` + `~/flexaidds_results/coarse_orient256_search_miss_*` |
| G4.2 niche Cartesian σ=2 | **FAIL** | NICHE_CARTESIAN; mean ΔBCR −0.441 (misses −0.5) | `G4_2_NICHE_CART.md` |
| G4.2 R=5 near-miss | **FAIL/null** | R=5 on 1N1M/1L7F | `~/flexaidds_results/g4_2_r5_near_miss_20260726_175237` |
| DUMP_POP full-pop ceiling | **MEASURED** | 0 sub-2 in full pop; 1N1M pop_best≈2.36 elect≈6.41 | `~/flexaidds_results/dump_pop_search_miss_20260726_172356` |
| 1J3J disease label | **SCORING_PULL** | CF prefers decoy vs crystal | pins in `tests/test_campaign_methodology_gates.py` |
| **G4.1 BOOM near-miss** | **FAIL (null magnitude); L4 PASS** | BOOM_FRAC {0.05,0.1,0.2} vs control; 1N1M+1L7F; R=2 | `~/flexaidds_results/g4_1_boom_near_miss_20260726_200953` + `workorders/G4_1_NEAR_MISS_POSTERIORI.md` |
| **ELECTION_V135** | **IN_FLIGHT** | `FLEXAIDDS_ELECTION_V135=1` vs control; R=5 | `~/flexaidds_results/election_v135_near_miss_20260726_225823` |

## G4.1 detail (CLOSED 2026-07-27)

- L4: control **0** `[BOOM]`; treatments **236** each (stderr-aware scanner `70ed4f51`)
- Magnitude: best mean_dBCR **−0.0192** (frac010); floor −0.5 or BCR&lt;2 → **null**
- Matrix 9dc9, NO_SEC=1, binary `a3fa78c1…`
- Flip: **election_fix_P0** (1N1M offline pop 2.36 / elect 6.41)
- Enforcer: `scripts/benchmark_self_eval.py` + `scripts/campaign_flip_order.py`

## Flip order (current)

1. election_fix_P0 ← **active (ELECTION_V135)**
2. G4.1_BOOM secondary (exhausted for magnitude)
3. G4.3 mutation later

## Full-85 block

Phase-4 sampling ACCEPT not yet (BOOM null; election pending).  
Residual path: `publication_residual_path` in SCRATCH / prior workorder notes.
