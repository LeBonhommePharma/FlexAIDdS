# A posteriori gate ledger (canonical)

Status enum: **PASS / FAIL / VOID / INVALID / MISSING_OUT / IN_FLIGHT**.  
Source of floors: `METHODOLOGY.md`, `workorders/PHASE4_GATES_ACTUALIZED.md`, `workorders/BENCHMARK_SELF_EVAL_CONTRACT.md`.

| Gate | Status | One variable / check | Evidence |
|------|--------|----------------------|----------|
| E10 + M2 triple | **PASS** | offline election vs scoring | `E10_election_vs_scoring.md` |
| Native–Elected CF inversion | **PASS** | SEARCH-MISS=5 / SCORING-LOCKED=3 | `INVERSION_MAP.md` |
| WAL_COERCIVE | **VOID** | structural no-op | `WALL_ORACLE.md` |
| pb_clash SEARCH-MISS | **VOID** | wrong panel | `PB_CLASH_ORACLE.md` |
| pb_clash SCORING-LOCKED 2b′ | **FAIL** | 0 sign flips | `PB_CLASH_SCORING_LOCKED.md` |
| Burial | **VOID/retired** | empty weight window | CAMPAIGN_GATE_SUMMARY |
| COARSE 64 vs 256 | **FAIL** | mean ΔBCR +3.44 | MATCHED_AB_GATE |
| G4.2 niche σ=2 | **FAIL** | mean ΔBCR −0.441 | G4_2_NICHE_CART |
| G4.2 R=5 near-miss | **FAIL/null** | R=5 | `g4_2_r5_near_miss_20260726_175237` |
| DUMP_POP ceiling | **MEASURED** | 1N1M pop≈2.36 elect≈6.41 | `dump_pop_search_miss_20260726_172356` |
| G4.1 BOOM near-miss | **FAIL (null mag); L4 PASS** | BOOM_FRAC 0.05/0.1/0.2 | `g4_1_boom_near_miss_20260726_200953` |
| **ELECTION_V135** | **FAIL (null mag)** | V135=1 vs control; R=5; elect Δ=0 | `election_v135_near_miss_20260726_225823` |

## ELECTION_V135 detail (CLOSED 2026-07-27)

- Protocol ON: `election_v135=true`, `election_score_tau=25` on treatment; control false
- 1N1M elect **6.3999** both arms; 1L7F **3.9907** both; BCR 4.04/3.92 identical
- Floor elect≤2.5 or gap shrink ≥1.0: **FAIL**
- Flip next: **G4.3 mutation** (or new selection architecture)

## Full-85 block

Phase-4 sampling ACCEPT still not met (BOOM null + election null).
