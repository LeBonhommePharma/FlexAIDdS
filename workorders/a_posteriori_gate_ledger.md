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
| **G4.3 MUTATION_GRANULAR** | **PASS_LIVENESS; null mag** | MUT_GRAN=1 vs control; R=2; mean_dBCR **+0.118** | `g4_3_mut_gran_near_miss_20260727_122215` |

## ELECTION_V135 detail (CLOSED 2026-07-27)

- Protocol ON: `election_v135=true`, `election_score_tau=25` on treatment; control false
- 1N1M elect **6.3999** both arms; 1L7F **3.9907** both; BCR 4.04/3.92 identical
- Floor elect≤2.5 or gap shrink ≥1.0: **FAIL**
- Flip next: **G4.3 mutation** (or new selection architecture)

## G4.3 detail (CLOSED 2026-07-27)

- One var: `FLEXAIDDS_MUTATION_GRANULAR=1` (+ gene_lim ±1-bin mutate; per-gen phenotype-dup clear)
- L4: control **0** / mut_gran **8** `[MUT-GRAN]` → **PASS**
- 1L7F elect 3.9907→**6.2458** (regression); BCR 3.9907→4.1128
- 1N1M elect 6.3999→6.0053; BCR 4.1954→4.3088
- mean_dBCR **+0.118** (null; floor ≤−0.5 or BCR&lt;2 or elect≤2.5)
- status **PASS_LIVENESS**; `ACCEPT_G4_3=False`
- Flip residual: **`new_search_arch`** / scoring-locked decoy work; full-85 still blocked
- Workorder: `G4_3_MUTATION_POSTERIORI.md`

## Full-85 block

Phase-4 sampling ACCEPT still not met (BOOM null + election null + mutation null).

## Publication freeze (offline)

Near-miss Phase-4 null stack table (SI/Methods ready):  
**[`PHASE4_NEAR_MISS_NULL_STACK.md`](PHASE4_NEAR_MISS_NULL_STACK.md)** — freezes G4.1 / ELECTION_V135 / G4.3 statuses, OUTs, and flip residual.

| Residual step | Artifact | Status |
|---------------|----------|--------|
| S1 null-stack freeze | `PHASE4_NEAR_MISS_NULL_STACK.md` | **done** |
| S2 pin schema (`accept.txt` + per-arm SHA) | `benchmark_self_eval.py validate-pins` + contract S2 | **done** |
| S3 SCORING-LOCKED SI | `SCORING_LOCKED_SI_PACKAGE.md` | **done** |
| S4 new_search_arch a priori | `NEW_SEARCH_ARCH_APRIORI.json` / `.md` | **APRIORI_ONLY_NO_LAUNCH** |
| S5 CF-proxy claim language | `CLAIM_LANGUAGE_FREEZE.md` | **done** |
| Full-85 | — | **blocked** |