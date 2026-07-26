# pb_clash burial oracle — **VOID** (ROADMAP_v2)

**Written:** 2026-07-25 (original) · **VOID labeled:** 2026-07-26  
**One variable:** `FLEXAIDDS_PB_CLASH_WEIGHT=1.0` (OFF arm = 0)  
**Panel (wrong):** SEARCH-MISS clean probes 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z  
**OUT:** `~/flexaidds_results/workorders/pb_clash_oracle_20260725_191818`

## Verdict: **VOID** — not a scientific PASS

Originally reported **PASS** (dCF moved toward native 5/5). **ROADMAP_v2** voids this:

| Flaw | Evidence |
|------|----------|
| **No magnitude floor** | ΔdCF 0.0001–0.0228 kcal on CF 27–138 (~1e-4 relative) = **noise** |
| **Wrong panel** | All 5 are **SEARCH-MISS** (inversion map): crystal already CF-min; scoring lever has nothing to fix |
| **Decoys not buried** | min lig–rec heavy dist 2.20–2.96 Å; **zero pairs &lt;2.0 Å**; pb_clash almost nothing to penalize |

### Historical micro-table (do not cite as steric PASS)

| PDB | dCF_off | dCF_on | ΔdCF |
|-----|--------:|-------:|-----:|
| 1J3J | -29.4515 | -29.4529 | -0.0014 |
| 1K3U | -120.6470 | -120.6564 | -0.0095 |
| 1L7F | -39.5676 | -39.5695 | -0.0019 |
| 1N1M | -106.1117 | -106.1345 | -0.0228 |
| 1M2Z | -27.4901 | -27.4902 | -0.0001 |

**MEMETIC_UNLOCK:** **false**  
**WALL_PILOT_PASS:** do not set  
**FLEXAIDDS_PB_CLASH_PHASE2_PASS:** do not set from this artifact

## Successor

See **`workorders/PB_CLASH_SCORING_LOCKED.md`** — revised Phase 2 on **1OQ5 / 1SQ5 / 1YGC** with elected decoys and magnitude-floor ACCEPT (≥1.0 kcal decrease + sign flip ≥2/3).
