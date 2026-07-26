# pb_clash SCORING-LOCKED oracle (ROADMAP_v2 Phase 2b′) — **FAIL**

**Written:** 2026-07-26  
**One variable:** `FLEXAIDDS_PB_CLASH_WEIGHT` (ladder 1.0 / 5.0 / 10.0; OFF=0)  
**Panel:** SCORING-LOCKED **1OQ5, 1SQ5, 1YGC** (class-matched)  
**Decoy:** each target’s **actual elected pose** from  
`~/flexaidds_results/pilot_w1_boom_interval_20260725_134740/{PDB}/elected_pose.pdb`  
**Config:** pilot `dock_config.json` (production LOCCLF; ops/gates missing these codes)  
**Binary sha256:** `7f05640a2a5723a18fa170ee58072897688a5e7659a29616df3bf68f1bc386ac`  
**Git @ run:** `fb70bf925bbdfae5552fd6556f24c5156acf1f29`  
**Script:** `scripts/pb_clash_burial_oracle.py --mode scoring-locked`

## ACCEPT (magnitude floor — can fail)

dCF decreases by **≥1.0 kcal** AND **flips sign** (OFF dCF&gt;0 → ON dCF&lt;0) on **≥2/3**,  
with **no SEARCH-MISS native CF-min regression** under the ON weight.

## Verdict: **FAIL** (all ladder arms)

| Weight | ge_floor ≥1.0 | sign flip | both | SM regress | OUT |
|-------:|--------------:|----------:|-----:|-----------:|-----|
| 1.0 | 0/3 | 0/3 | 0/3 | 1/5 (1J3J) | `pb_clash_scoring_locked_w1_20260726_000854` |
| 5.0 | 1/3 (1OQ5) | 0/3 | 0/3 | 1/5 (1J3J) | `pb_clash_scoring_locked_w5_20260726_001510` |
| 10.0 | 2/3 (1OQ5,1YGC) | 0/3 | 0/3 | 1/5 (1J3J) | `pb_clash_scoring_locked_w10_20260726_002104` |

### Weight 10.0 detail (strongest arm)

| PDB | dCF_off | dCF_on | decrease | sign_flip | clash_dec ON |
|-----|--------:|-------:|---------:|:---------:|-------------:|
| 1OQ5 | +17.89 | +13.34 | **+4.55** | N | 9.81 |
| 1SQ5 | +28.80 | +28.80 | 0.00 | N | 0.0 |
| 1YGC | +70.20 | +68.67 | **+1.53** | N | 1.53 |

**Interpretation:** pb_clash **can** tax some elected poses (1OQ5 clash rises with weight) but **never flips** the false-min inversion (elected still CF-better than crystal). **1SQ5** elected has **zero** clash under this term — lever is blind to that attractor. **1J3J** SEARCH-MISS recheck shows native not CF-min vs elected under ON weight (dCF_on ≈ +0.96) — regression gate also trips.

## Memetic

**Still locked.** Do **not** set `FLEXAIDDS_PB_CLASH_PHASE2_PASS=1`.  
Gate re-keyed: memetic needs `FLEXAIDDS_MEMETIC=1` **and** (`PB_CLASH_PHASE2_PASS` **or** legacy `WALL_PILOT_PASS`). Claim default remains OFF.

## Class-matched rule (durable)

| Class | Codes | Lever class |
|-------|-------|-------------|
| SEARCH-MISS | 1J3J 1K3U 1L7F 1N1M 1M2Z | **sampling only** |
| SCORING-LOCKED | 1OQ5 1SQ5 1YGC | **scoring only** |

Never judge sampling on SCORING-LOCKED or scoring on SEARCH-MISS.

## Next levers (after this FAIL)

- Stronger / different scoring term that can flip SCORING-LOCKED attractors (not weight thrash on same term)
- Or accept honest SCORING-LOCKED ceiling and focus sampling on SEARCH-MISS only
- Full-85 still blocked
