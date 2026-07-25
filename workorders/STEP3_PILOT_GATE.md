# STEP 3 W1 pilot gate — FAIL

**Phase:** STEP 3 W1 serial pilot  
**One variable:** `FLEXAIDDS_BOOM_INTERVAL=50` (SIGMA unset; memetic OFF; WALL_PILOT_PASS OFF)  
**OUT:** `/Users/lp.more/flexaidds_results/pilot_w1_boom_interval_20260725_134740`  
**Baseline source:** E10 offline extract of `v_autonomous_20260724_160919` (same campaign; thin result.csv leaf was lost mid-session)  
**Written:** 2026-07-25T21:53:46.665053+00:00  

**Verdict:** **FAIL** — clean-probe regression: [('1N1M', 'rmsd', 3.3796999999999997)]

## Panel metrics (genuine = seed_echo=0 AND rmsd_hungarian < 2)

| Metric | Pilot (BOOM=50) | Baseline (frozen) |
|--------|----------------:|------------------:|
| Genuine top-1 | 1/8 | 0/8 |
| BCR < 2 Å | 3/8 | 3/8 |
| Election gap (BCR<2, not genuine) | 2/8 | 3/8 |
| Mean Δ elected RMSD (pilot−base) | +0.303 | — |
| Mean Δ BCR | +0.479 | — |

## Per-target

| PDB | class | elect RMSD P/B | BCR P/B | genuine P/B | ΔRMSD | ΔBCR |
|-----|-------|---------------:|--------:|:-----------:|------:|-----:|
| 1J3J | clean | 62.22/62.31 | 22.96/21.35 | N/N | -0.09 | +1.60 |
| 1K3U | clean | 11.47/12.70 | 11.78/11.94 | N/N | -1.23 | -0.16 |
| 1L7F | clean | 3.92/4.22 | 3.96/3.98 | N/N | -0.29 | -0.02 |
| 1N1M | clean | 5.66/2.28 | 4.04/4.08 | N/N | +3.38 | -0.03 |
| 1M2Z | clean | 13.79/11.69 | 13.04/11.49 | N/N | +2.10 | +1.55 |
| 1OQ5 | gap | 3.95/3.78 | 1.65/1.06 | N/N | +0.16 | +0.59 |
| 1SQ5 | gap | 5.10/5.10 | 1.65/1.12 | N/N | -0.01 | +0.53 |
| 1YGC | gap | 1.75/3.34 | 1.01/1.24 | Y/N | -1.59 | -0.23 |

## Cadence

- Phase: STEP 3  
- One variable: BOOM_INTERVAL=50  
- Genuine: 1/8; BCR: 3/8; election gap: 2  
- **FAIL**  

Pilot only — not a full-85 claim. STEP 4a memetic still **blocked** by wall efficacy FAIL.
STEP 5 full-85 **not** authorized.
