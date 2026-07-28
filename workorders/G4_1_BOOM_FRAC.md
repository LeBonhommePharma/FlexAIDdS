# G4.1 BOOM_FRAC matched panel — **IN PROGRESS**

**One variable:** `FLEXAIDDS_BOOM_FRAC` unset (control) | 0.05 | 0.1 | 0.2  
**Interval:** claim JSON 100 (unchanged)  
**All arms:** `FLEXAIDDS_NO_SEC=1`, R=**2** (deviation from PHASE4 R=5: ~16 GB free + wall-clock), workers=2, pop=1000, gen=2000, matrix 9dc9  
**Panel:** SEARCH-MISS 1J3J 1K3U 1L7F 1N1M 1M2Z  
**OUT:** `~/flexaidds_results/g4_1_boom_frac_20260726_101238`

## ACCEPT (when complete)

mean ΔBCR (best treatment − control) ≤ **−0.5 Å** **or** ≥1 target BCR&lt;2 under treatment; no elect RMSD regression; L4 `[BOOM]` on treatments only; control zero inject; no CF≈0 wipeout.

## Status

Driver launched 2026-07-26; arm_control running. Metrics filled on completion via `analyze_g4_1.py`.
