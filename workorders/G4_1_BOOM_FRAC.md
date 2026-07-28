# G4.1 BOOM_FRAC matched panel — **FAIL (null magnitude) / PASS_LIVENESS**

**One variable:** `FLEXAIDDS_BOOM_FRAC` unset (control) | 0.05 | 0.1 | 0.2  
**Interval:** claim JSON 100 (unchanged)  
**All arms:** `FLEXAIDDS_NO_SEC=1`, R=**2** then R=**5** near-miss follow-up, workers=2, pop=1000, gen=2000, matrix 9dc9  
**Panel:** SEARCH-MISS 1J3J 1K3U 1L7F 1N1M 1M2Z (near-miss focus 1L7F+1N1M)  
**OUT (near-miss final):** `~/flexaidds_results/g4_1_boom_near_miss_20260726_200953`

## ACCEPT (a priori)

mean ΔBCR (best treatment − control) ≤ **−0.5 Å** **or** ≥1 target BCR&lt;2 under treatment; no elect RMSD regression; L4 `[BOOM]` on treatments only; control zero inject; no CF≈0 wipeout.

## Status — CLOSED

| gate | result |
|------|--------|
| L4 BOOM liveness | **PASS** — control 0 markers; treatments 236 each (stderr.log path) |
| Magnitude floor | **FAIL (null)** — best mean ΔBCR ≈ −0.019 Å (frac010) |
| `accept_g4_1` | **False** |

Primary a-posteriori: `workorders/G4_1_NEAR_MISS_POSTERIORI.md`  
Evidence: `workorders/g4_1_evidence/`  
Flip residual → `election_fix_P0` / G4.3 (see `workorders/a_posteriori_gate_ledger.md`).

**Note:** Prior 1N1M-only liveness is **not** claimed as G4.1 science success; full SEARCH-MISS panel + magnitude contract was required and failed magnitude.
