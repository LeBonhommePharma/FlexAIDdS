# Wall coercive oracle — STEP 2 (score-only)

**Date:** 2026-07-25
**Phase:** W2 wall oracle
**One variable:** `FLEXAIDDS_WAL_COERCIVE` OFF vs ON
**Panel:** 1J3J 1K3U 1L7F 1N1M 1M2Z (1G9V excluded)
**probe_cf / binary:** build/probe_cf · build/FlexAIDdS
**Raw dir:** `/Users/lp.more/flexaidds_results/workorders/wall_oracle_20260725_125933`

## Results

| PDB | dCF_off | native_wins_off | dCF_on | native_wins_on |
|-----|--------:|:---------------:|-------:|:--------------:|
| 1J3J | -0.4474440000000044 | True | -0.4474440000000044 | True |
| 1K3U | 0.38928700000002436 | False | 0.38928700000002436 | False |
| 1L7F | -1.2441829999999925 | True | -1.2441829999999925 | True |
| 1N1M | 0.08167899999999406 | False | 0.08167899999999406 | False |
| 1M2Z | -0.25623999999999114 | True | -0.25623999999999114 | True |

## Summary

| Metric | Value |
|--------|------:|
| n scored | 5 |
| native wins OFF | 3/5 |
| native wins ON | 3/5 |
| need (≥ ceil(7n/8)) | 5 |
| regressed already-min natives | 0 |
| **VERDICT** | **FAIL** |

## Cadence

- **Phase:** STEP 2 wall oracle
- **One variable:** FLEXAIDDS_WAL_COERCIVE
- **Genuine / BCR / gap:** n/a (score-only)
- **PASS/FAIL:** **FAIL**

## Gate consequence

**STOP before memetic / W3 sampling that depends on wall un-cap.** Re-diagnose wall term; do not set FLEXAIDDS_WALL_PILOT_PASS=1.

## Follow-up (rebuild retest)

Rebuilt `probe_cf` and re-scored 1K3U decoy OFF vs ON: **identical** `cf_total` and `cf_wal=32.148` (below `WAL_CONTACT_CAP=50`). So on these falsemin poses the wall is **not saturating**; un-capping cannot change dCF. Flag is wired in `vcfunction.cpp` but **ineffective on this panel**. Methodology STOP stands: do not enable memetic / WALL_PILOT_PASS until a panel with saturating wall contacts is designed or wall physics is re-diagnosed.

