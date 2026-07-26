# G4.2 Cartesian niche metric — matched A/B **IN PROGRESS**

**One variable:** `FLEXAIDDS_NICHE_CARTESIAN` OFF (A) vs ON (B)  
**B only:** `FLEXAIDDS_NICHE_SIGMA_ANG=2.0`  
**Both arms (G4.4 fix):** `FLEXAIDDS_NO_SEC=1` (full generation budget)  
**Panel:** SEARCH-MISS 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z  
**Protocol:** R=2, workers=2, OMP=1, pop=1000, gen=2000, BOOM unset, COARSE unset, matrix 9dc9  
**OUT:** `~/flexaidds_results/g4_2_niche_cart_ab_20260726_004752`

## L1–L4 (instrumentation)

| Check | Evidence |
|-------|----------|
| L1 | `gaboom.cpp` reads `FLEXAIDDS_NICHE_CARTESIAN` / `FLEXAIDDS_NICHE_SIGMA_ANG` |
| L2 | env-OFF default; claim JSON does not force cart niche |
| L3 | Cartesian path precomputes ligand heavy coords via `calc_rmsd_chrom`, pair RMSD in Å |
| L4 | B arm logs `[NICHE-CART] enabled` + periodic `n_lonely` / `mean_pshare` |

## ACCEPT (sampling magnitude floor)

- mean ΔBCR (B−A) ≤ **−0.5 Å** **or** ≥1 target crosses BCR&lt;2 Å under B  
- no SEARCH-MISS elected RMSD regression (B elect worse than A)  
- L4 cart fire on B; A must **not** print `[NICHE-CART] enabled`  
- n_lonely / niche diversity proxy must move (liveness of metric change)

## Status

Driver launched 2026-07-26; Arm A (gene RMSP) running. Metrics table and PASS/FAIL filled when both arms complete.

## Memetic

Option **(a)** — still locked; not tied to this lever.
