# G4.2 Cartesian niche matched A/B — **FAIL**

**Written:** 2026-07-26T07:21:03.720271+00:00  
**Root:** `/Users/lp.more/flexaidds_results/g4_2_niche_cart_ab_20260726_004752`  
**One variable:** `FLEXAIDDS_NICHE_CARTESIAN` OFF (A) vs ON (B); B `SIGMA_ANG=2.0`  
**Both:** `FLEXAIDDS_NO_SEC=1` (G4.4 full budget)

## L1–L4

| Arm | [NICHE-CART] enabled hits | NO_SEC |
|-----|--------------------------:|--------|
| A gene RMSP | 0 (must be 0) | yes |
| B Cartesian | 5 | yes |

L4 sample B: ['1J3J: enabled', '1J3J: has lines', '1K3U: enabled']

## M2

| Arm | Genuine | BCR&lt;2 |
|-----|--------:|--------:|
| A | 0/5 | 0/5 |
| B | 0/5 | 0/5 |

Mean ΔBCR (B−A): **-0.441** Å  
Mean Δ elect (B−A): **-0.251** Å  
Elect regressions: 0/5

## Per-target

| PDB | elect A/B | BCR A/B | ΔRMSD | ΔBCR | gen A/B |
|-----|----------:|--------:|------:|-----:|:-------:|
| 1J3J | 62.22/62.22 | 24.62/23.64 | +0.00 | -0.98 | N/N |
| 1K3U | 11.95/11.60 | 12.20/11.92 | -0.35 | -0.27 | N/N |
| 1L7F | 4.32/3.49 | 4.09/3.54 | -0.83 | -0.54 | N/N |
| 1N1M | 5.66/5.66 | 4.55/4.05 | +0.00 | -0.50 | N/N |
| 1M2Z | 13.79/13.71 | 13.04/13.13 | -0.07 | +0.09 | N/N |

## ACCEPT

| Check | Result |
|-------|--------|
| L4 cart on B only | True |
| Magnitude (mean dBCR≤−0.5 or BCR&lt;2 on B) | False |
| No elect regression | True |
| Niche log live | True |
| **Verdict** | **FAIL** |

## Cadence

- Phase 4 G4.2 SEARCH-MISS only  
- Full-85 still blocked  
- Memetic still locked (option a)  
- G4.1 BOOM_FRAC next if G4.2 fails without blocking  

## Follow-on

- **G4.1** BOOM_FRAC panel not run in this session (G4.2 took full-budget wall-clock; next).
- **G4.3** mutation granularity remains separate arm, not bundled.
- mean ΔBCR **−0.441 Å** narrowly misses −0.5 floor → honest **FAIL** (directional improvement without magnitude PASS).
