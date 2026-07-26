# COARSE_ORIENTATIONS=256 W1 pilot — FAIL

**One variable:** COARSE_ORIENTATIONS=256  

**Verdict:** **FAIL** — no directional BCR/RMSD improvement vs baseline pilot

| Metric | Pilot 256 | Baseline pilot |
|--------|----------:|---------------:|
| Genuine | 0/5 | 0/5 |
| BCR<2 | 0/5 | 0/5 |
| Mean Δ elect RMSD | +0.294 | — |
| Mean Δ BCR | +3.980 | — |
| Liveness 256 | 5/5 | — |

| PDB | elect P/B | BCR P/B | ΔRMSD | ΔBCR | gen P/B |
|-----|----------:|--------:|------:|-----:|:-------:|
| 1J3J | 62.22/62.22 | 42.85/22.96 | +0.00 | +19.90 | N/N |
| 1K3U | 12.55/11.47 | 12.01/11.78 | +1.08 | +0.23 | N/N |
| 1L7F | 4.31/3.92 | 3.98/3.96 | +0.38 | +0.01 | N/N |
| 1N1M | 5.66/5.66 | 3.79/4.04 | +0.00 | -0.25 | N/N |
| 1M2Z | 13.79/13.79 | 13.06/13.04 | +0.00 | +0.02 | N/N |

Pilot only — not full-85 claim.

