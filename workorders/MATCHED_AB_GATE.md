# Matched COARSE 64 vs 256 A/B — gate

**Written:** 2026-07-26T03:40:13.703171+00:00
**Root:** `/Users/lp.more/flexaidds_results/coarse_ab_matched_20260725_222652`
**One variable:** `FLEXAIDDS_COARSE_ORIENTATIONS` 64 (A) vs 256 (B)
**Verdict:** **FAIL**

## Skeptic fix

Prior COARSE=256 vs pilot_w1 was **VOID** (multi-var: different binary/git/R/BOOM).
This matched control holds binary, git, R=2, BOOM unset on both arms.

## Provenance (criterion-4)

| Field | Value |
|-------|-------|
| binary_sha256 | `7f05640a2a5723a18fa170ee58072897688a5e7659a29616df3bf68f1bc386ac` |
| git | `25b2121696fa8548035c9ef4619129bc5d4d13f8` |
| branch | `main` |
| matrix_md5 | `9dc93717dfed0698006d88dd6a9627bc` |
| restarts | 2 |
| workers | 2 · OMP=1 |
| ga | pop=1000 gen=2000 |
| BOOM_FRAC/INTERVAL | unset both arms |
| seed_elitism | 0 |
| panel | 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z (SEARCH-MISS only) |

## L1–L4 liveness

| Arm | expect n | L4 log hits matching | sample |
|-----|---------:|---------------------:|--------|
| A64 | 64 | 5 | `[COARSE-INIT] Scanning 2194 candidate grid points × 64 orientations` |
| B256 | 256 | 5 | `[COARSE-INIT] Scanning 2194 candidate grid points × 256 orientations` |

- L1 knob read: A=True B=True
- L2 env overrides JSON n_orientations=64 (claim default): documented + observed via log n
- L3 physically able: orientations 64/256 supported
- L4 log-observable: both arms need `[COARSE-INIT] … N orientations`

## M2 triple

| Arm | Genuine | BCR&lt;2 | Election gap |
|-----|--------:|--------:|-------------:|
| A64 | 0/5 | 0/5 | 0/5 |
| B256 | 0/5 | 0/5 | 0/5 |

Mean Δ elect RMSD (B−A): **+0.117** Å  
Mean Δ BCR (B−A): **+3.435** Å

## Per-target

| PDB | elect A/B | BCR A/B | ΔRMSD | ΔBCR | gen A/B |
|-----|----------:|--------:|------:|-----:|:-------:|
| 1J3J | 62.22/62.22 | 24.62/42.85 | +0.00 | +18.23 | N/N |
| 1K3U | 11.96/12.55 | 12.21/12.01 | +0.60 | -0.20 | N/N |
| 1L7F | 4.32/4.31 | 4.09/3.98 | -0.01 | -0.11 | N/N |
| 1N1M | 5.66/5.66 | 4.55/3.79 | +0.00 | -0.76 | N/N |
| 1M2Z | 13.79/13.79 | 13.04/13.06 | +0.00 | +0.02 | N/N |

## ACCEPT / FAIL

| Check | Result |
|-------|--------|
| Liveness L4 both arms | True |
| Directional BCR or elect gain | False |
| No genuine regression | True |
| One variable + BOOM unset | True |
| **Verdict** | **FAIL** |

## Blocks still in force

- Full-85 not unlocked by this pilot alone
- No memetic / WALL_PILOT_PASS
- Prior multi-var COARSE pilot remains VOID for attribution

Artifacts: `/Users/lp.more/flexaidds_results/coarse_ab_matched_20260725_222652/MATCHED_AB_GATE.json` · SCRATCH/MATCHED_AB_GATE.json · pilot_metrics.md

