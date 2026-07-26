# Matched COARSE A/B (64 vs 256) — FAIL

**Written:** 2026-07-26T03:30:08.864396+00:00  
**One variable:** `FLEXAIDDS_COARSE_ORIENTATIONS` **64 vs 256**  
**OUT:** `/Users/lp.more/flexaidds_results/coarse_ab_matched_20260725_222652`  

**Verdict:** **FAIL** — matched A/B: no directional BCR/RMSD improvement (64 vs 256)

## Matched protocol (skeptic fix)

| Field | Arm A | Arm B |
|-------|------:|------:|
| COARSE_ORIENTATIONS | **64** | **256** |
| Restarts | 2 | 2 |
| Workers | 2 | 2 |
| GA | 1000×2000 | 1000×2000 |
| BOOM_FRAC / INTERVAL | unset | unset |
| Binary sha256 | `7f05640a2a5723a1…` | same |
| Git | `25b2121696fa` | same |
| Matrix md5 | `9dc93717dfed0698006d88dd6a9627bc` | same |

## Liveness L1–L4

| Check | Evidence |
|-------|----------|
| L1 knob read | `config_parser.cpp` `FLEXAIDDS_COARSE_ORIENTATIONS` |
| L2 not stuck at JSON | dock_config still `n_orientations: 64`; engine log shows requested orient count |
| L3/L4 observable | A live-ish 5/5 · B ×256 5/5 |

### L4 log samples
```
A/1J3J: [COARSE-INIT] Scanning 2194 candidate grid points × 64 orientations
B/1J3J: [COARSE-INIT] Scanning 2194 candidate grid points × 256 orientations
A/1K3U: [COARSE-INIT] Scanning 16 candidate grid points × 64 orientations
B/1K3U: [COARSE-INIT] Scanning 16 candidate grid points × 256 orientations
A/1L7F: [COARSE-INIT] Scanning 16 candidate grid points × 64 orientations
B/1L7F: [COARSE-INIT] Scanning 16 candidate grid points × 256 orientations
A/1N1M: [COARSE-INIT] Scanning 28 candidate grid points × 64 orientations
B/1N1M: [COARSE-INIT] Scanning 28 candidate grid points × 256 orientations
A/1M2Z: [COARSE-INIT] Scanning 31 candidate grid points × 64 orientations
B/1M2Z: [COARSE-INIT] Scanning 31 candidate grid points × 256 orientations
```

## M2 triple (panel n=5, genuine metric only)

| Arm | Genuine | BCR&lt;2 | Election gap (BCR ok, not genuine) |
|-----|--------:|--------:|-----------------------------------:|
| A (64) | 0/5 | 0/5 | 0 |
| B (256) | 0/5 | 0/5 | 0 |
| Mean Δ (B−A) elect RMSD | +0.117 | | |
| Mean Δ (B−A) BCR | +3.435 | | |

## Per-target

| PDB | elect A/B | BCR A/B | ΔRMSD | ΔBCR | gen A/B |
|-----|----------:|--------:|------:|-----:|:-------:|
| 1J3J | 62.22/62.22 | 24.62/42.85 | +0.00 | +18.23 | N/N |
| 1K3U | 11.96/12.55 | 12.21/12.01 | +0.60 | -0.20 | N/N |
| 1L7F | 4.32/4.31 | 4.09/3.98 | -0.01 | -0.11 | N/N |
| 1N1M | 5.66/5.66 | 4.55/3.79 | +0.00 | -0.76 | N/N |
| 1M2Z | 13.79/13.79 | 13.04/13.06 | +0.00 | +0.02 | N/N |

## Provenance (criterion-4)

- binary_sha256: `7f05640a2a5723a18fa170ee58072897688a5e7659a29616df3bf68f1bc386ac`
- runner_sha256: `d1c3fb5564c9b2c42d7560a909ab0dd17e8625a387baffa2bac2be6ab99af14f`
- git: `25b2121696fa8548035c9ef4619129bc5d4d13f8`
- matrix_md5: `9dc93717dfed0698006d88dd6a9627bc`
- FLEXAIDDS_* per arm: see `armA_orient64/ARM_ENV.json`, `armB_orient256/ARM_ENV.json`
- dock_config sha256: per-target in MATCHED_AB_GATE.json rows

## Explicit non-claims

- Prior `coarse_orient256_search_miss_*` vs `pilot_w1_boom_*` comparison is **void** (multi-variable).
- Not full-85 · not memetic · not genuine 45% claim.

**Cadence:** Phase 4 matched A/B · one var COARSE_ORIENTATIONS · **FAIL**
