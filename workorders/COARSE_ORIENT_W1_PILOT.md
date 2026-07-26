# COARSE_ORIENTATIONS matched A/B — authoritative gate

**Supersedes** multi-var FAIL note that compared COARSE=256 to pilot_w1_boom_interval (different binary/git/R/BOOM).

**One variable:** `FLEXAIDDS_COARSE_ORIENTATIONS` = 64 vs 256  
**Verdict:** **FAIL**  
**Root:** `/Users/lp.more/flexaidds_results/coarse_ab_matched_20260725_222652`  
**Written:** 2026-07-26T03:40:13.703171+00:00

## Provenance (criterion-4)

| Field | Value |
|-------|-------|
| binary_sha256 | `7f05640a2a5723a18fa170ee58072897688a5e7659a29616df3bf68f1bc386ac` |
| git tip | `25b2121696fa8548035c9ef4619129bc5d4d13f8` |
| matrix_md5 | `9dc93717dfed0698006d88dd6a9627bc` |
| R / workers / OMP | 2 / 2 / 1 |
| ga | pop=1000 gen=2000 |
| FLEXAIDDS_BOOM_* | **unset** both arms |
| FLEXAIDDS_COARSE_ORIENTATIONS | A=64 · B=256 |
| seed_elitism | 0 |
| panel | SEARCH-MISS: 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z |

## L1–L4

| Check | A64 | B256 |
|-------|-----|------|
| L1 knob read | env + claim JSON coarse_init | env override |
| L2 not wrongly overridden | JSON stays n_orientations=64; env forces arm | env=256 |
| L3 physically able | 64 supported | 256 supported |
| L4 log-observable | 5/5 `[COARSE-INIT] … × 64 orientations` | 5/5 `… × 256 orientations` |

Dock config hash (1J3J Arm A): `31eeb168d978aedc226d920ef84476bf76835751e88e8bc1b3fec9982735456d`  
`ga.boom_inject_fraction: 0` in claim JSON; BOOM env unset both arms.

Full table: `workorders/MATCHED_AB_GATE.md`.

## M2 triple

| Arm | Genuine | BCR&lt;2 | Election gap |
|-----|--------:|--------:|-------------:|
| A64 | 0/5 | 0/5 | 0/5 |
| B256 | 0/5 | 0/5 | 0/5 |

Mean Δ elect (B−A): **+0.117** Å · Mean Δ BCR (B−A): **+3.435** Å (worse; 1J3J BCR 24.6→42.9)

| PDB | elect A/B | BCR A/B | ΔRMSD | ΔBCR | gen A/B |
|-----|----------:|--------:|------:|-----:|:-------:|
| 1J3J | 62.22/62.22 | 24.62/42.85 | +0.00 | +18.23 | N/N |
| 1K3U | 11.96/12.55 | 12.21/12.01 | +0.60 | -0.20 | N/N |
| 1L7F | 4.32/4.31 | 4.09/3.98 | -0.01 | -0.11 | N/N |
| 1N1M | 5.66/5.66 | 4.55/3.79 | +0.00 | -0.76 | N/N |
| 1M2Z | 13.79/13.79 | 13.04/13.06 | +0.00 | +0.02 | N/N |

## Prior VOID comparison (do not cite for science)

`coarse_orient256_search_miss_20260725_214924` vs `pilot_w1_boom_interval_20260725_134740` mixed binary/git/R/BOOM — **not one-variable**.

## Cadence

- Phase 4.1 matched control (skeptic fix) → **FAIL** honest
- Full-85 still blocked
- No memetic unlock
- Next one-var levers: niche Cartesian (code+flag), `NO_SEC` budget honesty, or strong burial decoys for SCORING-LOCKED
