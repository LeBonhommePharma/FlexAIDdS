# Wall coercive oracle — STEP 2

## Skeptic fix: production configs

`scripts/wall_coercive_oracle.py` uses **`ops/gates/configs/{PDB}_dock_config.json`** via `panel_manifest.tsv` (not `diagnostic/probe_config.json`).

**Spot-check:** 1M2Z native CF = **−117.74** (matches methodology production LOCCLF).

## A. Production falsemin decoys (arm-A false minima)

| Metric | Value |
|--------|------:|
| n | 5 |
| native wins OFF | 4/5 |
| native wins ON | 4/5 |
| OFF≡ON | yes |
| **VERDICT** | **FAIL** (cannot demonstrate wall un-cap; 2 targets still inverted) |

## B. Saturating burial decoys (redesign)

Falsemin translated toward receptor COM until `cf_wal` ≥ 45 (`scripts/wall_saturating_panel.py`).

| Metric | Value |
|--------|------:|
| n | 5 |
| native wins OFF | 5/5 |
| native wins ON | 5/5 |
| rescued | 0 |
| OFF≡ON identical CF | **yes (all 5)** |
| Scoring competitiveness (≥5/5 native CF-min vs buried decoy) | **PASS** |
| **WAL_COERCIVE efficacy** (OFF≠ON or rescued>0) | **FAIL** |

| PDB | dCF_off | dCF_on | identical |
|-----|--------:|-------:|:---------:|
| 1J3J | -29.451502 | -29.451502 | True |
| 1K3U | -17.693798 | -17.693798 | True |
| 1L7F | -39.567561 | -39.567561 | True |
| 1N1M | -106.111734 | -106.111734 | True |
| 1M2Z | -27.49011 | -27.49011 | True |

## Cadence / gate consequence

- **Do NOT** set `FLEXAIDDS_WALL_PILOT_PASS=1` (WAL_COERCIVE does not change CF even when `cf_wal`>CAP — probe path may clamp before env, or cap is not the active limiter on summed wal).
- **Memetic (STEP 4a) blocked.**
- **W1 non-memetic sampling knobs (BOOM/coarse) allowed** for STEP 3 under methodology (STOP is for memetic after wall fail of un-cap efficacy).
- Production CF policy: **PASS** (configs correct).

## Files

- Falsemin prod run: `~/flexaidds_results/workorders/wall_oracle_prod_*`
- Saturating panel: `~/flexaidds_results/workorders/wall_sat_panel_*`
- Saturating A/B: `~/flexaidds_results/workorders/wall_oracle_sat_*`
