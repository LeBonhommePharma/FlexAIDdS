# pb_clash burial oracle — PASS

**Written:** 2026-07-25T23:48:47.218943+00:00  
**One variable:** `FLEXAIDDS_PB_CLASH_WEIGHT=1.0` (OFF arm = 0)  
**Replaces:** WAL_COERCIVE STEP 2 (B3 structural no-op)  
**Configs:** `ops/gates/configs/{PDB}_dock_config.json` (production LOCCLF)  
**Panel:** 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z  
**OUT:** `/Users/lp.more/flexaidds_results/workorders/pb_clash_oracle_20260725_191818`  

**Verdict:** **PASS** — dCF moved toward native on 5/5; no native CF-min regression

## Metrics

| Metric | Value |
|--------|------:|
| Scored | 5 |
| dCF moved toward native (ON vs OFF) | 5/5 |
| Native CF-min under ON | 5/5 |
| Native regressions | 0 |
| OFF≡ON identical | 0/5 |

## Per-target

| PDB | dCF_off | dCF_on | ΔdCF | clash_dec OFF/ON | nat_min OFF/ON | moved |
|-----|--------:|-------:|-----:|-----------------:|:--------------:|:-----:|
| 1J3J | -29.4515 | -29.4529 | -0.0014 | 0.0/0.0014 | Y/Y | Y |
| 1K3U | -120.6470 | -120.6564 | -0.0095 | 0.0/0.0095 | Y/Y | Y |
| 1L7F | -39.5676 | -39.5695 | -0.0019 | 0.0/0.0019 | Y/Y | Y |
| 1N1M | -106.1117 | -106.1345 | -0.0228 | 0.0/0.0228 | Y/Y | Y |
| 1M2Z | -27.4901 | -27.4902 | -0.0001 | 0.0/0.0001 | Y/Y | Y |

## ACCEPT (revised STEP 2)

dCF moves toward/below 0 with weight ON vs OFF on ≥4/5; no clean native CF-min regression.

## Cadence

- Phase: STEP 2 replacement (pb_clash burial)  
- One variable: PB_CLASH_WEIGHT=1.0  
- **PASS**  
- Do **not** set WALL_PILOT_PASS from WAL_COERCIVE evidence  
- Memetic still blocked until a burial/steric oracle PASSes

## Scientific caveat (effect size)

Formal ACCEPT is met (5/5 dCF moved toward native; 0 native regressions; OFF≢ON).

However **|ΔdCF| is ≪ 1 CF unit** on every target and decoy `cf_clash` stays near 0 even after COM-burial translation. Deeper burial inflates `cf_total` via wall/com terms, not PoseBusters-ratio `pb_clash`. Treat this as:

1. **Env liveness + directional sign check for `FLEXAIDDS_PB_CLASH_WEIGHT`** — PASS  
2. **Not** a strong burial-opponent demonstration for re-keying memetic / `WALL_PILOT_PASS` without a stronger deep-interpenetration decoy panel  

Do **not** set `FLEXAIDDS_WALL_PILOT_PASS=1` from this micro-effect alone. Memetic remains blocked pending product re-key of the interlock and a panel with non-trivial `cf_clash`.

Live OUT: `/Users/lp.more/flexaidds_results/workorders/pb_clash_oracle_20260725_191818`  
Script: `scripts/pb_clash_burial_oracle.py`
