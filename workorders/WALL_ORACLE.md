# Wall coercive oracle (score-only, production configs)

Binary: `/Users/lp.more/Projects/FlexAIDdS/build/FlexAIDdS`
probe_cf: `/Users/lp.more/Projects/FlexAIDdS/build/probe_cf`
Manifest/config policy: **ops/gates/configs/{PDB}_dock_config.json** (not diagnostic/probe_config.json)
Panel scored ok: **5**
Native wins OFF: **4/5**
Native wins ON: **4/5**
Rescued (fail OFF → pass ON): **0**
Need ≥ ceil(7n/8) = **5**
Regressed already-min: **0**
**VERDICT: FAIL**

| PDB | config | cf_nat_off | cf_dec_off | dCF_off | dCF_on | win_off | win_on |
|-----|--------|-----------:|-----------:|-------:|------:|:-------:|:------:|
| 1J3J | `1J3J_dock_config.json` | -47.932918 | -47.66979 | -0.2631280000000018 | -0.2631280000000018 | True | True |
| 1K3U | `1K3U_dock_config.json` | -139.958755 | -140.635856 | 0.6771009999999933 | 0.6771009999999933 | False | False |
| 1L7F | `1L7F_dock_config.json` | -107.245341 | -105.54911 | -1.6962309999999974 | -1.6962309999999974 | True | True |
| 1N1M | `1N1M_dock_config.json` | -51.4955 | -51.461836 | -0.03366400000000169 | -0.03366400000000169 | True | True |
| 1M2Z | `1M2Z_dock_config.json` | -117.744882 | -117.136094 | -0.6087880000000041 | -0.6087880000000041 | True | True |

## Cadence

- Phase: STEP 2 wall oracle
- One variable: FLEXAIDDS_WAL_COERCIVE
- PASS/FAIL: **FAIL**

**STOP before memetic / WALL_PILOT_PASS.** Re-diagnose wall / panel.

## Spot-check (methodology)

| Check | Value | Expected |
|-------|------:|----------|
| 1M2Z native cf_total (OFF) | **-117.7449** | ~**-117.74** (production LOCCLF) |
| Non-production diagnostic CF | n/a | ~-187 was wrong |

Production CF match: **YES**

## Gate consequence

**VERDICT:** **FAIL**
Steps 3–5 blocked if FAIL. Do not set WALL_PILOT_PASS.
