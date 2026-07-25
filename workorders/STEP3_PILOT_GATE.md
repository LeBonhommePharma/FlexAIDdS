# STEP 3 W1 pilot gate — IN PROGRESS

**Phase:** STEP 3 W1 serial pilot  
**One variable:** `FLEXAIDDS_BOOM_INTERVAL=50` (SIGMA unset; memetic OFF; WALL_PILOT_PASS OFF)  
**Panel:** 1J3J,1K3U,1L7F,1N1M,1M2Z,1OQ5,1SQ5,1YGC  
**Workers:** 2 · OMP=1 · seed OFF · matrix **9dc9** · **restarts=5**  
**OUT:** `/Users/lp.more/flexaidds_results/pilot_w1_boom_interval_20260725_134740`  
**Git at launch:** `05e1fa21` · binary sha256 `a8204fb7…`  
**R5 addendum:** `…/STEP3_PROVENANCE_ADDENDUM.json` (records BOOM=50 from live worker env)

## Status

Detached `benchmark_datasets` PID active. Confirmed in logs:

- `[SEARCH-COVERAGE] boom_interval 100→50` on both 1J3J and 1K3U
- 1K3U finished first restart (pose PDBs written) and started **r1**
- 1J3J still on first restart (large multi-cleft grid)

**Not yet PASS/FAIL** — need all 8 `result.csv` then `aggregate_step3_pilot.py`.

Watcher: `STEP3_WATCHER.log` in OUT; on completion writes this file + `step3_pilot_gate.json`.

## ACCEPT (when complete)

1. No clean-probe regression vs frozen baseline  
2. Directional BCR or elected RMSD effect from BOOM alone  
3. Report genuine / BCR / election-gap on panel only (pilot ≠ claim)

## Next after gate

- **PASS or FAIL** for this knob — do **not** auto-launch full-85  
- STEP 4a memetic still **blocked** (wall un-cap efficacy FAIL)  
- If FAIL: next one-variable W1 knob (e.g. coarse-init / diversity), not dual knobs  
