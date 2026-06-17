# Astex Diverse Set Tier-1 Report

## Execution
- Date: 2026-06-07T10:10:10Z
- Targets: 1sq5 2hb1 1r1h 1t46 2c69
- Passed: 0, Failed: 0, Skipped: 5

## Baselines (Hartshorn et al. 2007 J Med Chem 50:726-741)

| Metric | Baseline | Tolerance | Measured | Status |
|:-------|:--------:|:---------:|:--------:|:------:|
| docking_power_top1 | 0.70 | ±0.05 | — | pending |
| docking_power_top3 | 0.85 | ±0.05 | — | pending |
| mean_rmsd | 2.30 Å | ±0.05 | — | pending |
| entropy_rescue_rate | 0.30 | ±0.05 | — | pending |

> Metrics are populated by DatasetRunner once pose output is available.
> Full 85-complex validation: `benchmark_datasets --benchmark astex`
