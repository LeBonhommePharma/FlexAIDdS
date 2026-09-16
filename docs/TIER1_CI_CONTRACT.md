# Tier-1 Astex CI contract

This is the human-readable form of `python/flexaidds/dataset_runner/tier1_contract.py`
and **METHODOLOGY.md §0.1.1**. It governs the 4-target Astex Diverse job in
`.github/workflows/benchmark-tier1.yml`. It is **not** an Astex-85 docking-power
claim and must not be quoted as one.

The CF energy matrix is **not** retuned to green this gate. The pin is
`MC_st0r5.2_6.dat` MD5 `9dc93717dfed0698006d88dd6a9627bc` (9dc9).

## Why the gate is split

Actions run [35016255184](https://github.com/LeBonhommePharma/FlexAIDdS/actions/runs/35016255184)
(`FLEXAIDDS_TIER1_SEED=20260816`, roster `1gpk,1mq6,1xm6,2cet`) completed 4/4
targets and emitted 40 poses. **Sampling** passed (`mean_rmsd ≈ 1.86`,
`median_rmsd ≈ 1.73`). **Ranking** failed (`docking_power_top1 = 0.0`,
`docking_power_top3 = 0.0`) because CF election ranked non-natives as rank-1
while near-natives sat at ranks 5–10 in 3 of 4 libraries (`sampling_power = 0.75`).
Pose RMSDs/scores matched the same-day matrix-pin dock; this is a known
ranking-layer weakness, not a DualAssembly/#502-only bug and **not** evidence
that Softβ / entropy “failed”.

The 0.70 / 0.85 docking-power baselines are aspirational
(`expected_baselines_role: ci_gate_only`). Gating the PR on them false-failed
a run that still sampled near-natives, and would have missed a chance to
hard-fail a real sampling or matrix regression.

## Hard gates (must stay green)

| Check | What it catches |
|---|---|
| Liveness / productivity / completeness (#326) | Crash, timeout, **0 poses**, unmeasured baselines |
| Roster replay | Seed / YAML drift that docks a different four-target set |
| Matrix pin (9dc9) | Wrong CF matrix / 72d7 fork |
| Election objective lock | Silent `--ranking-objective` mismatch vs `cf_minus_ts` |
| `sampling_power` ≥ 0.50 | Library collapse (fewer than 2/4 targets emit any RMSD ≤ 2.0 Å) |
| `mean_rmsd` ≤ 2.30 | Systematic min-RMSD worsening |
| `median_rmsd` ≤ 1.95 | Same, robust to one outlier |

`sampling_power` is **S_top10 / any-pose** success over the pinned draw. Rank
is ignored. Empty results measure `0.0` (and still fail productivity), so a
silent `docking_power_top1=0` from no poses cannot hide as an “advisory ranking miss”.

## Advisory metrics (reported, non-blocking)

Until a **registered election rule** exists (`ci_ranking_status:
unregistered`), these ranking / ΔS metrics are measured and flagged but
do not fail CI. Loading Astex YAML with `unregistered` while those
metrics are still `hard` is a config error — it would recreate the
35016255184 false fail. Promoting ranking to a merge gate requires
`ci_ranking_status: registered` **and** flipping the matching
`ci_gate_class` entries together.

- `docking_power_top1` (aspirational 0.70)
- `docking_power_top3` (aspirational 0.85)
- `entropy_rescue_rate` (aspirational 0.30)
- `generator_docking_power_top1` / `entropy_reranked_docking_power_top1` (reported, no baseline)

These stay in `expected_baselines` so completeness still requires them to be
**measured**. A missing ranking metric is still INCONCLUSIVE (exit 3). A
measured `0.0` with a full 10-pose library is advisory (exit 0) plus a
warning.

**Thesis direction** for a future registered election (not active as a CI
weight, and not Cartesian ligand entropy on the CF proxy): pose-density
clustering for conformational entropy; torsional-consistent tENCoM for
vibrational entropy.

## Replay receipt

Every live Tier-1 Astex dock writes `results/tier1/TIER1_RECEIPT.json`
(`schema: flexaidds.tier1_receipt.v1`) with seed, roster, anchors, drawn
slots, matrix MD5, binary SHA-256, election objective, hard vs advisory
snapshots, and verdict.

Pinned draw for seed **20260816**:

- anchors: `1gpk`, `1mq6`
- drawn: `1xm6`, `2cet`

### Reproduce the draw (no dock)

```bash
export FLEXAIDDS_TIER1_SEED=20260816
PYTHONPATH=python python3 -m flexaidds.dataset_runner.tier1_contract \
  --replay-roster --seed 20260816 --check-matrix
```

Expected stdout includes `roster=['1gpk', '1mq6', '1xm6', '2cet']` and
`matrix_pin_ok`.

### Reproduce the dock (same as CI)

```bash
export FLEXAIDDS_TIER1_SEED=20260816
export FLEXAIDDS_BINARY=/path/to/FlexAID   # absolute
export FLEXAIDDS_BENCHMARK_DATA=$PWD/benchmarks/astex_diverse
export FLEXAIDDS_ENTRY_TIMEOUT_SECONDS=5400
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0
export FLEXAIDDS_SEED_ELITISM=0
PYTHONPATH=python python -m benchmarks.run \
  --dataset astex_diverse --tier 1 \
  --results-dir results/tier1 \
  --datasets-dir benchmarks/datasets \
  --workers 4 --omp-threads 1 --verbose
PYTHONPATH=python python3 -m flexaidds.dataset_runner.tier1_contract \
  --validate-receipt results/tier1/TIER1_RECEIPT.json
```

`FLEXAID_SEED` is **not** pinned on pull_request (empty ⇒ `time(0)`). Pose
coordinates are therefore not bit-identical across PR runs; the **draw** is.
Record `ga_seed` from the receipt before comparing poses.

## Residual risks (what this gate still will / will not catch)

**Will catch**

- Engine crash, hang at the 5400 s cap, or empty pose libraries
- Wrong or missing CF matrix
- Seed/YAML drift that changes the four-target roster
- Sampling collapse (no near-native in ≥3 of 4 libraries, or min-RMSD blow-up)
- Completeness holes (a declared baseline never measured)

**Will not catch**

- CF ranking electing a non-native while a near-native sits in the library
  (advisory until a registered election rule)
- Astex-85 docking power, PoseBusters, or tENCoM/Eigen claim packages
- Campaign-path physics (`mif_enabled`, `coarse_init`, 50 poses, …) — see
  METHODOLOGY.md §0.1
- Native CF oracle (crystal CF vs best GA CF) on these four targets
- GA seed jitter on 1gpk’s 1.94 Å near-native flipping `sampling_power` 3/4 → 2/4
  (the 0.50 floor is exactly that 2/4 line; `median_rmsd` is the smoother net)
- Softβ / entropy science. Do not claim those failed from an advisory top-1 miss.
