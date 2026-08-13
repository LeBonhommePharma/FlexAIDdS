# E1 determinism campaign — state at stand-down (2026-08-06 23:53 EDT)

Stood down to hand the machine to Claude Science. **E1 completed in full before the
stand-down** — `driver.log` ends with `ALL_DONE Thu Aug  6 23:53:08 EDT 2026`.

## Status: COMPLETE — no cell was interrupted, nothing to discard

All 10 cells have a `DONE` file (`rc=1`) and 2 per-entry JSONs each.
**No partial cell exists, so no `DONE` file needs removing.** A resume would
correctly skip all 10 as complete.

`rc=1` is the benchmark *regression* gate firing on known-below-baseline metrics
(same "REGRESSION DETECTED" the original A/B logged). It is **not** a crash or a
failed run: every cell produced 10 poses per target with real RMSDs.
`TimeoutExpired` count across all cells: **0**.

## Results (best-RMSD per repeat, Astex tier 1)

| cell | 1gpk RMSD | 1gpk wall | 1mq6 RMSD | 1mq6 wall |
|---|---|---|---|---|
| t4_rep1 | 6.180 | 523 s | 8.246 | 850 s |
| t1_rep1 | 6.180 | 538 s | 8.464 | 894 s |
| t4_rep2 | 6.180 | 507 s | 7.977 | 822 s |
| t1_rep2 | 6.180 | 533 s | 8.464 | 990 s |
| t4_rep3 | 6.101 | 624 s | 8.246 | 893 s |
| t1_rep3 | 6.180 | 573 s | 8.464 | 1337 s |
| t4_rep4 | 6.180 | 825 s | 7.877 | 1248 s |
| t1_rep4 | 6.180 | 857 s | 8.464 | 1368 s |
| t4_rep5 | 6.180 | 772 s | 8.464 | 956 s |
| t1_rep5 | 6.180 | 529 s | 8.464 | 880 s |

### Arm summary (n=5)

| arm | values | spread |
|---|---|---|
| OMP=1 1gpk | 6.180 ×5 | **0.000** |
| OMP=1 1mq6 | 8.464 ×5 | **0.000** |
| OMP=4 1gpk | 6.180, 6.180, 6.101, 6.180, 6.180 | **0.079** |
| OMP=4 1mq6 | 8.246, 7.977, 8.246, 7.877, 8.464 | **0.587** |

The OMP=4 1gpk spread (0.079) reproduces the original A/B post_cpu cell exactly.

## Verdict

**Determinism returns completely at 1 thread, on both targets.** Nondeterminism is
thread-gated. This excludes address-dependent ordering (ASLR / allocation address)
as a primary cause, since that would survive at 1 thread.

## Exact pinning used (reproduce with this, do not use main-repo HEAD)

Main repo HEAD moved to `033aac16` mid-campaign (#400/#401 changed tier-1 from a
fixed 2-target set to a seeded random draw of 4). Everything below is pinned to the
A/B-era worktree instead, which is clean at `5f39203`:

- worktree:  `~/Projects/FlexAIDdS/ab_mac_20260806T133329/wt_post_cpu`  (git `5f39203`, clean)
- binary:    `<worktree>/build/FlexAID`
- sha256:    `88edbc85466a0f3803b5e9c0e3f50445e6ee842e68ec060f57c0407d39b78958`
- harness:   `<worktree>/benchmarks/run.py`, `PYTHONPATH=<worktree>/python`
- datasets:  `<worktree>/benchmarks/datasets`  (`tier1_subset_size: 2`, fixed -> 1gpk, 1mq6)
- data root: `FLEXAIDDS_BENCHMARK_DATA=<worktree>/benchmarks/astex_diverse`
  (canonical tree; the nested `benchmarks/astex_diverse/data/` copy is DEPRECATED and
  lacks `*_binding_site.pdb` — using it yields `rmsd = -1.0` sentinels and 1 pose)
- env: `FLEXAID_SEED=12345`, `FLEXAIDDS_PARALLEL_REPRODUCE=0`,
  `FLEXAIDDS_PARALLEL_RESTARTS=0`, `FLEXAIDDS_SEED_ELITISM=0`
- `--workers 1`, `--omp-threads {1,4}`
- **per-target timeout held at 3600 s, unchanged, for every cell.** Not raised.
  Worst observed cell was 1368 s = 38% of ceiling.

## Resume command

    ~/Projects/FlexAIDdS/determinism_e1/drive_e1.sh

Skips any cell with a `DONE` file. As of now that is all 10, so it is a no-op.
To re-run a cell, delete its directory (not just `DONE`).

## Contention log

- 20:44:32–20:46:20 — Claude Science `hdrfix_smoke_20260806_204432` (1N2J/1G9V,
  `--mode autonomous`, main-repo build at `033aac16`). Overlapped **t1_rep2** only.
- Earlier load (~19:20) was system daemons (`fileproviderd` 187%, `contactsd` 83%),
  not Claude Science.
- Wall times after ~21:00 are inflated by contention and must **not** be read as
  performance data. RMSD values are unaffected: the GA has no wall-clock budget, the
  `time(0)` seed fallback at `gaboom.cpp:348-360` is dead when `FLEXAID_SEED` is set,
  and there is no `omp_set_dynamic` so team size is load-independent.

## Outstanding — NOT run

Three interventions, n=3 each at OMP=4, harness ready at
`~/Projects/FlexAIDdS/interventions/drive_interventions.sh`:

1. `FLEXAIDDS_CLEFT_SORT=1`        — isolates cleft-grid merge order
2. `seeding.mif_enabled=false`     — via wrapper `interventions/FlexAID_mifoff.sh`
                                     (+ `interventions/mif_off.json`); decides whether
                                     #372's default flip routes cleft nondeterminism
                                     into gene 0
3. `FLEXAID_DETERMINISTIC=1`       — aggregate control

E5 bisect: deprioritised by agreement (the window's determinism-relevant change is a
one-line default flip in #372, which a bisect would only re-derive).
