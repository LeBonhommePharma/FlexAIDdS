# Validated reference config: COM floor + P1 steer (2026-08-01)

**Status:** validated on the RIGID target 1gpk; does NOT yet generalize to the flexible
target 1mq6 (see #356). Promoted as the reference config for the COM-floor docking regime,
NOT (yet) as a global engine default.

## The config

```
FLEXAIDDS_COM_FLOOR        = 130     # soft floor on CF.com (vcfunction.cpp), S=5 fixed
FLEXAIDDS_PB_POCKET_WEIGHT = 1.0     # P1 GA-time pocket steer (top.cpp)
FLEXAIDDS_PB_CLASH_WEIGHT  = 1.0     # P1 GA-time clash steer  (top.cpp)
# determinism (as swept):
FLEXAID_SEED = 42 · OMP_NUM_THREADS = 1 · FLEXAIDDS_PARALLEL_RESTARTS = 0 · FLEXAIDDS_SEED_ELITISM = 0
```

## Evidence (corrected in-place RMSD metric, PR #354; clean build @ b352ca54)

Tier-1 astex_diverse sweep, 2 targets (1gpk, 1mq6), run `30715200428`:

```
cell                docking_power_top1
baseline (OFF)          0.000
P1-only (floor off)     0.000     <- P1 alone insufficient
F = 130                 0.500     <- docks 1gpk (rank-1 1.60 A); 1mq6 misses
F = 20000 .. 100000     0.000     <- above threshold, floor reverts to baseline
```

Per target at F=130:

```
1gpk (0 rot bonds):  rank-1 1.60 A, best-of-10 0.77 A  -> PASS (search + score both work)
1mq6 (10 rot bonds): rank-1 11.03 A, best-of-10 ~4.5 A -> FAIL at SEARCH (see #356)
```

## Findings

1. The COM floor at low F makes the engine dock a rigid ligand end-to-end — the score ranks
   the near-native pose #1. First green docking_power of the effort.
2. The floor is necessary (P1 alone = 0) and sharply F-thresholded (works at 130, reverts by 2e4).
3. Flexible-ligand docking fails at the search stage, not scoring (#356).
4. `mean_rmsd` is byte-identical (2.278) across all cells — the runner's INI-pollution bug;
   `docking_power` is the trustworthy metric (it moved 0->0.5). Separate follow-up.

## Scope caveat for promotion

`docking_power = 0.5` is 1 of 2 targets, below the 0.70 Astex bar. Do NOT flip the global
engine getenv defaults to this config: it is validated only on rigid ligands and regresses
flexible-ligand docking. Promote as the documented reference for the rigid regime and the
Tier-1 investigation; gate/engine default changes wait on #356.
