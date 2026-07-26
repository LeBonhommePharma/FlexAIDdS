# Phase 4 W1 pilot — COARSE_ORIENTATIONS=256 (SEARCH-MISS only)

**Status:** LAUNCHED (2026-07-25)  
**OUT:** `~/flexaidds_results/coarse_orient256_search_miss_20260725_214924`  
**One variable:** `FLEXAIDDS_COARSE_ORIENTATIONS=256` (JSON still emits 64; env overrides)

## Liveness (L1–L4) — PRECHECK PASS

Smoke: `~/flexaidds_results/coarse_orient_liveness_20260725_214424`

| Check | Evidence |
|-------|----------|
| L1 read | `LIB/config_parser.cpp` `FLEXAIDDS_COARSE_ORIENTATIONS` |
| L2 not stuck at JSON | dock_config `n_orientations: 64` but log `× 256 orientations` |
| L3 acts | `[COARSE-INIT] Injected 25 pre-screened seeds into gen-0` |
| L4 log | `Scanning … × 256 orientations` |

## Panel (from inversion map SEARCH-MISS)

`1J3J 1K3U 1L7F 1N1M 1M2Z` — **exclude** SCORING-LOCKED `1OQ5 1SQ5 1YGC`

## Protocol

- autonomous · R=2 · workers=2 · OMP=1 · GA 1000×2000 · seed OFF  
- matrix **9dc9** · BOOM_FRAC unset  

## ACCEPT (when complete)

1. Log shows 256 orientations / coarse inject on each target  
2. No elect-RMSD regression vs `pilot_w1_boom_interval_20260725_134740` on same codes  
3. Directional BCR or elect RMSD improvement  
4. Panel M2-style counts only (pilot ≠ claim)

## Blocks

No full-85 · no memetic · no WAL re-panel · SCORING-LOCKED not in this pilot  
