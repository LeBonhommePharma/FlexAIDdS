# BOOM_FRAC env liveness — verified (pre A/B)

**Date:** 2026-07-25  
**Context:** Audit B1 + OPS caveat that claim path hardcodes `boom_inject_fraction: 0.0` deliberately (anti-collapse at frac=1.0 blind).

## Does env win over JSON?

**Yes.** Order in `LIB/config_parser.cpp`:

1. JSON `ga.boom_inject_fraction` (DatasetRunner emits **0.0** at `DatasetRunner.cpp:6058`)
2. `FLEXAIDDS_BOOM_FRAC` via `getenv` → sets `GB->boom_inject_fraction`
3. `ProtocolConfig.boom_frac` (also from `FLEXAIDDS_BOOM_FRAC` in `ProtocolConfig.cpp:116`) → same value again

JSON does **not** win. Env is required for any blind-mode BOOM test.

## Smoke evidence

**OUT:** `~/flexaidds_results/boom_liveness_smoke_20260725_183333`  
**Env:** `FLEXAIDDS_BOOM_FRAC=0.1` `FLEXAIDDS_BOOM_INTERVAL=50` `FLEXAIDDS_RESTARTS=1`  
**GA:** 150 gens × 1000 pop · 1N1M · autonomous

| Check | Result |
|-------|--------|
| `dock_config.json` `boom_inject_fraction` | **0** (claim hardcode) |
| `dock_config.json` `boom_inject_interval` | **100** (JSON; interval env applied in-engine only) |
| `[BOOM]` lines | **2** (gen 50, gen 100) |
| n_inject | **50/1000** = 0.1 × (pop/2) → **frac=0.1 live** |
| gens timed | **150/150** (full smoke budget) |

```
[BOOM] injection #1 at gen 50: re-randomized worst 50/1000 chromosomes ...
[BOOM] injection #2 at gen 100: re-randomized worst 50/1000 chromosomes ...
```

## STEP 3 invalidation (confirmed)

Prior W1 pilot set only `FLEXAIDDS_BOOM_INTERVAL=50` with claim JSON frac=0.0 → **zero** `[BOOM]` lines. Gate measured seed/early-stop noise, not BOOM.

## Design rules for BOOM experiments

1. **Always** assert liveness: JSON frac==0 AND `[BOOM]` present with n_inject matching intended fraction.
2. Use **small** fraction (0.05–0.2). **Never** 1.0 on blind/claim path (documented wipe every 100 gens → CF≈0 @ gen~300).
3. Acceptance includes **generation reach / no CF≈0 wipe**, not only RMSD.
4. One variable: e.g. FRAC only (leave interval at 100) unless interval is the single knob under test **and** frac>0.

## Follow-on A/B (in flight)

`~/flexaidds_results/boom_frac_ab_*`  
- Arm A: control (no BOOM_FRAC)  
- Arm B: `FLEXAIDDS_BOOM_FRAC=0.1` only · R=2 · gen=2000 · 1N1M  

See `workorders/BOOM_FRAC_AB.md` when complete.
