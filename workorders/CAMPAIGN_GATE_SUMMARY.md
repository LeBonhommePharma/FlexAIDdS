# Campaign methodology follow-through — gate summary

**Source:** `docs/implementation/CAMPAIGN_METHODOLOGY_for_Grok.md`  
**OPS correction:** `~/flexaidds_results/workorders/WALL_ORACLE_FAIL_EXPLAINED.md`  
**Audit:** `workorders/DOCKING_BUG_AUDIT_2026-07-25.md` (if present) / B1–B3  
**Updated:** 2026-07-25

## Critical invalidations (instrumentation — not docking-quality fails)

| Gate | Label | Why |
|------|--------|-----|
| **STEP 2 WAL_COERCIVE** | **Structurally unpassable** | **B3:** cap is per-pair and binds only for *o* > 1 Å; Voronoi wall loop never enumerates deep interpenetration (~23× undercount). OFF≡ON is expected forever. **Not** a missing-env bug. |
| **STEP 3 BOOM_INTERVAL=50** | **Scientifically invalid as a BOOM test** | **B1:** inject requires `interval>0 && fraction>0` (`gaboom.cpp:986`); claim path emits `boom_inject_fraction: 0.0` deliberately (`DatasetRunner.cpp:6058`). Interval-only → zero `[BOOM]`. **Not** a docking success/fail on BOOM. |

Neither invalidation authorizes full-85, memetic, or `WALL_PILOT_PASS=1`.

## Gate table

| Step | Phase | One variable | Result | PASS/FAIL |
|------|-------|--------------|--------|-----------|
| 0 | Merge methodology tooling on main | n/a | methodology + workorders on main | **PASS** |
| 1 | E10 offline | n/a | N=85; election-gap ~18.8%; sampling primary | **PASS** (continue) |
| 2a | Wall WAL_COERCIVE (prod LOCCLF) | `WAL_COERCIVE` | production CF OK (1M2Z=−117.74); OFF≡ON | **STRUCTURAL FAIL** (unpassable) |
| 2b | Saturating wal panel | same | OFF≡ON even when Σ cf_wal≥45 | confirms B3 |
| **2′** | **pb_clash burial oracle (replacement)** | `FLEXAIDDS_PB_CLASH_WEIGHT=1.0` | 5/5 dCF toward native; 0 regressions; **micro |ΔdCF|** — see caveat in PB_CLASH_ORACLE | **PASS formal**; **not** memetic unlock |
| 3 | W1 BOOM_INTERVAL=50 only | interval only | zero inject; clean 1N1M RMSD noise | **INVALID as BOOM** (instrumentation) |
| 3′ | Small-frac BOOM liveness A/B | `BOOM_FRAC=0.1` | JSON still 0; B has `[BOOM]` n_inject=50/1000; A silent; no CF≈0 wipe | **PASS (liveness)** — not a success-rate claim |
| 4a | Memetic | — | Blocked until burial/steric oracle PASSes | **NOT RUN** |
| 5 | Full-85 claim | — | Blocked until W1/W3 gates allow | **NOT RUN** |

## BOOM product caveat

Claim `boom_inject_fraction: 0.0` is a **deliberate** anti-collapse fix (frac=1.0 wiped blind GA every 100 gens → CF≈0 @ ~300). Env `FLEXAIDDS_BOOM_FRAC` **does** override JSON (verified). Use **0.05–0.2 only**; never 1.0 on blind path.

## Explicit blocks

- No dual full-85; WORKERS≤4; OMP=1/worker  
- No memetic / no `FLEXAIDDS_WALL_PILOT_PASS=1` from WAL-only evidence  
- No re-run of WAL_COERCIVE expecting OFF≠ON  
- No treating STEP 3 interval-only pilot as BOOM efficacy  
- Matrix **9dc9** (`md5 9dc93717dfed0698006d88dd6a9627bc`) for dock pilots; score-only oracles record binary sha + env  

## Artifacts

| Artifact | Path |
|----------|------|
| E10 | `workorders/E10_election_vs_scoring.md` |
| Wall (WAL) | `workorders/WALL_ORACLE.md` |
| STEP 3 invalid pilot | `workorders/STEP3_PILOT_GATE.md` |
| BOOM liveness | `workorders/BOOM_FRAC_LIVENESS.md` |
| BOOM A/B | `workorders/BOOM_FRAC_AB.md` |
| pb_clash STEP2′ | `workorders/PB_CLASH_ORACLE.md` |
| Script | `scripts/pb_clash_burial_oracle.py` |

## Next allowed

**Primary (recommended):** score-only **Native–Elected CF inversion map** — see [`NEXT_CAMPAIGN_STEP.md`](NEXT_CAMPAIGN_STEP.md). Classifies SCORING-LOCKED vs SEARCH-MISS before another GA pilot.

Then (one variable each, after map):
1. SEARCH-MISS heavy → `FLEXAIDDS_COARSE_ORIENTATIONS=256` W1 pilot (matrix 9dc9).  
2. SCORING-LOCKED heavy → stronger deep-interpenetration decoys for pb_clash (micro formal PASS is not memetic unlock).  
3. Full-85 only after remaining cheap gates pass.  
4. Product decision: re-key memetic interlock from WAL to a **strong** burial oracle — do not auto-enable.
