# Campaign methodology follow-through — gate summary

**Source:** `docs/implementation/CAMPAIGN_METHODOLOGY_for_Grok.md`  
**Main tip:** post-STEP3 docs  
**Updated:** 2026-07-25

| Step | Phase | One variable | Result | PASS/FAIL |
|------|-------|--------------|--------|-----------|
| 0 | Merge methodology + wall/E10 tooling on main | n/a | main pushed | **PASS** |
| 1 | E10 offline (frozen archive) | n/a | N=**85**; elected proxy 21/85; BCR 24/85; election-gap **16/85=18.8%** | **PASS** (continue) |
| 2a | Wall WAL_COERCIVE (prod LOCCLF configs) | `FLEXAIDDS_WAL_COERCIVE` | 1M2Z native CF=**−117.74**; falsemin 4/5 OFF≡ON | **FAIL** (efficacy) |
| 2b | Saturating burial panel (`cf_wal`≥45) | same | 5/5 native CF-min OFF≡ON; rescued=0 | scoring competitiveness **PASS**; wall un-cap efficacy **FAIL** |
| 3 | W1 serial pilot | `FLEXAIDDS_BOOM_INTERVAL=50` only | 8/8 done; genuine 1/8 vs 0/8 base; BCR 3/8=3/8; **1N1M clean elect RMSD 2.28→5.66** | **FAIL** |
| 4a | Memetic | — | Blocked: wall efficacy FAIL | **NOT RUN** |
| 4b/c | Niche / coarse-init | — | Next one-variable after STEP3 FAIL | **NOT RUN** |
| 5 | Full-85 claim | — | Blocked until Steps 1–4 gates | **NOT RUN** |

## STEP 3 detail

- **OUT:** `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740`
- **Workers 2 · R=5 · matrix 9dc9 · boom_interval 100→50** confirmed in logs
- **Genuine:** 1/8 (1YGC elect 1.75 Å) vs baseline 0/8 on panel
- **BCR&lt;2:** 3/8 (1OQ5, 1SQ5, 1YGC) — same count as baseline; gap targets BCR *worse* (1.06→1.65, 1.12→1.65)
- **Clean regression:** 1N1M elected 2.28→5.66 Å → **FAIL ACCEPT**
- Mean Δ elect RMSD **+0.30**; mean Δ BCR **+0.48** (worse overall)
- Full gate: `workorders/STEP3_PILOT_GATE.md`

## Explicit non-claims

- Baseline **25.3%** is pre-`free_energy_strict` — not election-fix proof
- Do **not** set `FLEXAIDDS_WALL_PILOT_PASS=1`
- Do **not** dual full-85; WORKERS≤4
- BOOM_INTERVAL=50 is **not** a validated sampling lever on this panel

## Next (methodology)

1. STEP 3 FAIL → try **another one-variable** W1 knob (e.g. coarse-init / diversity / pop scale), not memetic  
2. Wall un-cap still open diagnosis (OFF≡ON even when `cf_wal`>CAP)  
3. No full-85 until cheap gates pass  
