# Campaign methodology follow-through — gate summary

**Source:** `docs/implementation/CAMPAIGN_METHODOLOGY_for_Grok.md`  
**Main tip at launch:** `05e1fa21` (E10 N=85 + saturating wall panel)  
**Updated:** 2026-07-25

| Step | Phase | One variable | Result | PASS/FAIL |
|------|-------|--------------|--------|-----------|
| 0 | Merge methodology + wall/E10 tooling on main | n/a | main == origin/main @ 05e1fa21 | **PASS** |
| 1 | E10 offline (frozen archive) | n/a | N=**85**; elected proxy 21/85=24.7%; BCR 24/85=28.2%; election-gap **16/85=18.8%**; sampling still primary | **PASS** (continue) |
| 2a | Wall WAL_COERCIVE (prod LOCCLF configs) | `FLEXAIDDS_WAL_COERCIVE` | 1M2Z native CF=**−117.74**; falsemin 4/5 OFF≡ON | **FAIL** (efficacy) |
| 2b | Saturating burial panel (`cf_wal`≥45) | same | 5/5 native CF-min OFF≡ON; rescued=0 | scoring competitiveness **PASS**; wall un-cap efficacy **FAIL** |
| 3 | W1 serial pilot | `FLEXAIDDS_BOOM_INTERVAL=50` only | OUT `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740`; workers=2; matrix 9dc9; live log shows `boom_interval 100→50` | **IN PROGRESS** (not PASS/FAIL yet) |
| 4a | Memetic | — | Blocked: needs `WALL_PILOT_PASS` after Step 2 efficacy PASS | **NOT RUN** |
| 4b/c | Niche / coarse-init | — | After Step 3 PASS/FAIL | **NOT RUN** |
| 5 | Full-85 claim | — | Blocked until Steps 1–4 gates | **NOT RUN** |

## Explicit non-claims

- Baseline **25.3%** is pre-`free_energy_strict` — do not cite as election-fix proof.
- Do **not** set `FLEXAIDDS_WALL_PILOT_PASS=1` (OFF≡ON even when `cf_wal`>CAP).
- Do **not** dual full-85; WORKERS≤4; no build while pilot holds binary.

## Artifacts

| Artifact | Path |
|----------|------|
| E10 workorder | `workorders/E10_election_vs_scoring.md` |
| Wall workorder | `workorders/WALL_ORACLE.md` |
| Pilot OUT | `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740/` |
| Pilot R5 addendum | `…/STEP3_PROVENANCE_ADDENDUM.json` |
| Campaign leaf (E10 input) | `~/flexaidds_results/v_autonomous_20260724_160919/` (85× result.csv) |

## STEP 3 ACCEPT (when complete)

1. No clean-probe regression vs frozen baseline on elected RMSD / BCR for 1J3J,1K3U,1L7F,1N1M,1M2Z.  
2. Directional effect on BCR or elected RMSD from BOOM_INTERVAL=50 alone.  
3. Report genuine / BCR / election-gap on the 8-target panel only (pilot, not claim).

## Baseline panel snapshot (frozen v_autonomous_20260724_160919)

| PDB | rmsd_hungarian | best_cluster_rmsd | seed_echo | note |
|-----|---------------:|------------------:|----------:|------|
| 1J3J | 62.31 | 21.35 | 0 | clean probe; far fail |
| 1K3U | 12.70 | 11.94 | 0 | clean probe |
| 1L7F | 4.22 | 3.98 | 0 | clean probe |
| 1N1M | 2.28 | 4.08 | 0 | clean probe; near 2Å elect |
| 1M2Z | 11.69 | 11.49 | 0 | clean probe; elected_cf sentinel |
| 1OQ5 | 3.78 | **1.06** | 0 | election-gap near-miss |
| 1SQ5 | 5.10 | **1.12** | 0 | election-gap near-miss |
| 1YGC | 3.34 | **1.24** | 0 | election-gap near-miss |
