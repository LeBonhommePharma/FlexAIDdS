# Campaign methodology follow-through — gate summary

**Sources:**  
- Downloads `CAMPAIGN_METHODOLOGY_for_Grok.md` (order; STEP 2 wall **superseded**)  
- Downloads `BENCHMARKING_ROADMAP.md` (**authoritative** on liveness L1–L4 + pb_clash Phase 2)  
- On-main: `docs/implementation/CAMPAIGN_METHODOLOGY_for_Grok.md`  

**Updated:** 2026-07-26 · **main tip at write:** see git log  

## Regime (Phase 1 — fixed)

| Metric | Value | Source |
|--------|------:|--------|
| Genuine top-1 | **20/79 = 25.3%** | `v_autonomous_20260724_160919` pre-`free_energy_strict` |
| BCR (sampling ceiling) | **22/79 = 27.8%** | same |
| Election gap | **2** (baseline) / **16/85=18.8%** (E10 full) | baseline + `E10_election_vs_scoring.md` |
| **Regime** | **SAMPLING-LIMITED** | BCR≈genuine; election not the primary wall |

Do **not** cite 25.3% as proof election fix worked.

## Critical invalidations (instrumentation — not docking-quality fails)

| Gate | Label | Why |
|------|--------|-----|
| **STEP 2 WAL_COERCIVE** | **Structurally unpassable (B3)** | Per-pair cap; Voronoi never sees deep burial. OFF≡ON expected. |
| **STEP 3 BOOM_INTERVAL only** | **Scientifically invalid as BOOM (B1)** | Claim JSON `boom_inject_fraction: 0.0`; inject needs fraction>0. INVALID as BOOM efficacy. |

## Gate table (BENCHMARKING_ROADMAP phases)

| Phase | Step | One variable / check | Result | PASS/FAIL |
|-------|------|----------------------|--------|-----------|
| 0 | Foundation | CF gate --config+--ligand; binaries; matrix 9dc9 | gate scores n=5, 0 skip; md5 9dc9 | **PASS** |
| 1 | E10 + M2 triple | offline | sampling-limited | **PASS** |
| 1.5 | Native–Elected CF inversion | pose role (fixed LOCCLF) | 8/8; **SEARCH-MISS=5** clean; **SCORING-LOCKED=3** gap | **PASS** |
| 2a | WAL wall | WAL_COERCIVE | structural no-op | **VOID / withdrawn** |
| 2b | pb_clash SEARCH-MISS (legacy) | `PB_CLASH_WEIGHT=1.0` | ΔdCF 1e-4..0.02; wrong panel | **VOID** (ROADMAP_v2) |
| 2b′ | pb_clash SCORING-LOCKED | weight 1/5/10; elected decoys | 0 sign flips; max decrease 4.55 @w=10 | **FAIL** (magnitude floor) |
| 3′ | BOOM small frac | `BOOM_FRAC=0.1` | live inject; no wipe; same false-min elect | **PASS liveness** |
| 3 | BOOM interval pilot | interval only | void under L2 | **INVALID** |
| 4.1 | COARSE matched 64 vs 256 (SEARCH-MISS) | COARSE_ORIENTATIONS only | L4 5/5 both; genuine 0/5 both; mean ΔRMSD +0.12; mean ΔBCR +3.44 | **FAIL** (matched; no directional gain) |
| 5 | Full-85 claim | — | blocked | **NOT RUN** |

## Explicit blocks

- No dual full-85; WORKERS≤4; OMP=1/worker; no build while docking holds binary  
- No memetic / no `WALL_PILOT_PASS` from WAL-only or micro-pb_clash  
- No WAL_COERCIVE re-panels expecting OFF≠ON  
- No interval-only BOOM; no BOOM_FRAC=1.0 blind  
- Matrix **9dc9** for docks; genuine metric only for rates  

## Artifacts

| Topic | Path |
|-------|------|
| E10 | `workorders/E10_election_vs_scoring.md` |
| Wall WAL | `workorders/WALL_ORACLE.md` |
| Inversion map | `workorders/INVERSION_MAP.md` |
| Next-step design | `workorders/NEXT_CAMPAIGN_STEP.md` |
| BOOM | `BOOM_FRAC_LIVENESS.md`, `BOOM_FRAC_AB.md` |
| pb_clash VOID | `PB_CLASH_ORACLE.md` |
| pb_clash 2b′ | `PB_CLASH_SCORING_LOCKED.md`, `scripts/pb_clash_burial_oracle.py --mode scoring-locked` |
| ROADMAP_v2 | `ROADMAP_v2_PANEL_CORRECTION.md` |
| Audit B1–B3 | `DOCKING_BUG_AUDIT_2026-07-25.md` |
| Inversion script | `scripts/native_elected_cf_inversion_map.py` |

## Next allowed (after matched COARSE FAIL)

1. ~~COARSE_ORIENTATIONS=256~~ — **matched FAIL** (see MATCHED_AB_GATE). Do not re-run same lever.  
2. Next **one-var** sampling levers on SEARCH-MISS: niche Cartesian (flagged), or `FLEXAIDDS_NO_SEC` budget honesty.  
3. SCORING-LOCKED (1OQ5/1SQ5/1YGC): strong burial decoys / scoring — not BOOM thrash.  
4. Full-85 only after Phase 4 sampling gates pass.

## Post-COARSE status

Prior multi-var COARSE=256 vs pilot_w1 is **VOID** (skeptic). Matched control
`~/flexaidds_results/coarse_ab_matched_20260725_222652` (same binary `7f05640a…`,
git `25b21216`, R=2, BOOM unset): **FAIL** — genuine 0/5 both arms; mean Δ elect
RMSD **+0.12** Å; mean Δ BCR **+3.44** Å (worse under 256 on 1J3J BCR). L1–L4 PASS
both arms (`[COARSE-INIT] … 64/256 orientations` ×5).

**Workorders:** `MATCHED_AB_GATE.md`, expanded `COARSE_ORIENT_W1_PILOT.md`.

Next single levers (still one variable): niche Cartesian distance (code+flag),
`FLEXAIDDS_NO_SEC` budget honesty, or strong burial decoys for SCORING-LOCKED.
Full-85 still blocked.

## ROADMAP_v2 Phase 2 correction (2026-07-26)

- Prior SEARCH-MISS pb_clash **VOID** (wrong panel + no magnitude floor).
- Class-matched rule: sampling → SEARCH-MISS only; scoring → SCORING-LOCKED only.
- Magnitude floor: scoring-oracle PASS requires **≥1.0 kcal** dCF decrease **and** sign flip on **≥2/3** inverted targets.
- Revised 2b′ on 1OQ5/1SQ5/1YGC elected decoys: **FAIL** (ladder w=1/5/10; 0 sign flips).
- Memetic re-keyed to `FLEXAIDDS_PB_CLASH_PHASE2_PASS` (or legacy WALL); **still OFF** until 2b′ PASS.
- Full-85 still blocked.
