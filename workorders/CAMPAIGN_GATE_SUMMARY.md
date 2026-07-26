# Campaign methodology follow-through — gate summary

**Sources:**  
- Downloads `CAMPAIGN_METHODOLOGY_for_Grok.md` (order; STEP 2 wall **superseded**)  
- Downloads `BENCHMARKING_ROADMAP.md` (**authoritative** on liveness L1–L4 + pb_clash Phase 2)  
- Downloads `PHASE4_GATES_ACTUALIZED.md` / `ROADMAP_v2_PANEL_CORRECTION.md`  
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

## Gate table (BENCHMARKING_ROADMAP + PHASE4_GATES_ACTUALIZED)

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
| **G4.4** | Early-stop audit | offline gens-reached on pilot/matched/boom logs | 45/45 truncated vs 2000-gen budget | **PASS audit** — truncation **common**; Phase 4 docks require `FLEXAIDDS_NO_SEC=1` |
| **G4.2** | Cartesian niche matched A/B | `FLEXAIDDS_NICHE_CARTESIAN` OFF vs ON (σ=2.0 Å); both `NO_SEC=1` | L4 cart B only; genuine 0/5 both; mean ΔBCR **−0.441** Å; 0 elect reg | **FAIL** (misses ≤−0.5 floor; directional only) |
| G4.1 | BOOM_FRAC panel {0.05,0.1,0.2} | — | — | **DEFERRED** (not run this session) |
| G4.3 | Mutation granularity | — | — | **NOT RUN** (separate arm later) |
| 5 | Full-85 claim | — | blocked | **NOT RUN** |

### G4.2 OUT / provenance

- Root: `~/flexaidds_results/g4_2_niche_cart_ab_20260726_004752`  
- Workorder: `workorders/G4_2_NICHE_CART.md` · gate JSON under OUT  
- Ship: `LIB/niche_distance.h` + `gaboom.cpp` env gate (OFF default)

## Explicit blocks

- No dual full-85; WORKERS≤4; OMP=1/worker; no build while docking holds binary  
- No memetic / no `WALL_PILOT_PASS` / no `PB_CLASH_PHASE2_PASS` from burial or wall data  
- No WAL_COERCIVE re-panels expecting OFF≠ON  
- No interval-only BOOM; no BOOM_FRAC=1.0 blind  
- **Burial opponents retired** (SCORING-LOCKED false mins are clash-free)  
- Matrix **9dc9** for docks; genuine metric only for rates  
- Sampling levers only on SEARCH-MISS; scoring only on SCORING-LOCKED  

## Artifacts

| Topic | Path |
|-------|------|
| E10 | `workorders/E10_election_vs_scoring.md` |
| Wall WAL | `workorders/WALL_ORACLE.md` |
| Inversion map | `workorders/INVERSION_MAP.md` |
| Next-step design | `workorders/NEXT_CAMPAIGN_STEP.md` |
| BOOM | `BOOM_FRAC_LIVENESS.md`, `BOOM_FRAC_AB.md` |
| pb_clash VOID | `PB_CLASH_ORACLE.md` |
| pb_clash 2b′ | `PB_CLASH_SCORING_LOCKED.md` |
| ROADMAP_v2 | `ROADMAP_v2_PANEL_CORRECTION.md` |
| Phase 4 actualized | `PHASE4_GATES_ACTUALIZED.md` |
| G4.4 early-stop | `G4_4_EARLY_STOP.md` |
| G4.2 niche cart | `G4_2_NICHE_CART.md` |
| Memetic status | `MEMETIC_STATUS.md` (option **a** locked) |
| Matched COARSE | `MATCHED_AB_GATE.md` |
| Niche unit API | `LIB/niche_distance.h` |

## Next allowed (after G4.2 FAIL)

1. ~~COARSE 64 vs 256~~ — matched **FAIL** (`MATCHED_AB_GATE`). Do not re-run.  
2. ~~Burial / pb_clash weight ladder~~ — **retired** (empty weight window; clash-free SCORING-LOCKED). Do not re-open.  
3. ~~G4.2 Cartesian niche (σ=2.0)~~ — matched **FAIL** (mean dBCR −0.441). Optional retune σ / pairing as a **new** one-var experiment, not a silent re-run.  
4. **Next primary:** **G4.1** BOOM_FRAC ∈ {0.05, 0.1, 0.2} on SEARCH-MISS, matched control, `NO_SEC=1`, magnitude floor, wipeout auto-FAIL.  
5. Then **G4.3** mutation granularity as a **separate** arm.  
6. Full-85 only after Phase 4 sampling gates clear.  

## Post-COARSE status

Matched control `~/flexaidds_results/coarse_ab_matched_20260725_222652`: **FAIL** — genuine 0/5 both arms; mean Δ elect +0.12 Å; mean Δ BCR +3.44 Å.

## ROADMAP_v2 Phase 2 correction (2026-07-26)

- Prior SEARCH-MISS pb_clash **VOID**; 2b′ SCORING-LOCKED **FAIL** honest.  
- Class-matched + magnitude floor rules permanent.  
- Memetic: option **(a)** locked (`MEMETIC_STATUS.md`).  

## Phase 4 actualized (PHASE4_GATES_ACTUALIZED.md)

- **Burial retired** · **Memetic (a)** · Order G4.4 → G4.2 → G4.1 → G4.3  
- **G4.4 PASS audit** (truncation common → `NO_SEC=1` for Phase 4 docks)  
- **G4.2 FAIL** science (L4 OK; magnitude floor missed)  
- **G4.1 DEFERRED** — not completed  
