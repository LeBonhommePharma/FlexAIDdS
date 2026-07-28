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
| G4.1 | BOOM_FRAC near-miss {0.05,0.1,0.2} | SEARCH-MISS; NO_SEC=1; R=2 | L4 PASS; best mean_dBCR −0.019 | **FAIL (null mag)** |
| ELECTION_V135 | election_v135 + τ=25 | R=5 near-miss | elect identical 6.40/3.99 | **FAIL (null mag)** |
| G4.3 | MUTATION_GRANULAR ±1-bin | near-miss; NO_SEC=1; R=2 | L4 PASS (8); mean_dBCR **+0.118** | **PASS_LIVENESS / null mag** |
| 5 | Full-85 claim | — | Phase-4 sampling stack null | **BLOCKED** |
| S2 | Closed-gate pins | `validate-pins` accept.txt + per-arm SHA | live G4.1/election/G4.3 PINS_OK | **PASS tooling** |
| S3 | SCORING-LOCKED SI package | offline class split | `SCORING_LOCKED_SI_PACKAGE.md` | **PACKAGED** |
| S4 | new_search_arch A+B | env-gated code | `new_search_arch.h` + a priori | **CODE IN; NO DOCK** |
| S5 | Claim language freeze | CF proxy vs STRICT ΔG | `CLAIM_LANGUAGE_FREEZE.md` | **FREEZE** |
| S4 A pilot | PHENOTYPE_UNIQUE near-miss | mean_dBCR −0.057; L4 PASS | `S4_PHENOTYPE_UNIQUE_POSTERIORI.md` | **PASS_LIVENESS** |
| Pre-gate triage | P2 + A/B/C bins | P2 HOLD; A/B reconstruction | `PUB_PREGATE_TRIAGE.md` | **TRIAGED** |

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

## Multi-session dock coordination (Sol #9)

- **Hold:** `~/flexaidds_results/BENCHMARK_HOLD.json` — any presence refuses new docks.
- **Lock:** `~/flexaidds_results/BENCHMARK_DOCK_LOCK/` — atomic `mkdir`; one owner.
- **CLI:** `python3 scripts/benchmark_coord.py status|preflight|release`
- **Offline queue:** `workorders/OFFLINE_BENCHMARKS_QUEUE.md`
- **WORKERS≤4** hard refuse; disk floor **20 GiB** (override only with `FLEXAIDDS_DISK_FLOOR_OVERRIDE=1`).
- Live G4.1 may already hold the lock; do not steal ownership.


## G4.1 BOOM near-miss (2026-07-27 FINAL)

- OUT: `g4_1_boom_near_miss_20260726_200953`
- L4 BOOM: LIVE on treatments; control zero
- Magnitude: **NULL** (best mean_dBCR=−0.0192 at frac010; floor −0.5)
- accept_g4_1: **False**
- Flip: **election_fix_P0** (1N1M offline pool 2.36 / elect 6.41)
- Next: `election_v135_near_miss_20260726_225823` (R=5, V135 vs control)
- Evidence: `workorders/g4_1_evidence/`, `workorders/G4_1_NEAR_MISS_POSTERIORI.md`


## ELECTION_V135 near-miss (2026-07-27 FINAL)

- OUT: `election_v135_near_miss_20260726_225823`
- One var: `FLEXAIDDS_ELECTION_V135=1` (tau=25) vs control; R=5; matrix 9dc9; NO_SEC
- Result: **NULL** — elect identical (1N1M 6.40 / 1L7F 3.99 both arms)
- accept: **False**
- Next: **G4.3 mutation** a priori draft `workorders/G4_3_MUTATION_APRIORI.json`
- Evidence: `workorders/ELECTION_V135_POSTERIORI.md`
