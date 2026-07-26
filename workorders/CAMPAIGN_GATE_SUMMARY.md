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
| 2b | pb_clash burial | `PB_CLASH_WEIGHT=1.0` | 5/5 micro ΔdCF; cf_clash≈0 | **PASS formal**; **not** memetic unlock |
| 3′ | BOOM small frac | `BOOM_FRAC=0.1` | live inject; no wipe; same false-min elect | **PASS liveness** |
| 3 | BOOM interval pilot | interval only | void under L2 | **INVALID** |
| 4.1 | COARSE_ORIENTATIONS=256 SEARCH-MISS | COARSE_ORIENTATIONS | live 5/5; genuine 0/5; BCR 0/5; mean ΔRMSD +0.29 | **FAIL** (no directional gain) |
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
| pb_clash | `PB_CLASH_ORACLE.md`, `scripts/pb_clash_burial_oracle.py` |
| Audit B1–B3 | `DOCKING_BUG_AUDIT_2026-07-25.md` |
| Inversion script | `scripts/native_elected_cf_inversion_map.py` |

## Next allowed (after inversion map)

1. If **SEARCH-MISS** dominates clean probes → one-var `FLEXAIDDS_COARSE_ORIENTATIONS=256` W1 (matrix 9dc9, WORKERS≤2).  
2. If **SCORING-LOCKED** dominates gap targets → strong burial decoys / scoring (not BOOM thrash).  
3. Full-85 only after Phase 4 sampling gates.

## Post-COARSE status

COARSE_ORIENTATIONS=256 **FAIL**ed directional ACCEPT on SEARCH-MISS panel (liveness OK). Next single levers: niche Cartesian distance (code+flag), or `FLEXAIDDS_NO_SEC` budget honesty, or strong burial decoys for SCORING-LOCKED. Full-85 still blocked.
