# SI package: SCORING-LOCKED vs SEARCH-MISS (class-matched residual)

**Status:** OFFLINE packaging (no GA dock)  
**Date:** 2026-07-28  
**Sources:** `INVERSION_MAP.md` / `INVERSION_MAP.json`, `PB_CLASH_SCORING_LOCKED.md`,  
`PHASE4_NEAR_MISS_NULL_STACK.md`, `PHASE4_GATES_ACTUALIZED.md`  
**Purpose:** Publication SI decision tree so near-miss **sampling nulls** are not mixed with  
**scoring-locked false minima**.

---

## Class rule (non-negotiable)

| Class | Codes | Lever class only |
|-------|-------|------------------|
| **SEARCH-MISS** | 1J3J, 1K3U, 1L7F, 1N1M, 1M2Z | **Sampling** (GA inject, election, mutation, niche, coarse, …) |
| **SCORING-LOCKED** | 1OQ5, 1SQ5, 1YGC | **Scoring / landscape** (not more mid-run diversity on these codes) |

A lever judged on the wrong class is **VOID** (methodology).

---

## Inversion map (score-only, production LOCCLF)

**Verdict:** PASS (8/8 classified). ε = 0.5 CF units.  
**OUT:** `~/flexaidds_results/workorders/inversion_map_20260725_213932`  
**Pilot poses:** `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740`

### SCORING-LOCKED only (this package)

| PDB | class | CF native | CF elected | dCF(n−e) | elect RMSD | BCR |
|-----|-------|----------:|-----------:|---------:|-----------:|----:|
| 1OQ5 | SCORING-LOCKED | −107.870 | −125.757 | **+17.886** | 3.9457 | 1.6481 |
| 1SQ5 | SCORING-LOCKED | −135.211 | −164.592 | **+29.381** | 5.0963 | 1.6467 |
| 1YGC | SCORING-LOCKED | −146.591 | −216.077 | **+69.486** | 1.7516 | 1.0103 |

**Meaning:** Elected pose scores **better** on CF than crystal ligand coordinates while  
near-native BCR ≤2 Å exists in the pool → **false-min attractor on the scoring surface**,  
not “need more BOOM.”

### SEARCH-MISS (context only — sampling residual is frozen null)

See `PHASE4_NEAR_MISS_NULL_STACK.md` for G4.1 / ELECTION_V135 / G4.3 on near-miss  
**1N1M + 1L7F** (subset of SEARCH-MISS). Sampling stack: **null magnitude**; full-85 blocked.

---

## Burial / pb_clash residual (already closed — do not reopen)

From `PB_CLASH_SCORING_LOCKED.md` (weights 1/5/10 on SCORING-LOCKED elected decoys):

| Result | Detail |
|--------|--------|
| **FAIL** | 0/3 sign flips at every weight; empty usable weight window |
| Mechanism | Elected decoys are essentially **clash-free**; pb_clash inert on 1SQ5 (decrease 0) |
| Memetic | **Still locked** — do not set `PB_CLASH_PHASE2_PASS` or `WALL_PILOT_PASS` |

**Forbid in residual campaign:** burial re-panel, WAL_COERCIVE re-run expecting OFF≠ON,  
pb_clash weight ladder thrash, memetic unlock from these data.

---

## Decision tree for Methods/SI

```
inversion map (native CF vs elected CF, fixed LOCCLF)
        │
        ├─ SEARCH-MISS ──► sampling levers only
        │                    └─ Phase-4 near-miss stack: NULL (freeze table)
        │                    └─ residual: new_search_arch a priori (S4), not full-85
        │
        └─ SCORING-LOCKED ─► scoring / landscape only
                             └─ burial/pb_clash: FAIL / retired
                             └─ residual: non-burial scoring hypothesis OR
                                honest “landscape prefers decoy” claim
```

---

## Publication products

1. **SI table:** SCORING-LOCKED CF native vs elected (above).  
2. **SI figure (optional):** class split pie (5 SEARCH-MISS / 3 SCORING-LOCKED).  
3. **Text:** “Near-miss sampling levers (BOOM, election V135, mutation granular)  
   did not clear magnitude floors; SCORING-LOCKED targets require a scoring residual,  
   not further diversity inject.”  
4. **Claim language:** CF/contact-function **scoring proxy** only (`CLAIM_LANGUAGE_FREEZE.md`).

---

## Explicit non-claims

- Not genuine top-1 / PoseBusters rates.  
- Not true thermodynamic ΔG.  
- Not permission to re-run VOID levers.  
- Not a new GA dock.  
- Not full-85 unlock.
