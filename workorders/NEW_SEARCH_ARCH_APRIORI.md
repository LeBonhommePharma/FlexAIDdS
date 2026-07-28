# new_search_arch — a priori only (NO LAUNCH)

**Status:** `APRIORI_ONLY_NO_LAUNCH`  
**Machine JSON:** [`NEW_SEARCH_ARCH_APRIORI.json`](NEW_SEARCH_ARCH_APRIORI.json)  
**Trigger:** Flip residual after Phase-4 near-miss null stack  
(`G4.3_null_phase4_sampling_exhausted` → `new_search_arch`).

## What this is

A **pre-registration** of the residual sampling architecture experiment.  
It does **not** authorize a dock, claim recipe change, or full-85.

## Panel

- Class: **NEAR_MISS** only — `1N1M`, `1L7F`  
- Not SCORING-LOCKED (1OQ5/1SQ5/1YGC) — see `SCORING_LOCKED_SI_PACKAGE.md`

## One variable (choose before any future launch)

Exactly **one** env-gated change vs matched control (env unset):

| Option | Sketch | Env (example) |
|--------|--------|----------------|
| A | Phenotype uniqueness policy isolated from granular mut | `FLEXAIDDS_NEW_SEARCH_ARCH=phenotype_unique` |
| B | Basin-aware reinject on Cartesian ligand RMSD | `FLEXAIDDS_NEW_SEARCH_ARCH=basin_reinject` |
| C | Selection architecture change | `FLEXAIDDS_NEW_SEARCH_ARCH=selection_v2` |

Default in product builds: **OFF**. Never bundle BOOM + election + mutation.

## Floors

- L1–L4: env read; not JSON-stuck; acts; stderr marker `[NEW-SEARCH-ARCH]` control zero  
- Magnitude: mean ΔBCR ≤ −0.5 Å **or** ≥1 BCR&lt;2 **or** 1N1M elect ≤2.5; no wipeout  
- Pins: matrix **9dc9**, same binary SHA both arms, NO_SEC=1, Sol #9, `validate-pins`

## Explicit non-launch

- No OUT path claiming a new run  
- No full-85  
- No memetic / WAL / burial reopen  
- Implementation of engine code is a **separate** goal after design picks A/B/C
