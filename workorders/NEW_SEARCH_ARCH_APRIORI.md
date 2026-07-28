# new_search_arch — A+B implemented (NO DOCK LAUNCH yet)

**Status:** `IMPLEMENTED_A_B_NO_DOCK_YET`  
**Machine JSON:** [`NEW_SEARCH_ARCH_APRIORI.json`](NEW_SEARCH_ARCH_APRIORI.json)  
**Code:** `LIB/new_search_arch.h` + wiring in `LIB/gaboom.cpp`  
**Trigger:** Flip residual after Phase-4 near-miss null stack.

## Implemented gates (default OFF)

| Option | Env | Behavior | L4 |
|--------|-----|----------|-----|
| **A** | `FLEXAIDDS_PHENOTYPE_UNIQUE=1` or `FLEXAIDDS_NEW_SEARCH_ARCH=phenotype_unique` | Classic mutate: if bit-flips leave phenotype bins unchanged, force ±1-bin step; reproduce uniqueness uses phenotype-bin hash | `[NEW-SEARCH-ARCH] phenotype_unique=1 …` |
| **B** | `FLEXAIDDS_BASIN_REINJECT=1` or `…=basin_reinject` | On diversity collapse (first half gens), reinject worst fraction preferring Cartesian ligand RMSD &gt; `FLEXAIDDS_BASIN_SIGMA_ANG` (default 2.0) vs best | `[NEW-SEARCH-ARCH] basin_reinject=1 …` + `[BASIN-REINJECT]` |
| Convenience | `FLEXAIDDS_NEW_SEARCH_ARCH=1` | Enables **both** A and B (dev only; **not** one-var for claim pilot) | both markers |
| **C** | — | Not implemented | — |

## Matched pilot (future — one variable only)

- Panel: NEAR_MISS `1N1M,1L7F`; R=2; matrix 9dc9; NO_SEC=1; Sol #9  
- Treatment: **either** PHENOTYPE_UNIQUE=1 **or** BASIN_REINJECT=1 (not both for claim)  
- Control: both unset; **same binary SHA**  
- Floors: mean ΔBCR ≤ −0.5 or BCR&lt;2 or 1N1M elect≤2.5; L4 control zero  

## Explicit non-claims

- No dock has been launched under this workorder yet  
- No full-85; no claim-recipe default ON  
- Do not bundle with BOOM / election / MUTATION_GRANULAR in one arm  
