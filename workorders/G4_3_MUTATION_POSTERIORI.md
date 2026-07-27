# G4.3 MUTATION_GRANULAR near-miss — a posteriori (FINAL)

**OUT:** `/Users/lp.more/flexaidds_results/g4_3_mut_gran_near_miss_20260727_122215`  
**Completed:** 2026-07-27T19:08:51Z (`ALL_ARMS_DONE` after mut relaunch)  
**Contract:** `workorders/G4_3_MUTATION_APRIORI.json` + `BENCHMARK_SELF_EVAL_CONTRACT.md`  
**Binary (treatment relaunch):** `19f300d9798d4985423fb501697ee3b397cc57596040ccfe7c84a8f8165225f6`  
**Code tip:** `b8b19468` (per-gen phenotype-duplicate clear under granular)

## Protocol

| Item | Value |
|------|--------|
| One variable | `FLEXAIDDS_MUTATION_GRANULAR=1` vs unset |
| Panel | NEAR_MISS `1N1M,1L7F` |
| Restarts | R=2 |
| Matrix | 9dc9 |
| NO_SEC | 1 |
| Matched control | yes (same OUT tree, sequential arms) |

### Operational note (stall → fix → relaunch)

First mut arm stalled after ~gen 4100 at ~0.08 gen/s. Root cause: lifetime
`duplicates` map keys **phenotype `to_ic`**; ±1-bin local steps exhaust the
neighborhood and `reproduce()` rejection-samples forever. Control (classic
bit-flip) can still jump to unexplored phenotypes via high-bit flips.

**Fix (treatment only, control already complete):** clear `duplicates` each
generation when `FLEXAIDDS_MUTATION_GRANULAR` is on (`LIB/gaboom.cpp`, commit
`b8b19468`). L4 also logs
`[MUT-GRAN] per-generation phenotype-duplicate clear…`. Mut arm wiped and
relaunched with restamped binary; control `result.csv` preserved.

Evidence: `OUT/evidence/mut_stall_diagnosis.txt`, `stalled_mut_partial/`.

## L4 instrumentation

| arm | n `[MUT-GRAN]` markers |
|-----|------------------------:|
| control | **0** |
| mut_gran | **8** |

L4 **PASS** (treatment live, control zero).

## Magnitude (best_cluster_rmsd / elect)

| arm | 1L7F elect | 1L7F BCR | 1N1M elect | 1N1M BCR |
|-----|-----------:|---------:|-----------:|---------:|
| control | 3.9907 | 3.9907 | 6.3999 | 4.1954 |
| mut_gran | 6.2458 | 4.1128 | 6.0053 | 4.3088 |
| **Δ** | **+2.255** | **+0.122** | **−0.395** | **+0.113** |

- mean_dBCR = **+0.118 Å** (treatment slightly **worse**)
- Floor mean_dBCR ≤ −0.5 or BCR&lt;2 or elect≤2.5 → **FAIL (null magnitude)**
- 1L7F elect **regressed** (3.99 → 6.25)

## Decision

| Layer | Result |
|-------|--------|
| L4 | **PASS** |
| Magnitude | **FAIL (null)** |
| status | **PASS_LIVENESS** |
| `ACCEPT_G4_3` | **False** |

### Flip residual

```
rule: G4.3_null_phase4_sampling_exhausted
action: Phase4_sampling_stack_null; residual_new_search_arch_or_scoring_locked_work
priority_order: [new_search_arch, scoring_locked_decoy_work, full85_still_blocked]
```

Phase-4 sampling levers on this near-miss panel are now exhausted nulls:

1. G4.1 BOOM_FRAC — null mag, L4 PASS  
2. ELECTION_V135 — null mag (identical elect)  
3. G4.3 MUTATION_GRANULAR — null mag (slight BCR regression), L4 PASS  

**Full-85 remains blocked** until a Phase-4 ACCEPT.

## Non-claims

- Not genuine top-1 / PoseBusters success rates.  
- CF/contact-function scoring proxy only.  
- Not permission to enable `MUTATION_GRANULAR` in claim recipe by default.
