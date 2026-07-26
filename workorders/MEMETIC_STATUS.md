# Memetic interlock — Phase 4 decision **(a)**

**Written:** 2026-07-26  
**Decision:** **Leave memetic locked** (PHASE4_GATES_ACTUALIZED option a — recommended).

## Why both unlock paths are closed

| Key | Status |
|-----|--------|
| `FLEXAIDDS_WALL_PILOT_PASS` | Structurally unpassable (B3) |
| `FLEXAIDDS_PB_CLASH_PHASE2_PASS` | SCORING-LOCKED pb_clash **FAIL** — empty weight window; clash-free false mins |

Neither failure is evidence about memetic refinement risk. Do **not** set either flag from current burial/wall data.

## Shipped gate (still default-OFF)

`LIB/memetic_gate.h`: requires `FLEXAIDDS_MEMETIC=1` **and** (PHASE2_PASS **or** WALL_PASS).  
Claim/default path leaves all unset → `use_memetic=0`.

## Option (b) not implemented this phase

Refinement-direction oracle (near-native start, mean dRMSD ≤ 0 on ≥4/5 SEARCH-MISS) is the correct future memetic gate if needed. Not built in this Phase 4 sampling-first pass.
