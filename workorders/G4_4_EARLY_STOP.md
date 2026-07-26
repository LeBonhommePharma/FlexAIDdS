# G4.4 Early-stop audit — truncation is COMMON

**Written:** 2026-07-26T04:44:56.518264+00:00
**Nominal max_generations:** 2000
**Logs scored:** 45
**Truncated (&lt;90% budget or early reason):** 45
**Truncation common:** **True**

## Finding

SEARCH-MISS pilot/matched logs almost never reach 2000 gens. Dominant reasons:
- **fitness_stagnation** (CF stagnant 300–400 gens + gene-space collapsed) — e.g. 1N1M, 1K3U, 1M2Z
- **entropy_convergence / H plateau** — e.g. 1J3J, 1L7F

## Per-log (gens_reached approx)

| Leaf | PDB | tag | gens | reason | truncated |
|------|-----|-----|-----:|--------|:---------:|
| pilot_w1 | 1J3J | main | 390 | entropy_convergence | Y |
| pilot_w1 | 1J3J | r1 | — | entropy_convergence | Y |
| pilot_w1 | 1J3J | r2 | 490 | fitness_stagnation | Y |
| pilot_w1 | 1J3J | r3 | 140 | fitness_stagnation | Y |
| pilot_w1 | 1J3J | r4 | — | entropy_convergence | Y |
| pilot_w1 | 1K3U | main | 1140 | fitness_stagnation | Y |
| pilot_w1 | 1K3U | r1 | 1340 | h_plateau | Y |
| pilot_w1 | 1K3U | r2 | 1340 | h_plateau | Y |
| pilot_w1 | 1K3U | r3 | 1000 | fitness_stagnation | Y |
| pilot_w1 | 1K3U | r4 | 1300 | entropy_convergence | Y |
| pilot_w1 | 1L7F | main | 1200 | entropy_convergence | Y |
| pilot_w1 | 1L7F | r1 | 290 | entropy_convergence | Y |
| pilot_w1 | 1L7F | r2 | 160 | entropy_convergence | Y |
| pilot_w1 | 1L7F | r3 | 670 | fitness_stagnation | Y |
| pilot_w1 | 1L7F | r4 | 900 | entropy_convergence | Y |
| pilot_w1 | 1N1M | main | 500 | fitness_stagnation | Y |
| pilot_w1 | 1N1M | r1 | 500 | fitness_stagnation | Y |
| pilot_w1 | 1N1M | r2 | 300 | entropy_convergence | Y |
| pilot_w1 | 1N1M | r3 | 280 | entropy_convergence | Y |
| pilot_w1 | 1N1M | r4 | — | fitness_stagnation | Y |
| pilot_w1 | 1M2Z | main | 500 | fitness_stagnation | Y |
| pilot_w1 | 1M2Z | r1 | — | fitness_stagnation | Y |
| pilot_w1 | 1M2Z | r2 | 470 | fitness_stagnation | Y |
| pilot_w1 | 1M2Z | r3 | 600 | fitness_stagnation | Y |
| pilot_w1 | 1M2Z | r4 | 340 | fitness_stagnation | Y |
| coarse_A64 | 1J3J | main | 390 | entropy_convergence | Y |
| coarse_A64 | 1J3J | r1 | — | entropy_convergence | Y |
| coarse_A64 | 1K3U | main | 1140 | fitness_stagnation | Y |
| coarse_A64 | 1K3U | r1 | 1340 | h_plateau | Y |
| coarse_A64 | 1L7F | main | 1200 | entropy_convergence | Y |
| coarse_A64 | 1L7F | r1 | 290 | entropy_convergence | Y |
| coarse_A64 | 1N1M | main | 500 | fitness_stagnation | Y |
| coarse_A64 | 1N1M | r1 | 500 | fitness_stagnation | Y |
| coarse_A64 | 1M2Z | main | 500 | fitness_stagnation | Y |
| coarse_A64 | 1M2Z | r1 | — | fitness_stagnation | Y |
| coarse_B256 | 1J3J | main | 310 | entropy_convergence | Y |
| coarse_B256 | 1J3J | r1 | 520 | entropy_convergence | Y |
| coarse_B256 | 1K3U | main | 780 | entropy_convergence | Y |
| coarse_B256 | 1K3U | r1 | 1200 | entropy_convergence | Y |
| coarse_B256 | 1L7F | main | 390 | entropy_convergence | Y |
| coarse_B256 | 1L7F | r1 | 860 | entropy_convergence | Y |
| coarse_B256 | 1N1M | main | 480 | fitness_stagnation | Y |
| coarse_B256 | 1N1M | r1 | 400 | entropy_convergence | Y |
| coarse_B256 | 1M2Z | main | 500 | fitness_stagnation | Y |
| coarse_B256 | 1M2Z | r1 | 270 | fitness_stagnation | Y |

## Gate consequence (PHASE4_GATES_ACTUALIZED)

Truncation **is common** → fix measurement conditions **before** G4.1/G4.2 science claims.

**Fix (shipped path, not a new default product change):** Phase 4 sampling docks MUST set
`FLEXAIDDS_NO_SEC=1` (and preferably `FLEXAIDDS_BENCHMARK=1`). In `gaboom.cpp`, `no_sec`
disables H-plateau/SEC early exits and converts stagnation into exploration boost instead of break.
Documented on G4.2/G4.1 drivers. Interactive/default claim path unchanged.

Do **not** interpret prior COARSE/BOOM pilots as full-budget results without this audit.

