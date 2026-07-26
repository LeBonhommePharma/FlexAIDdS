# PHASE 4 GATES — actualized from Phase 2/3 results
OPS 2026-07-26. Supersedes BENCHMARKING_ROADMAP.md §6 Phase 4. Read before the next experiment.

===============================================================================
## WHAT PHASE 2/3 SETTLED (all verdicts verified against the data, not the label)
===============================================================================
P2 pb_clash on SCORING-LOCKED (1OQ5/1SQ5/1YGC), weights 1/5/10 -> **FAIL, definitively**.
  w=10: 1OQ5 +17.89->+13.34 (-4.55) | 1SQ5 +28.80->+28.80 (0.0000) | 1YGC +70.20->+68.67 (-1.53)
  ge_floor 2/3 but sign flips 0/3, and SEARCH-MISS regressions 1/5.
  This is not "needs more weight". Linear extrapolation of the observed slopes:
    1OQ5 sign flip at weight ~39 · 1YGC at weight ~460 · 1SQ5 NEVER (slope exactly 0).
  And at w=10 it ALREADY inverted a target it was supposed to protect:
    1J3J dCF -2.35 (native better) -> +0.96 (decoy better), native_min N, REGRESSED.
  So the usable weight window is empty: below 10 nothing flips, at 10 it breaks clean targets.

**THE MECHANISM (measured, and it retires a long-running thesis).**
I measured the clash geometry of each SCORING-LOCKED decoy (= that target's ACTUAL elected pose):
    1OQ5 min 1.99 A, 2 pairs <2.0 A · 1SQ5 min 2.38 A, ZERO <2.0 A · 1YGC min 2.09 A, ZERO <2.0 A
pb_clash penalizes overlap under ~0.75*(vdw_i+vdw_j) ~ 2.4-2.6 A. These decoys are essentially
CLASH-FREE. 1SQ5's exact 0.0000 at every weight is not a bug — pb_clash is structurally inert on
a pose with nothing to penalize.
=> **The false minima are NOT over-burial artifacts.** They are clash-free poses that the
ATTRACTIVE term scores better than the crystal. Every burial-opponent lever (wall un-cap,
pb_clash, com cap) is aimed at a mechanism that is not operating on these targets.
Two independent oracles have now failed for two different structural reasons. Stop building
burial opponents.

P3 sampling pilots -> BOOM_INTERVAL alone VOID (unwired); BOOM_FRAC=0.1 liveness PASS;
   COARSE 64 vs 256 matched A/B FAIL (genuine 0/5 both arms, mean dBCR +3.44).

===============================================================================
## MEMETIC: BOTH UNLOCK PATHS ARE NOW CLOSED — a decision is required
===============================================================================
config_parser.cpp now accepts FLEXAIDDS_PB_CLASH_PHASE2_PASS **or** the legacy
FLEXAIDDS_WALL_PILOT_PASS. Wall oracle: structurally unpassable. pb_clash oracle: FAILED with an
empty weight window. So memetic is blocked with NO remaining path, and neither failure is
evidence about memetic itself — they are evidence that BURIAL is the wrong opponent.
DO NOT set either flag from the current data. Instead choose ONE:
  (a) Leave memetic locked and spend Phase 4 entirely on sampling (RECOMMENDED — see below).
  (b) Re-key the interlock to a criterion that actually tests the memetic risk. The risk memetic
      poses is "refining against a defective objective walks away from the crystal". The correct
      gate for that is a REFINEMENT-DIRECTION oracle, not a burial oracle:
      start from a near-native pose (BCR pose, <2 A), run refinement ON the intended objective,
      measure whether RMSD-to-crystal DECREASES. ACCEPT: mean dRMSD <= 0 on >=4/5 SEARCH-MISS
      probes. That directly measures the failure mode the interlock exists to prevent.
  Do not unlock memetic by weakening the gate.

===============================================================================
## PHASE 4 GATES (actualized) — sampling only, SEARCH-MISS panel only
===============================================================================
PANEL RULE (from INVERSION_MAP, non-negotiable):
  sampling levers -> SEARCH-MISS ONLY: 1J3J 1K3U 1L7F 1N1M 1M2Z
  scoring levers  -> SCORING-LOCKED ONLY: 1OQ5 1SQ5 1YGC
  A lever judged on the wrong class is void. This error produced 2 of the 4 void gates.

EVERY Phase 4 experiment carries, in order:
  1. LIVENESS L1-L4 (knob read · not overridden by the DatasetRunner-emitted JSON · can
     physically act on these inputs · fires observably in the log).
  2. MAGNITUDE FLOOR — direction alone is not evidence. For sampling: dBCR <= -0.5 A mean OR
     >=1 target crossing BCR<2 A. For scoring: |dCF| >= 1.0 kcal AND a sign flip.
  3. MATCHED CONTROL — same binary, same R, all other knobs unset. Grok's matched COARSE re-run
     is the template; an unmatched comparison is void regardless of outcome.
  4. CLEAN-PROBE NON-REGRESSION — no SEARCH-MISS probe's elected RMSD worsens.

### G4.1 BOOM_FRAC small-fraction panel  [liveness already PASS]
  One variable: FLEXAIDDS_BOOM_FRAC in {0.05, 0.1, 0.2}, interval at claim JSON 100.
  Panel: the 5 SEARCH-MISS. R=5. Matched control = same binary, BOOM unset.
  EXTRA ACCEPT (the documented catastrophe): every restart must reach full generations.
  Early termination at gen ~300 with CF~0 is the known wipeout signature -> auto-FAIL.
  ACCEPT: mean dBCR <= -0.5 A or >=1 target crosses BCR<2 A, no clean-probe regression.

### G4.2 Niche metric — calc_rmsp  [HIGHEST-CONFIDENCE REMAINING LEVER]
  gaboom.cpp:4111-4123 is not a distance in pose space: an unweighted RMS over a gene vector
  whose gene 0 is a cleft-grid ORDINAL (0..15131) mixed with genes 1-9 in degrees. With
  sigma_share=204.1945: a 0.375 A step in z exits the niche while 7.9 A in y stays inside, and
  flipping ALL NINE angles 180 deg gives rmsp=170.8 < 204.19 (two fully reoriented poses = one
  niche). This is a structural defect, not a tuning parameter, and it directly starves the
  basin diversity that SEARCH-MISS needs.
  Change: replace with Cartesian RMSD over ligand heavy atoms; re-calibrate sigma_share in
  Angstroms (start 2.0 A).
  Gate it env-OFF by default; A/B on the 5 SEARCH-MISS with a matched control.
  ACCEPT: as above, plus n_niches occupied must increase (a liveness proxy for the fix acting).

### G4.3 Mutation granularity
  Only ~71 gene bits are live (gene0 14, each orientation gene 7, each torsion ~5) because a bit
  is dead unless weight >= 2^31/nbin, so no operator can make a SMALL move and a near-native
  chromosome is unimprovable. Fix jointly with G4.2 — but as SEPARATE A/B arms, never one run.

### G4.4 Early-stop audit  [do this FIRST — it is free and may explain the others]
  Restarts terminating at ~300-1000 of 2000 generations was observed on 1N1M. If the GA is
  being truncated, every sampling lever above is being measured on a shortened search and the
  results are misleading. Read the existing pilot logs; no new docking required.
  ACCEPT: report per-target generations-reached distribution. If truncation is common, fix the
  stopping rule BEFORE running G4.1/G4.2.

### ORDER: G4.4 (free, offline) -> G4.2 (highest confidence) -> G4.1 -> G4.3.
Rationale: G4.4 may invalidate the others' measurement conditions; G4.2 is a proven structural
defect while G4.1 is re-enabling something deliberately disabled.

===============================================================================
## PHASE 5 UNCHANGED
===============================================================================
Full 85 · R=10 · WORKERS<=4 · seed OFF · matrix 9dc9 · detached · full provenance ·
ONE config delta vs the 25.3% genuine / 27.8% BCR baseline. Report the triple. Never the raw column.

===============================================================================
## VOID-GATE SCOREBOARD — 4, all wrong-instrument, none wrong-execution
===============================================================================
1. Wall oracle — term structurally blind to the triggering values.
2. BOOM_INTERVAL — knob not wired (fraction hardcoded 0.0).
3. pb_clash on clean probes — wrong panel class + no magnitude floor (OPS error).
4. pb_clash on SCORING-LOCKED — right panel, but the decoys are clash-free so the lever is
   inert; empty weight window. This one FAILED HONESTLY and is a real negative result.
The liveness precheck + magnitude floor + class-matched panel now cover all four failure modes.
