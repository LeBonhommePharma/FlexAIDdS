# ROADMAP CORRECTION v2 — the panel was wrong, and Grok's inversion map fixes it
2026-07-25, OPS. Corrects BENCHMARKING_ROADMAP.md Phase 2. Read before the next experiment.

## MY ERROR
Roadmap Phase 2 specified the burial oracle on the 5 "clean probes" (1J3J 1K3U 1L7F 1N1M 1M2Z).
Grok's INVERSION_MAP proves those 5 are ALL **SEARCH-MISS**: the crystal already scores BETTER
CF than the elected pose (dCF -2.4 to -236.2). On those targets the SCORING IS ALREADY CORRECT —
there is nothing for a scoring lever to fix. Running a burial oracle there measures nothing.

## WHY THE pb_clash "PASS" IS NOT A PASS
Reported: dCF moved toward native 5/5 -> PASS. Actual magnitudes:
  1J3J -0.0014 | 1K3U -0.0095 | 1L7F -0.0019 | 1N1M -0.0228 | 1M2Z -0.0001  (kcal)
That is 1e-4..2e-2 on CF values of 27-138 — i.e. ~1e-4 relative. It is numerical epsilon, not a
physical effect, and "moved 5/5" passes on the SIGN of noise.
ROOT CAUSE, measured directly on the decoys the oracle used (min ligand-receptor heavy-atom dist):
  1J3J 2.20 A | 1K3U 2.25 A | 1L7F 2.56 A | 1M2Z 2.96 A | 1N1M 2.24 A
  pairs <2.0 A: ZERO on all five.
pb_clash fires on interpenetration (ratio 0.75 x sum-vdW ~ 2.4-2.6 A and closer). The "buried"
decoys are barely in contact, so pb_clash has almost nothing to penalize. The decoys are not
buried. My own contingency clause applied and was not triggered: "IF ZERO EFFECT: the panel
lacks deep-clash decoys — BUILD decoys by deliberate burial." An epsilon effect IS zero effect.

## THE ACCEPTANCE TEST WAS ALSO WRONG
"dCF moves toward 0 on >=4/5" has no magnitude floor, so noise passes it. Any future scoring
oracle MUST require a magnitude: |dCF| shift >= 1.0 kcal AND a sign flip on at least one
genuinely inverted target. Direction alone is not evidence.

## THE CORRECT PANEL (Grok derived this; adopt it)
| class | targets | dCF(nat-elect) | meaning | correct lever |
|---|---|---|---|---|
| SEARCH-MISS | 1J3J 1K3U 1L7F 1N1M 1M2Z | -2.4 to -236.2 | native already CF-min; GA never found it | SAMPLING only |
| SCORING-LOCKED | 1OQ5 1SQ5 1YGC | +17.9 / +29.4 / +69.5 | elected beats crystal while BCR<=2 A | SCORING only |
The SCORING-LOCKED three are REAL inversions with near-native poses already in the pool
(BCR 1.65 / 1.65 / 1.01). They are the only valid substrate for a burial/scoring oracle, and
their elected poses are REAL false-minimum attractors — no synthetic decoy construction needed.

## REVISED PHASE 2 (run this instead)
Panel: 1OQ5 1SQ5 1YGC. Decoy = each target's ACTUAL elected pose (a true false-min attractor).
Native = crystal. Arms: FLEXAIDDS_PB_CLASH_WEIGHT ladder, one value per arm, --config + --ligand.
ACCEPT: dCF (currently +17.9/+29.4/+69.5) DECREASES by >= 1.0 kcal and FLIPS SIGN on >= 2/3,
with no SEARCH-MISS probe regressing (native must stay CF-min on all 5).
This is a real test: it can fail, and passing it means the objective stopped preferring a
false minimum over the crystal on targets where the right answer was already sampled.

## PHASE 4 SPLIT (Grok's implication 1-3; adopt verbatim)
- Sampling levers (COARSE, BOOM, niche, mutation granularity) -> SEARCH-MISS codes ONLY.
- Scoring levers (pb_clash, burial) -> SCORING-LOCKED codes ONLY.
- Never judge a sampling lever on a SCORING-LOCKED target or vice versa; that is the
  wrong-substrate error that produced two void gates.

## MEMETIC INTERLOCK — STILL DEADLOCKED (unchanged, still blocking)
config_parser.cpp:386-401 still requires FLEXAIDDS_WALL_PILOT_PASS=1, gated on the wall oracle
that cannot succeed (L3). Grok correctly refuses to set it. Re-key memetic to the REVISED
Phase 2 above (a real burial oracle on SCORING-LOCKED), or memetic is permanently unevaluable.

## WHAT GROK GOT RIGHT (do not re-litigate)
- BOOM_FRAC liveness: env DOES beat the claim-JSON 0.0, [BOOM] fires at frac 0.1, no CF~0
  wipeout signature. My L2 concern is resolved; the small-fraction caution stands.
- COARSE 64 vs 256: VOIDed its own earlier multi-variable result and re-ran a matched control
  (same binary, R=2, BOOM unset) -> FAIL, genuine 0/5 both arms. Exemplary discipline.
- CF gate now scores n=5 with 0 skips using --config + --ligand. Phase 0 tool gate PASS.

## SCOREBOARD OF VOID GATES (3 so far — all wrong-instrument, none wrong-execution)
1. Wall oracle: term structurally blind to the triggering values.
2. BOOM_INTERVAL pilot: knob not wired (fraction hardcoded 0.0).
3. pb_clash burial oracle: right lever, WRONG PANEL + no magnitude floor.
The liveness precheck caught #1 and #2. Add the magnitude floor and the class-matched panel
rule so it also catches #3's failure mode.
