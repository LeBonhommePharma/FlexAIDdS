# LANE D — WASTE AND SELECTION (smallest, most CERTAIN, fully offline)
BRANCH: `lane/d-selection-waste`      SEAT: Cursor Pro
DEPENDS ON: nothing.   MERGES FIRST.   MACHINE TIME: zero.

## HARD BOUND ON THIS LANE — read before you start
The ENTIRE remaining selection envelope is **15 targets**: those where a sub-2 A pose exists
in the pool but min-CF does not pick it. min-CF is 26/84 = 31.0%; the ceiling is 41/84 = 48.8%.
**If your combined number implies more than +15 targets, it is wrong.** And per R6 you must
report ONE combined re-election number, never a sum of per-defect deltas — min-CF over all
poses already contains the void-recovery and representative effects.

## THE DEFECTS, MEASURED
D1 — SENTINEL POSES. 4,284 of 35,873 emitted poses (11.9%) carry a sentinel CF (1e4-1e5) with
     all REMARK terms zeroed, from reconstruction failure (cf_bad.wal = 1e12) and grid escape
     (cf_oob.com = 99999) being EMITTED rather than dropped. The nearest-native sentinel sits
     at 1.59 A — a near-correct geometry, unelectable by construction.
D2 — CLUSTER REPRESENTATIVE. For 7 targets a sub-2 A pose is in the pool but the cluster's
     chosen representative is > 2 A: 1IGJ 2.04, 1K3U 2.98, 1N1M 2.64, 1N2J 2.04, 1SG0 2.07,
     1W1P 2.05, 2BYS 2.44. Fix: representative = min-CF member, not centroid/head.
D3 — ANY-FAILURE-WINS VOIDING. LIB/DatasetRunner.cpp:6466
       if (rp.ri == 0) ret = ri_ret; else if (ri_ret != 0) ret = ri_ret;
     One failed restart voids the whole target even when the other nine succeeded. 8/84 voided.
     Fix: elect over whichever restarts completed; void only if ALL fail.
D4 — THE ONE ANTI-WRONG-BASIN GUARD IS DISABLED IN ORACLE MODE. DatasetRunner.cpp:1266 gates
     cluster-spread demotion on `proto.cluster_spread_max > 0 && !oracle_mode`. The campaign
     ran oracle mode, so the guard that exists for exactly the big-but-wrong-basin failure was
     inert. Fix: decouple from oracle_mode.

## OWNERSHIP
YOU OWN, in LIB/DatasetRunner.cpp, THESE REGIONS ONLY:
    ~1005-1024   the election site (use_shannon_G, soft_T)
    ~1266   the cluster-spread guard
    ~6433-6466   the voiding propagation
    the pose-emission path that writes sentinel poses
FORBIDDEN, in the SAME FILE:
    ~2008, ~3249, ~3312  cofactor_blacklist / receptor construction  -> LANE B
    ~5949-6320           .inp emission + receptor_rotamer_prep       -> LANE E
FORBIDDEN OUTRIGHT:
    LIB/vcfunction.cpp (A/B), LIB/gaboom.cpp (C), LIB/read_input.cpp (E),
    LIB/SdfReader.cpp (A), LIB/config_parser.cpp (A/B)
You merge FIRST because you are smallest and fully offline-verifiable. Lane E rebases onto you.

## PRE-EXISTING — F5 NaN rank guard (do not re-implement)
`FLEXAIDDS_NAN_RANK_GUARD` in `LIB/gaboom.h` (QS_ASC / QS_DSC only) is already
in tree, DEFAULT OFF, parsed by `flexaids::env_bool`. PR #420 landed it gated.
You own the *question* of whether enabling it changes the 15-target selection
envelope — report that from `score_canonical.py`. Do not turn it on from other
lanes. Do not revert the `EnvFlags.h` include. `LIB/gaboom.cpp` stays FORBIDDEN
for D1–D4; this is not a license to edit gaboom.cpp.

## THE WORK
Fix D1-D4, each behind its own env gate, all default OFF:
    FLEXAIDDS_DROP_SENTINEL_POSES, FLEXAIDDS_REP_MINCF,
    FLEXAIDDS_VOID_ONLY_IF_ALL_FAIL, FLEXAIDDS_SPREAD_GUARD_ALWAYS

## ACCEPTANCE GATES (all offline, against the frozen benchmark)
  G1. Report the top-1 rate on N=84 BEFORE and AFTER, from score_canonical.py, spyRMSD.
  G2. ONE combined number from a SINGLE re-election pass. Not a sum. (R6)
  G3. Sanity: combined gain <= +15 targets. If not, your metric or denominator is wrong.
  G4. With all four gates OFF, score_canonical.py --frozen still prints 31.0% / 48.8%.
  G5. For D1, report how many of the 4,284 sentinels would be dropped and confirm no
      legitimate pose is dropped with them (the 1.59 A one is a sentinel, not a real pose —
      verify which side of your filter it falls on and say so).

## WHAT TO AVOID
  * Do NOT sum D1+D2+D3+D4 deltas. This is the single most repeated error in this project.
  * Do NOT re-diagnose the soft-beta T=300 election; it is settled and the fix is an env var
    (min-CF), which is already the scorer's baseline.
  * Do NOT enable FLEXAIDDS_CLEFT_SORT — it moves single-thread results (settled, PR #403
    superseded it).
  * Do NOT touch the receptor or .inp regions of DatasetRunner.cpp; they belong to B and E.
  * Do NOT "improve" clustering wholesale. Consensus re-ranking was measured WORSE (25.0%).
