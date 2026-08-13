# FlexAIDdS SWARM — COMMON PREAMBLE  (v2, 2026-08-13)
# Paste this at the TOP of every agent's task, before its lane file.

## 1. The numbers you are trying to move
Astex Diverse, N=84, spyRMSD graph-automorphism, campaign astex84_dG_20260809_141245
(35,873 poses; 83 targets with poses).

    as-run soft-beta T=300 election    15/84 = 17.9%
    min-CF election (free, env var)    26/84 = 31.0%
    linear-reweight cap (LOTO-CV)              32.5%   <- proven ceiling of RESCORING
    pool ceiling at 10 restarts        41/84 = 48.8%   <- cap of ANY selection fix
    PUBLISHED BAR                              45.2%   (FlexAID-2015 top-1; top-10 66.7%)

Consequence: **ranking alone cannot reach the bar.** Every lane below targets either the
scoring INPUTS, the search, or waste — never the election alone (that is lane D, and it is
bounded at 15 targets).

## 2. NEW (measured 2026-08-13) — restarts are a real ceiling lever
30 restarts vs 10, on the 20 cheapest of the 43 ceiling-miss targets:
  * 4/19 converted to a sub-2 A ceiling (1HQ2 4.66->1.55, 1T40 2.23->0.92,
    1T46 2.39->1.19, 1S3V 2.90->1.74); 3 of those are also picked by min-CF
  * ceiling delta: mean +0.48 A, MEDIAN +0.02 A, max +3.11 A
  * 11/19 targets gained NOTHING (delta <= 0.05 A)
  * 1UNL went from ZERO poses to 1,400 (ceiling 4.38 A — recovered, still a failure)
The effect is CONCENTRATED, not diffuse, and costs 3x machine time. Do not treat it as a
general fix; it is the reason the ceiling has headroom at all.

## 3. THE FROZEN BENCHMARK IS THE UNIT OF EXPERIMENT
ASTEX84_FROZEN_POSE_BENCHMARK.csv
  artifact 3cc422aa-6ba6-4691-a060-3f705fd63c14  version 559a2075-6a6f-4429-b871-f1ac86ec0192
  local copy: /Users/lp.more/flexaidds_results/workorders/ASTEX84_FROZEN_POSE_BENCHMARK.csv
One row per pose: target, pose_path, restart, rmsd_spyrmsd, n_heavy, and every CF term
(CF_total, com, sas, wal, hbond, metal, elec, gist, con). Testing a scoring or election
hypothesis against it takes SECONDS. Use it before you touch C++.

## 4. THE ONE SCORER — no exceptions
    python3 /Users/lp.more/flexaidds_results/workorders/score_canonical.py \
        --frozen /Users/lp.more/flexaidds_results/workorders/ASTEX84_FROZEN_POSE_BENCHMARK.csv
    python3 .../score_canonical.py --run /Users/lp.more/flexaidds_results/<run_dir>
It must print 26/84 = 31.0% and 41/84 = 48.8% on the frozen file. If it does not, STOP —
your checkout or the file is wrong. A number that did not come out of this script does not
enter any report, PR description, or message.

## 5. HARD RULES — breaking any of these makes your work unusable
  R1. DO NOT LAUNCH A DOCKING CAMPAIGN. Claude Science owns the box and holds the live run.
      A second campaign on this Mac corrupts both — it has already happened twice.
      Single-target probes for a gate are fine. A dataset run is not.
  R2. ONE METRIC: spyRMSD graph automorphism from the crystal SDF bond block, no
      superposition. NEVER use result.csv "rmsd_to_crystal" (ordered, over-strict) or
      "rmsd_hungarian" (element-blocked, over-permissive — it inflated a measured ceiling
      from 48.8% to 57.8%). They disagree by up to 5.9 A on one target.
  R3. ONE DENOMINATOR: N=84. 2HR7 is excluded (its ligand is PEG / CCD P33; the deposit is
      an apo insulin-receptor ectodomain with no cognate ligand). Report N=85 alongside if
      you must, never instead.
  R4. result.csv "success" means "docking ran", NOT "docking succeeded". Never gate on it.
  R5. SHIP EVERY FEATURE ENV-GATED **OFF**. After your change, the scorer must still print
      31.0% / 48.8% on the frozen set with your flag unset. A PR that moves a DEFAULT is an
      automatic reject — a knob documented as OFF being ON is exactly what cost 13 points.
  R6. NEVER SUM GAINS ACROSS FIXES. min-CF over all poses ALREADY contains void-recovery and
      cluster-representative effects. Re-elect ONCE, report ONE combined number. This
      double-count has been made three times in this project, twice in one night.
  R7. MEASURE, DO NOT INFER. Reading source and asserting runtime behaviour is the single
      most expensive mistake made here. Every claim needs a command you ran and a number you
      saw. "The code says X" is not evidence that X happened.
  R8. STAY IN YOUR LANE'S FILES. Your lane file lists OWNED and FORBIDDEN paths. If you need
      a forbidden file, write the required diff into your PR description and say so — do not
      edit it.
  R9. One branch per lane, named in your lane file. Do NOT merge. Do not rebase another
      lane's branch.
 R10. If a gate cannot be met, say so plainly and hand back. A lane that reports "could not
      measure" is useful; a lane that reports an unmeasured number is worse than silence.

## 6. Paths
  source   : /Users/lp.more/Projects/FlexAIDdS        (main @ aa15464e)
  campaign : /Users/lp.more/flexaidds_results/astex84_dG_20260809_141245/
  deep-restart : /Users/lp.more/flexaidds_results/deeprestart_20260813_000719/
  inputs   : /Users/lp.more/flexaidds_results/cache_v2/astex_diverse/<PDB>/
  sites    : /Users/lp.more/flexaidds_results/astex85_sites_clean/<PDB>/
  workorders (shared drop box): /Users/lp.more/flexaidds_results/workorders/
  score-only tool: build/probe_cf   (--receptor --pose --ligand --config --mode --pdb)
  cmake is NOT on PATH: use /opt/homebrew/bin/cmake

  probe_cf TRAP: without --config it runs the DEFAULT engine config, the LOCCLF sphere-prune
  never applies, the WHOLE receptor enters optres, and CF.com inflates ~200x (measured:
  -24,560 without --config vs -117.7 with, same pose, same binary). ALWAYS pass --config.

## 7. Already REFUTED — do not re-derive or propose these
  * GA budget / early stopping is not a lever (90% of restarts stop early; success-median
    719 vs fail 649 generations, p=0.053).
  * Search box volume is not inflated (volume vs success r=-0.10, p=0.37).
  * No RNG race; no REMARK truncation (66-line/2138-byte blocks vs a 5000 limit); no NaN in
    the partition function.
  * Consensus re-ranking is WORSE than min-CF (25.0%).
  * Ligand or receptor size does not predict failure (p=0.24 / p=0.11).
  * Linear reweighting of the four live CF terms cannot pass 32.5% (leave-one-target-out CV).
  * A rescoring bisect to the 45.2% engine is IMPOSSIBLE: repo history starts 2026-05-21;
    the FlexAID-2015 codebase is not in it.
  * Cleft-grid nondeterminism is FIXED (PR #403, bucket merge). Do not re-diagnose it.
  * FLEXAIDDS_CLEFT_SORT is NOT the fix — it imposes a third canonical ordering and moves
    single-thread results. Do not enable it.

## 8. How to report
Post a drop into /Users/lp.more/flexaidds_results/workorders/ named <LANE>_RESULT.md with:
  1. what you measured, and the exact command
  2. the scorer's output block, verbatim
  3. root cause with file:line AND the runtime observation that confirms it
  4. gain or cost in points on N=84, or "could not measure" and why
  5. what you did NOT do and why
Then open the PR. Claude Science re-verifies every number against the frozen benchmark
before it becomes a project number.
