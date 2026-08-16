# CAMPAIGN METHODOLOGY — FlexAIDdS Astex-85 genuine success-rate campaign

> **Landed on `main` 2026-07-25** from OPS handoff `CAMPAIGN_METHODOLOGY_for_Grok.md`.
> Companion hub: [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md).
> **Critical correction:** the 25.3% OPS session record **does not** prove `free_energy_strict` worked — that run predated the measured product default; see STEP 1 E10. **Not a published FlexAIDdS success rate** — unverified / no METHODOLOGY.md §0 receipt.

For Grok Build 4.5 (/loop + /goal). OPS-authored 2026-07-25. Follow in order; do not skip.
Repo: $FLEXAIDDS_ROOT (or git root)

## GOAL (single sentence)
Raise GENUINE Astex-85 top-1 (seed-echo excluded) from the measured 25.3% baseline toward the
FlexAID-2015 anchor of 45.2%, changing ONE variable at a time and never citing the raw success column.

## THE MEASURED STARTING POINT (do not re-derive; this is the reference)
Run v_autonomous_20260724_160919, autonomous blind, Workers=4, PRE PR#300/#301:
  GENUINE top-1        = 20/79 = 25.3%
  best-cluster <2A     = 22/79 = 27.8%   <-- SAMPLING CEILING
  election gap         = 2 targets only
  seed-echo            = 0 (clean)
  median rmsd_hungarian= 3.99 A
Archive: <iCloud>/FlexAIDdS_benchmarks/archived_from_ssd/archive_batch_20260725T095624Z/
         flexaidds_results/v_autonomous_20260724_160919/
WHAT THIS RUN DOES AND DOES NOT SHOW (read carefully):
  DOES: on THIS build the ceiling is 27.8% and elected is 25.3% -> only ~2 targets are lost to
  election, so the DOMINANT limiter is SAMPLING. Work sampling.
  DOES NOT: this run PREDATES commit 7b1aeb4a (acf -> free_energy_strict). A branch scan at launch
  time showed free_energy_strict_in_cluster=0 on every branch, so this binary ran the LEGACY acf
  election. Do NOT cite this run as evidence the election fix worked — the fix is still UNMEASURED.
  A post-7b1aeb4a run is required to quantify its effect (and it may move the 2-target gap either
  way). Treat 25.3% as a pre-merge reference point, not a post-fix result.

## NON-NEGOTIABLE RULES (violating any one invalidates the campaign)
R1  GENUINE METRIC ONLY. success = rank-0 elected pose, rmsd_hungarian < 2.0 A, seed_echo
    EXCLUDED, sentinel-guarded (rmsd >= 0). NEVER the raw success_rmsd column. (A regressed build
    once read 45.2% raw that was 100% seed-echo = zero genuine docking.)
R2  BOX: M3 Pro 18 GB. ONE science owner at a time. WORKERS=4 MAX (6 OOM-killed a campaign:
    7x "Killed: 9"). OMP_NUM_THREADS=1 per worker. NO dual full-85. NEVER `cmake --build` while a
    run holds the binary (has invalidated runs twice).
R3  ONE VARIABLE PER RUN. If two knobs change, the result is uninterpretable and must be discarded.
R4  RELAUNCH AFTER ANY SCORING-SURFACE MERGE. A run launched before a merge dockes pre-merge code
    and cannot baseline the next change.
R5  RECORD THE FULL ENV in provenance (every FLEXAIDDS_* var + matrix md5 + binary sha256 + git
    commit). A com-cap run was rendered UNCITABLE because its cap value was never written down.
R6  DETACH: launch under `setsid` (macOS has no setsid(1) -> use python double-fork os.setsid()).
    A driver died with the parent shell and lost a campaign.
R7  NO full-85 until the cheap gates pass (W0/W1). Six baselines were re-fired and none completed —
    the single largest compute waste in this effort.

## CLAIM PROTOCOL (only these count as citable rates)
matrix 9dc9 (md5 9dc93717dfed0698006d88dd6a9627bc) · GA 1000x2000 · R=10 restarts when citing ·
seed OFF · label defined-cleft vs autonomous explicitly. Anything else is a pilot, not a claim.

## ENVIRONMENT GOTCHAS (pre-solved; do not rediscover)
- System python3 has NO numpy -> both comparative CLIs traceback.
  Use: PY=~/.claude-science/conda/envs/python/bin/python
- That env has NO pytest -> use ~/.claude-science/conda/envs/cpp-python-core/bin/python
- Always: export PYTHONPATH=$FLEXAIDDS_ROOT (or git root)/python
- run_comparative_phases.py flag is --pipeline-dry (NOT --dry-run).
- probe_cf REQUIRES --config or it scores the WHOLE receptor (no LOCCLF prune) and inflates CF ~200x.
  PROVEN 1M2Z native: without --config cf_total=-24560.96 ; with --config cf_total=-117.74.
  probe_cf also needs --ligand <topology.sdf> when the pose is a PDB, else rc=256 / no [NATIVE_CF].
- Recursive `find` over the iCloud tree or ~/Documents/PhD HANGS on evicted stubs. Use bounded globs.
- ~/.claude and ~/Library/Application Support/Claude were moved to iCloud -> app/agent context hangs.
  Restore before a long campaign (rsync commands in MOVED_TO_ICLOUD_20260725T095624Z.txt).

## PHASE SEQUENCE (the repo's own pipeline enforces this; gates are fail-closed)
Driver:  PYTHONPATH=<repo>/python $PY scripts/run_comparative_phases.py --pipeline-dry
Gate:    PYTHONPATH=<repo>/python $PY scripts/comparative_phase_gate.py --dry-run
Tests:   PYTHONPATH=<repo>/python <cpp-python-core python> -m pytest python/tests/test_comparative_phases.py -q
Current gate state: P0 pass, P1 pass, P2 HOLD (live blocker), P3/P4 pending, next_allowed=P2.
P3 requires P2=pass; P4 requires P2=pass AND P3=pass. Do NOT force-pass (--force-p2-pass exists;
using it voids the claim).

### STEP 1 — W0.1 / E10: confirm election-vs-sampling on the frozen archive (OFFLINE, no dock)
  script: scripts/e10_election_vs_scoring.py  (currently stranded on branch 6ca30452 — MERGE IT FIRST)
  run:    $PY scripts/e10_election_vs_scoring.py --target-dir <run>/<TGT>/ [...] \
                --out-json e10.json --out-csv e10.csv --out-md e10.md
  NOTE: --target-dir is append; --campaign-dir is NOT (a second one overwrites the first).
  Input must have result.csv or *_0.pdb per target dir.
  DELIVERABLE: workorders/E10_election_vs_scoring.md
  ACCEPT: reports fraction of targets where an independent scorer ranks the near-native head above
  the elected decoy. EXPECTED (from the baseline): small — election is already fixed. If E10
  disagrees and shows a large election gap, STOP and reconcile before touching sampling.

### STEP 2 — W2 wall oracle: decide the steric-wall un-cap (SCORE-ONLY, no dock)

> **2026-07-25 correction (B3):** `FLEXAIDDS_WAL_COERCIVE` is **structurally unpassable**
> as a Voronoi-wall un-cap (per-pair cap; deep clashes not Voronoi-visible). Do **not**
> re-run WAL_COERCIVE expecting OFF≠ON. **Replacement burial/steric oracle:** one-variable
> `FLEXAIDDS_PB_CLASH_WEIGHT` via `scripts/pb_clash_burial_oracle.py` (production LOCCLF
> configs). See `workorders/PB_CLASH_ORACLE.md` and `WALL_ORACLE_FAIL_EXPLAINED.md`.
> Memetic / `WALL_PILOT_PASS` remain blocked until a **strong** burial oracle (non-trivial
> `cf_clash`) PASSes — not micro-ΔdCF alone.
>
> **STEP 3 BOOM note (B1):** claim path emits `boom_inject_fraction: 0.0` deliberately;
> inject requires `interval>0 && fraction>0`. Interval-only pilots are invalid BOOM tests.
> Env `FLEXAIDDS_BOOM_FRAC` overrides JSON (use 0.05–0.2 only; never 1.0 blind).

  WHY FIRST: WAL_CONTACT_CAP=50.0 (soft_wall.h:13) makes the wall SATURATE at exactly 1.000 A
  overlap while CF.com is unbounded below (per-atom floor -198.9) -> past 1 A, deeper burial is FREE.
  Memetic refinement is INTERLOCKED behind this (FLEXAIDDS_MEMETIC=1 is a no-op without
  FLEXAIDDS_WALL_PILOT_PASS=1) precisely to prevent burial walk-away.
  run: build/probe_cf --receptor <apo> --pose <native.sdf> --config <dock_config.json>   [flag OFF]
       then again with FLEXAIDDS_WAL_COERCIVE=1                                          [flag ON]
       over the 5 clean probes: 1J3J 1K3U 1L7F 1N1M 1M2Z  (NOT 1G9V — distal-heme/oracle-offset
       confound; hold it as chain-control only)
  ACCEPT: dCF = cf_native - cf_bestdecoy <= 0 on >= 7/8 of the panel with the flag ON where it fails
  OFF, AND no clean probe whose native is already CF-min gets worse.
  IF IT FAILS: STOP. Do not proceed to memetic. Re-diagnose.
  ALSO FIX FIRST: ops/gates/cf_gate_probe_cf.sh passes neither --config nor --ligand, so the gate
  currently SKIPS every PDB decoy (silent no-op) and would measure the 200x artifact. Patch both.

### STEP 3 — W1 serial pilot (6-9 targets, WORKERS=2)
  One knob vs the 25.3% baseline. Matrix 9dc9, seed OFF, genuine metric.
  Panel: the 5 clean probes + 2-4 of the election-gap near-misses
         (1OQ5 bcr=1.06, 1SQ5 1.12, 1YGC 1.24, 1YVF 1.70).
  ACCEPT: no clean probe regresses; the knob shows a directional effect on BCR or elected RMSD.

### STEP 4 — W3 sampling (the actual bottleneck; BCR raisers)
  Target: raise best_cluster<2A above 27.8%. Levers, one at a time:
   (a) memetic local refinement — ONLY after Step 2 PASS (trust-region <=1.5 A, refine against the
       steric/wall term, NOT full CF; refining against a defective objective makes top-1 WORSE)
   (b) niche metric: calc_rmsp (gaboom.cpp:4111-4123) is an unweighted RMS over a gene vector whose
       gene 0 is a cleft-grid ORDINAL (0..15131) mixed with degrees; sigma_share=204.19 means a
       0.375 A step in z leaves the niche while 7.9 A in y stays inside, and flipping ALL NINE
       angles 180 deg still counts as the SAME niche. Replace with a Cartesian pose distance.
   (c) coarse-init / restart pooling / GA diversity
  ACCEPT per lever: BCR up, no clean-probe regression, genuine top-1 not down.

### STEP 5 — W4 full-85 claim (ONLY after Steps 1-4 gates pass)
  Full 85, R=10, WORKERS=4, seed OFF, matrix 9dc9, detached via setsid, full env in provenance.
  ACCEPT: genuine >= 45% OR a documented, honest fail. Report genuine AND best-cluster AND the
  election gap. Never the raw column.

## EXPLICIT DO-NOT (each already cost real time)
- Do NOT re-fire the com-cap full-85 or softbeta@72d7 autonomous as "entropy validation".
- Keep COM_BURIAL_CAP OFF main (a per-optres cap does not bound TOTAL com: observed -1144 at 9x a
  -130 cap; and the -130 value was tuned to one target's native com = single-target over-fit).
- Do NOT elect by raw cluster frequency (the freq-1666 cluster on 1G9V is 10.4 A — picks WORSE).
- Do NOT trust any absolute CF from probe_cf run WITHOUT --config.
- Do NOT add memetic before the wall oracle passes.
- Do NOT re-verify the two-lever diagnosis; it is confirmed five times over.

## REPORTING CADENCE
At each gate: state the phase, the ONE variable changed, the genuine number, the BCR, the election
gap, and PASS/FAIL vs the stated acceptance test. If a run dies, say so and give the kill reason —
a dead run is not a partial result.