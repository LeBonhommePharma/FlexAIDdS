# FlexAIDdS Reviewer Benchmark Analysis: 1G9V & 1GM8
Date: 2026-06-24
Location: /Users/lp.more/FlexAIDdS_reviewer_benchmark/
Source of run: reproduce_astex85.sh with FLEXAIDDS_* exports (RESTARTS=7, THERMO=1, T_EFF=0.596, NATIVE_SEED_FRAC=0.90, RECEPTOR_ROTAMER_PREP=1, CONSENSUS=1)
Current success observed in this data: task states 43/84 (genuine low vs claimed 80/85)

## Dock Config Excerpt (common, matches methodology)
From 1G9V/dock_config.json and 1GM8/:
{
  "flexibility": { ..., "receptor_rotamer_prep": true },
  "reference_ligand": { "seed_fraction": 0.9, ... },
  "thermo_engine": {"enabled": true, "T_eff": 0.596, "tencom_scale": 1}
}
See reproduce_astex85.sh provenance section for full env: RESTARTS=7, PARALLEL=1, NATIVE_SEED_FRAC=0.90 etc.
benchmark_astex_native_85.json used (oracle self-dock pairs).

## result.csv for targets (reported vs best_cluster)
1G9V/result.csv:
pdb_id,best_score,rmsd_to_crystal,rmsd_hungarian,...,success,cf_native,best_cluster_rmsd,best_cluster_idx,...,pose_source,...,G_bind,...
1G9V,-95.9867,3.1234,2.9474,...,0,-51.6580,0.7210,1,...,ga_cluster,...,-38.4754,...

1GM8/result.csv:
1GM8,-98.1847,2.5596,2.4121,...,0,-59.5929,0.4521,1,0,ga_cluster,...,-39.2942,...

Observations:
- reported rmsd_hungarian >2.0 (fail, success=0)
- best_cluster_rmsd <<2.0 (0.72 and 0.45) -- poses capable of success WERE generated
- G_bind column matches the *top-level* (r0) [THERMO] not min-G
- pose_source=ga_cluster (not ini)
- best_cluster_idx=1

This indicates two-stage selector (if active) did not elect a low-RMSD pose for the reported metrics.

## [THERMO] G_bind per restart (from **/stdout.log)
For 1G9V (min G = most negative = best):
- r1/stdout.log: G_bind=-26.856604
- r2/stdout.log: G_bind=-40.133560
- r4/stdout.log: G_bind=-57.049618   <--- MIN G (best by far)
- r5/stdout.log: G_bind=-32.540318
- r6/stdout.log: G_bind=-40.499645
- top stdout.log (used for csv): G_bind=-38.475361

r3/stdout.log: NO [THERMO] line (grep returned 0 matches; file contains \0 bytes at offset ~186k -- possible truncated/resume corruption?)

For 1GM8:
- top: -39.294205
- r1: -48.386234
- r2: -53.029716   <--- MIN
- r3: -47.356712
- r4: -51.371246
- r5: -44.632332
- r6: -48.614330

Pattern identical: csv G_bind = r0/top's value, not the min-G restart.

## Selector Code Evidence (from source)
LIB/ReportedPoseSelector.cpp implements:
- build_cross_restart_pool(prefixes): for each ri, parse G from (prefix + "/stdout.log" or parent/stdout.log fallback), attach same g_bind to all its PoseCandidates (_N.pdb + INI).
- elect_reported_pose(pool, thermo_on): 
  if (thermo_on):
    // TWO-STAGE (v88/2015 JCIM):
    // (a) min finite G_bind across pool -> chosen_ri
    // (b) within chosen_ri only, max boltzmann_composite (Z+H * pop^w)
    // Then return that path for rmsd_hungarian computation.
  Note: no console output of chosen_ri / elected path / min_g value.
  parse_g_bind_from_log: rfind("[THERMO]") + find G_bind=

In LIB/DatasetRunner.cpp:
- all_prefixes built as ri=0: out_prefix (top), ri>=1: out_dir/r{ri}/<pdb>
- thermo_on = result.has_thermo || env FLEXAIDDS_THERMO
- best_pose_pdb = elect...(pool, thermo_on)  [or fallback]
- THEN: rmsd_hungarian computed on best_pose_pdb
- SEPARATELY: best_cluster_rmsd = min over ALL pool candidates' hungarian_rmsd
- [THERMO] parse for result.thermo_G_bind (written to csv G_bind column) happens on *main stdout_path* (top level log) ONLY -- last [THERMO] seen. Not from elected restart.
- !thermo_on branch does consensus re-rank; for thermo, two-stage is "terminal"

Csv write (approx 6821+):
G_bind column <- result.thermo_G_bind (always r0's value)

## r*/ and top _N.pdb inspection
- Top level and each rN/ contain their own set of <pdb>_<0-9>.pdb + _INI.pdb + .mcf + stdout/stderr + dock_config
- REMARK CF=... present; some have "Frequency: N" and "Cluster X: ... Frequency: N"
- Example from 1G9V top/1G9V_1.pdb (high freq candidate?): Frequency:108 , CF~-73
- r4 (minG) 1G9V_1.pdb: CF~-77 , Frequency:1
- No "elected path" or "selected restart" logged in stdout (no prints in elect fn or thermo path).
- No _N.pdb special "reported" copy visible; selection is transient for RMSD calc + csv row.

## Why reported RMSD high?
- If two-stage active + parse succeeded for r4's G=-57, should have restricted to r4 poses and picked highest Z+H composite within r4.
- But reported rmsd=~2.95 matches likely a top/r0 elected pose (G~-38).
- Likely causes (no direct elected-path log to confirm):
  1. G_bind parse in build_cross... yielded non-finite for the best restarts (e.g. r3 had no [THERMO] at all; resume may have mangled logs; rfind logic picked wrong last line).
  2. chosen_ri remained -1 -> fallback to full-pool Z+H composite (non-thermo path) which picked a bad-RMSD high-freq? pose.
  3. Even within min-G restart, the freq-gated Z+H representative was *not* the low-RMSD one (best_cluster=0.72 was a different member/pose; selector prefers thermo-basin over oracle RMSD).
  4. G_bind column in result.csv decoupled from elected restart (always r0).
  5. Resume artifacts (resume_driver.log, before_resume json, binary-ish logs) may have caused inconsistent logs/poses at selection time.
- Note best_cluster_rmsd scans *same* all_prefixes pool as elect, so "search found good pose" but selector did not report it.

## Recommendations (do not claim done)
- Add debug prints in elect_reported_pose / build: e.g. fprintf for chosen_ri, min_g, elected_path, per-restart g_binds seen.
- Ensure G_bind column in per-target result.csv reflects the *elected restart's* [THERMO] G (tie to best_pose_pdb's restart).
- Verify parse_g_bind_from_log robustness on resumed/truncated logs (handle multiple [THERMO], NaN, last vs first).
- Consider logging the elected restart's full thermo line when thermo_on.
- Cross-check during run: after elect, parse G from the *specific elected log* and store in result for csv.
- Validate selector on these exact cases by unit test replay (see tests/test_statmech.cpp for existing ReportedPoseSelector tests).
- Reproduce with sequential (no PARALLEL) + no resume to rule out log interleaving.
- Confirm v88 G_bind formula details (see thermo_engine, T_EFF=0.596 usage in statmech/BindingMode) vs what is emitted in [THERMO].
- After fixes, re-run targeted on 1G9V/1GM8 (and known failures) and verify min-G restart's elected pose has sub-2.0 hungarian + csv G_bind matches it.
- Check emit_aggregate_from_run_trees.py and validate_benchmark for whether they override reported RMSD.

## Evidence files / excerpts saved
- This report: reviewer_benchmark_findings_1g9v_1gm8.md (relative to worktree root)
- Raw:
  /Users/lp.more/FlexAIDdS_reviewer_benchmark/1G9V/result.csv
  /Users/lp.more/FlexAIDdS_reviewer_benchmark/1G9V/stdout.log (G_bind=-38.47)
  /Users/lp.more/FlexAIDdS_reviewer_benchmark/1G9V/r4/stdout.log (G_bind=-57.05)
  Similar for 1GM8, r*/1G9V_*.pdb (REMARK CF/Frequency)
  reproduce_astex85.sh (exports + 7 restarts)
  LIB/ReportedPoseSelector.{h,cpp} (two-stage impl)
  LIB/DatasetRunner.cpp (call site, all_prefixes= rN, G_bind csv from main log only, best vs reported split)
- No dedicated "SCRATCH" dir found at time of analysis; this md + source reads serve as captured evidence. Consider mkdir -p ~/FlexAIDdS_reviewer_benchmark/SCRATCH or worktree/scratch/ for future.

Do not claim completion or high success restored. Analysis only via file inspection (list_dir / read_file / grep on absolute paths to reviewer_bench + relative in worktree for source).
