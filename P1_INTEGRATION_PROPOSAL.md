# P1 Core C++ Integration — Survey + Chunked Implementation Proposal

**Agent**: Grok Build subagent (PARALLEL TRACK, read-only survey + proposal only)
**Date**: 2026-07-08
**Scope**: Accelerate P1 by producing precise, chunked proposal. NO edits to sources, NO builds/tests executed (respects P0 build-p0-gpf/ isolation + AGENTS.md read-only constraints for this track). All info from targeted `read_file` + `grep` + `list_dir` on workspace root.
**Worktree**: Operated in `/Users/lp.more/Projects/flexaidds-gpf-grandcanon` (current cwd).
**Inputs read (minimal exploration)**: GPF_IMPLEMENTATION_PLAN.md (P1 sections), AGENTS.md (core rules), LIB/top.cpp, LIB/ParallelDock.{h,cpp}, LIB/ParallelCampaign.{h,cpp}, LIB/DatasetRunner.{h,cpp}, LIB/BindingMode.{h,cpp}, LIB/GrandPartitionFunction.{h,cpp}, LIB/TargetServer.{h,cpp}, LIB/statmech.h (key), LIB/gaboom.h (cluster decls), LIB/FastOPTICS_cluster.cpp, benchmarks/datasets/competition_example.yaml, various targeted greps for usage/conc/CCBM.
**Cross-refs validated**:
- Existing GPF/TargetServer usage ONLY in DatasetRunner.cpp (partial, approx-based).
- `BindingPopulation::get_global_ensemble()` (returns `statmech::StatMechEngine`) at BindingMode.cpp:132 (impl), BindingMode.h:194 (decl); used also in MultiModelDock.cpp:186, tests.
- Ensemble builds from `pose.total_energy()` (CCBM-aware: `CF + receptor_strain`).
- `StatMechEngine::compute().log_Z` (Thermodynamics struct: statmech.h:38).
- GPF `add_ligand(name, log_Z, conc_M=1.0)` + engine overload (GrandPartitionFunction.h:69,73); `add_or_overwrite` (used by TargetServer).
- TargetServer: `create_session(name) -> DockingSession`, `register_result(sess)` (which does `grand_xi_.add_or_overwrite(name, sess.log_Z)` ignoring conc for now; DockingSession lacks conc/best_center/conformer_pop yet).
- No usage yet in top/Parallel*/Binding output paths.
- Conc sources: defaults=1.0M (GPF.h:42 `c_standard`, add_* defaults); competition_example.yaml (per-ligand `conc_M` under `competition_sets[*].ligands[*]`); NO support yet in config_parser.cpp/read_input.cpp/ga.inp/JSON/TargetConfig/DatasetEntry. Temp from TargetConfig/FA->temperature.
- CCBM note: MUST use `total_energy()` (Pose.h:49, BindingMode.cpp:270/392, BindingPopulation.cpp:144) for log_Z ensembles — NOT raw CF or predicted_dG. receptor_strain from multi-model. get_global_ensemble explicitly avoids double-weighting.
- DatasetRunner approx (lines ~6913-6915): `sess.log_Z = -result.predicted_dG / (kB * T)` (from dG proxy, not ensemble); TODOs for best_center/conformer_pop from BindingMode.
- Post-GA: top.cpp ~2504 (post_engine on chrom_snapshot), ~2570 clustering, ~2589 (FastOPTICS/DP/CF).
- Cluster entry: gaboom.h:272-274; impls create local BindingPopulation (FastOPTICS_cluster.cpp:123+); output_Population inside.
- Parallel paths: ParallelDock has own get_global_engine() (but plan specifies BindingPop hook); ParallelCampaign per-ligand loop (surrogate dG now).

**Todo tracking (internal for this survey, per AGENTS)**:
- [x] Read plan + AGENTS (P1 sections prioritized)
- [x] Targeted reads/greps on must-read files (top, Parallel*, DatasetRunner, Binding*)
- [x] Cross-ref GPF usage in DatasetRunner + BindingPop get_global
- [x] Identify conc sources + CCBM total_energy
- [x] Locate exact hook points + non-breaking accessor spots
- [x] Produce chunked proposal + locations + sketches + risks + tests
- [x] Write ONLY P1_INTEGRATION_PROPOSAL.md (small, conventional)

**P1 Gate (from plan)**: multi-ligand synthetic test produces correct pA/pB/selectivity; ctest green; single-ligand outputs unchanged.

---

## Chunked Implementation Proposal (P1 only; post-P0 gate)

Follow AGENTS.md strictly on real track (todo_write one-at-a-time, fresh build+ctest, commit immediate per logical chunk with `Fix:/Add:`, zero failures, no batch, hygiene, inspect-first). This proposal is read-only output to enable fast wiring.

**High-level strategy**:
- Add *minimal non-breaking accessors* first (get_log_Z, get_canonical_engine, perhaps get_best_center-ish, get_conformer_populations).
- Wire TargetServer (or lightweight GPF) *conditionally* ("if active") in post-GA/post-cluster paths + Parallel* result aggregation loops.
- Use `BindingPopulation::get_global_ensemble().compute().log_Z` (or BindingMode per-mode if needed) + engine overload.
- Default conc=1.0M for P1 (context/wiring of real concs = P3).
- Harden DatasetRunner to use real (when surfaced) + populate TODO fields (may require result extensions or REMARK parsing).
- Optional REMARK augmentation only when GPF active (non-breaking).
- Preserve CF/GA ranking, use total_energy for CCBM, log-space, single-ligand unchanged.
- Test gates: synthetic 2-ligand (known Z ratio + c) via direct GPF or small harness; ctest -R 'grand|Target|binding_mode_statmech'; single-ligand regression.
- Chunks sized for immediate commit+verify (no monolithic).

**Dependencies on P0**: P0 adds CMake FLEXAIDS_GRAND_CANONICAL (default ON) + ensures compilation + test hardening. P1 assumes Grand/Target/MultiSite always available in test builds.

---

## Exact Files + Change Locations + Pseudocode/Diff Sketches

1. **LIB/BindingMode.h** (add minimal non-breaking public accessors to BindingMode + BindingPopulation; CCBM note in comments)
   - Location: After existing statmech APIs (~line 90-110 for BindingMode; ~194 for Population).
   - Sketch (add):
     ```cpp
     // In class BindingMode (public):
     double get_log_Z() const;                    // convenience: get_thermodynamics().log_Z (post-corrections)
     const statmech::StatMechEngine& get_canonical_engine() const; // or return copy? view via ref (const& to engine_)
     // Optional for P1 harden:
     // std::array<float,3> get_best_center() const; // from Rep pose or compute
     // std::vector<double> get_conformer_populations() const; // if CCBM

     // In class BindingPopulation (public, after get_global_ensemble):
     double get_log_Z() const { auto e = get_global_ensemble(); return e.compute().log_Z; }
     statmech::StatMechEngine get_canonical_engine() const { return get_global_ensemble(); } // or ref if cached
     ```
   - Non-breaking: new methods only; no sig changes to existing.
   - Impl in .cpp (see below).

2. **LIB/BindingMode.cpp** (impls + ensure total_energy usage documented)
   - Locations:
     - BindingMode::get_log_Z (new, after get_heat_capacity ~351): `return get_thermodynamics().log_Z;`
     - BindingMode::get_canonical_engine (new): `rebuild_engine(); return engine_;` (or const ref)
     - In get_global_ensemble (already correct @132): comment "Uses total_energy() for CCBM (CF + strain)".
     - Optional: augment output_BindingMode (~518) for grand if context:
       ```cpp
       // After existing REMARKs, optional:
       // if (Population && Population->target_server) {
       //   auto& gpf = ...; snprintf(..., "REMARK Grand Xi=%.4g p_bound=%.4f\n", ...);
       // }
       ```
   - CCBM: all paths (rebuild_engine:270, get_thermo etc.) already use `pose.total_energy()` — preserve.

3. **LIB/gaboom.h** (extend cluster decls non-breaking for TargetServer hook)
   - Location: lines 272-274.
   - Sketch (add defaulted param at end for backward):
     ```cpp
     void cluster(..., char* gainp, target::TargetServer* ts = nullptr, const std::string& ligand_name = "");
     // same for DensityPeak_cluster, FastOPTICS_cluster
     ```
   - Or (preferred minimal): keep sigs, add internal context via FA_Global extension or separate wiring. (Alternative: surface BindingPopulation* from clusters via out-param.)

4. **LIB/FastOPTICS_cluster.cpp** (and cluster.cpp, DensityPeak_Cluster.cpp) — inside post-populate
   - Locations: After `BindingPopulation Population1(...)`; after `Population1.output_Population(...)` (~170).
   - Sketch (if ts param added):
     ```cpp
     if (ts && !ligand_name.empty()) {
         auto sess = ts->create_session(ligand_name);
         sess.completed = true;
         sess.n_poses = ...; // from pop
         auto eng = Population1.get_global_ensemble();  // or best mode
         sess.log_Z = eng.compute().log_Z;
         // populate best_center from Rep, conformer_pop if CCBM
         ts->register_result(sess);
     }
     ```
   - For 3 populations (current code), decide policy (e.g. use Population1 primary, or aggregate?); document. Use global_ensemble per plan.
   - Call sites (top.cpp) will pass nullptr for now (non-breaking).

5. **LIB/top.cpp** (core post-GA / post-clustering + Parallel paths)
   - Locations:
     - Includes: add `#include "TargetServer.h"` after statmech.h (~11).
     - ParallelDock block (~2335): after `auto global_thermo = pdm.aggregate();` (or after get_best), if ts active: `double lz = ...; ts->...` (note: ParallelDock path may need BindingPop construction or use its engine; plan calls for BindingPop hook — may require minor bridge or post-cluster).
     - Campaign (~2428): in result loop or after run_campaign, feed per-ligand.
     - Standard GA + post (~2504-2510 post_engine; after clustering ~2570-2589):
       ```cpp
       target::TargetServer* ts = ...; // from context (see below)
       std::string lig_name = derive_from_inputs(dockinp); // or FA
       if (ts) {
           // post cluster: but pop local inside cluster fn — rely on cluster hook or
           // construct temp BindingPop here from chrom_snapshot? (less ideal)
           // Preferred: after cluster, if results written, or thread ts down.
       }
       // Post-GA thermo already builds engine — can feed:
       // if (ts) ts->...add_or_overwrite(lig_name, post_thermo.log_Z); but prefer full BindingPop for consistency.
       ```
     - Context creation (new, minimal):
       ```cpp
       // Near FA init or before GA (single-receptor for main top path):
       std::unique_ptr<target::TargetServer> local_ts;
       if (/* multi-ligand or env or flag */) {
           target::TargetConfig tcfg{ .temperature_K = FA->temperature };
           local_ts = std::make_unique<target::TargetServer>(tcfg);
           // validate etc.
       }
       target::TargetServer* active_ts = local_ts.get();
       // Pass active_ts + ligand_name down to cluster calls and Parallel*.
       ```
     - After cluster: optional REMARK or post-process for Xi if active_ts.
   - For campaign: group by receptor (like DatasetRunner).

6. **LIB/ParallelDock.cpp / .h**
   - Locations: aggregate() ~387, get_global_engine ~373; after run in top.
   - Sketch: expose or use BindingPop? (ParallelDock currently bypasses full BindingPop clustering for speed; proposal: after aggregate, or add optional ts to ParallelDockManager ctor/run, then in result collection: `if (ts) { auto e = ...; ts->grand... or create/register with log_Z = e.compute().log_Z; }`
   - Add to RegionResult or manager: pass-through for ligand context.
   - Non-breaking: new optional fields/params.

7. **LIB/ParallelCampaign.cpp / .h**
   - Locations: per-ligand loop ~240 (`for (int li=0; ...)`), after consensus ~315; result struct.
   - Sketch (in loop, after success):
     ```cpp
     if (shared_ts /* per-receptor */) {
         auto sess = shared_ts->create_session(lr.name);
         sess.completed = true;
         // For now (P1): if full GA later, use real; currently use from dG_consensus approx or skip until P1 full
         // sess.log_Z = -lr.dG_consensus / (kB*T); // temp, replace with ensemble
         // Later: lr.log_Z = eng.compute().log_Z from per-ligand BindingPop
         shared_ts->register_result(sess);
     }
     ```
   - CampaignSummary + LigandResult: add optional `double log_Z = 0; double p_bound=0;` etc (P2 surface).
   - Create shared TargetServer in run_campaign if receptor shared (config has single receptor).

8. **LIB/DatasetRunner.cpp / .h** (harden existing partial usage)
   - Locations:
     - Creation ~5076: already has tcfg (add conc later in P3).
     - Per-ligand create_session ~5199.
     - Register ~6905-6920: **replace approx**:
       ```cpp
       // BEFORE (approx):
       // sess.log_Z = -static_cast<double>(result.predicted_dG) / (statmech::kB_kcal * T);
       // AFTER (P1 harden, when real log_Z available in result or parsed):
       if (have_real_log_Z_from_ensemble) {
           sess.log_Z = result.log_Z_from_binding_pop;  // new field or parse
       } else {
           // keep approx temporarily or compute from F if F== -kT logZ + corrections
       }
       // Populate:
       sess.best_center[0] = ...; // from parsed best pose CoM or REMARK
       if (!result.conformer_populations.empty()) sess.conformer_populations = ...;
       ```
     - Cross analysis ~7037: already uses rank etc — good.
     - Add to DockingResult (h) + parsing: `double ensemble_log_Z = 0.0; std::vector<double> conformer_pops; float best_cx,cy,cz;`
     - Parsing sites: after CSV/stdout/pose parse (~6350+), surface from child output (P1 will make binary emit "REMARK Ensemble log_Z=..." or CSV field when ts active).
     - Also populate in autonomous vs subprocess paths.
   - Cross-ref: keep TargetServer grouping by receptor_path.

9. **LIB/TargetServer.h / .cpp** (minor non-breaking if needed for conc/best_center)
   - Locations: DockingSession struct (~40): add `double concentration_M = 1.0;`
   - In register_result (~50): pass conc to add_or_overwrite if set:
     ```cpp
     grand_xi_.add_or_overwrite(session.ligand_name, session.log_Z, session.concentration_M);
     ```
   - create_session can take/ set default conc.
   - Optional: `void set_default_concentration(double c);` or per-session.
   - In ctor of ts: temp from config.
   - (Keep thread-safe.)

10. **Other supporting (small)**:
    - LIB/statmech.h (if needed): nothing; or add `double log_Z() const { return compute().log_Z; }` convenience (non-breaking).
    - cmake / CMakeLists.txt: nothing for P1 (P0 handles); ensure no new sources.
    - tests/test_grand_integration.cpp (new? or extend test_grand_partition.cpp + test_target_server.cpp + test_binding_mode_statmech.cpp): synthetic 2-ligand using direct BindingPop construction + TargetServer.
    - No change to gaboom.cpp / Vcontacts / core GA.
    - For output: optional in BindingMode::output_BindingMode or top post-write.

**Pseudocode for main hook (post-cluster in cluster path or top)**:
```cpp
// After BindingPopulation pop = ...; pop.output_Population(...);
if (active_ts && !ligand_name.empty()) {
    auto sess = active_ts->create_session(ligand_name);
    sess.completed = true;
    sess.n_poses = pop.get_Population_size(); // or count
    auto eng = pop.get_global_ensemble();  // CCBM total_energy inside
    auto td = eng.compute();
    sess.log_Z = td.log_Z;
    // best_center, conformers from pop.get_binding_mode(0) or Rep
    active_ts->register_result(sess);
    // optional: double p = active_ts->binding_probability(lig_name);
    // emit REMARK "Grand log_Xi=... p_bound=..." to remark buf before write_pdb
}
```

**Concentrations in P1**: Always default 1.0 unless context passed (future: from FA or YAML entry). Warn on extremes per GPF.

---

## Risks + Mitigations (from survey + AGENTS guardrails)

- **Signature changes to cluster funcs**: Could affect other callers (FastOPTICS etc). Mit: default=nullptr overloads; keep existing  sigs + internal dispatch.
- **Subprocess in DatasetRunner**: Real log_Z/BindingPop data not directly available to parent. Mit: P1 core emits extra in stdout/CSV/REMARK (e.g. "Ensemble log_Z=XX"); parser hardens to prefer it. Or in-process path for benchmarks.
- **Multiple Populations (FO code creates 3)**: Which log_Z? Mit: primary (Pop1), or merge via StatMechEngine::merge; document choice. Single-ligand unchanged.
- **ParallelDock bypasses BindingPop**: Uses own engine on snapshots. Mit: either construct BindingPop post or use Parallel's engine for log_Z (plan prefers BindingPop — propose bridge or note).
- **CCBM + total_energy mismatch**: Old approx used predicted_dG (may ignore strain). Mit: proposal insists on get_global_ensemble(); add asserts in tests.
- **Numerical/empty ensemble**: GPF already guards (size>0); BindingPop returns engine even if empty. Mit: use existing checks.
- **Concurrency in Parallel/Campaign/Dataset**: TargetServer create/register are thread-safe (per header). Mit: verify in P1 tests.
- **No conc yet**: P1 uses 1M; p(empty)/selectivity only meaningful for >=2 ligands. Single unchanged.
- **REMARK schema drift**: Optional only; keep under feature or env.
- **Build overlap**: Proposal only — real impl uses separate worktree/branch per plan.
- **Overclaim thermo**: Stick to "ensemble log_Z", "GPF-derived p_bind" etc. per scientific guardrails.

---

## Test Ideas (for real P1 impl; synthetic first)

- Extend `tests/test_grand_partition.cpp` or new `tests/test_grand_integration.cpp`: 
  - Construct minimal BindingPopulation (or mock poses with total_energy), get log_Z, feed to TargetServer + GPF for 2 ligands + known concs; assert pA/pB, selectivity, log_Xi match hand logsumexp.
  - CCBM case: poses with + receptor_strain; verify != pure-CF log_Z.
- In `test_binding_mode_statmech.cpp` / `test_target_server.cpp`: roundtrip log_Z from pop -> ts -> gpf.binding_probability.
- Synthetic multi-ligand (use benchmarks/grand_synthetic/*.json fixtures for Z+c + expected).
- Regression: single-ligand run (no ts or 1 ligand) produces identical PDB/CSV/RRD/stdout (diff before/after).
- DatasetRunner harden: feed known log_Z via result, check cross_ligand_results.ranked_ligands + _cross_ligand.md.
- Edge: empty ensemble (log_Z=0 or -inf handling), c=0, identical Z, extreme ratios (log-space).
- HW parity note (from plan P4): feed same Z (synthetic) to GPF under scalar/Metal builds; identical grand obs.
- ctest: `-R 'grand|Grand|Target|binding_mode_statmech|statmech|multi_site' --output-on-failure`
- Manual: small 2-ligand same-receptor; inspect Ts->grand_partition().rank() + p_bound.
- Future (P2+): Python roundtrips.

**Verification commands (to run in real P1 after each chunk)**:
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release ...
cmake --build build -j
ctest --test-dir build --output-on-failure -R 'grand|Target|binding_mode'
python3 scripts/check_repo_hygiene.py
# (no full benchmark runs in P1 gate)
```

---

## Summary Bullet List (Files + Locations)

- **BindingMode.h:194** (Population decls) + ~90 (Mode): add get_log_Z(), get_canonical_engine().
- **BindingMode.cpp:132** (global_ensemble), ~299 (get_thermo), ~444 (output_BindingMode): impl + optional REMARK.
- **gaboom.h:272-274**: optional defaulted ts/ligand_name params on cluster* decls.
- **FastOPTICS_cluster.cpp:123,170** (and siblings): post-pop hook to register if ts.
- **top.cpp:11** (include), ~2335 (PD), ~2428 (campaign), ~2504 (post-GA), ~2570 (cluster call site), ~2600 (context creation).
- **ParallelDock.{h,cpp}:~373 (get_global), ~387 (aggregate), ctor/run**: conditional ts feed.
- **ParallelCampaign.{h,cpp}:~240 (ligand loop), ~315 (consensus), LigandResult**: shared ts + register.
- **DatasetRunner.{h,cpp}:~455 (target_servers_), ~5199 (create), ~6913 (register: replace approx + populate best/conc), DockingResult struct + parse sites (~5350,6360), entry creation**.
- **TargetServer.{h,cpp}:~40 (DockingSession add conc), ~55 (register use conc)**.
- **New/minor**: tests/test_grand_integration.cpp (or extend existing); update competition yaml usage later (P3).
- No changes to GA/scoring/statmech core, CMake sources list (P0), Python (P2).

**Next after proposal**: Real impl on feat branch + worktree (per plan), spawn subagents per track, todo_write, verify execution before claims.

This proposal is self-contained for fast handoff post-P0. All facts validated against files in this session.

(End of P1 survey proposal. Only this .md was written.)
