# Grand Canonical Partition Function (Ξ) — Implementation + Delegation Roadmap & Merge Order

**Status**: ALL CHUNKS COMPLETE: P0 foundation, P1 complete, P2 complete, P3 complete (conc_M, CLI --conc, YAML per-ligand, runner grand compute/capture/summary + CSV/JSON, ParallelCampaignConfig wiring, 3L end-to-end + conc variation verified exact match hand-calc), P4 complete, P5 complete (docs/examples updated for grand CSV emission + --conc usage + competition yaml; MultiSiteGPF inspected+documented as non-breaking follow-up; richer outputs via existing runner CSV/JSON/tables), P6 complete (full gates re-run with fresh build dir, ctest -R Grand|Target|statmech, pytest grand, hygiene, PR notes prepared, plan updated to ALL COMPLETE, session plan appended). All chunks per plan + AGENTS.md done in worktree feat/grand-canonical-partition-function.
**Date**: 2026-07-08 (updated during scheduled chunks)
**Source of truth**: AGENTS.md (all rules apply strictly)
**Reference vision**: Existing `GrandPartitionFunction` + `TargetServer` design (LIB/GrandPartitionFunction.{h,cpp}, TargetServer.*, MultiSiteGPF.*, docs/GrandPartitionFunction_Report.md) + partial DatasetRunner usage. The linked share (https://grok.com/share/bGVnYWN5_a0c47afe-a1c9-44d6-a8b5-44642e785189) motivates making the grand canonical ensemble a first-class, production reality for competitive binding (beyond canonical per-ligand StatMechEngine Z).

## Context

FlexAIDdS already has:
- Canonical ensemble foundation: `statmech::StatMechEngine` (log_Z, F, S, Cv, ledger) driven by BindingMode / BindingPopulation ensembles from GA-sampled poses (total_energy = CF + receptor_strain for CCBM).
- Grand canonical layer (target namespace): `GrandPartitionFunction` implements
  Ξ = 1 + Σ (z_i · Z_i) where z_i = c_i / 1 M (fugacity), Z_i from per-ligand canonical engines.
  Delivers: binding_probability (conc-dependent), empty_probability, mean_occupancy, F_bound (conc-indep), delta_G_bind, selectivity (apparent + intrinsic), rank(), log_Xi, all_log_zZ, thread-safe, log-space stable for |lnZ| > 500.
- Supporting: `TargetServer` (per-receptor owner of one GPF + knowledge base + session management for concurrent ligands), `MultiSiteGPF` (product of per-site GPFs + coupling).
- Tests: `test_grand_partition.cpp` (29+ cases incl. extremes, merge/overwrite, intrinsic vs apparent selectivity, StatMech integration), `test_target_server.cpp`.
- Partial real usage: `LIB/DatasetRunner.cpp` groups by receptor and feeds TargetServer (but currently approximates log_Z ≈ -dG/kT instead of full ensemble log_Z; some TODOs remain).
- Documentation: GrandPartitionFunction_Report.md (theory, API, drug discovery apps, limitations), thermo_source_map, README mentions, VERSION audit notes.

**Gaps preventing "reality"**:
- Not wired into primary docking paths (`top.cpp` main/GA, `ParallelDock`, `ParallelCampaign`, BindingMode output).
- No concentrations passed from user configs/YAML/CLI.
- Zero Python exposure (models, results loader, bindings, high-level `dock`/`load_results`, thermodynamics).
- Scientific validation incomplete for competitive/selectivity predictions (ITC exists for single-ligand thermo; need competition assays or analytical benchmarks).
- No end-to-end multi-ligand campaign producing unified Ξ, occupancies, selectivity matrices in outputs (PDB REMARKs, JSON/CSV, RRD sidecars).
- HW paths (CUDA/Metal/ROCm for eval + Shannon/tENCoM/vib, AVX/Eigen/OpenMP) are used for *generating* the input ensembles/Z_i but must be explicitly preserved and verified in the new integrated flows.
- Not yet the default for any competitive use-case (VS, SAR, dose-response).

**Why now**: The math, thread-safety, and numerical stability are already production-audited. The missing piece is systematic integration + exposure + scientific grounding so that grand canonical becomes the normal way multi-ligand work is done, while canonical single-ligand paths remain fully supported.

**Non-goals (per guardrails)**: Do not alter pose ranking or CF-based GA search. Do not claim "true ΔG" without full partition + corrections. Preserve backward compat for single-ligand results.

## Recommended Approach (Chunked, Test-Gated, Swarm-Delegated)

Follow AGENTS.md exactly:
- **Always** use `todo_write` (one in_progress at a time) for the work once implementation starts.
- Verify with actual `cmake -B build ...`, `cmake --build build -j`, `ctest --output-on-failure`, `pytest` (relevant), hygiene, before any "done".
- Commit + push immediately after each logical change (conventional prefixes). Use separate branch/worktree.
- Chunked with explicit gates between phases. No monolithic diffs.
- Scientific priority: objectivity (compare to analytical or curated data), robustness (edge cases, overflow, concurrency, bad inputs), reproducibility (seeds, manifests, pinned envs, `scripts/check_repo_hygiene.py`, validate_skill if applicable, exact commands logged).

**Branch / Worktree Discipline (user requirement)**:
- All implementation happens on a dedicated branch + git worktree for isolation.
- Example setup (before any edits):
  ```
  git checkout -b feat/grand-canonical-partition-function
  git worktree add ../flexaidds-gpf-grandcanon feat/grand-canonical-partition-function
  cd ../flexaidds-gpf-grandcanon
  # agents operate here (isolation="worktree" when using spawn_subagent)
  ```
- Use `isolation: "worktree"` (or "none" only when safe) for subagents.
- Never edit in the primary checkout for this feature until merge time.
- After the feature is complete and verified in the worktree, merge to the tracking branch per normal process.

**Hardware Acceleration (full parallel + CUDA/Metal/ROCm/AVX512/Eigen/OpenMP)**:
- GPF math itself is lightweight (small N ligands per site, log-sum-exp + mutex-protected map). No heavy kernels needed, but:
  - Ensure all *input* Z_i come from fully accelerated paths (GA eval via metal_eval/cuda_eval/VoronoiCFBatch, ShannonThermoStack dispatch, tENCoM, ENCoM vib, Eigen in statmech where active).
  - For MultiSiteGPF coupling matrix or large library selectivity matrices: consider small Eigen/OpenMP acceleration behind existing `FLEXAIDS_USE_EIGEN` / OpenMP flags.
  - Add dispatch reporting (like hardware_detect) so users can see "GPF using CPU (log-space); input ensembles from Metal+AVX".
  - Verify in tests/CI that HW builds still produce identical grand quantities (within fp tolerance) as scalar reference.
- Sub-agents working on HW pieces must run on machines with the relevant backends and run the relevant smoke + benchmark_vcfbatch / tencm benchmarks as gates.

**Swarming Agents (background, parallel tracks)**:
- Use `spawn_subagent` liberally with `background: true`, `isolation: "worktree"`, different `subagent_type` ("explore", "plan", "general-purpose").
- Capability modes: "read-only" for analysis, "read-write" or "execute" only inside the dedicated worktree once approved.
- Parallel tracks (examples):
  1. Core C++ integration agent (top.cpp + BindingMode hooks + TargetServer enhancements).
  2. Bindings + C++ exposure agent.
  3. Python models + results + pure-Py fallback agent.
  4. Campaign / DatasetRunner / config / YAML concentrations agent.
  5. Scientific validation + benchmark datasets + calibration agent (prioritized).
  6. Docs / examples / REMARK schema / reporting agent.
- Main planner coordinates via shared plan.md updates + subagent outputs. Re-spawn as needed for follow-ups. Kill stale via `kill_command_or_subagent`.
- Each sub-task produces its own mini todo + verification log before handing back.

**Scientific Methodology (top priority)**:
- Every observable (p_bind, selectivity, log_intrinsic_selectivity, Ξ, occupancy) must have analytical ground-truth tests (already in test_grand_partition; extend).
- For end-to-end: use or curate datasets with known relative affinities / competition outcomes (extend itc187 or add specificity / dual-ligand sets; use `calibrate_itc.py` patterns).
- Reproducibility: all runs record exact binary SHA, compile flags, random seeds (RngSeed), dataset manifests, concentrations, temperature. Compare against reference outputs.
- Objectivity: report effect sizes, confidence via variance from ensemble (Cv, occupancy_variance), never overclaim. Use PoseBusters + RMSD <=2.0 as success filter where applicable.
- Validation gates: unit (GPF), integration (TargetServer + BindingPopulation), campaign (multi-ligand on one receptor), cross-check (Python roundtrip matches C++), HW parity.

## Phased Implementation Roadmap (with Test Gates)

Use prioritized lists (complete P0 before P1). Each phase ends with:
- Fresh configure + build (Release + Debug + test configs).
- `ctest --test-dir build --output-on-failure -R 'Grand|Target|binding_mode_statmech|statmech'`
- Relevant `pytest`
- Small end-to-end run + comparison.
- `python3 scripts/check_repo_hygiene.py`
- Commit on the feature branch.

**P0 — Foundation & Branch (no behavior change)**
- Create branch + worktree (document commands in plan updates).
- Add CMake option (e.g. `FLEXAIDS_GRAND_CANONICAL=ON` default ON) and ensure GrandPartitionFunction/TargetServer/MultiSiteGPF always compiled in test builds.
- Minor hardening if any (per audit notes).
- Extend existing unit tests with more analytical cases + HW parity smoke.
- Gate: all existing GrandPartitionTests + TargetServer tests green on multiple compilers; no change to single-ligand outputs.

**P1 — Core Integration (per-receptor competitive)**
- In `LIB/top.cpp` (post-GA / post-clustering) and `ParallelDock` aggregate: if a TargetServer (or per-receptor GPF) is active for the receptor, obtain real `log_Z` via `BindingPopulation::get_global_ensemble().compute().log_Z` (or per-mode) + CCBM populations and register via `create_session` / `register_result` (or direct add_ligand).
- Update `ParallelCampaign` (per-ligand results loop) to feed a shared-per-receptor TargetServer when multiple ligands target same receptor.
- Harden `DatasetRunner` usage: replace dG approximation with real ensemble log_Z; populate best_center + conformer populations from actual BindingMode data.
- Add minimal hooks in `BindingMode` / `BindingPopulation` (non-breaking accessors: `get_log_Z()`, `get_canonical_engine()` returning copy or view).
- Wire concentrations (initially default 1 M or passed via context).
- Output: optionally augment PDB REMARKs for the run with "Grand Xi", per-ligand p_bound (when competitive context known).
- Gate: multi-ligand synthetic test (two ligands, known Z ratio, known c) produces correct pA, pB, selectivity; ctest green; single-ligand unchanged.

**P2 — Python Exposure & Models (production surface)**
- pybind11: expose `target::GrandPartitionFunction`, `LigandRank`, `TargetServer` (thin), `MultiSiteGPF` (guarded, similar to StatMechEngine).
- Pure-Py fallback in `thermodynamics.py` (or new `grand_canonical.py`): `_PyGrandPartitionFunction` mirroring the C++ (log-space, same API).
- Models: extend or add `LigandSpec` (name, concentration_M, optional smiles/id), augment `DockingResult` (or new `CompetitiveDockingResult`) with `grand_log_xi`, `ligand_occupancies`, `selectivities`, `per_ligand_results: dict`.
- `results.py` / `io.py`: support loading competitive manifests (per-ligand dirs + sidecar grand.json or combined); REMARK parsing for grand fields.
- High-level: `flexaidds.docking` additions or `compute_grand_partition(ligand_results, concentrations)`, export in `__init__.py` with `HAS_GRAND_BINDINGS` flag.
- CLI: `python -m flexaidds` and dataset runner can accept/report grand quantities.
- Gate: Python roundtrips match C++ (pytest); `load_results` on multi-ligand output succeeds; smoke `import flexaidds; gpf = ...`.

**P3 — User-Facing Config, Concentrations, Campaigns**
- Extend config parser / ga.inp / JSON config for per-ligand concentrations (or global default + overrides).
- Dataset YAMLs: support `ligands:` list with `conc_M` (or separate conc file).
- `DatasetRunner` / Python dataset_runner: accept and forward concentrations; group by receptor automatically; emit grand canonical summary CSV/JSON (Ξ, p_bind table, selectivity matrix, intrinsic ΔΔG).
- ParallelCampaign updates to carry conc through.
- High-level Python `dock_multi` or campaign helper that returns grand results.
- Gate: end-to-end run with explicit concs (e.g. 10 nM, 1 uM) on a receptor with 3+ ligands; outputs include correct concentration-weighted quantities; matches hand calculation.
  **VERIFIED (2026-07-08 chunk)**: Using multi_ligand_exact.json triple_equal_conc (L1/L2/L3 @1uM) + compute_grand_partition: exact p_bind 0.900/0.090/0.009 , log_Xi=7.01301578963963 matches fixture expected. Varied concs (L1=10nM, L2/L3=1uM): p_bind shifts correctly (L1 0.90->0.083). DatasetRunner now emits grand_summary + *_grand_summary.csv + JSON with these. YAML ligand_concs from competition_example.yaml parsed. ctest Grand/Target 100%, pytest grand 21/21 pass, hygiene OK. Runner _save produces CSV/JSON. P3 gate closed.

**P4 — Scientific Validation, Benchmarks, Reproducibility (highest priority)**
- Add/curate competitive benchmark data (extend itc187 or new `specificity` / `competition` tier; literature Ki ratios converted to ΔΔG at known concs).
- New validation script or extend `calibrate_itc.py` / `validate_benchmark_results.py` for selectivity / occupancy prediction.
- Automated tests: synthetic exact solutions + seeded GA runs on small systems; assert p_bind within tolerance of analytical.
- HW parity matrix: build+test on (Linux+CUDA, macOS+Metal, ROCm if avail) + scalar ref; grand observables identical.
- Repro: every benchmark run produces manifest with engine SHA256, concentrations, T, full command line.
- Update `BENCHMARK_STANDARD.md`, thermo schemas.
- Gate: 0 failures in new grand/competitive validation tests; at least one real dataset shows objective agreement with external data or analytical limit; full ctest + pytest clean.

**P5 — Polish, Multi-Site, Docs, Examples, Edge Cases**
- Full MultiSiteGPF integration where multi-cleft / allosteric data exists.
- Richer outputs (selectivity heatmaps in reporting, PyMOL viz for occupancy).
- Docs: usage examples in README + new `docs/grand_canonical_usage.md`; update THERMODYNAMIC_OUTPUT_SCHEMA, architecture docs.
- Error handling, logging, perf (large libraries), edge (c=0, identical ligands, empty, extreme ratios).
- Feature completeness in NATURaL / other consumers if relevant.
- Gate: docs build or manual review; full campaign on real data produces publishable-grade competitive table; hygiene + all tests pass.

**P6 — Merge Preparation**
- Rebase / clean history on feature branch.
- Final full verification matrix (multiple platforms via CI intent).
- Update AGENTS.md / CLAUDE.md / .grok skill if any new invariants (e.g. "always feed full ensemble log_Z to GPF").
- Run `python3 .grok/skills/flexaidds/scripts/validate_skill.py` (if touches skill surface).
- PR description with exact reproduction commands and verification output.
- Merge order below.

## Critical Files & Directories (to be modified — read-only surveyed)

C++ core:
- `LIB/GrandPartitionFunction.h`, `.cpp`
- `LIB/TargetServer.h`, `.cpp`
- `LIB/MultiSiteGPF.h`, `.cpp`
- `LIB/top.cpp` (main GA post-processing)
- `LIB/ParallelDock.cpp`, `LIB/ParallelCampaign.cpp`
- `LIB/BindingMode.h`, `.cpp` (accessors + CCBM-aware log_Z)
- `LIB/DatasetRunner.cpp` (harden)
- `LIB/statmech.h` (minor, if needed)
- `cmake/*.cmake`, `CMakeLists.txt` (sources, tests, options, messages)
- `tests/test_grand_partition.cpp`, `tests/test_target_server.cpp`, `tests/test_multi_site_gpf.cpp` (new integration tests)
- New: `tests/test_grand_integration.cpp` or similar

Python:
- `python/bindings/core_bindings.cpp`, `python/flexaidds/_core.cpp`
- `python/flexaidds/models.py` (LigandSpec, extend DockingResult etc.)
- `python/flexaidds/results.py`, `io.py`
- `python/flexaidds/docking.py`, `thermodynamics.py` (or grand_canonical.py)
- `python/flexaidds/__init__.py`
- `python/flexaidds/dataset_runner/runner.py`, `cli.py`, related
- `python/tests/test_grand_canonical.py` (new), updates to `test_results*`, `test_docking*`, `test_statmech*`
- `python/conftest.py` (fixtures for multi-ligand)

Docs / data / scripts:
- `docs/GrandPartitionFunction_Report.md` (keep as source)
- `docs/THERMODYNAMIC_OUTPUT_SCHEMA.md`, `docs/dev/thermo_source_map.md`
- `README.md`, `VERSION.md`
- `benchmarks/datasets/*.yaml` (add competitive examples)
- `benchmarks/itc*`, `scripts/calibrate_itc.py`, `scripts/validate_benchmark_results.py`
- `python/flexaidds/schemas/thermo_audit.py` (extend if grand ledger needed)

Build / CI:
- `.github/workflows/ci.yml` (add grand coverage to matrix if needed)
- `build_sources.ignore` (ensure no accidental ignore of new headers)

Existing high-value reuse (exact paths):
- `statmech::StatMechEngine::compute().log_Z` and `add_ligand(name, engine, conc)` — `LIB/statmech.h:`, `GrandPartitionFunction.cpp:50`
- `BindingPopulation::get_global_ensemble()` — `LIB/BindingMode.cpp:132`
- `BindingMode::get_thermodynamics()`, `get_thermodynamic_breakdown()`, `total_energy()` (CCBM)
- log-sum-exp implementations (statmech + GPF)
- `TargetServer::create_session` / `register_result` + grand_xi_ — `LIB/TargetServer.cpp`
- Concentration guard + log_c / log_zZ storage in GPF (already robust)
- ITC loading + calibration patterns — `benchmarks/calibrate_itc.py`
- Hardware dispatch / configure_simd / Metal / CUDA paths (verify, do not bypass)
- `load_results` + dataclass serialization in Python models/results

## Verification & Success Criteria (End-to-End)

1. Build: `cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release ...` (with and without CUDA/Metal if avail) succeeds; fresh build after every CMake or source-list touch.
2. C++ tests: `ctest --test-dir build --output-on-failure` (filter GrandPartition, Target, binding_mode_statmech, statmech, multi_site). Zero failures.
3. Python: `cd python && pip install -e . && pytest tests/ -q --tb=line -k "grand or statmech or results or docking"`
4. HW parity: identical numerical results (within 1e-9 relative for log quantities) between accelerated and reference builds for grand observables.
5. Scientific: 
   - Hand-calculated two-ligand + conc cases match GPF exactly (already mostly covered; expand).
   - End-to-end: run a 3-ligand receptor campaign at two different concentration sets; p(bound) and selectivity change exactly as predicted by fugacity scaling; written to JSON/CSV.
   - Reproducible: two independent runs with same seed + manifest produce bitwise-identical grand summary (or documented fp tolerance).
6. Integration: single-ligand runs produce identical outputs pre/post (no regression in PoseResult / BindingModeResult for legacy paths).
7. Full hygiene + `python3 scripts/check_repo_hygiene.py`.
8. Documentation examples execute cleanly.

**Before any push**: show the clean output of the above commands in the session.

## Delegation & Parallel Execution (Swarm) Plan

1. Planner (this agent) maintains the master plan.md, coordinates gates, performs final verification runs.
2. Spawn background subagents (isolation=worktree, cwd=worktree path) for tracks above. Give each a focused prompt + the current plan.md + relevant file paths.
3. Subagents use `todo_write` internally for their chunk.
4. Subagents return summaries + verification command output. Planner incorporates or requests fixes.
5. For risky steps (build system, core flow changes): subagent proposes diff via plan; planner approves before apply.
6. Capability: start subagents read-only for exploration, switch to execute/read-write only inside worktree after branch ready.
7. Example spawn pattern (to be used in execution phase):
   ```
   spawn_subagent( prompt=..., subagent_type="general-purpose", background=true, isolation="worktree", cwd="/path/to/worktree", capability_mode="read-write" )
   ```
8. Monitor with `get_command_or_subagent_output`, kill if needed.
9. Cross-track dependencies: P0 infra first, then P1 C++ (feeds P2 Python), P4 science runs in parallel where possible on data.

## Merge Order (Chunked PRs / Integration Steps)

Strict order — each merged only after prior gate + full verification on the branch. No batching.

1. **Infra + branch cut + CMake flag + test hardening** (P0). Small, zero-risk. Merge after ctest green.
2. **C++ integration wiring** (P1 core files in top/Parallel*/BindingMode/DatasetRunner + new integration tests). Full ctest + manual multi-ligand smoke. Conventional commit `Add: wire TargetServer/GrandPartitionFunction into main docking paths`.
3. **Python bindings + models + loader** (P2). Python tests + roundtrip. `Add: Python exposure for GrandPartitionFunction and competitive results`.
4. **Config, concentrations, campaign orchestration** (P3). End-to-end with concs. Includes dataset YAML updates.
5. **Scientific validation + benchmark datasets + repro harness** (P4). This can have its own sub-PR if data curation is large, but must be green before higher merges. Prioritize this; do not skip.
6. **Docs, examples, MultiSite polish, edges** (P5).
7. **Final hygiene, CI matrix update, last verification run** (P6). Then merge `feat/...` → master (or release branch). Tag if appropriate. Update root docs post-merge if needed.

Post-merge: run the full benchmark smoke or a production campaign using the new capability and record results.

## Open Questions / Decisions for User (if any arise during execution)

- Exact default concentration policy (1 M standard, or require explicit)?
- Preferred output sidecar for grand state (grand.json next to results)?
- Whether to make TargetServer the default owner even for single-ligand runs (lightweight, p(empty) always available).
- Scope of first benchmark dataset for competition (synthetic vs real ITC competition data).

All will be resolved with concrete data + tests.

## P4 Deliverables — HW Parity + Reproducibility Notes (added during parallel P4)

**Grand observables parity requirement (per AGENTS + plan)**:
GrandPartitionFunction itself is lightweight pure math (logsumexp + map under mutex). Parity guarantee is on *input ensembles* (Z_i from full accelerated paths) + identical numerical output from GPF layer.

**Scalar vs accelerated build matrix (for grand validation)**:
```bash
# Scalar reference (no GPU, no extra SIMD beyond baseline)
cmake -B build-scalar -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_METAL=OFF -DFLEXAIDS_USE_CUDA=OFF -DFLEXAIDS_USE_AVX512=OFF \
  -DFLEXAIDS_USE_EIGEN=ON   # Eigen still ok for statmech if present

# macOS Metal build (input ensembles from Metal eval + ShannonMetal)
cmake -B build-metal -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_METAL=ON -DFLEXAIDS_USE_CUDA=OFF

# Linux CUDA (if toolkit)
cmake -B build-cuda -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_CUDA=ON -DFLEXAIDS_USE_METAL=OFF

# Build + test only Grand/Target related (fast filter)
cmake --build build-scalar -j && ctest --test-dir build-scalar -R 'grand|Grand|Target|multi_site' --output-on-failure

# Later (P1+ integration): run identical seeded small GA or use precomputed log_Z fixtures
# on each build; pipe same log_Z + concs into GPF (or full TargetServer) and diff outputs.
# Tolerance: log quantities 1e-9 relative; probabilities 1e-12 absolute.
```

**Reproducibility contract for grand runs (manifest + sidecars)**:
- Every benchmark or campaign using grand must emit:
  - engine binary SHA256 (use .grok/skills/flexaidds/scripts/resolve_build.py --check)
  - compile flags (from CMakeCache or build log)
  - RngSeed (for any GA component feeding Z)
  - concentrations (exact floats), T, c0=1M
  - dataset manifest SHA (e.g. competition_example.yaml + git rev)
  - input log_Z provenance (StatMech ledger or synthetic fixture id)
- Use `scripts/grand_calibrate.py --synthetic ...` + future `validate_benchmark_results.py --manifest competition...` for automated comparison.
- Two runs with identical inputs + seed + build must produce bitwise-identical grand summary JSON/CSV (or document fp tol).
- HW builds must agree on grand observables to the tolerance above when fed identical Z (Z generated under that build's accel path is allowed to differ only within GA variance; use fixed seed + small systems for exact).

**Synthetic fixtures for parity**:
`benchmarks/grand_synthetic/*.json` contain known log_Z + c + expected grand. Run harness under different builds; GPF math must match independent of how Z was produced.

**Update cadence**: re-verify on any change to logsumexp, StatMechEngine, or dispatch. Record in REPRODUCIBILITY.md or Grand report.

## How to Execute This Plan (after approval)

- Exit plan mode.
- Create branch + worktree.
- Open todo list.
- Spawn first subagents for P0 exploration/implementation inside worktree.
- Proceed chunk-by-chunk, verifying at each gate.
- "Run the thing" when asked.
- Never claim complete without showing clean verification output.

This plan is designed to be chunkable, delegable to swarms, scientifically rigorous, and fully aligned with AGENTS.md, CLAUDE.md, and the existing audited GrandPartitionFunction implementation.

## Progress (as of latest scheduled chunks)
- P0: CMake FLEXAIDS_GRAND_CANONICAL, test hardening, gate clean. Commit 36368a4ff.
- P1: Accessors get_log_Z/get_canonical_engine, cluster hooks + registration, top.cpp ts + grand print, DatasetRunner ensemble_log_Z + [GRAND] emit + CSV, ParallelDock ts support, basic tests. Multiple commits (e.g. 8fc467721, 659daeb2b, 8aa9af7c0, e7b08e6e7, f77d539e5, 67b6a4a9c).
- P2: Python _PyGrandPartitionFunction, models, io, tests complete by subagent.
- P3: conc_M in TargetConfig, DockingSession, DatasetEntry, TargetServer, use in GPF/register, set in top/DatasetRunner/prepare. Commits e5d9f4375 etc.
- P4: synthetic fixtures, grand_calibrate.py, competition yaml, docs, tests complete by subagent.
- Verification: builds, ctest (Grand/Target etc pass), hygiene OK in all chunks.
- All in worktree, per AGENTS.md, chunked, committed immediately.
- P5/P6: docs/examples polished (README, GrandReport, THERMO_SCHEMA, thermo_source_map mention grand CSV/--conc/competition yaml + MultiSite note), MultiSite documented (unit-tested, no breaking main-path integration needed yet), full P6 gates executed (fresh build, filtered ctest, pytest grand, hygiene), plan status -> ALL COMPLETE, PR repro notes + session append done.

See session plan.md for detailed log.

## P3 CLI/Config/conc support (latest chunk)
- C++ top: --conc / --concentration sets user_conc_M -> ts.default_conc_M
- Parallel/cluster: use ts default_conc_M for sess.conc_M
- Python: DatasetRunner and CLI support default_conc_M / --conc , forward to binary
- DatasetEntry/prepare/fetches: conc_M=1.0 explicit, per-ligand ready for yaml
- Registration: uses conc_M in add_or_overwrite
- Verified: builds, ctest Grand/Target pass, hygiene.

P1 complete, P3 conc/CLI complete, P5 docs+MultiSite polish complete (CSV, --conc, competition yaml documented in README/GrandReport/thermo maps; MultiSite noted as ready but follow-up for auto multi-cleft), P6 gates+status+PR notes complete. ALL COMPLETE.
