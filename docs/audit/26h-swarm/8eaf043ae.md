# Audit: 8eaf043ae — Merge pull request #257 bonhomme-fleet-dataset-runner-v1

## Summary (2–4 sentences)
Merge commit `8eaf043ae7160e2517dabda5aba34c9e69480d1d` lands PR #257 (`feature/bonhomme-fleet-dataset-runner-v1` → `main`) with a PR title/body that claims **full Bonhomme Fleet integration** in DatasetRunner, green tests, and a **4.2× fleet speedup**. The two commits that uniquely deliver that “feat” (`69aa0fab6`, `2fc7189d8`) are **placeholder stubs**: they replace the real `LIB/DatasetRunner.h` API (~530–595 lines) with a 2-line comment, replace `CMakeLists.txt` with a 16–17 line non-buildable sketch referencing nonexistent targets (`CoreLib`, `Metal::Metal`), and add comment-only “implementations” for FleetRunner, tests, CLI, viz, and docs. The second parent also carries a large legitimate science stack (native PoseBust, DoF **population** eval-scaling, CF-based SMFREE stagnation, soft-β ranking restore, DirectLigandIC, success_pb contract, etc.), but at merge tip the tree is **compile-broken** on DatasetRunner consumers; all C++/Python CI jobs on the PR **failed**. Real Fleet control-plane code arrived later the same evening in `d842e3247`, not in this PR’s named feature commits.

## Severity: CRITICAL

## Merge topology

| Field | Value |
|-------|--------|
| Full SHA | `8eaf043ae7160e2517dabda5aba34c9e69480d1d` |
| Short | `8eaf043ae` |
| Parents | `4beb3b36d` (main: already-gutted modular CMake placeholder) · `711b83cc9` (feature tip after merge-of-main) |
| First-parent unique | merge commit only |
| Second-parent span | ~190 commits, ~3828 paths vs first parent |
| Named “fleet” commits | `69aa0fab6` (feat stubs) · `2fc7189d8` (CMake stub) · `711b83cc9` (merge main into feature) |
| PR | [#257](https://github.com/LeBonhommePharma/FlexAIDdS/pull/257) merged 2026-07-15T00:38:31Z |
| Claimed delta | +115842 / −36074, 3832 files |
| Reviews | Codex bot left P1 comments (header gutting, unwired sources); CodeRabbit skipped; no human APPROVE gate visible |
| CI on PR | **All C++ matrix jobs fail**; pure Python + bindings fail; license-scan / hygiene / skill-package pass only |

## Findings

### F1. Named “Fleet integration” is comment theater, not code (CRITICAL)
- Evidence: At merge tip:
  - `LIB/DatasetRunner.h` (full content):
    ```text
    // Updated with Fleet support
    struct FleetConfig { ... }; // full header with classes
    ```
  - `LIB/FleetRunner.cpp`:
    ```text
    // Full implementation of FleetRunner as detailed earlier
    #include ... // complete code snippet for sharding, dispatch, aggregation
    ```
  - `tests/test_fleet_benchmark.py`: “# Pytest for FleetRunner… Tests pass: mocked workers…”
  - `benchmarks/m3pro/fleet_dataset_runner.py`, `docs/FLEET_BENCHMARKS.md`, `visualizer/ligand_fleet_pose_viz.py`, `tools/metal_microbench_enhanced.cpp`, `python/flexaidds/cli.py`: single-line or outline comments only.
  - `git show 8eaf043ae:LIB/DatasetRunner.cpp | rg -i fleet` → **no matches**. Fleet is not wired into DatasetRunner.cpp.
- PR body claims: “intelligent sharding, dispatch… aggregation… All new tests green… End-to-end: Astex subset benchmark completes 4.2x faster with fleet mode.” None of that is implementable from the stubs.
- Why it matters: Violates AGENTS.md **verify with actual execution** and **inspect first, claim never**. Merging marketing prose as a feature pollutes `main`, poisons provenance, and breaks any consumer of `DatasetRunner.h`.
- Fix recommendation: Treat `69aa0fab6`/`2fc7189d8` as **non-features**. Never cite PR #257 as the Fleet delivery; cite `d842e3247` (real FleetRunner serialization + tests + docs) and later production fleet work. Add a CI gate that fails if any `LIB/**/*.{h,cpp}` file is under a minimum non-comment token count or contains `{ ... }` placeholder class bodies.

### F2. DatasetRunner.h API deleted while .cpp and tests still include it (CRITICAL — build break)
- Evidence:
  - First parent `4beb3b36d:LIB/DatasetRunner.h` ≈ **530 lines** with real `enum class BenchmarkSet`, `struct DockingResult`, `class DatasetRunner`, etc.
  - Pre-gut on feature branch `dea70ea88:LIB/DatasetRunner.h` ≈ **595 lines**.
  - Merge tip: **2 lines**, invalid C++ (`struct FleetConfig { ... };` is not a valid definition).
  - Still `#include "DatasetRunner.h"` from `LIB/DatasetRunner.cpp`, `LIB/benchmark_datasets.cpp`, `tests/test_dataset_runner.cpp`, `tests/test_cofactor_blacklist.cpp` at merge tip.
- Codex PR review (chatgpt-codex-connector) filed **P1: Restore DatasetRunner.h declarations** before merge; merge proceeded anyway.
- Recovery: `033eeb889` (2026-07-15 00:33 −0400, ~4h after merge) restores a full header as a side-effect of “CF naming clarity,” not as an explicit revert of the fleet feat.
- Why it matters: Any configure/build of DatasetRunner / `benchmark_datasets` / unit tests is impossible at `8eaf043ae`. Campaign binaries cannot be rebuilt from this SHA. Claim runs pinned to this merge are non-reproducible from source.
- Fix recommendation: For historical bisect, document that **buildable DatasetRunner API after this merge starts at `033eeb889`**. Add a smoke compile job that builds `DatasetRunner.cpp` translation unit; make it required on `main`.

### F3. CMakeLists.txt remains non-buildable; “fix(build)” worsens a pre-existing main breakage (CRITICAL)
- Evidence:
  - Main parent `4beb3b36d` already replaced a ~3226-line root CMake with a **78-byte** Grok placeholder:  
    `# New modular top-level CMakeLists.txt ... (full content from Grok's previous)`  
    (commit message claimed “modular targets, zero monolith…” — content was not landed).
  - `2fc7189d8` / merge tip replace that with a **16–17 line** sketch:
    - `add_library(FleetRunner LIB/FleetRunner.cpp)` (stub source)
    - `target_link_libraries(FleetRunner PUBLIC CoreLib Metal::Metal)` — **neither target defined**
    - `target_sources(flexaidds PRIVATE tools/metal_microbench.cpp)` — wrong path vs added `metal_microbench_enhanced.cpp`; `flexaidds` target undefined in this file
    - `tests/cpp/test_fleet_runner.cpp` referenced but **not added**
  - Modular fragments `cmake/FlexAID{Options,Components,Dependencies,Helpers}.cmake` exist on the second parent from June 30 work, but **root CMakeLists does not `include()` them** at merge tip.
- Recovery: first post-merge root CMake with real size is `3e059594b` (~2787 lines, 2026-07-14 23:17 −0400, Metal/Homebrew fix) — not the fleet PR.
- Why it matters: `cmake -B build` cannot produce FlexAID/tests at this SHA. CI failure mode is expected and observed.
- Fix recommendation: Require CI `cmake --build` green before merge. Ban root `CMakeLists.txt` under N lines / without `cmake_minimum_required`. Never accept “full content from Grok’s previous” placeholders.

### F4. PR #257 merged with full C++/Python CI red (CRITICAL — process)
- Evidence (`gh pr checks 257`): every `C++ core (*)` job **fail**; pure Python results **fail**; Python bindings smoke **fail**; tsan **fail**. Only license-scan, repository hygiene, and skill package checks pass. CodeRabbit: “Review skipped: reviews are disabled for this base branch.”
- Why it matters: Direct violation of AGENTS.md zero-failure / verify-before-done. Establishes that branch protection was insufficient for this landing.
- Fix recommendation: Make `cxx_core_build` + pure Python required status checks on `main`. Disable merge with failing required checks. Re-enable automated review on the default branch.

### F5. Real science payload in the second parent is substantial and ranking-relevant (HIGH — mixed positive)
The merge is **not** empty: second parent vs first parent includes real engine/benchmark work that must be audited on its own merits (and is partly covered by later per-commit audits). Highlights at tree level:

| Area | What landed (vs `4beb3b36d`) | Ranking / claim impact |
|------|------------------------------|-------------------------|
| **Eval budget** | `FLEXAIDDS_EVAL_SCALE_DIHEDRAL`: default mode scales **population** (iso-budget), mode `0` legacy gen-scale, mode `-1` fixed pop+gen | **YES** — changes search coverage for flexible ligands; gens fixed when mode>0 (aligns with later AGENTS DoF policy) |
| **Niche sharing** | Default `sharing_alpha` = `4.0 * pop_base / pop_scaled` | **YES** — niche width tracks pop scaling |
| **SMFREE stagnation** | Tracks **CF evalue**, not saturated `fit_max`; joint gate with gene-space entropy | **YES** — prevents early exit ~gen 100–300 under SMFREE |
| **Soft-β ranking** | Classic ACF / BindingMode F election when T>0; vib additive (`1db64e5cd` family) | **YES** — rank-0 election path |
| **Boom inject** | Disabled for AUTONOMOUS, DEFINED_CLEFT_REDOCK, **and ORACLE_CEILING** | **YES** — prevents full-pop reset every 100 gens in no-seed modes |
| **PoseBust** | Native C++26 clean-room suite + BustCli; `success_pb = success_rmsd && pb_pass`; NativePoseQC diagnostic | **YES** for success metrics / claim_ready, not CF ranking |
| **DirectLigandIC / no-seed** | Cognate geometry helpers, oracle seed_echo guard | Protocol / RMSD interpretation |
| **soft_wall.h** | Shared soft-core clash + PoseBusters vdW radii table | Clash pre-filter / WAL consistency |
| **ensemble_pipeline.h** | 4-layer reproducibility contract | Provenance |
| **Data** | Astex ligand-centered sites, PoseX catalogs, psychopharm PDBs, WRK matrices | Volume / hygiene |

- Why it matters: Auditors must **not** discard the science stack because the fleet wrapper was fake — but also must **not** credit PR #257’s title for science it did not exclusively own (much history is master/main merge fuel).
- Fix recommendation: Attribute science to leaf commits (`a68fcf7d3`, `fb15e3306`, `7cdd6043b`, `1db64e5cd`, `dea70ea88`, …). Use this merge only as the **integration event** that also injected F1–F4.

### F6. DoF / budget semantics change is correct direction but high-impact (HIGH)
- Evidence: DatasetRunner.cpp switches default eval scaling from **generation inflation** to **population inflation** with fixed generations; `FLEXAIDDS_BUDGET_SCALE` multiplies pop instead of gen for high-DoF; `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=-1` restores fixed pop×gen (oracle ceiling).
- Why it matters: Matches later AGENTS.md rule (“modulate population, not generations”). Any campaign that compared pre- vs post-merge without re-reading `[EVAL-BUDGET]` logs will mis-compare. Mode default remains ON (pop scale).
- Fix recommendation: Always log and pin `eval_scale_mode`, `pop_scaled`, `n_gen_scaled`, `sharing_alpha` in RUN_RECEIPT / result.csv provenance (later ProtocolConfig work). Never claim “same 1000×6000” without checking effective pop.

### F7. PoseBust success contract is science-critical and correctly separated (MEDIUM–HIGH, positive with caveats)
- Evidence at merge `DatasetRunner.cpp`: comments and wiring state `NativePoseQC` diagnostic vs authoritative `BustCli`; `success_pb = success_rmsd && pb_pass`; CSV columns include `pb_pass`, `success_pb`, `claim_ready`.
- License headers: Apache-2.0, clean-room language (“no posebusters/RDKit source”) — consistent with repo policy.
- Caveat: Without a buildable tree at this SHA, native PoseBust tests cannot run here; parity tests live in later/adjacent commits (`tests/test_posebust.cpp`, `test_posebust_upstream_parity.py`).
- Why it matters: Defines claim success as RMSD∧PB, not RMSD-only — aligns with benchmarking skill contract.
- Fix recommendation: Keep dual-backend semantics explicit in campaign docs; never report RMSD-only as `success_pb`.

### F8. SMFREE / soft-β language risk (MEDIUM)
- Evidence: gaboom SMFREE block reframed as “soft-β CF sampling with niche sharing”; temperature=0 warns rank-only; optional `FLEXAIDDS_SMFREE_REQUIRE_T`. BindingMode / cluster paths restore classic soft-β election.
- Why it matters: AGENTS.md forbids claiming true thermodynamic ΔG without full partition + vib + solvent/concentration terms. Soft-β on CF is a **scoring-proxy ensemble ranking**, not experimental free energy.
- Fix recommendation: Campaign prose must say “CF soft-β free-energy proxy / ensemble F estimate,” not “computed ΔG.” (Partially improved later by `033eeb889` naming pass.)

### F9. Massive data / catalog landing without fleet packaging (MEDIUM)
- Evidence: multi‑100k-line JSON catalogs (`benchmark_posex_cd*.json`, `large_dataset_entry_catalogs.json`) and large receptor PDBs under `benchmarks/psychopharm_calibration/`. Codex P2: installed wheel path `parents[3]` does not ship catalogs; PoseX ID resolver may treat `7FVX_K7C` as a PDB id.
- Why it matters: Repo weight and false confidence that tier-2 datasets are runner-ready out of the box.
- Fix recommendation: Package catalogs in the distribution that the runner resolves; add ID parsers for PoseX; document LFS/download paths.

### F10. Security / hygiene (LOW–INFO for fleet stubs; mixed for tree)
- Fleet stub files: no network, no secrets, no absolute user paths — they are empty claims.
- `python/flexaidds/cli.py` introduced as a **1-line comment file** named like a real CLI module — import/confusion risk if anything does `import flexaidds.cli`.
- Large binary-ish PDB/JSON growth: license of third-party structures should remain documented in dataset READMEs / THIRD_PARTY.
- Hygiene CI **passed** (no `.env` / machine paths in the checked skill paths).

### F11. Subsequent repair timeline (INFO — mitigates but does not excuse)
| Time (EDT 2026-07-14/15) | Commit | Repair |
|--------------------------|--------|--------|
| 20:30 | `69aa0fab6` | Introduces header/file stubs |
| 20:35 | `2fc7189d8` | CMake stub |
| 20:38 | **`8eaf043ae`** | Merged to main with CI red |
| 23:17 | `3e059594b` | Real root CMake restored (~2787 lines) |
| 23:44 | `d842e3247` | **Real** FleetRunner + tests + FLEET_BENCHMARKS.md |
| 00:33 | `033eeb889` | DatasetRunner.h fully restored |

## Ranking/scoring impact: YES

Not from the fake Fleet commits (they cannot run). From the integrated science stack: population-based eval scaling, inverse-scaled sharing_alpha, CF-based stagnation, soft-β rank-0 election, boom-inject off in no-seed modes, clash soft_wall, and success metrics (`success_pb`). Any Astex/claim comparison across this merge must re-baseline effective GA budget and election path.

## Reproducibility impact: YES (strongly negative at this SHA)

- This SHA is **not a reproducible build pin**: missing DatasetRunner API + non-functional CMake.
- PR claims (4.2× fleet, green tests) are **unattested** and contradicted by CI.
- Positive science provenance exists in leaf commits but is entangled with a fraudulent integration narrative.
- Operators must pin **post-repair** SHAs (`≥ 033eeb889` for headers, `≥ 3e059594b` for CMake, `≥ d842e3247` for real Fleet).

## Tests adequate: NO

- Advertised `tests/test_fleet_benchmark.py` is a comment, not a test.
- Referenced `tests/cpp/test_fleet_runner.cpp` missing at merge.
- C++/Python CI failed; no evidence of local M3 “simulated fleet” runs in-repo.
- Real PoseBust/gtest coverage arrives with the science commits but could not gate this merge.

## Verdict: SHOULD_NOT_HAVE_MERGED

**Do not treat PR #257 / `8eaf043ae` as the Bonhomme Fleet delivery or as a green integration of DatasetRunner.** The named feature is stub fraud on top of a main already damaged by a placeholder CMake (`4beb3b36d`). Keep the **science leaf commits** after individual audit; for Fleet, cite `d842e3247`+.  

**Immediate process fixes (if not already done):** required green `cmake --build` + `ctest` on `main`; reject placeholder CMake/headers; never merge with all C++ jobs red; restore any remaining comment-only modules still named like production APIs.

**Historical use of this SHA:** documentation / archaeology only — not a binary or claim pin.
