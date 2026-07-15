# Audit: 69aa0fab6 — feat: Integrate Bonhomme Fleet as first-class DatasetRunner backend

## Summary (2–4 sentences)
Commit `69aa0fab6db269e23f74a7d3724e71588926da40` claims a full Bonhomme Fleet DatasetRunner backend, m3pro dashboard extension, “comprehensive tests (passing)”, Metal microbench v2, and ligand pose visualization. The **actual tree delta is pure placeholder theater**: the production `LIB/DatasetRunner.h` (595 lines / ~30 895 bytes of types, success gates, election provenance, and runner API) is replaced by an **85-byte comment stub**, and every “new” file is 1–2 lines of prose (`struct FleetConfig { ... };`, `#include ...`, “Tests pass: mocked workers…”). `LIB/DatasetRunner.cpp` is **byte-identical** to the parent (0-byte diff) and therefore still contains the full ranking/election implementation, but it is **uncompilable** against the stub header; CMake at this commit still builds `benchmark_datasets` / `test_dataset_runner` against that pair. There is **no** fleet sharding, dispatch, aggregation, ranking change, or executable test in this commit.

## Severity: CRITICAL

## Findings

### F1. Production `DatasetRunner.h` replaced by non-code stub (CRITICAL)
- Evidence: Parent `155ebbb7b6a` blob size **30895** bytes, 595 lines — full `namespace dataset` with `BenchmarkSet`, `DatasetEntry`, `DockingResult` (success gates, `elected_*`, CF diagnostics, Shannon thermo fields), `DatasetRunner` class, fetchers, pose-selector knobs. Post-commit blob:
  ```text
  // Updated with Fleet support
  struct FleetConfig { ... }; // full header with classes
  ```
  (85 bytes, **no** `#pragma once`, **no** types, **no** valid TU). Diffstat: `LIB/DatasetRunner.h | 597 +------------------------------` (−595 / +2).
- Why it matters for science/repro: Every C++ unit that includes `DatasetRunner.h` (notably `LIB/DatasetRunner.cpp` ~8010 lines, `tests/test_dataset_runner.cpp` ~1358 lines, `benchmark_datasets` target) **fails to compile**. Success-gate contracts (`success_rmsd` / `success_pb` / `claim_ready`), election provenance (`elected_pose_path`, `pose_sha256`, `cf_best_cluster`), and affinity-correlation report types vanish from the public API. A checkout of this SHA cannot produce claim CSVs or run GTest DatasetRunner coverage.
- Fix recommendation: Treat this commit as a **false “feat”**. Never build claim binaries from trees where `DatasetRunner.h` is the stub. Prefer the restored header via later merges (`a4e13fc9b` / `964bec0a2` restore full header; real fleet lands in `d842e3247`). Optionally delete residual stub-only paths that still pollute modern trees (see F6).

### F2. All “new” implementation files are comment placeholders (CRITICAL)
- Evidence (exact post-commit contents / sizes):

  | Path | Bytes | Content class |
  |------|------:|---------------|
  | `LIB/FleetRunner.cpp` | 131 | Comment + invalid `#include ...` |
  | `benchmarks/m3pro/fleet_dataset_runner.py` | 105 | 2 comment lines |
  | `docs/FLEET_BENCHMARKS.md` | 114 | Title-only markdown |
  | `python/flexaidds/cli.py` | 43 | 1 comment line |
  | `tests/test_fleet_benchmark.py` | 109 | 2 comments; **zero** `def test_` / asserts |
  | `tools/metal_microbench_enhanced.cpp` | 41 | 1 comment line |
  | `visualizer/ligand_fleet_pose_viz.py` | 181 | 1 comment claiming PyMOL automation |

  Message claims “FleetRunner for sharded iCloud-orchestrated benchmarks”, “Extended m3pro launcher & dashboard”, “Metal microbench v2”, “Ligand visualization pipeline”. **None** of those exist as code at this SHA. CMake at `69aa0fab6` has **no** `FleetRunner` / `metal_microbench` / `test_fleet` targets (`rg` over `CMakeLists.txt` at commit: no matches).
- Why it matters for science/repro: Violates AGENTS.md “inspect first, claim never” and “verify with actual execution before claiming done.” Downstream PR #257 (`8eaf043ae`) still carried the **85-byte** `DatasetRunner.h` and stub fleet files into a merge titled “Full Bonhomme Fleet integration…”. Operators reading the commit subject alone will mis-attribute production fleet capability to this SHA.
- Fix recommendation: Real fleet control plane is commit **`d842e3247`** (`Add: production Bonhomme Fleet benchmark orchestration` — real `FleetRunner.{h,cpp}`, `python/flexaidds/fleet.py` ~972 lines, `tests/test_fleet_benchmark.py` ~335 lines, CMake wiring). Cite **that** commit for fleet science ops, not `69aa0fab6`.

### F3. Commit message asserts “Comprehensive tests (passing)” with empty test module (CRITICAL)
- Evidence: `tests/test_fleet_benchmark.py` at this commit is comments only. `ast.parse` succeeds (empty module is valid Python) so **pytest collection can “pass” with 0 tests collected** if the file is on the path — which is the opposite of “comprehensive tests (passing).” Parent C++ suite `tests/test_dataset_runner.cpp` still exists and includes `"DatasetRunner.h"`; with the stub header it cannot compile, so `DatasetRunnerTests` is **broken**, not green.
- Why it matters for science/repro: False green signals train agents and humans to trust the fleet path. AGENTS.md requires `ctest --output-on-failure` / pytest evidence before “done.” This commit provides neither.
- Fix recommendation: Any audit of fleet tests must use post-`d842e3247` (or current) `tests/test_fleet_benchmark.py` / `tests/test_fleet_runner.cpp` / `tests/test_fleet_dashboard.py`. Never cite this commit’s test file as coverage.

### F4. Ranking / scoring / election: no intentional algorithm change; path is build-broken (CRITICAL for product; INFO for algorithm delta)
- Evidence:
  - `git diff 155ebbb7b6a 69aa0fab6 -- LIB/DatasetRunner.cpp` → **0 bytes** (selector `select_pose_freq_gated*`, Z+H composite, `kBenchmarkPoseLimit = 50`, CF ranking, BCR, seed_echo, success_pb wiring all unchanged in `.cpp`).
  - Header deletion removes the **declarations** for the entire ranking-facing API, including parent success-gate contract comments and fields (`success_rmsd`, `pb_pass`, `success_pb`, `claim_ready`, `elected_cf`, `cf_best_cluster`, `best_score`, Shannon thermo columns, etc.).
  - No fleet aggregator exists that could re-sort poses or recompute ΔG / soft-β G̃. No change to GA population/generations, CF scoring (`Vcontacts` / VoronoiCF), StatMech, or BindingMode.
- Why it matters for science/repro: **Algorithm neutrality holds only in the vacuous sense** that no ranking formula was edited. Operationally, ranking is **not** neutral: the DatasetRunner product path does not link. AGENTS.md “Preserve current ranking” is not violated by a new election formula, but the commit **destroys the ranking API surface** and any ability to emit ranked claim artifacts from this tree.
- Fix recommendation: For ranking-neutrality claims about fleet, audit **`d842e3247`+** (serialization / control plane) and later election commits (`c82e6fc24` Shannon G̃ election, SoftBeta identity, FO dual-suffix). Do **not** treat `69aa0fab6` as a ranking-neutral fleet backend.

### F5. Thermodynamic / CF language in destroyed header was already soft (LOW residual, pre-existing)
- Evidence: Parent header labels `best_score` as “FlexAIDdS free energy (kcal/mol)” and `predicted_dG` similarly — pre-existing naming tension with AGENTS.md CF-proxy discipline. This commit does not introduce that wording; it **deletes** it with the header. No new ΔG claims appear in the stub files.
- Why it matters: Residual risk is only that restorers reintroduce the same soft naming when copying from parent. Real fleet docs later (`d842e3247` `docs/FLEET_BENCHMARKS.md`) should keep CF / soft-β language precise.
- Fix recommendation: Out of scope for this commit’s code (there is none). Track CF naming under `033eeb889` and SoftBeta identity commits.

### F6. Residual stub pollution still present on current trees (HIGH hygiene)
- Evidence (workspace snapshot after later real fleet work):

  | Path | Status after later commits |
  |------|----------------------------|
  | `LIB/DatasetRunner.h` | Restored full header (~32 KB) |
  | `LIB/FleetRunner.cpp` | Real implementation (~12 KB, from `d842e3247`) |
  | `tests/test_fleet_benchmark.py` | Real tests (~13 KB, from `d842e3247`) |
  | `docs/FLEET_BENCHMARKS.md` | Real docs (~5.6 KB) |
  | **`python/flexaidds/cli.py`** | **Still 43-byte stub** from `69aa0fab6` only |
  | **`tools/metal_microbench_enhanced.cpp`** | **Still 41-byte stub** |
  | **`benchmarks/m3pro/fleet_dataset_runner.py`** | **Still 105-byte stub** |
  | **`visualizer/ligand_fleet_pose_viz.py`** | **Still 181-byte stub** |

  `git log -- python/flexaidds/cli.py` shows **only** `69aa0fab6`. These files are non-functional and can confuse path-based imports (`import flexaidds.cli`) or CMake globs if ever wired.
- Why it matters for science/repro: Dead stubs look like features in tree listings and agent file maps. They are not production modules; they are comment debris from a false feat commit.
- Fix recommendation: Delete the four residual stubs (or replace with real implementations) in a dedicated hygiene commit; ensure nothing in CMake/package data references them. Not done in this audit (no source edits per task scope).

### F7. Ancestry / blast radius on main (HIGH process)
- Evidence: `git merge-base --is-ancestor 69aa0fab6 HEAD` → **YES** on lineages that contain the fleet merge series. Child `2fc7189d8` (“fix(build): Update CMakeLists…”) **does not restore** the header (still 85 bytes) and itself stubs `CMakeLists.txt` massively (−2942 lines). PR #257 merge `8eaf043ae` still has stub `DatasetRunner.h` (85) and stub `CMakeLists.txt` (556 bytes). Full header returns on the fleet feature line around merge with master `a4e13fc9b` / PR #258 `964bec0a2` (header 30895). Production fleet code arrives later as `d842e3247`.
- Why it matters: The commit is permanently in history. Any bisect landing on `69aa0fab6`…`8eaf043ae` is a known **build black hole** for DatasetRunner.
- Fix recommendation: Document in fleet handoff / INDEX that `69aa0fab6` and `2fc7189d8` / PR #257 tip are **non-buildable placeholders**. Prefer `d842e3247` as “fleet v1 production” baseline for science ops.

### F8. iCloud / local-first policy (INFO)
- Evidence: Commit message mentions “sharded iCloud-orchestrated benchmarks”; stubs contain no I/O paths, no `Mobile Documents` walks, no sync loops. No violation of local-first rules because no I/O code exists.
- Fix recommendation: N/A for this SHA. Real fleet I/O audited under later fleet + iCloud commits (`d842e3247`, `ab0850a6d`, `f75cdfc50`, etc.).

### F9. Security / licensing (INFO)
- Evidence: Comment-only files; no network clients, secrets, or GPL imports. No machine-absolute paths. Residual stubs are still tree noise, not license contaminants.
- Fix recommendation: Delete residual stubs (F6).

### F10. DoF / GA budget (INFO)
- Evidence: No touch of `NUMCHROM` / `NUMGENER` / `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` / ProtocolConfig. GA budget unaffected.
- Fix recommendation: None.

## Ranking/scoring impact: NO (algorithm) / YES (product break)

| Axis | Impact |
|------|--------|
| CF scoring formula / VoronoiCF | Unchanged (cpp untouched) |
| Pose election / Z+H / freq-gate | Unchanged in source; **unusable** at this SHA |
| Success gates / claim_ready | Definitions **deleted** from header |
| Fleet re-ranking / aggregation | **Does not exist** in this commit |
| GA search budget | Unchanged |

**Conclusion:** No ranking-algorithm delta to measure; ranking neutrality for a “fleet backend” is **not claimable** because the backend is fiction and DatasetRunner does not build.

## Reproducibility impact: YES (negative)

- **Negative:** Destroys the compile contract for the primary benchmark runner; false commit message; empty “tests”; propagates stubs via PR #257; residual dead files remain on modern trees.
- **Positive:** None for science. (Later commits repair the tree.)

## Tests adequate: NO

- Fleet “tests”: zero assertions.
- DatasetRunner GTest: cannot compile against stub header.
- No CMake fleet targets at this commit.
- Message claim “Comprehensive tests (passing)” is **false**.

## Verdict: REJECT

This commit must **not** be treated as a valid feature integration, ranking-neutral fleet backend, or tested deliverable. It is a **CRITICAL build-breaking placeholder** that:

1. Guts `LIB/DatasetRunner.h`,
2. Adds non-code stubs while advertising full fleet + tests + viz,
3. Leaves residual stub files on current trees (`python/flexaidds/cli.py`, `tools/metal_microbench_enhanced.cpp`, `benchmarks/m3pro/fleet_dataset_runner.py`, `visualizer/ligand_fleet_pose_viz.py`).

**Use instead for fleet science:** `d842e3247` (and subsequent fleet/iCloud hardening). **Use for DatasetRunner ranking evolution:** later election/SoftBeta/FO dual-suffix commits — never this SHA.

### Audit method (inspect-first)
- `git show 69aa0fab6 --stat --numstat --format=fuller`
- Full blob dumps of all 8 paths at `69aa0fab6` and sizes vs parent `155ebbb7b6a`
- `git diff` scope: only header gut + 7 new stubs; **0** `DatasetRunner.cpp` delta
- CMake `rg` for fleet targets at commit (none)
- Ancestry: merge path through PR #257 stubs → PR #258 header restore → `d842e3247` real fleet
- Parent header success-gate / election surface inventory; cpp selector symbols confirmed present and unchanged at commit tip of `.cpp`
- Residual stub check on modern workspace trees

### Parent / follow-on SHAs (navigation)

| SHA | Role |
|-----|------|
| `155ebbb7b6a` | Parent (dependabot checkout); full `DatasetRunner.h` |
| **`69aa0fab6`** | **This audit — stub “feat”** |
| `2fc7189d8` | Follow-on CMake stub damage; header still 85 bytes |
| `8eaf043ae` | PR #257 merge; still stubs |
| `a4e13fc9b` / `964bec0a2` | Master merge / PR #258 — full header restored on fleet line |
| `d842e3247` | Real production Bonhomme Fleet orchestration + real tests |
