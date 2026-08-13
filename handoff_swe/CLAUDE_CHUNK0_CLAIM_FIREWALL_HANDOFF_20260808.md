# Claude handoff: Chunk 0 statistical-mechanics claim firewall

Prepared: 2026-08-08T21:52:11Z  
Shannon task: `flexaidds_chunk0_claim_firewall_20260808`  
Recommended executor: Claude Code (`implement_only`)  
Alternative: Claude Dispatch coordinating disjoint ownership lanes

## Mission and stop condition

Finish the already-started **Chunk 0 scientific-provenance and claim-firewall
remediation** in the live FlexAIDdS checkout. Chunk 0 is metadata,
serialization, claim presentation, and fail-closed integration work. It must
not silently change statistical-mechanics numeric kernels, GA ranking,
clustering, or output ordering.

Chunk 0 is complete only when:

- proxy/contact-function ensembles cannot be presented as calibrated canonical
  thermodynamics or physical binding affinity in any supported consumer;
- physical claims require schema-v2 provenance, strict structured evidence,
  `available == true`, and the appropriate reference state;
- C++, Python, Swift, TypeScript, PDB/JSON loaders, UI/intelligence surfaces,
  docs, and the validator agree on that contract;
- known invalid vibrational/eigenvalue wiring is disabled rather than
  fabricated;
- producer-to-consumer fixtures prove wire interoperability;
- all relevant fresh build/test gates pass after the final edit;
- only Chunk 0 changes are committed and pushed, with user-owned work excluded.

Do not begin Chunk 1 scientific-model corrections in the same commit. In
particular, raw finite-sample partition sums, CF-to-physical calibration,
standard-state binding cycles, GPF normalization, mutual-information sign,
multi-site coupling, CCBM coordinate scoring, and Shannon entropy composition
belong to later chunks unless a minimal fail-closed guard is required here.

## Live snapshot

- Checkout: `/Users/lp.more/Projects/FlexAIDdS`
- Branch: `fix/statmech-claim-firewall`
- HEAD at freeze: `dfb99308c7d6074d54418caeb63a9bb9d0c1f4e3`
- `origin/main` at freeze: `ea044869636529d157356444e6aa6463446603d0`
- Divergence: local HEAD is 0 ahead / 1 behind `origin/main`.
- The feature branch has no configured upstream and
  `origin/fix/statmech-claim-firewall` did not exist at freeze.
- Index: empty; nothing staged.
- Worktree: heavily dirty by design; all implementation is uncommitted.
- `git diff --check`: passed at freeze.
- `python3 -m py_compile python/flexaidds/models.py`: passed at freeze.
- All Codex workers were interrupted and returned no-edit final receipts.
- Full tracked binary-diff SHA-256:
  `f5e218232946ab27a8330df258edbbf365a0f584a02745faa847cfcaaff27513`.
- `git status --short` SHA-256:
  `23945d7c0680ae6c341db226af0d3014bf29da38f337ffc70ad5f89683bba945`.
- Task-owned untracked content-set SHA-256 (handoff excluded):
  `870d0d713c6cf3186a5a2fd4d3a7bfca41aa4a30ab9ca9cb5734f0d3171a32db`.

Re-inspect every live file before editing. This document is a transfer receipt,
not permission to trust stale line numbers.

### Shannon receipt and lifecycle

The dry-run plan returned `ok: true`, `pair_mode: implement_only`, agent
`claude_code`, role `implement`, slice `full`, with both fabricated-code and
fabricated-review fields null. No live agent was attached.

At the final check, `gate-status` reported `/tmp/shannon.sock`, but `monitor`
returned connection refused. Treat the socket as stale. Before registering a
live receiver:

```bash
cd /Users/lp.more/Projects/Shannon
./scripts/shannon gate
```

In another terminal, require a successful monitor before spawn:

```bash
cd /Users/lp.more/Projects/Shannon
export PYTHONPATH="$PWD/hub${PYTHONPATH:+:$PYTHONPATH}"
python3 -m agent_manager monitor
python3 -m agent_manager spawn claude_code \
  --task flexaidds_chunk0_claim_firewall_20260808 \
  --reason "Resume verified Chunk 0 handoff"
python3 -m agent_manager control claude_code \
  "handoff loaded; live reinspection started" \
  --task flexaidds_chunk0_claim_firewall_20260808
```

During work, send real phase updates with `control`; submit only real evidence
with `result`. On completion or abort, detach with:

```bash
python3 -m agent_manager kill claude_code \
  --task flexaidds_chunk0_claim_firewall_20260808 \
  --reason "handoff session ended"
```

## Authority and safety boundary

Read `AGENTS.md`, then `METHODOLOGY.md`, before editing. Their rules are
binding. In particular:

- do not use `git reset --hard`, `git checkout --`, rebase, merge, force-push,
  branch deletion, or history rewriting;
- do not stash or overwrite the dirty worktree;
- preserve default pose ranking, clustering, and output order;
- keep new thermodynamic behavior fail-closed and behind tests/flags;
- do not claim physical `Delta G` from CF/contact-function proxy data;
- use `apply_patch` for edits;
- run full relevant tests after the last edit and before commit/push;
- before push, verify `gh api user --jq .login` is exactly
  `LeBonhommePharma`;
- do not commit secrets, absolute machine paths, benchmark output directories,
  or ignored local fixtures.

A request to prepare this handoff did not authorize a merge/rebase or any
destructive cleanup.

## User-owned pre-existing changes: preserve exactly

These were present before Chunk 0 and are not part of the claim-firewall commit:

- `LIB/CleftDetector.cpp`: user-owned cleft-dump work; do not stage it.
- `LIB/DatasetRunner.cpp` has mixed ownership:
  - user hunk near the run-command assembly adds per-target
    `FLEXAIDDS_CLEFT_DUMP` forwarding;
  - user hunk near report CSV output adds
    `election_mode,consensus_count,rank0_demoted` columns;
  - Chunk 0 hunks include `DatasetThermoLog.h` and current/legacy thermo-label
    parsing. Use partial staging.
- `LIB/gaboom.cpp` has mixed ownership:
  - the large `FLEXAIDDS_GENTRACE` block near the beginning of `GA()` is
    user-owned;
  - later proxy/provenance labels are Chunk 0. Use partial staging.
- Untracked result/work directories are user-owned and must remain untracked:
  - `ab_mac_20260806T133329/`
  - `determinism_e1/`
  - `determinism_e1_INVALID_wrongdataroot/`
  - `interventions/`
- A local Git-excluded Python test,
  `python/tests/test_validate_benchmark_results.py`, is known to fail because
  it expects a missing `docking_mode`. It is not part of Chunk 0; do not edit
  or stage it. Run tracked Python tests explicitly when validating the package.
- This handoff file is a local operational artifact and contains the exact
  machine checkout path. Leave it untracked for the Chunk 0 commit unless LP
  explicitly requests a sanitized, repository-portable handoff document.

### Exhaustive ownership manifest at freeze

The only user-only tracked file is `LIB/CleftDetector.cpp`. The only
mixed-ownership tracked files are `LIB/DatasetRunner.cpp` and `LIB/gaboom.cpp`,
with hunk ownership defined above. Every other modified tracked path below is
Chunk 0 task-owned:

```text
.github/workflows/ci.yml
LIB/BindingMode.cpp
LIB/MultiModelDock.cpp
LIB/ParallelDock.cpp
LIB/cluster.cpp
LIB/statmech.cpp
LIB/statmech.h
LIB/top.cpp
README.md
docs/SCORING.md
docs/USERGUIDE.md
docs/dev/thermo_invariants.md
docs/dev/thermo_source_map.md
docs/entropy-help/DOMAIN_SETUP.md
docs/entropy-help/MANIFESTO.md
docs/entropy-help/README.md
docs/entropy-help/THERMODYNAMIC_OUTPUT_SCHEMA.md
docs/entropy-help/audit-report-example.json
docs/entropy-help/audit-report-template.md
docs/entropy-help/audits/audits.json
python/bindings/core_bindings.cpp
python/flexaidds/__init__.py
python/flexaidds/_core.cpp
python/flexaidds/models.py
python/flexaidds/thermodynamics.py
python/tests/test_import_fallback.py
python/tests/test_thermodynamics.py
python/tests/test_thermodynamics_dataclass.py
site/entropy-help/index.html
site/entropy-help/ledger.html
site/entropy-help/request.html
swift/Package.swift
swift/Sources/FleetScheduler/DeviceCapability.swift
swift/Sources/FleetScheduler/FleetScheduler.swift
swift/Sources/FlexAIDCore/FXGA.mm
swift/Sources/FlexAIDCore/FXStatMechEngine.mm
swift/Sources/FlexAIDCore/include/FXStatMechEngine.h
swift/Sources/FlexAIDCore/include/FXTypes.h
swift/Sources/FlexAIDdS/DockingRunner.swift
swift/Sources/FlexAIDdS/FlexAIDRunner.swift
swift/Sources/FlexAIDdS/Models.swift
swift/Sources/FlexAIDdS/ThermodynamicResult.swift
swift/Sources/HealthIntegration/BindingEntropyScore.swift
swift/Sources/Intelligence/BindingModeNarrator.swift
swift/Sources/Intelligence/CampaignJournalist.swift
swift/Sources/Intelligence/CleftAssessor.swift
swift/Sources/Intelligence/ConvergenceCoach.swift
swift/Sources/Intelligence/FleetExplainer.swift
swift/Sources/Intelligence/HealthEntropyInsight.swift
swift/Sources/Intelligence/IntelligenceOracle.swift
swift/Sources/Intelligence/LigandFitCritic.swift
swift/Sources/Intelligence/SelectivityAnalyst.swift
swift/Sources/Intelligence/ThermoReferee.swift
swift/Sources/Intelligence/ThermoRefereeTools.swift
swift/Sources/Intelligence/VibrationalInterpreter.swift
swift/Tests/FlexAIDdSTests/FleetSchedulerTests.swift
swift/Tests/FlexAIDdSTests/IntelligenceFeatureTests.swift
swift/Tests/FlexAIDdSTests/StatMechEngineTests.swift
swift/Tests/FlexAIDdSTests/ThermoRefereeTests.swift
tests/test_binding_mode_statmech.cpp
tests/test_dataset_runner.cpp
tests/test_statmech.cpp
typescript/apps/viewer/package.json
typescript/apps/viewer/src/App.tsx
typescript/apps/viewer/src/FleetDashboard.tsx
typescript/apps/viewer/src/IntelligenceEngine.ts
typescript/apps/viewer/src/IntelligencePanel.tsx
typescript/apps/viewer/src/MolstarViewer.tsx
typescript/apps/viewer/src/RefereePanel.tsx
typescript/apps/viewer/src/__tests__/IntelligenceAnalyzers.test.ts
typescript/apps/viewer/src/__tests__/IntelligenceEngineReferee.test.ts
typescript/apps/viewer/src/intelligence/BindingModeAnalyzer.ts
typescript/apps/viewer/src/intelligence/CleftAnalyzer.ts
typescript/apps/viewer/src/intelligence/ConvergenceAnalyzer.ts
typescript/apps/viewer/src/intelligence/PoseQualityAnalyzer.ts
typescript/apps/viewer/src/intelligence/SelectivityAnalyzer.ts
typescript/packages/flexaidds/src/StatMechEngine.ts
typescript/packages/flexaidds/src/index.ts
typescript/packages/flexaidds/src/resultLoader.ts
typescript/packages/flexaidds/src/types.ts
typescript/packages/shared/package.json
typescript/packages/shared/src/BindingPopulation.ts
typescript/packages/shared/src/PoseQualityContext.ts
typescript/packages/shared/src/RefereeVerdict.ts
typescript/packages/shared/src/SelectivityAnalysis.ts
typescript/packages/shared/src/SelectivityContext.ts
typescript/packages/shared/src/index.ts
```

Task-owned untracked paths included in the content-set checksum are:

```text
LIB/DatasetThermoLog.h
scripts/validate_thermo_claims.py
swift/Sources/FlexAIDCore/FXStatMechBridgeInternal.hpp
swift/Sources/FlexAIDdS/ScientificProvenance.swift
swift/Sources/FlexAIDdS/ShannonEntropyDecomposition.swift
swift/Tests/FlexAIDdSTests/ScientificProvenanceTests.swift
tests/test_thermo_claim_firewall.py
typescript/apps/viewer/src/__tests__/ClaimPresentation.test.ts
typescript/apps/viewer/src/claimPresentation.ts
typescript/apps/viewer/tsconfig.json
typescript/apps/viewer/vitest.config.mjs
typescript/package-lock.json
typescript/package.json
typescript/packages/flexaidds/src/scientificClaims.test.ts
typescript/packages/shared/scientificProvenance.test.mjs
```

If the live checksums differ before Claude edits, stop and reclassify the drift;
do not assume a new path belongs to Chunk 0.

## Implemented work already present

Treat this as implemented-but-not-fully-verified because later edits make some
earlier test receipts stale.

### C++ provenance core and producers

- `LIB/statmech.h/.cpp` define schema-v2 `ScientificProvenance`, exact enum
  vocabulary, derived fail-closed claim validity, strict SHA-256 receipt syntax,
  result/ledger propagation, and conservative merge downgrade.
- Contact-function producers use explicit proxy provenance in BindingMode,
  ParallelDock, MultiModelDock, top-level post-scoring, gaboom, and classic
  clustering.
- C++/PDB/stdout labels were changed from physical-looking thermodynamic claims
  to proxy/diagnostic language while retaining legacy numeric fields for
  migration.
- Latest interrupted C++ patch added `provenance_for_breakdown()` in
  `LIB/statmech.cpp`: any enabled or nonzero unreceipted correction downgrades
  breakdown provenance to proxy-only.
- Static `make_breakdown()` no longer returns an authorized all-zero ledger for
  an empty engine; it now reaches `compute()` and should throw consistently.
- `LIB/DatasetThermoLog.h` and `tests/test_dataset_runner.cpp` now parse exact
  current proxy labels plus legacy labels. `LIB/DatasetRunner.cpp` consumes the
  helper. These edits have not been compiled.

### Python facade

- `python/flexaidds/thermodynamics.py`, both native-binding surfaces, and public
  imports carry the schema-v2 contract while keeping a stable Python facade.
- Serialization emits nested snake-case `scientific_provenance`, writes the
  corrected `heat_capacity_kcal_mol_K` key, reads the legacy `_K2` key, ignores
  serialized `claim_validity`, and rejects hostile schema/evidence types.
- The current partial edit to `python/flexaidds/models.py` adds provenance,
  `proxy_free_energy`, and `soft_beta_G` fields and begins serialization/ranking
  propagation. It compiles syntactically but is incomplete across constructors,
  `io.py`, and `results.py`.

### Swift

- Schema-v2 provenance, bridge metadata, fail-closed claim helpers, and proxy
  presentation gates are present across FlexAIDdS, FleetScheduler,
  HealthIntegration, ThermoReferee, and IntelligenceOracle.
- The latest Swift pass also gated `BindingModeNarrator`,
  `CampaignJournalist`, `SelectivityAnalyst`, and `VibrationalInterpreter`,
  including deterministic/FoundationModels/follow-up paths, and added hostile
  tests in `IntelligenceFeatureTests.swift`.
- Swift package header paths were repaired. The actual C++ core is still not
  linked into SwiftPM tests.

### TypeScript and viewer

- SDK/shared provenance types, helpers, viewer presentation gates, tests,
  workspace package metadata/lockfile, and CI wiring exist.
- Existing TypeScript tests are green, but they encode known fail-open and
  wire-shape gaps listed below. No assigned blocker-fix edit landed before the
  worker was interrupted.

### Public docs/site and validator

- entropy.help hard-coded “published audits” were downgraded to
  `PLANNED_UNVERIFIED`/`EXAMPLE_UNVERIFIED`; fake hashes/signatures were removed.
- `scripts/validate_thermo_claims.py` and
  `tests/test_thermo_claim_firewall.py` exist.
- Major scoring/user/developer docs were relabeled, but README and validator
  coverage still contain blockers below.

## Interrupted or incomplete edits

1. `LIB/statmech.cpp/.h`: the latest empty-engine and correction-downgrade
   changes are structurally complete but have no tests and have not compiled.
2. `LIB/BindingMode.cpp`: the assigned remediation did not land. It still:
   - adds vibrational/NATURaL corrections only to F, violating `F = H - TS`;
   - reads `atoms[0].eigen[m][0]` as an eigenvalue even though production stores
     XYZ eigenvectors on real atoms and leaves sentinel atom 0 null.
3. `python/flexaidds/models.py`: partial field/serialization/ranking work is
   present. Complete every construction/round-trip path before testing.
4. TypeScript blocker fixes did not start; current green tests are insufficient.
5. Latest Swift source build is stale by one final Selectivity wording edit.

## Remaining blockers and acceptance tests

### P0: C++ ledger fail-closed behavior

- Add tests proving empty physical `make_breakdown()` throws.
- Add tests proving any unreceipted correction flag or nonzero correction value
  downgrades total provenance while leaving all numeric formulas unchanged.
- Preserve the base `Thermodynamics` numeric parity test.

### P0: disable invalid BindingMode vibration path

- Make `compute_vibrational_correction()` return zero/fail closed until a real
  eigenvalue channel exists; do not reinterpret eigenvector X components.
- Stop adding vibration/NATURaL only to scalar `Thermodynamics`. Keep the base
  CF-proxy ledger internally coherent. If corrections remain visible, expose
  them only as explicit proxy breakdown diagnostics.
- Rewrite the fabricated atom-0 sentinel test to prove the invalid layout is
  rejected/ignored. Ranking must remain unchanged.

### P0: finish Python PDB/JSON integration

- Parse `proxy_free_energy`, `soft_beta_G`, schema/domain/measure/reference, and
  claim validity/evidence fields from PDB REMARKs.
- Build `ScientificProvenance` from source fields; absent receipts must remain
  proxy-only even if a serialized claim says otherwise.
- Propagate the new fields through every `PoseResult`, `BindingModeResult`,
  `DockingResult.from_dict/from_json/to_json/to_records`, and loader constructor.
- Rank new outputs by emitted `soft_beta_G`; preserve legacy free-energy order
  only when the election field is absent.
- Add a two-mode PDB fixture whose soft-beta order is opposite legacy F and
  assert the emitted election order wins.
- Rewrite model/results examples and docstrings so units are conditional on
  provenance, not automatically kcal/mol or physical Delta G.

### P0: strict TypeScript availability and wire normalization

- Both SDK and shared predicates must require literal `available === true`.
  Missing, false, string, numeric, or malformed values fail closed.
- Create one runtime normalization boundary accepting:
  - SDK camel-case `bindingModes`;
  - shared/viewer `modes`;
  - Python `binding_modes`, unit-suffixed thermodynamic keys, and nested
    snake-case provenance.
- Every normalized accepted record must set availability deliberately.
- Add producer-to-normalizer-to-viewer golden tests, including Python output
  and hostile payloads.
- Remove the viewer's unsupported `.rrd` advertisement unless a real parser is
  implemented.
- For proxy selectivity, return public driver `inconclusive` and suppress
  untrusted `explanation`/`designSuggestion` affinity or potency prose.
- Inspect `resultLoader.ts` against HEAD and restore unrelated numeric changes
  (unknown temperature and cluster parsing/count behavior) unless explicitly
  justified and separately tested. Chunk 0 must preserve numerics.

### P0: public claims and validator coverage

- Fix README claims that present `G_bind`/`dG_eff` as physical Delta G or an
  active ranking criterion.
- Remove or receipt the `78/85 = 91.8%` reproducibility claim.
- Correct the feature-flag description and the nonexistent Python
  `temperature=`/`add_energies` example.
- Expand validator coverage beyond `site/entropy-help` and
  `docs/entropy-help`, including README and relevant scoring docs.
- Bind each quantitative/completion claim to its own claim ID/evidence; the
  presence of one unrelated valid ID must not authorize a page.
- Reject signature-like keys such as `detached_signature`, not only exact
  `signature`.
- Add hostile tests for both loopholes.

### P0/P1: Swift completion

- Run `swift build --target Intelligence` after the last wording edit.
- Re-run hostile checks for the four latest Intelligence surfaces.
- Decide whether the retained numeric Selectivity `deltaG` compatibility field
  needs explicit availability metadata; it must not authorize a physical claim.
- Fix SwiftPM to link the real repository C++ implementation before calling
  XCTest green. Do not ship temporary `-undefined dynamic_lookup`, stubs, or
  fake symbols. Current XCTest link failure includes `StatMechEngine`,
  `BoltzmannLUT`, BindingMode/Population, ENCoM/tENCoM, Shannon stack, `GA`,
  `read_input`, and `ic2cf`.

### P1: inspect remaining physical-looking consumers

Audit and either gate or explicitly relabel:

- `LIB/ReferenceEntropy.cpp/.h`;
- `LIB/ParallelCampaign.cpp`;
- `LIB/benchmark_datasets.cpp`;
- DatasetRunner conversion of the proxy prediction to pKd/affinity metrics.

Do not invent calibration. If no calibrated energy/reference provenance exists,
the safe Chunk 0 behavior is unavailable/proxy diagnostic, not affinity.

### Explicit scope note: WASM

Public WASM parity is deferred. Do not advertise WASM as provenance-safe in
Chunk 0. If it remains publicly supported, `typescript/wasm/bindings.cpp` is a
blocker because it exposes physical-looking numeric fields without provenance.

## Verification ledger at freeze

These are historical receipts, not blanket current verification:

| Surface | Last observed result | Current? |
|---|---|---|
| Full CTest | 87/87 passed before later C++ edits | **Stale** |
| Direct `test_statmech` | 101/101 passed before latest breakdown edits | **Stale** |
| Python tracked suite | 1178 passed, 69 skipped before `models.py` partial edit | **Stale** |
| Python focused provenance | 63 passed, 17 skipped before `models.py` partial edit | **Stale** |
| TypeScript typecheck/tests | SDK 8, shared 5, viewer 80 passed | Current for pre-fix code, but coverage is insufficient |
| Swift `FlexAIDdS` target | built in an earlier pass | Stale relative to all later Swift edits |
| Swift `Intelligence` target | built before final one-line wording edit | **Stale** |
| Swift XCTest | sources compile; final link fails on real C++ symbols | **Blocked/failing** |
| Claim validator | passed current documents before loophole fixes | Green but insufficient |
| Firewall pytest | 8 passed | Green but insufficient |
| `git diff --check` | passed at freeze | Current |
| Python `models.py` syntax | `py_compile` passed at freeze | Current syntax only |

## Required execution order

Keep exactly one plan item in progress.

1. Re-read instructions, inspect status/diffs, and run cheap syntax/build probes
   on interrupted C++/Python/Swift files.
2. Finish C++ statmech/BindingMode/DatasetRunner blockers and focused tests.
3. Finish Python model/PDB/JSON propagation and focused round-trip/ranking tests.
4. Finish TypeScript strict availability, normalization, presentation, and
   producer-to-viewer golden tests.
5. Finish Swift build/link/presentation gates.
6. Fix README/public-doc claims and validator loopholes.
7. Audit the remaining physical-looking consumers listed above.
8. Run all fresh combined gates.
9. Inspect the complete diff and partial-stage mixed-ownership files.
10. Commit/push Chunk 0 only. Stop before Chunk 1 and report exact evidence.

## Fresh final gates

Use a fresh C++ build because sources/build configuration changed. Defer to
`METHODOLOGY.md §0` exactly: Homebrew CMake/CTest and `-j4`, not PATH-dependent
tools or all logical CPUs.

```bash
export PATH="/opt/homebrew/bin:$PATH"
claim_build_dir="$(mktemp -d /tmp/flexaidds-claim-firewall.XXXXXX)"
/opt/homebrew/bin/cmake -S . -B "$claim_build_dir" \
  -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
/opt/homebrew/bin/cmake --build "$claim_build_dir" -j4
/opt/homebrew/bin/ctest --test-dir "$claim_build_dir" --output-on-failure
```

Run Python focused tests first, then all tracked package tests so the local
excluded fixture is not collected accidentally:

```bash
python3 -m pytest \
  python/tests/test_results_io.py \
  python/tests/test_results_loader_models.py \
  python/tests/test_thermodynamics.py \
  python/tests/test_thermodynamics_dataclass.py \
  python/tests/test_thermo_breakdown.py -q

git ls-files 'python/tests/test_*.py' -z | \
  xargs -0 python3 -m pytest -q
```

Run TypeScript and Swift gates:

```bash
(cd typescript && npm ci && npm run typecheck && npm test)
(cd swift && swift build --target FlexAIDdS && swift build --target Intelligence)
(cd swift && swift test)
```

Run claim/hygiene checks:

```bash
python3 scripts/validate_thermo_claims.py
python3 -m pytest tests/test_thermo_claim_firewall.py -q
python3 scripts/check_repo_hygiene.py
git diff --check
```

Do not waive a red tracked test. If Swift remains unable to link the real core,
Chunk 0 is not fully verified; report the blocker rather than using stubs.

### Implementation gates versus OPS merge gates

Claude owns the implementation/build/unit/integration gates listed above and may
push the verified feature branch. Claude does **not** launch benchmark,
determinism, parity, or Astex campaigns under this handoff. OPS owns the
canonical pre-merge gates in `METHODOLOGY.md §§1–4 and §6` (including default-path
parity and any applicable determinism/science gates). Therefore a successful
feature-branch push is not, by itself, a merge-ready or scientific-validation
claim. Do not invent or abbreviate those procedures in this document.

## Commit and push protocol

1. Inspect `git diff --stat`, `git diff`, and every untracked candidate.
2. Leave `LIB/CleftDetector.cpp` and the four result directories unstaged.
3. Leave this local handoff file unstaged unless explicitly requested otherwise.
4. Use partial staging for `LIB/DatasetRunner.cpp` and `LIB/gaboom.cpp` so user
   hunks are excluded.
5. Inspect `git diff --cached --check` and the full cached patch.
6. Verify GitHub identity:

   ```bash
   gh auth status
   gh api user --jq .login
   ```

   The login must be `LeBonhommePharma`.
7. Use a conventional commit such as:
   `Fix: add fail-closed thermodynamic claim provenance`.
8. The branch has no upstream and no same-named remote branch. Push it without
   force using:

   ```bash
   git push --set-upstream origin HEAD
   ```

   Do not rebase/merge the one upstream commit without explicit user approval.

## Ready-to-paste Claude Code prompt

```text
Resume the in-flight FlexAIDdS Chunk 0 claim-firewall remediation from:
/Users/lp.more/Projects/FlexAIDdS/handoff_swe/CLAUDE_CHUNK0_CLAIM_FIREWALL_HANDOFF_20260808.md

You are the sole implementer under Shannon task
flexaidds_chunk0_claim_firewall_20260808. Read AGENTS.md, METHODOLOGY.md, and the
entire handoff before editing. Re-inspect the live dirty checkout; preserve all
listed user-owned changes and never reset/checkout/stash them. Use apply_patch.
Keep exactly one plan item in progress.

Finish every P0/P1 item in the handoff chunk by chunk with focused tests between
chunks. Preserve statistical-mechanics numerics, GA ranking, clustering, and
output order unless the handoff explicitly requires loader parity with the
already-emitted soft_beta_G election field. Fail closed whenever calibrated
energy, ensemble-measure, availability, evidence receipts, or matched reference
state is missing. Do not invent physical Delta G or affinity.

After the final edit, run the fresh full C++/CTest, tracked Python, TypeScript,
Swift, validator, and hygiene gates specified in the handoff. Do not use Swift
stubs or dynamic lookup to fake a passing test. Partial-stage only Chunk 0,
verify LeBonhommePharma identity, commit/push the feature branch without force,
then stop before Chunk 1 and report exact commands, counts, and blockers.
```

## Ready-to-paste Claude Dispatch prompt

```text
Coordinate completion of the existing FlexAIDdS Chunk 0 claim firewall using:
/Users/lp.more/Projects/FlexAIDdS/handoff_swe/CLAUDE_CHUNK0_CLAIM_FIREWALL_HANDOFF_20260808.md

Read AGENTS.md, METHODOLOGY.md, and the entire handoff first. This is a dirty
shared checkout. Preserve the user-owned hunks/directories listed there. Do not
reset, stash, merge, rebase, or force-push.

Delegate disjoint ownership lanes only:
1. C++ core: statmech/BindingMode/DatasetThermoLog and their tests.
2. Python: thermodynamics/models/io/results/bindings and Python tests.
3. TypeScript: SDK/shared/viewer normalization, availability, presentation,
   and golden tests.
4. Swift: bridge/provenance/Intelligence and real SwiftPM linking/tests.
5. Docs/validator: README, scoring/entropy-help surfaces, validator/tests.

Tell every worker: you are not alone in the codebase; stay in owned files, do
not revert others, do not commit, and return changed files, exact tests/counts,
and blockers. Keep one integration owner for shared files, staging, combined
tests, and cross-language contract review. Do not accept package-local green as
producer-to-consumer evidence. Finish all handoff blockers, run every final gate,
partial-stage around user-owned DatasetRunner/gaboom hunks, verify GitHub login
LeBonhommePharma, commit/push Chunk 0 without force, and stop before Chunk 1.
```

## Claude Fable five-session pack

This pack supersedes the solo-executor instructions only for the five Fable
sessions below. It is **prepared but not launched**. The five Shannon delegate
plans were dry runs only; `ok: true` does not prove the live `science` identity
or machine resources are free.

### C++ is the scientific truth

For this pack, the canonical scientific contract is the live C++ implementation
and its native tests:

- `LIB/statmech.h`
- `LIB/statmech.cpp`
- `tests/test_statmech.cpp`

`BindingMode`, transport adapters, Python, TypeScript, Swift, docs, and UI
consumers must conform to that C++ contract. They may not redefine or weaken
its vocabulary, receipt syntax, claim predicates, reference-state requirements,
units, or numeric invariants for convenience. If a contract change is required,
Session 1 changes and tests C++ first, then emits a truth receipt and a
mandatory language-neutral golden fixture at
`tests/fixtures/scientific_provenance_v2_golden.json`. Session 1 publishes the
fixture's SHA-256 in its receipt. Downstream sessions must verify that exact hash
before consuming the fixture; a missing or changed fixture blocks them.

Do not infer physical truth from field names such as `free_energy`, `deltaG`,
`enthalpy`, or `kcal`. C++ provenance predicates decide what claims are
authorized. Current docking CF/contact-function ensembles remain proxy-only.

`available` is a separate record-integrity/transport gate, not a replacement
for C++ `ScientificProvenance`. A physical presentation is allowed only when:

1. the record carries literal JSON/host-language boolean `available == true`;
2. the C++-truth provenance predicate authorizes the requested canonical or
   binding claim.

Absent, null, false, numeric `1`, strings such as `"true"`/`"false"`, arrays,
and objects all fail closed. The mandatory golden fixture must freeze these
transport cases separately from the provenance cases. C++ computation already
fails by throwing when no result exists; adapters may set availability true only
after an actual result is computed and transported.

### Dependency graph and write ownership

```text
Session 1: C++ truth
        |
        +--> Session 2: Python mirror
        +--> Session 3: TypeScript mirror
        +--> Session 4: Swift mirror
                         |
Sessions 1-4 receipts ---+--> Session 5: docs/validator/integration
```

Run Session 1 first. Sessions 2–4 may begin only after its receipt exists.
Session 5 begins only after Sessions 1–4 have returned receipts. No session may
edit another session's lane. None may stage, commit, push, merge, rebase, stash,
or reset. One later Claude Code/Dispatch integration owner performs staging and
git operations after all gates.

Durable receipt root:
`handoff_swe/fable_chunk0_receipts_20260808/`

Exact receipt paths:

```text
handoff_swe/fable_chunk0_receipts_20260808/s1_cpp_truth.json
handoff_swe/fable_chunk0_receipts_20260808/s2_python.json
handoff_swe/fable_chunk0_receipts_20260808/s3_typescript.json
handoff_swe/fable_chunk0_receipts_20260808/s4_swift.json
handoff_swe/fable_chunk0_receipts_20260808/s5_integration.json
```

These are authorized local operational artifacts and remain untracked. Each
session writes only its own receipt, computes `/usr/bin/shasum -a 256`, and
reports the path plus hash in its final response. A receipt cannot contain its
own hash; the coordinator injects the exact external path/hash pair into every
dependent session prompt. The dependent session verifies the hash before work
and records every consumed path/hash in its own receipt. Session 5 must consume
and verify all four upstream receipts. Any other new untracked path is drift and
requires reclassification.

### Existing Claude Science experiment fence

- Do not register, control, detach, kill, or reuse any existing Shannon
  `science` or `dataset_runner` identity/task.
- Do not launch these five sessions by any route—Shannon, Fable UI, Dispatch,
  CLI, or another orchestrator—while existing Science/DatasetRunner experiments
  are active or their state is ambiguous. The canonical `science` ID represents
  one live agent, not five parallel sessions.
- Do not touch experiment processes, binaries, builds, outputs, receipts,
  result roots, iCloud/CloudDocs, `ab_mac_*`, `determinism_*`, or
  `interventions/`.
- Do not launch docking, DatasetRunner, parity, determinism, Astex,
  tENCoM/Eigen, or any scientific campaign.
- Do not use broad `kill`, `pkill`, reapers, or cleanup commands.
- Use only uniquely named `mktemp -d` build directories. Never rebuild into or
  delete an experiment's build directory.
- Before any compile/test with meaningful CPU, memory, disk, or I/O load, verify
  experiment safety. If state is ambiguous or Science is active, defer the
  command and record it in the receipt.

### Dry-run task receipts

```text
S1 flexaidds_c0_fable_s1_cpp_truth_20260808       ok=true, live_spawn=false
S2 flexaidds_c0_fable_s2_python_20260808          ok=true, live_spawn=false
S3 flexaidds_c0_fable_s3_typescript_20260808      ok=true, live_spawn=false
S4 flexaidds_c0_fable_s4_swift_20260808           ok=true, live_spawn=false
S5 flexaidds_c0_fable_s5_integration_20260808     ok=true, live_spawn=false
```

### Required receipt from every Fable session

```text
session_id:
receipt_path:
consumed_receipts: exact path + externally supplied SHA-256
files_changed:
contract_observed:
golden_fixture_path:
golden_fixture_sha256:
availability_transport_gate_verified: yes/no + cases
numeric_behavior_changed: yes/no + exact evidence
tests_run: exact commands and pass/fail counts
tests_deferred_for_experiment_safety:
unresolved_blockers:
out_of_lane_files_touched: none/list
commit_created: no
```

### Session 1 prompt: C++ truth

```text
You are Fable Session 1, the C++ scientific-truth owner for FlexAIDdS Chunk 0.
Task ID: flexaidds_c0_fable_s1_cpp_truth_20260808.

Read AGENTS.md, METHODOLOGY.md, and the complete handoff at:
handoff_swe/CLAUDE_CHUNK0_CLAIM_FIREWALL_HANDOFF_20260808.md

The live dirty checkout is authoritative. You are not alone in the codebase.
Do not revert, stage, commit, push, stash, reset, or overwrite another session's
work. Existing Claude Science experiments are untouchable; obey the experiment
fence and defer resource-heavy tests when necessary.

Your exclusive write ownership:
- LIB/statmech.h
- LIB/statmech.cpp
- LIB/BindingMode.cpp and its header only if strictly required
- LIB/MultiModelDock.cpp
- LIB/ParallelDock.cpp
- LIB/cluster.cpp
- the Chunk 0 claim/provenance hunks in LIB/gaboom.cpp (mixed-ownership file)
- LIB/top.cpp
- LIB/DatasetThermoLog.h
- the Chunk 0 parser hunks in LIB/DatasetRunner.cpp (mixed-ownership file)
- LIB/ReferenceEntropy.cpp and LIB/ReferenceEntropy.h
- LIB/ParallelCampaign.cpp and its header only if required
- LIB/benchmark_datasets.cpp
- tests/test_statmech.cpp
- tests/test_binding_mode_statmech.cpp
- tests/test_binding_mode_vibrational.cpp
- tests/test_dataset_runner.cpp
- tests/fixtures/scientific_provenance_v2_golden.json (mandatory)

Do not edit LIB/CleftDetector.cpp, the user GENTRACE hunk in LIB/gaboom.cpp, the
user cleft/CSV hunks in LIB/DatasetRunner.cpp, any Python/TypeScript/Swift/docs/
site/CI file, or any result directory.

Finish the native P0 items in the handoff:
1. Verify empty make_breakdown fails closed.
2. Downgrade provenance for every unreceipted correction while preserving the
   exact numeric formulas.
3. Disable the invalid atoms[0].eigen-as-eigenvalue vibrational path.
4. Keep scalar Thermodynamics internally coherent; do not add corrections only
   to F. Keep corrections as explicit proxy diagnostics if retained.
5. Verify DatasetRunner parses current proxy labels and legacy labels without
   changing numeric CSV behavior.
6. Freeze C++ schema-v2 vocabulary, SHA receipt predicates, known rejected
   digests, and claim validity into native tests and a downstream-consumable
   truth receipt plus the mandatory golden fixture at the exact path above.
   Include separate transport cases proving only literal available=true passes
   the record gate; provenance predicates still independently control claims.
7. Gate or explicitly relabel the remaining C++ physical-looking consumers:
   ReferenceEntropy, ParallelCampaign, benchmark_datasets, and DatasetRunner's
   proxy-to-pKd/affinity conversion. Never invent calibration.
8. Review every existing Chunk 0 C++ producer/output hunk for conformance to the
   frozen contract, including gaboom, top, classic clustering, ParallelDock,
   and MultiModelDock.

Ranking, clustering, and statistical-mechanics numeric kernels must not change.
Use METHODOLOGY.md section 0 tool paths and -j4 for any safe build. Run focused
native tests only if they cannot contend with Science; otherwise defer them.

Write your receipt to
handoff_swe/fable_chunk0_receipts_20260808/s1_cpp_truth.json, compute and report
its external SHA-256, and include the mandatory golden-fixture SHA-256. State
the exact C++ truth files Sessions 2-4 must consume. Do not commit.
```

### Session 2 prompt: Python mirror

```text
You are Fable Session 2, the Python mirror owner for FlexAIDdS Chunk 0.
Task ID: flexaidds_c0_fable_s2_python_20260808.

Start only after the coordinator provides the exact path and SHA-256 for
Session 1's C++ truth receipt. Verify that hash, then verify the golden fixture
hash recorded in it. Read AGENTS.md, METHODOLOGY.md, the complete handoff,
Session 1's receipt, and the fixture. C++ statmech is authoritative. Do not
change, reinterpret, or weaken the C++ contract.

You are not alone in the codebase. Your exclusive write ownership is python/**,
except the local Git-excluded python/tests/test_validate_benchmark_results.py.
Do not edit LIB/**, tests/** outside python/tests, TypeScript, Swift, docs, site,
CI, handoff files, or experiment outputs. Do not stage or commit.

Complete the Python P0 mirror:
1. Match C++ enums, strict schema version, SHA receipt rules, rejected digests,
   claim predicates, and correction downgrade exactly.
2. Finish the partial models.py edit across every constructor and round trip.
3. Parse PDB REMARK proxy_free_energy, soft_beta_G, provenance domain/measure/
   reference/evidence, and fail-closed availability.
4. Rank new outputs by emitted soft_beta_G; preserve legacy ordering only when
   that field is absent.
5. Keep unit labels conditional on provenance and remove physical Delta G claims
   from proxy paths.
6. Add C++-golden-corpus parity, hostile deserialization, PDB, JSON, and
   opposite-order election fixtures.
7. Preserve all existing numeric results except the explicitly required loader
   election parity.

Run focused and tracked Python tests only when experiment-safe. Never collect or
edit the excluded local fixture. Write your receipt to
handoff_swe/fable_chunk0_receipts_20260808/s2_python.json, include the consumed
Session 1 path/hash and fixture hash, then report your receipt's external
SHA-256. Do not commit.
```

### Session 3 prompt: TypeScript mirror

```text
You are Fable Session 3, the TypeScript SDK/shared/viewer mirror owner for
FlexAIDdS Chunk 0.
Task ID: flexaidds_c0_fable_s3_typescript_20260808.

Start only after the coordinator provides the exact path and SHA-256 for
Session 1's C++ truth receipt. Verify that hash and the golden-fixture hash it
records. Read AGENTS.md, METHODOLOGY.md, the complete handoff, Session 1's
receipt, and the fixture. C++ statmech is authoritative. Do not invent a
TypeScript-specific scientific contract.

You are not alone in the codebase. Your exclusive write ownership is
typescript/**. Do not edit C++, Python, Swift, docs/site/CI, handoff files, or
experiment outputs. Do not stage or commit.

Complete the TypeScript P0/P1 mirror:
1. Require literal available === true; missing, false, string, numeric, and
   malformed values fail closed.
2. Normalize SDK bindingModes, shared modes, and Python binding_modes plus
   unit-suffixed snake thermodynamics and nested snake provenance at one runtime
   boundary.
3. Consume the C++ golden corpus and add producer-to-normalizer-to-viewer tests.
4. Suppress forged proxy affinity/potency explanation and design suggestion.
5. Return public proxy selectivity driver as inconclusive.
6. Remove unsupported .rrd advertisement unless a real parser exists.
7. Restore unrelated numeric loader changes from the Chunk 0 diff; preserve
   legacy temperature and cluster/count behavior.
8. Keep public WASM parity explicitly unsupported unless provenance is actually
   wired and tested.

Run npm ci/typecheck/tests only when experiment-safe. Write your receipt to
handoff_swe/fable_chunk0_receipts_20260808/s3_typescript.json, include the
consumed Session 1 path/hash and fixture hash, then report your receipt's
external SHA-256. Do not commit.
```

### Session 4 prompt: Swift mirror

```text
You are Fable Session 4, the Swift bridge and Intelligence mirror owner for
FlexAIDdS Chunk 0.
Task ID: flexaidds_c0_fable_s4_swift_20260808.

Start only after the coordinator provides the exact path and SHA-256 for
Session 1's C++ truth receipt. Verify that hash and the golden-fixture hash it
records. Read AGENTS.md, METHODOLOGY.md, the complete handoff, Session 1's
receipt, and the fixture. C++ statmech is authoritative. Swift must transport
and enforce that truth without inventing provenance.

You are not alone in the codebase. Your exclusive write ownership is swift/**.
Do not edit C++, Python, TypeScript, docs/site/CI, handoff files, or experiment
outputs. Do not stage or commit.

Complete the Swift P0/P1 mirror:
1. Match C++ vocabulary, strict receipts, rejected digests, predicates, and
   correction validity using the golden corpus.
2. Preserve real C++ provenance across the Objective-C++ bridge; do not invent
   proxy metadata for arbitrary/calibrated inputs.
3. Rebuild the latest Intelligence surfaces and hostile proxy fallbacks.
4. Make Selectivity's retained numeric deltaG compatibility field explicitly
   unavailable/diagnostic so it cannot authorize a physical claim.
5. Link SwiftPM tests to the real C++ implementation. Do not use stubs,
   fabricated symbols, or -undefined dynamic_lookup.
6. Keep every affinity, optimization, lead-design, and unit-bearing claim gated
   on binding_physical provenance.

Run target builds/tests only when experiment-safe. If real-core linking remains
blocked, report the exact linker evidence instead of claiming green. Write your
receipt to handoff_swe/fable_chunk0_receipts_20260808/s4_swift.json, include the
consumed Session 1 path/hash and fixture hash, then report your receipt's
external SHA-256. Do not commit.
```

### Session 5 prompt: docs, validator, CI, and integration

```text
You are Fable Session 5, the final docs/validator/CI and integration owner for
FlexAIDdS Chunk 0.
Task ID: flexaidds_c0_fable_s5_integration_20260808.

Start only after the coordinator provides exact paths and SHA-256 values for
Sessions 1-4. Verify all four hashes and the golden-fixture hash before work.
Read AGENTS.md, METHODOLOGY.md, the full handoff, all four receipts, and the C++
golden fixture. C++ statmech remains the sole scientific truth. Package-local
green tests are not producer-to-consumer evidence.

You are not alone in the codebase. Your exclusive write ownership:
- README.md
- docs/**
- site/entropy-help/**
- scripts/validate_thermo_claims.py
- tests/test_thermo_claim_firewall.py
- .github/workflows/ci.yml

You may inspect every other lane read-only. Do not repair another lane yourself;
return a concrete blocker to its owner. Do not edit user-owned CleftDetector,
DatasetRunner/gaboom hunks, experiment outputs, or the local handoff. Do not
stage or commit.

Complete final integration:
1. Remove or receipt README/public physical Delta G, ranking, affinity, 78/85,
   reproducibility, feature-flag, and nonexistent API claims.
2. Make validator claim evidence local to each claim and reject every
   signature-like placeholder key.
3. Verify Session 1 gated/relabelled ReferenceEntropy, ParallelCampaign,
   benchmark_datasets, and DatasetRunner pKd conversion, and verify Session 3
   handled WASM advertising. Return any defect to the owning session; never
   invent calibration or edit another lane.
4. Verify docs/schema/CI use the C++ vocabulary and golden corpus rather than a
   separate scientific definition.
5. Review cross-language wire fixtures and all receipts for numeric
   noninterference and fail-closed behavior.
6. Run combined implementation gates only if they cannot interfere with Claude
   Science. OPS retains METHODOLOGY.md sections 1-4 and 6 merge/science gates.
7. Produce a final ready/not-ready receipt and an exact task-only staging
   manifest. Leave actual staging/commit/push to Claude Code/Dispatch.

Write the final receipt to
handoff_swe/fable_chunk0_receipts_20260808/s5_integration.json, record all four
consumed path/hash pairs plus the golden-fixture hash, and report your receipt's
external SHA-256 with every deferred gate and unresolved blocker. Do not commit.
```
