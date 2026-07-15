# Audit: a4e13fc9b — Merge branch master into feature/bonhomme-fleet-dataset-runner

## Summary (2-4 sentences)
Commit `a4e13fc9b0abbe1186d82dd0c8e6c32d488952b5` merges **master** (`155ebbb7b`, tip after Dependabot PR #255) into **feature/bonhomme-fleet-dataset-runner** (`4beb3b36d`). The feature tip’s only unique commit was a **catastrophic CMakeLists stub** (78-byte placeholder replacing 3226 lines of the real build system). Conflict resolution correctly took **master’s** `CMakeLists.txt`; the resulting tree is **byte-identical to master** (`tree f6eb92ab…`). This merge is therefore a pure master-alignment / stub-discard with **no residual feature content**, no fleet DatasetRunner integration, and no new ranking logic of its own — and it made subsequent PR #258 (`964bec0a2`) a **tree-no-op** onto master.

## Severity: LOW

(Process/meta issues MEDIUM; merge resolution itself is correct and safe.)

## Merge topology

| Role | SHA (short) | Subject |
|------|-------------|---------|
| Merge | `a4e13fc9b` | Merge branch 'master' into feature/bonhomme-fleet-dataset-runner |
| Parent 1 (ours / feature) | `4beb3b36d` | CMake: entropy collapse — modular targets, zero monolith, PGO-first, Eigen3-centric, hardware acceleration paradise |
| Parent 2 (theirs / master) | `155ebbb7b` | Merge pull request #255 (actions/checkout-7.0.0) |
| Merge-base | `d8a422f93` | Chore: ignore .tools/ and .venv*/ for local tool installs |
| Feature-only commits | **1** | Only `4beb3b36d` |
| Master-only commits | **187** | Full science + packaging window (PoseBust, soft-β, eval-scale, …) |
| Both-modified paths | **1** | `CMakeLists.txt` only |
| Result tree | `f6eb92ab…` | **Identical** to parent2 (master) tree |

Author: Louis-Philippe Morency (LeBonhommePharma). Committer: GitHub merge UI. Date: 2026-07-14 20:37:20 -0400. Signed-off-by present. No conflict-marker residue in source (only decorative `# ====…` section banners).

## Findings

### F1. Feature tip was a build-system destruction stub — correctly discarded (HIGH if landed; FIXED by merge)
- Evidence: `4beb3b36d` rewrote `CMakeLists.txt` as a single line (78 bytes, no newline):
  ```text
  # New modular top-level CMakeLists.txt ... (full content from Grok's previous)
  ```
  Diffstat vs merge-base: `3227 +-------------------------------------------------------` (1 insertion, 3226 deletions). Blob `b27f97fc4`. No modular targets, no `add_executable`, no PGO, no Eigen wiring — only a comment claiming Grok would fill content later.
- Why it matters for science/repro: Configuring that tree cannot produce FlexAID / DatasetRunner / ctest. Any claim run or CI on the pre-merge feature tip would fail closed at configure. The commit **subject** oversells a completed “entropy collapse / hardware acceleration paradise.”
- Fix recommendation (already done in this merge): Keep master’s full `CMakeLists.txt` (result == master). Do **not** resurrect `4beb3b36d`’s blob. Any real modular CMake rewrite must land as a complete, building diff with `cmake --build` + `ctest` evidence (AGENTS.md).

### F2. Pure master tree takeover — merge is content-empty vs master (INFO / GOOD)
- Evidence: `git rev-parse 155ebbb7b^{tree} a4e13fc9b^{tree}` both yield `f6eb92abfb952b67d09b67073c345d7c3f980542`. `git diff 155ebbb7b a4e13fc9b` is empty. All inspected science paths (`LIB/DatasetRunner.cpp`, `gaboom.cpp`, `BindingMode.cpp`, `top.cpp`, `Vcontacts.cpp`, `cluster.cpp`, `FastOPTICS_cluster.cpp`) match master.
- Why it matters: This is the **correct** resolution when the only both-modified file is a stub vs a living 2948-line CMake. No hybrid merge bugs, no silent half-application of “entropy collapse.” Ranking/scoring code is exactly master’s production stack at 2026-07-14.
- Fix recommendation: None for resolution. For hygiene, prefer `git merge -s ours` only when intentionally discarding; here default merge already produced the desired tree because the feature side had no non-conflicting unique paths.

### F3. Branch name / PR title overclaim fleet work that is not present (MEDIUM process)
- Evidence: Branch `feature/bonhomme-fleet-dataset-runner`. Unique pre-merge content is CMake stub only. Merge result has no first-class FleetRunner CMake target; fleet hits in tree are pre-existing experimental `swift/Sources/FleetScheduler/*`, `typescript/.../FleetDashboard.tsx`, and `benchmarks/m3pro/dashboard/fleet_*` — not a DatasetRunner backend integration introduced here. Immediate downstream PR #258 (`964bec0a2`, Merge: `155ebbb7b` + `a4e13fc9b`) also has **identical tree** to master and still titles the PR after the discarded CMake subject.
- Why it matters for science/repro: Operators may believe PR #258 / this merge landed “Bonhomme Fleet DatasetRunner” capability. History claims do not match tree. Violates “inspect first, claim never.”
- Fix recommendation: Rename/retitle closed PR notes to “Align feature branch with master; discard CMake stub.” Land real fleet DatasetRunner work on a new branch with tests and a truthful subject. Do not cite this SHA as fleet delivery.

### F4. Master delta is large and ranking-relevant relative to pre-merge feature tip (INFO — inherited, not invented)
- Evidence: Combined diff vs parent1 is ~3821 files, +118214 / −34933 (dominated by `swift/` experimental tree + benchmarks + LIB). Science commits now on the branch include (non-exhaustive):
  - Soft-β / classic entropy ranking restore (`1db64e5cd`, ACF-vs-CF ablation `4c11e5a6c`)
  - GA CF stagnation vs SMFREE proxy (`6eb170a70`)
  - Eval budget gen→pop scale + `FLEXAIDDS_EVAL_SCALE_DIHEDRAL=-1` fixed budget (`fb15e3306`, `577f2f66a`)
  - Native PoseBust + `success_pb = success_rmsd && pb_pass` (`a68fcf7d3`, `d5cef3041`, …)
  - No-seed cognate stack / DirectLigandIC (`dea70ea88`), oracle-ceiling seed fix (`26cb99276`)
  - Skill rename `.grok/skills/flexaid-docking` → `flexaidds` + `resolve_build.py`
  - Dependabot Actions bumps at tip (#251–#255)
- Why it matters: Anyone who had checked out the feature tip *before* this merge and rebuilds after gets master’s full protocol stack — including ranking and PB gates — not a fleet-only delta. Relative to **master**, ranking is unchanged.
- Fix recommendation: Document branch rebuilds as “now == master @ 155ebbb7b.” Claim campaigns must pin binary SHA from this tree (or later) and use pop-scale DoF budget rules per AGENTS.md.

### F5. Absolute machine paths and iCloud paths enter via master docs/scripts (LOW–MEDIUM hygiene)
- Evidence: Merge tree greps (from master) include hardcoded `/Users/lp.more/Projects/FlexAIDdS` and iCloud Mobile Documents paths in `benchmarks/astex_entropy/README.md`, `benchmarks/BENCHMARK_STANDARD.md`, `benchmarks/astex_diverse/janitor.sh`, etc. Not introduced by feature-side edits (feature had none), but now present on the branch tip after merge.
- Why it matters: AGENTS.md forbids machine-specific absolute paths in shared agent skills/scripts; docs are softer but still break portability and encourage CloudDocs walks. Janitor log path is operator-hostile on other machines.
- Fix recommendation: Out of band (not this merge’s job): replace with `$FLEXAIDDS_ROOT` / `git rev-parse --show-toplevel` and local-first benchmark roots. Run `python3 scripts/check_repo_hygiene.py` on future skill edits.

### F6. No conflict-marker residue; CMake conflict was clean master-take (GOOD)
- Evidence: Only both-modified file was `CMakeLists.txt`. Merge result matches master exactly (2948 lines of real build). False-positive grep hits on `====` are section banners in Python/CMake comments, not `<<<<<<<` markers.
- Why it matters: Avoids classic “shipped conflict markers” silent compile failures.
- Fix recommendation: N/A.

### F7. Licensing surface inherited cleanly; PoseBust is in-tree Apache-oriented native code (INFO)
- Evidence: Master brings full `LIB/PoseBust/*` (Loaders, Checks*, Engine, BustCli) wired in `CMakeLists.txt` and `test_posebust`. `THIRD_PARTY_LICENSES.md` still documents NRGRank as scientific inspiration only (GPL isolation). No new GPL dependency introduced by the merge resolution itself.
- Why it matters: Campaign success contract uses NativePoseQC / PB extract gates — licensing and science both require the clean-room path.
- Fix recommendation: Keep using native PoseBust for claim gates; do not vendor GPL PoseBusters Python as a link dependency.

### F8. Security / CI surface (INFO)
- Evidence: Master tip is series of Dependabot Actions version bumps + existing `pypi-release.yml` (OIDC trusted publishing, pinned action SHAs in the workflow body sampled). No secrets committed by this merge. Feature stub introduced no network code.
- Fix recommendation: Continue Dependabot pin discipline; verify `gh api user` is `LeBonhommePharma` before any push (AGENTS.md).

### F9. Tests: inherit master suite; no merge-specific test (LOW)
- Evidence: Master delta adds/updates many tests (`test_posebust.cpp`, `test_classic_entropy_ranking.cpp`, `test_ensemble_pipeline.cpp`, `test_soft_wall.cpp`, `test_resolve_build.py`, …). This merge commit itself adds nothing unique and does not re-run CI as part of the commit object. Feature tip had **no** ability to build tests (CMake stub).
- Why it matters: Correct resolution restores testability; empty PR #258 still doesn’t prove CI green without an Actions run on the merge SHA.
- Fix recommendation: Treat green CI on `155ebbb7b` / `a4e13fc9b` (same tree) as the gate. Any future non-empty fleet work needs DatasetRunner + ctest coverage.

### F10. Immediate PR #258 is a null landing (MEDIUM process)
- Evidence: `964bec0a2` parents = `155ebbb7b` + `a4e13fc9b`; all three trees share `f6eb92ab…`. `git diff 155ebbb7b 964bec0a2` empty. PR subject still references the discarded CMake entropy-collapse work.
- Why it matters: GitHub history shows a “merged feature” that changed zero files — confuses bisect, release notes, and scientific provenance of “what fleet shipped when.”
- Fix recommendation: Annotate release notes: PR #258 empty after stub discard. Real fleet commits must be separate, non-empty PRs.

## Ranking/scoring impact: NO (vs master) / YES (vs pre-merge feature tip)

| Comparison | Impact |
|------------|--------|
| Merge result vs master `155ebbb7b` | **None** — identical tree; same CF proxy, soft-β mode ranking, PB success contract, eval-scale knobs. |
| Merge result vs feature tip `4beb3b36d` | **Restores entire master stack** — feature tip could not even link a scorer; after merge, production ranking + PB gates apply. |

This merge does **not** invent ranking changes; it refuses to ship a broken CMake and realigns the feature branch to master’s ranking code of record.

## Reproducibility impact: YES (positive net)

- Positive: Discards non-buildable CMake stub; restores full configure/build/test surface; aligns branch with master’s reproducible protocol (RUN_RECEIPT-era paths on master, PoseBust gates, eval-scale env contract).
- Negative: Process noise (branch/PR names, empty PR #258) can mis-document what was shipped; absolute paths in inherited benchmark docs harm multi-machine repro until cleaned.

## Tests adequate: YES (inherited) / N/A (merge-unique)

No merge-unique code. Adequacy equals master’s suite at `155ebbb7b`. Feature tip had zero testability.

## Verdict: SAFE

**Accept as-is.** Correct conflict resolution (master CMake wins), pure tree equality to master, no residual hybrid bugs, no ranking/scoring regression vs production master. Do not treat this commit or PR #258 as delivery of Bonhomme Fleet DatasetRunner or of a modular CMake rewrite — both claims are false against the tree. Any future fleet or CMake modularization must be a non-empty, verified change set with build+test evidence.
