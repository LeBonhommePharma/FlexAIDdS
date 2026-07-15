# Audit: `43a5dbfa3` — Merge PR #265 audit-datasetrunner-split-p0

| Field | Value |
|-------|--------|
| **Commit** | `43a5dbfa3a33b2edfc7fd6dd44888ed35a49e64c` |
| **Short** | `43a5dbfa3` |
| **Subject** | Merge pull request #265 from LeBonhommePharma/fix/audit-datasetrunner-split-p0 |
| **PR** | [#265](https://github.com/LeBonhommePharma/FlexAIDdS/pull/265) — *Refactor: DatasetRunner P0 stats leaf + split plan* |
| **Parents** | `599a66216` (main) · `e6032dc39` (branch tip, sole feature commit) |
| **Author date** | 2026-07-15 00:34:09 -0400 |
| **Merged at** | 2026-07-15T04:34:10Z |
| **Audit focus** | Split correctness · behavioral parity · CMake wiring · ranking safety · test gate · CI |
| **Scope** | Read-only audit (no source edits) |
| **Auditor** | Grok Build (26h-swarm #50) |
| **Verdict** | **PASS** — clean mechanical leaf extract; ranking/gates untouched; residual issues are pre-existing or out-of-band |

---

## 1. Executive summary

This merge lands **P0** of the documented DatasetRunner monolith split: pure correlation and serial-RMSD free functions move from `LIB/DatasetRunner.cpp` into a freestanding leaf `LIB/DatasetRunnerStats.{h,cpp}`, plus a multi-PR split plan under `docs/implementation/datasetrunner-split-plan.md`.

The functional delta is **zero by design**. Normalized function bodies for `compute_pearson_r`, `compute_ranks`, and `compute_rmsd` are **byte-identical** to the pre-extract versions; `compute_spearman_rho` / `compute_kendall_tau` differ only by removed comments. Public symbols remain in `namespace dataset` and stay reachable through `#include "DatasetRunner.h"` (which now includes `DatasetRunnerStats.h`).

**No ranking, Fix B / BCR election, PoseBusters, claim_ready, GA budget, or success-gate logic was touched.** Aggregate affinity correlations and `benchmark_datasets` thermo tables continue to call the same symbols.

Local verification at audit time: `./build/test_dataset_runner --gtest_filter='StatisticalMetrics.*:RMSDComputation.*'` → **24/24 PASS**.

---

## 2. Files touched (merge tree = PR tip)

| Path | Δ | Role |
|------|---|------|
| `LIB/DatasetRunnerStats.h` | +31 (new) | Public leaf API (Pearson / Spearman / Kendall / serial RMSD) |
| `LIB/DatasetRunnerStats.cpp` | +132 (new) | Implementations + internal `compute_ranks` |
| `LIB/DatasetRunner.cpp` | −117 / +1 | Implementations removed; pointer comment left |
| `LIB/DatasetRunner.h` | −16 / +3 | Declarations removed; `#include "DatasetRunnerStats.h"` |
| `CMakeLists.txt` | +3 | Source on `benchmark_datasets`, `test_dataset_runner`, `test_cofactor_blacklist` |
| `build_sources.ignore` | +1 | `LIB/DatasetRunnerStats.h` (header ignore pattern) |
| `docs/implementation/datasetrunner-split-plan.md` | +139 (new) | Region map + P0–P7+ chunk plan + test gates |
| `tests/test_dataset_runner.cpp` | +17 | Constant-series + mismatched-length leaf contracts |

**8 files, +327 / −133.** Single feature commit on the branch: `e6032dc39`. Merge commit adds no extra content beyond that tree (`git diff 43a5dbfa3^1 43a5dbfa3` matches the PR).

---

## 3. Split correctness

### 3.1 What was extracted (intended leaf)

| Symbol | Visibility | Role |
|--------|------------|------|
| `compute_pearson_r` | public | Linear correlation |
| `compute_spearman_rho` | public | Pearson-of-average-ranks |
| `compute_kendall_tau` | public | Kendall τ-b (pair enumeration) |
| `compute_rmsd` | public | Flat 3N float serial RMSD (no Kabsch) |
| `compute_ranks` | `static` TU-local | Average ranks for ties |

**Leaf properties (good):** no I/O, no process control, no `DatasetRunner` state, no ranking selectors, no PB/tENCoM, Apache-2.0 file headers, no machine-absolute paths, no secrets.

### 3.2 Behavioral parity evidence

Compared `e6032dc39^:LIB/DatasetRunner.cpp` function bodies vs `e6032dc39:LIB/DatasetRunnerStats.cpp` (whitespace-normalized):

| Function | Identical? | Notes |
|----------|------------|--------|
| `compute_pearson_r` | **Yes** | Same means, cov, denom &lt; 1e-15 → 0.0 |
| `compute_ranks` | **Yes** | Average 1-based ranks for exact ties |
| `compute_rmsd` | **Yes** | Sentinel −1.0 on empty / size mismatch / non-multiple-of-3 |
| `compute_spearman_rho` | Logic yes | Only dropped comment “Spearman ρ = Pearson r of ranks” |
| `compute_kendall_tau` | Logic yes | Only dropped “Kendall tau-b: handles ties” comment |

`nm` on the local object confirms exported symbols under `dataset::` for all four public functions — no accidental static/anonymous-namespace demotion.

### 3.3 API stability

- `DatasetRunner.h` still brings the stats API into any TU that only includes the facade.
- Call sites need no code changes:
  - `DatasetRunner::run()` aggregate block: unqualified `compute_*` (ADL / same namespace).
  - `LIB/benchmark_datasets.cpp`: `dataset::compute_pearson_r` / spearman / kendall on pred vs exp ΔG/ΔH/TΔS vectors.
- Tests additionally `#include "DatasetRunnerStats.h"` to assert standalone compilability of the leaf.

### 3.4 CMake / build hygiene

All three TUs that compile `DatasetRunner.cpp` also compile `DatasetRunnerStats.cpp` at this merge:

1. `benchmark_datasets`
2. `test_dataset_runner`
3. `test_cofactor_blacklist`

No other CMake target listed `DatasetRunner.cpp` without Stats at `43a5dbfa3`. Header listed in `build_sources.ignore` (same pattern as `DatasetRunner.h`).

**Note:** GitHub Actions on PR head `e6032dc39` failed **Configure** with STRICT source guard orphans — but the three listed orphans were **`LIB/FleetRunner.cpp`**, **`tools/metal_microbench_enhanced.cpp`**, and **`LIB/ProtocolConfig.h`**, **not** `DatasetRunnerStats.*`. Local PR claim “CMake source validator clean” is **overstated for STRICT CI**; the Stats leaf itself is correctly referenced and is not an orphan.

---

## 4. Findings

### F1. Pure extract — ranking / success gates unchanged (INFO → PASS criterion)

- Evidence: Diff does not touch `select_pose_freq_gated_pooled`, `compute_pose_ligand_rmsd`, Hungarian path, `success_rmsd` / `success_pb` / `claim_ready`, GA pop/gen, CF scoring, or Fix B election. Stats calls remain post-dock aggregates only (`DatasetRunner.cpp` ~7799–7802 and `benchmark_datasets.cpp` affinity tables).
- Why it matters: AGENTS.md forbids silent ranking changes. This PR is a textbook leaf-first split.
- Fix recommendation: None for this merge. Keep P7+ execute/validate extractions behind full DatasetRunnerTests + posebust suite as the plan states.

### F2. `compute_rmsd` remains production-dead after extract (LOW)

- Evidence: Production crystal/pose RMSD uses `compute_pose_ligand_rmsd` / `hungarian_rmsd` / `pose_pose_rmsd` inside `DatasetRunner.cpp`. `rg 'compute_rmsd\('` outside `DatasetRunnerStats.*` hits **only** `tests/test_dataset_runner.cpp`. Pre-existing methodology audit (`flexaidds_methodology_audit.md`) already flagged a “dead duplicate” flat `compute_rmsd`.
- Why it matters for science/repro: Leaf extract **preserves** the dual-implementation risk — future edits to the serial helper will still not affect benchmark success RMSD, and editors may assume the leaf is the success metric. Not a regression vs parent, but the split documents the helper as part of DatasetRunner’s public surface.
- Fix recommendation (follow-up, not merge-blocker): Either (a) wire serial `compute_rmsd` only as a unit-test utility in `tests/`, or (b) mark header docs “not used for success_rmsd; see hungarian path,” or (c) delete when P4 extracts real RMSD helpers and consolidates tests.

### F3. Correlation sentinel `0.0` collides with true zero correlation (LOW, pre-existing)

- Evidence: `compute_pearson_r` / spearman / kendall return `0.0` for n&lt;2, size mismatch, constant series, or degenerate denom. True r/ρ/τ of 0 is indistinguishable. Header now documents constant-series → 0.0 (documentation improvement). Call sites only compute when `pred_affinities.size() >= 3` (DatasetRunner) or `exp_*.size() >= 3` (benchmark_datasets), which reduces n&lt;2 issues but not constant-series or mismatch.
- Why it matters: Affinity-panel reports can print “Pearson 0.000” when the series is constant or ill-conditioned, looking like “uncorrelated” rather than “undefined.” Same behavior as before the extract.
- Fix recommendation (follow-up): Return `std::optional<double>` / NaN with `isfinite` gates (report already filters finite correlations in Markdown) or a parallel `ok` flag. Out of scope for a parity-preserving leaf move.

### F4. NaN / Inf inputs undefined (LOW, pre-existing)

- Evidence: No `isfinite` filtering in the leaf. NaN means poison ranks and correlations. Ties in ranks use exact `==`; Kendall ties use `abs(dx) < 1e-12` — inconsistent tie policy between Spearman and Kendall (pre-existing).
- Why it matters: Contaminated predicted_dG vectors could yield silent NaN report fields; write_report already skips non-finite Pearson/Spearman in one Markdown path but still writes raw CSV fields.
- Fix recommendation (follow-up): Reject non-finite pairs at call site or inside the leaf; align Kendall/Spearman tie epsilons if panels use near-ties.

### F5. O(n²) Kendall and affinity semantics are report-only (INFO)

- Evidence: Kendall double-loop is fine for Astex-scale n (~85) and small ITC panels; not a ranking hot path. `run()` correlates **pKd** (`-predicted_dG/1.3636` vs experimental affinity); `benchmark_datasets` prints labels **“ΔG (kcal/mol)”** correlating `predicted_dG` vs `-affinity*1.3636`. Unit conventions differ between the two consumers — **pre-existing**, not introduced here.
- Why it matters for science language: AGENTS.md — do not over-claim true thermodynamic ΔG. These metrics score **whatever `predicted_dG` is** (often CF-proxy-derived), not a validated partition-function ΔG. Extract does not change that.
- Fix recommendation: When P3 extracts `write_report`, rename columns to “pred vs exp affinity correlation” / “CF-proxy score correlation” where appropriate. Not a P0 defect.

### F6. Predicted affinity filter excludes true zero `predicted_dG` (LOW, pre-existing, call-site)

- Evidence: `if (entries[i].has_affinity() && r.predicted_dG != 0.0f)` before pushing correlation pairs; similar `!= 0.0f` in `benchmark_datasets.cpp`. Zero is used as “missing” sentinel.
- Why it matters: A legitimate zero prediction would be dropped from r/ρ/τ. Unchanged by this PR.
- Fix recommendation: Follow-up — use an explicit `has_predicted_dG` flag rather than zero-as-missing.

### F7. Split plan quality (INFO / positive) with minor doc drift

- Evidence: Plan correctly prioritizes pure leaves, defers execute/validate (P7+), states invariants (`success_pb := success_rmsd && pb_pass`, claim_ready pose SHA, no RMSD reseeding ranking, Apache-2.0). Region map line numbers are an explicit pre-P0 snapshot (8010 lines → 7896 after P0).
- PR body says “P0–P10”; plan body defines **P0–P7+** — harmless marketing drift.
- Current tree’s plan file has been updated by later P1 work (`✅ DONE` markers); **as merged in #265** it was the initial P0 plan only.
- Fix recommendation: Keep one plan doc; avoid renumbering mid-flight.

### F8. Tests for the leaf (PASS with small gaps)

- Evidence: Pre-existing suite covers perfect ±Pearson/Spearman/Kendall, known Kendall 0.2, ties, empty, single-element, RMSD identity / known / shift / empty / mismatch. P0 adds **ConstantSeriesReturnsZero** and **MismatchedLengthReturnsZero**. PR claimed 24/24 filtered + 63/63 full binary locally; auditor re-ran filtered suite → 24/24 PASS on current `build/test_dataset_runner`.
- Gaps: no NaN test; no large-n Kendall stress; no assertion that production RMSD path still does **not** call `compute_rmsd` (optional architectural test).
- Fix recommendation: Optional NaN tests when tightening sentinels (F3/F4).

### F9. CI red on PR is real but orthogonal (LOW for this change)

- Evidence: `gh pr checks 265` — multiple C++ core / macOS / pure-Python jobs **fail at Configure or early steps**; source guard orphans = FleetRunner / metal_microbench / ProtocolConfig (see §3.4). License scan, CodeQL analyze (partial), Metal shader smoke, Python bindings (linux), repo hygiene, skill DatasetRunner package job **passed**. CodeRabbit release notes incorrectly frame the extract as “New Features” / “Bug Fixes” for invalid-data handling — behavior was already present.
- Why it matters: Merge landed with red checks; AGENTS.md wants green before push for *this* change’s own faults. Strictly, the Stats extract did not create the orphans; main/branch pollution did.
- Fix recommendation: Treat STRICT orphans as a separate hygiene PR; do not blame P0 Stats for FleetRunner wiring.

### F10. Security / licensing (INFO)

- Evidence: No new network, shell, or path construction. Apache-2.0 headers on new files. No GPL. No secrets. Hygiene pattern for headers follows existing DatasetRunner ignore list.
- Fix recommendation: N/A.

---

## 5. Scientific / product impact matrix

| Axis | Impact of this merge |
|------|----------------------|
| Pose ranking / election | **None** |
| CF / contact-function scoring | **None** |
| `success_rmsd` / `success_pb` / `claim_ready` | **None** |
| GA budget (pop × gen / DoF scale) | **None** |
| Aggregate affinity correlations | **Parity** (same formulas, same call sites) |
| Serial `compute_rmsd` success metric | **Still unused** (same as parent) |
| Reproducibility of docks | **None** (no dock path change) |
| Maintainability | **Positive** — first real cut of ~8k-line monolith; plan gates future chunks |
| Binary size / link | + one small TU on three targets only |

**Ranking/scoring impact: NO**

**Reproducibility impact: NEUTRAL** (parity-preserving refactor; docs improve intent)

**Tests adequate: YES** for a pure leaf move (with known pre-existing metric-semantics gaps F2–F6)

---

## 6. Verification performed this audit

```text
# Parity
python3: normalized body compare e6032dc39^ DatasetRunner.cpp vs DatasetRunnerStats.cpp
  → pearson/ranks/rmsd identical; spearman/kendall comment-only delta

# Link symbols
nm …/DatasetRunnerStats.cpp.o | c++filt → dataset::compute_{pearson_r,spearman_rho,kendall_tau,rmsd}

# Call graph
rg compute_rmsd → tests only (production uses hungarian / pose-ligand path)
rg compute_pearson_r → DatasetRunner aggregates + benchmark_datasets tables + tests

# Tests (local binary present)
./build/test_dataset_runner --gtest_filter='StatisticalMetrics.*:RMSDComputation.*'
  → [  PASSED  ] 24 tests.

# CI (PR head e6032dc39)
Configure STRICT orphans: FleetRunner.cpp, metal_microbench_enhanced.cpp, ProtocolConfig.h
  (DatasetRunnerStats not orphaned)
```

---

## 7. Verdict: **PASS**

Safe, correctly scoped **P0** of the DatasetRunner split:

1. Behavioral parity of statistical helpers is verified.
2. CMake wires every `DatasetRunner.cpp` consumer.
3. Public API preserved via facade include.
4. No scientific ranking, claim, or budget surface changed.
5. Plan document correctly defers high-risk execute/validate extractions.

**Not merge-blockers (track separately):** dead serial `compute_rmsd` (F2), 0.0/NaN sentinel design (F3–F4), affinity zero-as-missing (F6), pre-existing STRICT orphans that reddened PR CI (F9), CodeRabbit “new feature” noise.

**Do not use this merge as evidence of improved affinity ranking science** — it only relocates report metrics. Continue to describe outputs as CF/contact-function scoring proxy and ensemble-derived estimates per AGENTS.md, not validated thermodynamic ΔG, unless the full StatMech + vib + solvent path is active and validated.

### Follow-ups (optional, other PRs)

1. P1+ per plan (provenance already landed later as PR #272 track).
2. Mark or quarantine dead `compute_rmsd` before P4 RMSD extraction.
3. Clear STRICT orphans (FleetRunner / metal microbench / ProtocolConfig ignore or wire).
4. When touching correlation call sites, replace zero-sentinels and align ΔG vs pKd labeling with scientific guardrails.
