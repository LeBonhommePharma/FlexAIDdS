# DatasetRunner.cpp Split Plan

**Status:** P0 + P1 landed (stats + provenance leaves). Full split is multi-PR; do not attempt monolithically.
**Source of truth for workflow:** `AGENTS.md`.
**Audit trigger:** `LIB/DatasetRunner.cpp` ≈ 7816 lines after P1 (was ~8010 at audit; P0 stats + P1 provenance extracted).

## Goals

1. Keep behavioural parity for ranking, success gates (`success_rmsd`, `success_pb`, `claim_ready`), and provenance fields.
2. Extract **leaf modules first** (pure logic, no process control, no GA wiring).
3. One PR chunk at a time with an explicit **test gate** after each chunk.
4. Never change pose ranking, clustering, or publication gates unless the user asks (see AGENTS.md scientific guardrails).

## Region map (pre-P0 line numbers, 8010-line snapshot)

Line numbers refer to `LIB/DatasetRunner.cpp` at the audit snapshot before subsequent
extractions. After each extract, ranges for remaining regions shift — re-map from section banners.

| Region | Approx. lines | Role | Risk |
|--------|---------------|------|------|
| Banner / includes | 1–70 | Headers, OS deps | — |
| Seeds + `SubprocessGuard` | 71–262 | Deterministic GA seeds; fork/wait/kill lifecycle | Medium (signals) |
| RMSD helpers (Hungarian / pose–pose) | 263–618 | Munkres, serial + Hungarian ligand RMSD, pose loaders | Medium (success path) |
| HVIB + frequency-gated pose selection | 619–1155 | Eigen/Hvib columns; Fix B / pooled Fix B selector | **High** (ranking) |
| **Statistical free functions** | **1156–1272** | Pearson / Spearman / Kendall / serial RMSD | **Low (P0)** |
| Site prep helpers | 1273–1711 | Centroids, chain prune, ligand-centered site, blinding | Medium |
| Residue sets | 1712–1875 | Water/ions, cofactor blacklist, glycans | Low–medium |
| Ctor / path / exec / download | 1876–2181 | Cache dirs, curl, PDB/CIF cache validation | Low–medium |
| Structure parse + ligand extract | 2182–3710 | mmCIF/PDB HETATM, bonds, SDF extract, receptor strip | Medium |
| `prepare_pdb_entry` | 3710–3807 | Single-complex prep | Medium |
| Dataset fetchers + `prepare*` API | 3808–4926 | Astex, HAP2, CASF, DUD-E, PoseBusters set, SAMPL, PDBbind, DOI | Medium (I/O) |
| **`run()` — binary locate + provenance** | **4931–5104** | FlexAID binary discovery; `provenance.json` | Low (provenance leaf) |
| `run()` — infra | 5105–5280 | Signals, TargetServer sketch, async I/O, schedule | High |
| `run()` — **execute** (per-entry dock) | 5275–6235 | Config JSON, restarts, blind/seed, fork_exec, multi-cleft | **High** |
| `run()` — parse results / thermo tags | 6236–6570 | stdout/stderr parse, CF, Shannon, native CF | Medium |
| `run()` — RMSD + pose election ledger | 6572–7234 | Crystal RMSD, consensus, BCR-gate, pose SHA ledger | **High** |
| `run()` — **validate / PB / tENCoM** | **7236–7493** | Upstream `bust`, NativePoseQC diagnostic, claim_ready | **High** (gates) |
| `run()` — session I/O + aggregates | 7494–7774 | Per-complex CSV, TargetServer cross-ligand, correlations | Medium |
| **`write_report`** | **7776–8008** | Markdown + aggregate CSV | Low–medium |

### Logical buckets (for future file layout)

```
LIB/
  DatasetRunner.h / .cpp          # facade: prepare(), run(), class state
  DatasetRunnerStats.h / .cpp     # P0 — pure metrics (done)
  DatasetRunnerProvenance.h / .cpp # P1 — provenance.json + hash helpers (done)
  DatasetRunnerRmsd.*             # Hungarian / pose RMSD (not ranking selector)
  DatasetRunnerPrep.*             # download, extract_ligand, residue sets
  DatasetRunnerFetch.*            # fetch_astex, fetch_casf, code lists
  DatasetRunnerReport.*           # write_report Markdown/CSV
  # Later only, behind heavy tests:
  # DatasetRunnerExecute.*        # dock loop / config emission
  # DatasetRunnerValidate.*       # PB + tENCoM claim path
```

## PR chunks and test gates

### P0 — Pure stats leaf ✅ DONE

**Extract:** `compute_pearson_r`, `compute_spearman_rho`, `compute_kendall_tau`, `compute_rmsd` (+ internal `compute_ranks`).

**Files:**
- `LIB/DatasetRunnerStats.h` / `.cpp`
- Wire into `benchmark_datasets`, `test_dataset_runner`, `test_cofactor_blacklist`
- `DatasetRunner.h` includes `DatasetRunnerStats.h` (stable API for existing includes)

**Why first:** Zero I/O, zero ranking, zero process control. Covered by `StatisticalMetrics` / `RMSDComputation` in `tests/test_dataset_runner.cpp`.

**Test gate:**
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(sysctl -n hw.ncpu) --target test_dataset_runner
./build/test_dataset_runner --gtest_filter='StatisticalMetrics.*:RMSDComputation.*'
```

### P1 — Provenance JSON writer (leaf) ✅ DONE

**Extract:** `run()` block (`cmd_token`, hashes, matrix path, `provenance.json`).

**Files:**
- `LIB/DatasetRunnerProvenance.h` / `.cpp`
- Wired into `benchmark_datasets`, `test_dataset_runner`, `test_cofactor_blacklist`
- `DatasetRunner.h` includes `DatasetRunnerProvenance.h`
- Call site in `DatasetRunner::run()` delegates to `write_dataset_run_provenance(...)`

**Why next:** Self-contained I/O + hash helpers; no ranking, no claim_ready, no GA.

**Test gate:**
```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(sysctl -n hw.ncpu) --target test_dataset_runner
./build/test_dataset_runner --gtest_filter='ProvenanceJson.*:StatisticalMetrics.*:RMSDComputation.*'
ctest --test-dir build -R DatasetRunnerTests --output-on-failure
```

### P2 — Residue sets + path utilities

**Extract:** `excluded_residues()`, `cofactor_blacklist()`, `glycan_residues()`, maybe path helpers.

**Test gate:** `test_cofactor_blacklist` + code-list tests.

### P3 — Report generation

**Extract:** `write_report` (~7776–8008).

**Test gate:** Synthetic `BenchmarkReport` → CSV headers/summary; no docking.

### P4 — Hungarian / serial pose RMSD helpers

**Extract:** munkres / hungarian / pose-ligand RMSD loaders (not Fix B ranking selector).

**Test gate:** Existing RMSD tests + symmetry fixture if missing.

### P5 — Dataset code lists + fetch dispatch

**Extract:** hardcoded code lists and thin `fetch_*` / `prepare*` (not `run()`).

**Test gate:** `DatasetRunnerCodes.*`.

### P6 — Ligand extraction / structure parse (chunked)

**Extract carefully:** CIF parsers, `extract_ligand`, receptor strip.

**Test gate:** blacklist + ligand cache tests; more fixtures before large moves.

### P7+ — Execute / validate (defer)

Do not extract until P0–P6 are green. Includes config/restart loop, Fix B/BCR ranking, PoseBusters/tENCoM claim path, SubprocessGuard/signals. Each needs its own design note and full DatasetRunnerTests + `tests/test_posebust.cpp`.

## Non-goals

- getenv / ProtocolConfig (other agent)
- CF naming (other agent)
- CI formatting, Astex dedup
- Changing success semantics or default ranking

## Invariants

1. `success_pb := success_rmsd && pb_pass`
2. `claim_ready` requires official PB + tENCoM/Eigen on the same elected pose SHA256
3. Crystal coords evaluation-only after blinding; never re-seed ranking from RMSD
4. Apache-2.0 only
5. GitHub identity `LeBonhommePharma`

## Verification

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(sysctl -n hw.ncpu) --target test_dataset_runner
./build/test_dataset_runner --gtest_filter='StatisticalMetrics.*:RMSDComputation.*:ProvenanceJson.*'
ctest --test-dir build -R DatasetRunnerTests --output-on-failure
```
