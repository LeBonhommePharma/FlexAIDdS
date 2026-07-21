# Orchestrator summary — CMA-ES validation dispatch

**UTC:** 2026-07-21T19:40:00Z  
**Protocol:** `VALIDATION.md` + `KICKOFF.md`  
**Workers:** `flexaidds-build-ab` (A–F), `flexaidds-harness` (G local half) — ran in parallel  

## Provenance of this run

| Field | Value |
|-------|--------|
| Branch at worker start | `perf/pb-clash-grid-hoist` |
| HEAD used for build + GA dock | `9723f9de41da1d185725384c6b774eca87fcbc4e` |
| HEAD subject | `wire CMA-ES search backend (opt-in FLEXAIDDS_SEARCH=cmaes)` |
| HEAD actual tree content | WebGPU + ThermoWhiteboard only — **no CMA-ES sources** |
| Also checked | `cmaes-ab` @ `9a5e8c8d3` — also lacks `LIB/cmaes_search.*` |
| Toolchain | `g++-16` Homebrew GCC 16.1.0; CMake 4.4.0 |
| Host | `LPmore.local` Darwin arm64 |

## CLOSED / OPEN (A–G)

| # | Status | One-line evidence |
|---|--------|-------------------|
| A | **OPEN** | No CMA-ES TUs; base cmake configure exit 0 (`validation_evidence/build_ab/A_cmake_configure.exit`) |
| B | **OPEN** | Base FAST binary linked sha256 `967d1946…` but adapter not in link (`B_binary.sha256`) |
| C | **OPEN** | Zero `FLEXAIDDS_SEARCH`/cmaes symbols (`C_doctor.txt`) |
| D | **OPEN** | GA short live 1G9V: elected CF=`-5.08119`, RMSD≈`13.0879` Å; CMA-ES arm N/A (`D_ga_result.txt`, `D_cmaes_result.txt`) |
| E | **OPEN** | Eval-matched 2e6 impossible without CMA-ES arm (`E_ab_summary.txt`) |
| F | **OPEN** | No CMA-ES entropy trace; GA-only partial only (`F_entropy_trace.txt`) |
| G | **OPEN** (local half) | No apptainer, no `.sif` recipe, no `collapse_fingerprint.py` (`harness/ITEM_G.md`) |

**CLOSED: 0 / OPEN: 7**

## Root blockers (real, not speculative)

1. **CMA-ES adapter not in tree:** `LIB/cmaes_search.cpp`, `LIB/cmaes_search.h`, `apply_integration.sh`, `CMAES_INTEGRATION.md` all **absent**. `git grep FLEXAIDDS_SEARCH` → 0 hits.
2. **Mislabelled commit:** HEAD message claims CMA-ES wiring; diff is WebGPU shaders + `ThermoWhiteboard.h` only (`baseline_inventory.txt`).
3. **Harness tooling:** `apptainer`/`singularity`/`sbatch` missing on macOS host (`G_tooling.txt`).
4. **No fingerprint tooling:** `analysis/collapse_fingerprint.py` absent (orchestrator re-confirmed).

## Orchestrator: collapse_fingerprint on build-ab CMA-ES trace

```
script=analysis/collapse_fingerprint.py
script_present=no
cmaes_trace_usable=no
```

Log: `validation_evidence/orchestrator/collapse_fingerprint_attempt.log`  
No CMA-ES entropy series was produced (only `D_cmaes_result.txt` blocker file). **Not fabricated.**

## Real numbers that *did* land on disk (GA-only, short budget)

From `validation_evidence/build_ab/D_ga_result.txt` (100×50, complex **1G9V**):

| Metric | Value |
|--------|--------|
| elected CF | `-5.08119` |
| best CF among ranks | `-16.70442` |
| elected ordered heavy RMSD vs crystal SDF | `13.0879` Å (n=25/25) |
| binary sha256 | `967d194698b048513edabf4038e105ff50b5ccd9eea83ceeb975fe22aaef2f6a` |
| GA entropy H_final (log) | `2.300787` |

These do **not** close D/E (both arms + eval-matched 2e6 required).

## Intended Narval submit (NOT EXECUTED)

See `validation_evidence/harness/G_narval_submit_command.txt` — host is not a login node.

## What must land before any A–G can CLOSE

1. Real CMA-ES integration sources + CMake wiring + `FLEXAIDDS_SEARCH=cmaes` gate (the P1–P7 package).
2. `analysis/collapse_fingerprint.py`.
3. Apptainer recipe + Linux host for G local half; Narval login node for G remote half.
NOTE: validation_evidence/build_ab/build_fast/FlexAIDdS is gitignored by pattern build* (see .gitignore:75). Provenance retained via B_binary.sha256 + local path on host.
