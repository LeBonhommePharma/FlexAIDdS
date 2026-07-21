# VALIDATION — FlexAIDdS CMA-ES search backend

Base commit: `2a60f65` (github.com/lmorency/FlexAIDdS). Verified in a sandbox with
**g++ 13.3, no cmake** — so everything in PROVEN was checked WITHOUT linking the full
engine. Do not read "proven" as "works as a dock." It means exactly what each line says.

## FOR THE WORKER — acceptance protocol (read before touching anything)
You (a Claude Code subagent) close the OPEN items below. Rules, non-negotiable:
- An item flips to **CLOSED** only with an **on-disk artifact** that proves it. No artifact,
  no close.
- In the Status/Evidence columns, write the artifact path + a one-line proof
  (e.g. `build/FlexAIDdS sha256=… ; cmake configure exit 0`). Numbers must be read from a
  file, never narrated from memory.
- An item you could not reach stays **OPEN** with the **real blocker** in Evidence
  (e.g. `OPEN — no GCC>=14 on host`). A truthful OPEN is a valid result; a fabricated CLOSED
  is a failure.
- Do not edit the PROVEN table (P1–P7) — those are already established. You only fill A–G.
- Never mark D or E closed without a real RMSD on a real complex. That is the whole point.

## PROVEN (live, real numbers — do not modify)

| # | Claim | How verified |
|---|---|---|
| P1 | Adapter compiles against the **real** engine headers (C++23) | `g++ -std=c++23 -Wall -Wextra -c LIB/cmaes_search.cpp` → exit 0, **0 warn / 0 err** |
| P2 | Coupling to engine = **exactly 5 seam fns**, correct signatures | `nm -uC` footprint: `set_gene_lim`, `set_bins`, `eval_chromosome`, `get_cf_evalue`, `get_apparent_cf_evalue` — nothing else engine-side |
| P3 | Those 5 symbols are **defined** in the engine (so it will link) | gaboom.cpp:3651 / :3450 / :2887; ic2cf.cpp:576 / :568 |
| P4 | Integration surface is minimal + additive | `top.cpp` patched (opt-in ternary); CMake `target_sources(flexaid_core …)` after `add_subdirectory(LIB)` (156→166); new files only; **`ic2cf.cpp` / `gaboom.cpp` untouched** |
| P5 | Search + seam plumbing functional vs a **mock** objective | seed 12345: best CF `2.477e-08`, pose max rel err `1.02e-04`, snapshot `5120 = 64×80` exact |
| P6 | Entropy trace emits and collapses as designed | `H_search` +0.95→−63.5 nats, `F`→−2.4641; `H_energy`→ln(64)=4.1589 (single-basin nuance — see CMAES_INTEGRATION.md) |
| P7 | `apply_integration.sh` on-ramp works from a clean tree | tested: correct CMake insertion point, idempotent re-run, smoke gate PASS |

## OPEN / CLOSED — host validation dispatch (2026-07-21)

**Workers:** CMA-ES swarm chunks 1–6 + orchestrator merge.  
**Integration branch:** `feat/cmaes-search-backend`.  
**Orchestrator build dir:** `.swarm/cmaes/orchestrator/build_fast`.  
**Toolchain:** `g++-16` (Homebrew GCC 16.1.0), CMake 4.4.0, host Darwin arm64.  
**Orchestrator summary:** `.swarm/cmaes/orchestrator/SUMMARY.md` (also `validation_evidence/orchestrator/SUMMARY.md`).

| # | Acceptance item | Owner | Status | Evidence (path + one-line proof) |
|---|---|---|---|---|
| A | `cmake configure` accepts the wiring AND `cmake/ValidateSources.cmake` passes with the new TUs | orchestrator | **CLOSED** | exit 0 (`validation_evidence/build_ab/A_cmake_configure.exit`); `cmaes_search.cpp` in `LIB/CMakeLists.txt` FLEXAID_CORE_SOURCES; log `validation_evidence/build_ab/A_cmake_configure.log`; inventory `A_inventory.txt` |
| B | Full engine + adapter **link** into one `FlexAIDdS` binary (`BUILD_FLEXAIDDS_FAST`) | orchestrator | **CLOSED** | link OK; binary `.swarm/cmaes/orchestrator/build_fast/FlexAIDdS` sha256 `404b3ccddc22c12bf3cfaced9b0eaf996faa16d5caa501d7ff691282d66f9eb1` (`B_binary.sha256`); build log `B_build.log` shows `cmaes_search.cpp.o` + `Built target FlexAIDdS` |
| C | Adapter runs against the **real `ic2cf`** (not the mock) with no exception; snapshot fills for a real ligand's DOF count | orchestrator | **OPEN** | Doctor proves symbols linked (`C_doctor.txt`: `cmaes_run_dock`, `FLEXAIDDS_SEARCH`, `[SEARCH] backend=cmaes`); mock path 4/4 PASS (`ctest_cmaes.log` via clang++ standalone). **No live dock yet** — real-ic2cf snapshot fill not exercised |
| D | Live dock on one Astex complex: elected-pose **RMSD** + best **CF**, both GA and CMA-ES arms | orchestrator | **OPEN** | Prior GA short 1G9V only (`D_ga_result.txt`); CMA-ES arm not run (`D_cmaes_result.txt`); no fabricated RMSD |
| E | GA-vs-CMA-ES **A/B**, eval-matched at 2e6 | orchestrator | **OPEN** | Manifest ready (`scripts/cmaes_ab_manifest.json` budget 1000×2000=2e6 validated); both arms not executed — ΔRMSD/ΔCF = N/A |
| F | Real entropy trace on the **rugged** real surface | orchestrator | **OPEN** (mock CLOSED) | Mock collapse fingerprint PASS → `validation_evidence/build_ab/F_fingerprint_mock.json` sha256 `53b0d3ed…db72b7`; real CMA-ES dock trace still absent |
| G | Locked-arch `.sif` builds; in-container dock; collapse fingerprint INVARIANT | orchestrator | **OPEN** (local half advanced) | Recipe + scripts landed (`containers/flexaidds_locked_x86_64.def`, `scripts/narval_cmaes_array.sh`, `G_harness_local.txt`); apptainer **MISSING** on Darwin; no `.sif`; mock fingerprint INVARIANT as F |

## Bottom line
Host integration **landed and linked**: `LIB/cmaes_search.{cpp,h}`, `FLEXAIDDS_SEARCH=cmaes` branch in `top.cpp`, `apply_integration.sh`, fingerprint tool, harness recipes, mock unit tests (4/4). Closed **A, B** with on-disk proof. **C** has binary/doctor proof but needs a live dock for snapshot fill. **D/E/F(real)/G(.sif)** remain OPEN — no fabricated RMSD/CF. Mock tests pass; `ctest` under g++-16 + Homebrew GTest is ABI-blocked (use clang++ standalone or rebuild GTest with g++-16).

### Unblock checklist (remaining)
1. ~~Land the CMA-ES package~~ (done on `feat/cmaes-search-backend`).
2. ~~Land `analysis/collapse_fingerprint.py`~~ (done).
3. Live CMA-ES short dock on 1G9V → close C/D/F; full 2e6 A/B → close E.
4. Build `.sif` on Linux x86_64 (Apptainer); in-container smoke + fingerprint → close G.
