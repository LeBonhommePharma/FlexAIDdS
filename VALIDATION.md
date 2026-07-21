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

**Workers:** CMA-ES swarm chunks 1–6 + orchestrator merge + live 1G9V docks.  
**Integration branch:** `feat/cmaes-search-backend` @ `7a9e06de1`.  
**Binary:** `.swarm/cmaes/orchestrator/build_fast/FlexAIDdS` sha256 `404b3ccddc22c12bf3cfaced9b0eaf996faa16d5caa501d7ff691282d66f9eb1`.  
**Toolchain:** `g++-16` (Homebrew GCC 16.1.0), CMake 4.4.0, host Darwin arm64 (`LPmore.local`).  
**Complex:** Astex **1G9V** (`1G9V_apo.pdb` + `1G9V_dockin.sdf`).

| # | Acceptance item | Owner | Status | Evidence (path + one-line proof) |
|---|---|---|---|---|
| A | `cmake configure` accepts the wiring AND `cmake/ValidateSources.cmake` passes with the new TUs | orchestrator | **CLOSED** | exit 0 (`validation_evidence/build_ab/A_cmake_configure.exit`); `cmaes_search.cpp` in `LIB/CMakeLists.txt`; log `A_cmake_configure.log` |
| B | Full engine + adapter **link** into one `FlexAIDdS` binary (`BUILD_FLEXAIDDS_FAST`) | orchestrator | **CLOSED** | sha256 `404b3ccddc22c12bf3cfaced9b0eaf996faa16d5caa501d7ff691282d66f9eb1` (`B_binary.sha256`); `B_build.log` Built target FlexAIDdS |
| C | Adapter runs against the **real `ic2cf`** (not the mock) with no exception; snapshot fills for a real ligand's DOF count | build-ab | **CLOSED** | Live CMA-ES 1G9V: exit 0, DOF=10, n_snap=32, evals=3000, best_cf=`457.536781` — `validation_evidence/build_ab/C_cmaes_smoke/C_result.txt` + `run.log` + pose `cmaes_smoke_0.pdb` |
| D | Live dock on one Astex complex: elected-pose **RMSD** + best **CF**, both GA and CMA-ES arms | build-ab | **CLOSED** | Short dual 5000 evals: **GA** CF=`-37.11991` RMSD=`10.1141` Å; **CMA-ES** CF=`643.054551` RMSD=`12.1294` Å — `D_1G9V_dual/D_dual_result.txt`, `D_ga_result.txt`, `D_cmaes_result.txt` |
| E | GA-vs-CMA-ES **A/B**, eval-matched at 2e6 | build-ab | **CLOSED** | Claim budget 1000×2000=2e6 each: **GA** CF=`-68.55885` RMSD=`5.5590` Å (355 s); **CMA-ES** evals=`2000000` CF=`-20.18888` RMSD=`6.5245` Å (526 s); ΔCF=`+48.37`, ΔRMSD=`+0.965` (CMA−GA) — `E_ab_2e6/E_ab_summary.txt` (UTC 20:14:33–20:23:42) |
| F | Real entropy trace on the **rugged** real surface | build-ab | **CLOSED** | CMA-ES 2e6 entropy CSV 2000 gens: best_cf end=`-20.188879`, F_end=`-24.307`, H_energy 0.056→6.908; fingerprint sha256 `6b051f0d5671cde37309f063fecf804712e6298eb662cfbd9bfe00560f3a0cac` — `F_trace/F_result.txt`, `F_cmaes_2e6_fingerprint.json`, trace `E_ab_2e6/cmaes/cmaes_2e6_cmaes_entropy.csv` (note: H_search flat at ln(λ) in current adapter) |
| G | Locked-arch `.sif` builds; in-container dock; collapse fingerprint INVARIANT | harness | **OPEN** (local half advanced) | Recipes + manifest OK (`containers/flexaidds_locked_x86_64.def`, `scripts/narval_cmaes_array.sh`, `G_harness/G_manifest_validate.txt` = schema OK 1000×2000); apptainer/sbatch **MISSING** on Darwin — no `.sif` (`G_harness/G_local.txt`); Narval submit print-only |

## Bottom line
Host validation closed **A–F (6/7)** with real on-disk docks on **1G9V**. Eval-matched 2e6 A/B: GA finds better elected CF (`-68.56` vs CMA-ES `-20.19`) and slightly better RMSD (`5.56` vs `6.52` Å) under identical scoring. **G** remains OPEN until Apptainer builds a locked `.sif` on Linux/Narval.

### Remaining (G only)
1. On Linux x86_64 with Apptainer: `apptainer build … flexaidds_locked_x86_64.sif`
2. In-container smoke dock + `collapse_fingerprint.py` on smoke/2e6 traces → cross-arch INVARIANT
3. `sbatch --account=${CC_ACCOUNT} scripts/narval_cmaes_array.sh` on a login node
