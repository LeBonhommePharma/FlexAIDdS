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

**Workers:** `flexaidds-build-ab` + `flexaidds-harness` (parallel).  
**Build/GA provenance HEAD:** `9723f9de41da1d185725384c6b774eca87fcbc4e` (subject claimed CMA-ES; tree content was WebGPU/ThermoWhiteboard only).  
**Re-check on this workspace:** integration sources still **absent** (`validation_evidence/orchestrator/baseline_inventory_main_2a60f65.txt`).  
**Toolchain:** `g++-16` (Homebrew GCC 16.1.0), CMake 4.4.0, host Darwin arm64 (`LPmore.local`).  
**Evidence pack commit:** `docs/cmaes-validation-host-ag` @ `251aa4ae9`.  
**Orchestrator summary:** `validation_evidence/orchestrator/SUMMARY.md`

| # | Acceptance item | Owner | Status | Evidence (path + one-line proof) |
|---|---|---|---|---|
| A | `cmake configure` accepts the wiring AND `cmake/ValidateSources.cmake` passes with the new TUs | build-ab | **OPEN** | Base FAST cmake exit 0 (`validation_evidence/build_ab/A_cmake_configure.exit`); **no CMA-ES TUs** — `A_cmaes_absence_proof.txt`, `A_inventory.txt` |
| B | Full engine + adapter **link** into one `FlexAIDdS` binary (`BUILD_FLEXAIDDS_FAST`) | build-ab | **OPEN** | Base binary sha256 `967d194698b048513edabf4038e105ff50b5ccd9eea83ceeb975fe22aaef2f6a` (`B_binary.sha256`); adapter not in link |
| C | Adapter runs against the **real `ic2cf`** (not the mock) with no exception; snapshot fills for a real ligand's DOF count | build-ab | **OPEN** | `C_doctor.txt` — zero cmaes/`FLEXAIDDS_SEARCH` binary symbols; code hits outside validation_evidence = 0 |
| D | Live dock on one Astex complex: elected-pose **RMSD** + best **CF**, both GA and CMA-ES arms | build-ab | **OPEN** | GA short 1G9V: elected CF=`-5.08119`, RMSD=`13.0879` Å, best CF=`-16.70442` (`D_ga_result.txt`); CMA-ES arm not run (`D_cmaes_result.txt`) |
| E | GA-vs-CMA-ES **A/B**, eval-matched at 2e6 | build-ab | **OPEN** | `E_ab_summary.txt` — no CMA-ES arm; ΔRMSD/ΔCF = N/A |
| F | Real entropy trace on the **rugged** real surface | build-ab | **OPEN** | `F_entropy_trace.txt` — no CMA-ES `H_search`/`H_energy`; GA partial only |
| G | Locked-arch `.sif` builds; in-container dock; collapse fingerprint INVARIANT | harness | **OPEN** (local half) | `harness/ITEM_G.md` — apptainer MISSING (`G_tooling.txt`); no `.sif`/recipe; no `analysis/collapse_fingerprint.py`; Narval sbatch print-only (`G_narval_submit_command.txt`) |

## Bottom line
Correct and functional **to the compiler gate** in the **sandbox that produced P1–P7** — those claims are not re-proven in this host tree because **`LIB/cmaes_search.cpp` and `apply_integration.sh` are not present**. Host dispatch closed **0 / 7**. Real partial progress only: FAST base link + short GA dock on **1G9V**. Nobody calls this "docking" validation until D and E have real RMSD on **both** arms.

### Unblock checklist (required before re-dispatch)
1. Land the CMA-ES package: `LIB/cmaes_search.{cpp,h}`, top/CMake wiring, `apply_integration.sh`, `CMAES_INTEGRATION.md`.
2. Land `analysis/collapse_fingerprint.py`.
3. Land Apptainer recipe + build `.sif` on Linux; fingerprint smoke trace.
4. Re-run KICKOFF workers; fill A–G only from new on-disk artifacts.
