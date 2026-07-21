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
**Build/GA HEAD:** `9723f9de41da1d185725384c6b774eca87fcbc4e` on `perf/pb-clash-grid-hoist`.  
**Toolchain:** `g++-16` (Homebrew GCC 16.1.0), CMake 4.4.0, host `LPmore.local` Darwin arm64.  
**Orchestrator summary:** `validation_evidence/orchestrator/SUMMARY.md`  
**Protocol note:** Commit subject claims CMA-ES wiring; tree content is WebGPU/ThermoWhiteboard only. Adapter sources (`LIB/cmaes_search.*`, `apply_integration.sh`, `analysis/collapse_fingerprint.py`) are **absent** (also absent on `cmaes-ab` @ `9a5e8c8d3`). Truthful OPEN with real logs is the valid result.

| # | Acceptance item | Owner | Status | Evidence (path + one-line proof) |
|---|---|---|---|---|
| A | `cmake configure` accepts the wiring AND `cmake/ValidateSources.cmake` passes with the new TUs | build-ab | **OPEN** | `validation_evidence/build_ab/A_cmake_configure.exit`=`0` for base FAST configure (`-DFLEXAIDS_USE_METAL=OFF`, `g++-16`); **no CMA-ES TUs** present — `A_cmaes_absence_proof.txt`, `A_inventory.txt` |
| B | Full engine + adapter **link** into one `FlexAIDdS` binary (`BUILD_FLEXAIDDS_FAST`) | build-ab | **OPEN** | Base binary linked: `validation_evidence/build_ab/B_binary.sha256`=`967d194698b048513edabf4038e105ff50b5ccd9eea83ceeb975fe22aaef2f6a` (`B_build.log` 100% FlexAIDdS); adapter TUs not in link → cannot close |
| C | Adapter runs against the **real `ic2cf`** (not the mock) with no exception; snapshot fills for a real ligand's DOF count | build-ab | **OPEN** | `validation_evidence/build_ab/C_doctor.txt` — zero `cmaes` / `FLEXAIDDS_SEARCH` symbols; source hits for `FLEXAIDDS_SEARCH` in `LIB/`/`src/` = 0 |
| D | Live dock on one Astex complex: elected-pose **RMSD** + best **CF**, both GA and CMA-ES arms | build-ab | **OPEN** | GA short (100×50) **1G9V**: elected CF=`-5.08119`, ordered heavy RMSD=`13.0879` Å, best CF=`-16.70442` (`D_ga_result.txt`, pose `D_1G9V/ga_short_0.pdb`); CMA-ES arm **not run** (`D_cmaes_result.txt` blocker) |
| E | GA-vs-CMA-ES **A/B**, eval-matched at 2e6 (λ=1000×2000gen ≡ pop-1000×2000gen); ΔRMSD/ΔCF attributable to the **search operator alone** (scoring identical per P4) | build-ab | **OPEN** | `validation_evidence/build_ab/E_ab_summary.txt` + `E_dock_config_2e6.json` prepared; no CMA-ES arm → ΔRMSD/ΔCF = N/A; full 2e6 GA not launched |
| F | Real entropy trace on the **rugged** real surface — expect `H_energy` to be informative here, unlike the single-basin mock (P6) | build-ab | **OPEN** | `validation_evidence/build_ab/F_entropy_trace.txt` — no CMA-ES `H_search`/`H_energy`; GA-only partial `H_final=2.300787` / F=`-16.7049` is not the required CMA-ES rugged-surface trace |
| G | Locked-arch `.sif` builds; in-container dock runs; cross-arch collapse fingerprint is **INVARIANT** within tol | harness | **OPEN** (local half) | `validation_evidence/harness/ITEM_G.md` — apptainer/singularity MISSING (`G_tooling.txt`); no container recipe / `.sif` (`G_sif_build.log`); no `analysis/collapse_fingerprint.py` (`G_collapse_fingerprint.out`); Narval `sbatch` **print-only** (`G_narval_submit_command.txt`) — host not login node |

## Worker deliverables (absolute)

- Build-ab table: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/ITEMS_A_F.md`
- Build-ab DONE: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/DONE.txt`
- Harness G: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/ITEM_G.md`
- Orchestrator summary: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/orchestrator/SUMMARY.md`
- Collapse attempt: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/orchestrator/collapse_fingerprint_attempt.log`
- Baseline inventory: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/orchestrator/baseline_inventory.txt`

## Orchestrator post-steps (executed)

1. `python3 analysis/collapse_fingerprint.py` on build-ab CMA-ES trace → **FAIL** (script absent; no CMA-ES entropy series). Log: `validation_evidence/orchestrator/collapse_fingerprint_attempt.log`.
2. Stitched summary: `validation_evidence/orchestrator/SUMMARY.md`.
3. This file rewritten with A–G **OPEN** + evidence paths (no fabricated CLOSED).

## Bottom line
PROVEN table (P1–P7) is unchanged and was **not re-verified in this host tree** — the adapter sources those claims name are missing from this checkout. Host dispatch closed **0 / 7** items. Real partial progress: FAST base binary link + short GA live dock on **1G9V** with on-disk CF/RMSD. Nobody calls this a validated CMA-ES dock until D and E carry a real RMSD on **both** arms at eval-matched budget, and G has a real `.sif` + fingerprint.
