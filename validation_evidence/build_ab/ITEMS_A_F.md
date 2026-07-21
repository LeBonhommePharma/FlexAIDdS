# ITEMS A–F — flexaidds-build-ab evidence

**Host tree:** `/Users/lp.more/Projects/FlexAIDdS`  
**Branch (at build):** `perf/pb-clash-grid-hoist`  
**HEAD (at build):** `9723f9de41da1d185725384c6b774eca87fcbc4e`  
**Also checked:** branch `cmaes-ab` (`9a5e8c8d3`) — no `LIB/cmaes_search.*` there either  

**Toolchain:** `g++-16` (Homebrew GCC 16.1.0), CMake 4.4.0  

Do not treat PROVEN P1–P7 as present in this tree. Those claims assume adapter sources that are **not on disk** here (see inventory).

| # | Status | Evidence |
|---|--------|----------|
| A | **OPEN** | Wiring missing: `LIB/cmaes_search.cpp` / `.h`, `apply_integration.sh` absent (`A_inventory.txt`, `A_cmaes_absence_proof.txt`, `A_external_search.txt`). CMake configure of base tree with `-DFLEXAIDS_USE_METAL=OFF` exit 0 (`A_cmake_configure.exit`=`0`; log `A_cmake_configure.log` ends “Build files have been written to …/build_fast”). ValidateSources ran with non-fatal orphans; **no CMA-ES TUs** to accept. First configure with Metal ON failed (metallib missing). Summary: `A_summary.txt`. |
| B | **OPEN** | Base FAST binary **did link** without adapter: `B_build.log` ends `[100%] Built target FlexAIDdS`; `B_binary.sha256`=`967d194698b048513edabf4038e105ff50b5ccd9eea83ceeb975fe22aaef2f6a`; `B_binary.ls` size 2663624 arm64 Mach-O. Acceptance requires **engine + adapter**; adapter TUs absent → cannot claim B CLOSED. |
| C | **OPEN** | Doctor: `C_doctor.txt` — `FLEXAIDDS_SEARCH` unset; binary help has **no** `--search`; `nm`/strings: **zero** cmaes / FLEXAIDDS_SEARCH symbols; source hits for `FLEXAIDDS_SEARCH` in `LIB/`/`src/` = **0**. Item needs CMA-ES adapter on real `ic2cf`; not present. (GA path did run real `ic2cf` under item D, but that does not close C.) |
| D | **OPEN** | GA arm live on Astex **1G9V** (short budget 100×50, not claim budget): `D_ga_run.log` exit 0, `D_ga_result.txt` — elected CF=`-5.08119` (pose `D_1G9V/ga_short_0.pdb`), best CF among ranks=`-16.70442`, ordered heavy RMSD vs crystal SDF ≈ `13.0879` Å (n=25/25). CMA-ES arm: **not run** — `D_cmaes_result.txt` blocker: no backend/symbols. Both arms required → OPEN. |
| E | **OPEN** | Config prepared for 1000×2000 (`E_dock_config_2e6.json`); intended commands in `E_ab_summary.txt`. Primary blocker: **no CMA-ES arm**. Full 2e6 GA not launched (session wall-clock; would not close E without CMA-ES). `ΔRMSD`/`ΔCF` = N/A. |
| F | **OPEN** | No CMA-ES entropy trace machinery / no CMA-ES run. Partial GA-only lines in `F_entropy_trace.txt` / `F_entropy_trace_ga_partial.log` (`H_final=2.300787`, Helmholtz F=`-16.7049`) — **not** the required CMA-ES `H_search`/`H_energy` rugged-surface trace. |

## CLOSED count: **0** / OPEN count: **6**

## Root blocker (shared)

CMA-ES search backend is **not integrated in this repository checkout**:

- `absent`: `LIB/cmaes_search.cpp`, `LIB/cmaes_search.h`, `apply_integration.sh`, `CMAES_INTEGRATION.md`, `analysis/collapse_fingerprint.py`
- `rg FLEXAIDDS_SEARCH` hits (code, excl. validation_evidence): **0**
- `cmaes-ab` branch tip also lacks `LIB/cmaes_search.cpp`
- HEAD subject says “wire CMA-ES…” but tree content is unrelated (WebGPU / other)

Until adapter sources land in tree + CMake `target_sources` + `FLEXAIDDS_SEARCH=cmaes` gate, items **A–F cannot close**.

## Key artifact paths (absolute)

- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/A_inventory.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/A_cmaes_absence_proof.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/A_cmake_configure.log`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/A_cmake_configure.exit`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/B_build.log`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/B_binary.sha256`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/B_binary.ls`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/build_fast/FlexAIDdS`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/C_doctor.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/D_ga_run.log`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/D_ga_result.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/D_cmaes_result.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/D_1G9V/`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/E_ab_summary.txt`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/E_dock_config_2e6.json`
- `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/build_ab/F_entropy_trace.txt`
