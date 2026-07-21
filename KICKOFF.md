# KICKOFF — CMA-ES validation dispatch (local host)

**Branch expectation:** integration applied + committed; base science commit `2a60f65`.
**Host:** macOS with GCC≥14 (Homebrew) + CMake≥3.28 for FAST build; Apptainer for harness `.sif`.

## Parallel workers

### 1) `flexaidds-build-ab` (items A–F)

1. Build the **FAST** target (`BUILD_FLEXAIDDS_FAST=ON`) with this machine's **GCC≥14** toolchain.
2. Doctor-check the binary / wiring (env, symbols, search backend gate).
3. Run the **eval-matched GA-vs-CMA-ES A/B** at **2e6 evals each** on **one Astex complex**
   (`λ=1000×2000 gen ≡ pop-1000×2000 gen`; scoring path identical — only the search operator differs).
4. Close VALIDATION items **A–F**, each with an **on-disk artifact** under
   `validation_evidence/build_ab/` (logs, sha256, RMSD/CF tables, entropy trace).

### 2) `flexaidds-harness` (item G local half)

1. Build the **locked-arch Apptainer `.sif`**.
2. Smoke-dock **one complex in-container**.
3. Validate the **manifest + array script**.
4. Run the **collapse fingerprint** on the smoke trace.
5. **Do NOT `sbatch` to Narval** unless this host is a Compute Canada login node —
   print the submit command instead.
6. Close VALIDATION item **G**'s **local half** with artifacts under
   `validation_evidence/harness/`.

## Orchestrator (after both return)

1. Run `analysis/collapse_fingerprint.py` on the **build-ab CMA-ES trace**.
2. Stitch one summary (paths + numbers from disk only).
3. Rewrite `VALIDATION.md` marking A–G **CLOSED/OPEN** with evidence paths per protocol.

## Non-negotiable

- Real on-disk numbers and provenance only.
- A gate failure with real logs is a **valid** result.
- Never fabricate a build, a pose, or a fingerprint.
