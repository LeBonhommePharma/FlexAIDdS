# ITEM G — Locked-arch harness (local half)

**Owner:** flexaidds-harness  
**Host:** LPmore.local (Darwin arm64) — **not** a Compute Canada / Narval login node  
**Branch/head:** `perf/pb-clash-grid-hoist` @ `9723f9de41da1d185725384c6b774eca87fcbc4e`  
**Timestamp (UTC):** 2026-07-21T19:33:00Z  

## Status (G local half): **OPEN**

CLOSED requires all three on disk: (1) locked-arch `.sif` built, (2) in-container smoke dock, (3) collapse fingerprint of the smoke trace. **None of the three exist.** Truthful OPEN with real blockers is the valid result.

## Acceptance vs evidence

| Gate | Result | Blocker (real) | Evidence path |
|------|--------|----------------|---------------|
| Tooling: apptainer/singularity | FAIL | both MISSING (exit 1) | `validation_evidence/harness/G_tooling.txt` |
| Narval login node? | NO | macOS laptop; no cvmfs/slurm/sbatch | `G_tooling.txt` |
| Container definition file | ABSENT | only FlexAID `AMINO.def`/`NUCLEOTIDES.def` (atom types) | `G_inventory.txt` |
| Locked-arch `.sif` build | SKIPPED | no tooling + no recipe | `G_sif_build.log`, `G_sif_build_command.txt` |
| `G_sif.sha256` | ABSENT | no binary to hash (by design) | — |
| In-container smoke dock | SKIPPED | no `.sif` | `G_smoke_dock_command.txt` |
| Manifest + array script | OPEN | no locked-arch manifest; no `#SBATCH` CMA-ES array | `G_manifest_validate.txt` |
| Collapse fingerprint | OPEN | no `analysis/collapse_fingerprint.py`; no smoke trace | `G_collapse_fingerprint.out` |
| Narval `sbatch` | NOT EXECUTED | host is not login node (rule) | `G_narval_submit_command.txt` |

## One-line proofs (paths)

- Inventory: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_inventory.txt` — no container `.def`, no `.sif`, no `analysis/collapse_fingerprint.py`, no Narval array script.
- Tooling: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_tooling.txt` — `apptainer`/`singularity`/`sbatch`/`sinfo` exit 1; hostname `LPmore.local`.
- SIF build log: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_sif_build.log` — `status=OPEN_BLOCKER`, skip reason recorded.
- Intended SIF build: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_sif_build_command.txt` — `STATUS=NOT_EXECUTED`.
- Smoke dock: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_smoke_dock_command.txt` — `STATUS=NOT_EXECUTED`.
- Manifest validate: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_manifest_validate.txt` — `STATUS=OPEN`.
- Fingerprint: `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_collapse_fingerprint.out` — script + trace both absent.
- Narval submit (print-only): `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/harness/G_narval_submit_command.txt`.
- Related tree gap (CMA-ES sources): `/Users/lp.more/Projects/FlexAIDdS/validation_evidence/orchestrator/baseline_inventory.txt` — `absent=LIB/cmaes_search.*`, `absent=analysis/collapse_fingerprint.py`, `apptainer=MISSING`.

## Full intended Narval submit command (NOT EXECUTED)

```bash
sbatch --job-name=flexaidds-g-smoke --account=${CC_ACCOUNT} --time=02:00:00 --cpus-per-task=8 --mem=16G --output=validation_evidence/harness/narval_smoke_%j.out --wrap='module load apptainer; apptainer exec --bind $SCRATCH/FlexAIDdS:/work $SCRATCH/FlexAIDdS/validation_evidence/harness/flexaidds_locked_x86_64.sif /opt/flexaidds/bin/FlexAIDdS --search cmaes --complex 1G9V --ga-population 1000 --ga-generations 10 --out /work/validation_evidence/harness/smoke_1G9V_narval'
```

Caveat: best-effort only — repo has no committed `scripts/narval_cmaes_array.sh` or locked `.sif`/`.def`; `${CC_ACCOUNT}` and image path must be real on Narval before submit.

## What would close G local half

1. Add an Apptainer recipe (locked arch/compiler) and build on Linux → write `.sif` + `G_sif.sha256`.
2. `apptainer exec …` smoke dock on one complex → logs with real RMSD/CF.
3. Land `analysis/collapse_fingerprint.py` and run it on the smoke entropy/trace → invariant fingerprint on disk.
4. (Remote half) Run the printed `sbatch` on a real Narval login node and compare fingerprints cross-arch within tol.

## VALIDATION.md

**Not rewritten** (orchestrator owns status table). This file is the harness deliverable for item G local half.
