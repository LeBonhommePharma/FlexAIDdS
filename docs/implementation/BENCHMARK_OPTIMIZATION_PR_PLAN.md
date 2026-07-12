# Benchmark Optimization PR Plan (DAG)

**Source:** Astex entropy benchmark post-mortem + Metal dispatch review (2026-07-07).  
**Target repos:** `Projects/FlexAIDdS` (core), `Documents/.../current_master_code` (orchestrator configs).

## DAG (topological order)

```
PR-01 ──┬──> PR-03 (fitness flag uses fixed sphere ingest)
        └──> PR-06 (benchmark configs)

PR-02 (grid cache) ──> PR-06

PR-04 (Shannon threshold) ── independent

PR-05 (resolve_build profiles) ──> PR-06
```

## PR stack

| ID | Branch | Scope | Test gate |
|----|--------|-------|-----------|
| PR-01 | `opt/read-spheres-hetatm` | `read_spheres.cpp` accepts HETATM; unit test | `ctest -R read_spheres` or new test |
| PR-02 | `opt/grid-cache-datasetrunner` | `DatasetRunner` sets `FLEXAIDDS_GRID_CACHE_DIR` | build `benchmark_datasets` |
| PR-03 | `opt/fitness-backend-env` | `FLEXAIDDS_FITNESS_BACKEND` env (default `cpu`) | `test_hardware_dispatch` |
| PR-04 | `opt/shannon-gpu-threshold` | Env-tunable Shannon Metal min-N | `test_hardware_dispatch` |
| PR-05 | `opt/resolve-build-profiles` | `resolve_build.py --profile metal\|lto` | `validate_skill.py` smoke |
| PR-06 | `opt/astex-config-profiles` | YAML `build_profile`, grid cache paths | orchestrate `--dry-run` |

## Out of scope (future milestones)

- Metal FITNESS_EVAL parity validation suite (needs PR-03 + dedicated benchmark branch)
- Orchestrator tool pipelining (Boltz while GA runs)
- ARM64 NEON VoronoiCFBatch

## RAM contract (18 GB)

All PRs must preserve: `threads=1`, `parallel_restarts=0`, `omp_threads≤4`.