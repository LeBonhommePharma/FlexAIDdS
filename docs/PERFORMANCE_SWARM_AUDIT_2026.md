# FlexAIDdS Performance Swarm Audit — June 2026

**Scope:** Cross-platform performance audit (Wave 1 agents) + P0/P1 implementation  
**Constraint:** Strict Astex/CASF accuracy — no GA budget or ranking changes  
**Baselines:** `results/perf_swarm/baseline_macos_metal.json` (182 harvested dock timings)

---

## Executive Summary

Simulation wall-clock is dominated by **steady-state GA offspring scoring** (`reproduce()` → serial `ic2cf`), not post-GA clustering. Wave 1 identified that GPU batch paths in `calculate_fitness()` are largely unused under default BOOM because offspring are scored inline with `status='n'`.

**P0 (implemented, `27e68e51`):** Parallel `reproduce()` eval, Vcontacts `inv_d12`, Metal SAS 512 + async drain, `flexaid_core` native flags.

**P1 (this wave):** CUDA wired into `FlexAIDdS`, orchestration budget logging, Python I/O fix, SIMD RMSD path, clustering timing probe, `perf.yml` baseline hook.

---

## Time Budget Model (Measured + Structural)

| Subsystem | Est. % of dock | Evidence |
|-----------|----------------|----------|
| GA scoring (`ic2cf` / Vcontacts) | **82–92%** | Wave 1 ga_hotpath agent |
| PSHARE niche sharing | 5–14% | `gaboom.cpp` O(pop²) |
| Post-GA clustering (CF/DP) | 0.3–3% | clustering_entropy agent |
| Post-GA clustering (FO) | 5–20% | When `clustering_algorithm=FO` |
| Orchestration / I/O | 1–5% | C++ path low; Python path higher |

**Harvested timings (macOS):** typical smoke ~**140 ms/gen**; outlier **1N1M ~1299 ms/gen** (v128 repro).

---

## Cross-Platform Impact Matrix

| Recommendation | macOS Metal | Linux AVX2 | Linux CUDA | Status |
|----------------|-------------|------------|------------|--------|
| Parallel `reproduce()` eval | High | High | High | **P0 done** |
| Metal async CF + SAS 512 | High | N/A | N/A | **P0 done** |
| `flexaid_core` `-mcpu=native` | High | High | High | **P0 done** |
| Vcontacts `inv_d12` | Medium | Medium | Medium | **P0 done** |
| CUDA on `FlexAIDdS` | N/A | N/A | High | **P1 done** |
| OMP/restart budget logging | Medium | Medium | Medium | **P1 done** |
| Python `DEVNULL` subprocess | Low | Low | Low | **P1 done** |
| `simd::rmsd` → `sum_sq_distances_f` | Medium | High | Medium | **P1 done** |
| Shannon Metal buffer reuse | Medium | N/A | N/A | P2 pending |
| Parallel coord cache | Blocked | Blocked | Blocked | `calc_rmsd_chrom` mutates `atoms[]` |
| FOPTICS OpenMP | Low–Med | Low–Med | Low–Med | P2 (FO path only) |
| SoA distances default ON | Medium | High | Low | P2 + tier-1 gate |

---

## Implementation DAG (Remaining)

### P2 — Medium effort, tier-1 per PR

1. Shannon Metal persistent buffers (`ShannonMetalBridge.mm`)
2. `FLEXAIDS_USE_SOA_DISTANCES` validation on real `ic2cf` paths
3. FOPTICS distance OpenMP / Metal precompute when FO enabled
4. `calc_rmsd_chrom` thread-local atom copies → parallel coord cache

### P3 — Structural

1. SoA atom layout (`flexaid.h`)
2. Runtime CPU dispatch (AVX2/AVX-512 fat binary)
3. PGO build option
4. Cross-generation Metal GPU pipelining in `gaboom.cpp`

---

## Accuracy Gate Checklist (per PR)

1. `ctest --test-dir build --output-on-failure` (relevant targets)
2. Tier-1 Astex subset: `python -m flexaidds.dataset_runner --dataset astex_diverse --tier 1`
3. Compare `regression_flags` in report JSON vs YAML baselines
4. Microbench: `benchmark_vcfbatch` speedup ≥ 1.0; `LookupPerformance` < 50 ns

---

## CI Integration

- **Tier-1/2 workflows:** accuracy gates (unchanged)
- **`perf.yml`:** manual dispatch; use `scripts/compare_perf_baseline.py` with `baseline_sha` artifact
- **Harvest script:** `scripts/harvest_perf_baselines.py` (low-I/O, queue-safe)

---

## Agent Reports

Full Wave 1 output: `results/perf_swarm/reports/`