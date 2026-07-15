# ProtocolConfig — typed protocol knobs (getenv consolidation)

**Status:** Chunk 1 (foundation)  
**Owner path:** `LIB/ProtocolConfig.{h,cpp}`  
**Related audit:** ~64 `getenv` / `std::getenv` call sites in
`LIB/DatasetRunner.cpp`, `LIB/gaboom.cpp`, `LIB/top.cpp` (repo-wide getenv
surface is larger; this doc tracks the audit trio first).

## Goal

Replace ad-hoc `FLEXAIDDS_*` environment reads with **one typed, serializable
`flexaids::ProtocolConfig`**. Environment variables remain a **compatibility
adapter only** (`ProtocolConfig::from_env()`). Future chunks can load JSON /
CLI / skill configs into the same struct without touching call sites again.

## What landed in chunk 1

| Piece | Location |
|-------|----------|
| Typed struct + `defaults()` / `from_env()` / `to_json()` / `from_json()` | `LIB/ProtocolConfig.h`, `LIB/ProtocolConfig.cpp` |
| Unit tests (defaults, env overrides, JSON round-trip) | `tests/test_protocol_config.cpp` |
| DatasetRunner constructor + dock-config emission path | `LIB/DatasetRunner.cpp` / `.h` |
| Engine init: `FLEXAIDDS_DATA_DIR`, `FLEXAIDDS_HBOND_WEIGHT` | `LIB/top.cpp` |

### Migrated env vars (go through `ProtocolConfig`)

| Env var | Field | Default | Migrated call sites |
|---------|-------|---------|---------------------|
| `FLEXAIDDS_SEED_BASE` | `seed_base` | `0` | DatasetRunner `deterministic_ga_seed` |
| `FLEXAIDDS_RESTARTS` | `restarts` | `5` | DatasetRunner multi-restart loop |
| `FLEXAIDDS_PARALLEL_RESTARTS` | `parallel_restarts` | `restarts > 1` | DatasetRunner restart launcher |
| `FLEXAIDDS_VCT_R0` | `vct_r0` | `7.0` | DatasetRunner dock_config scoring |
| `FLEXAIDDS_VCT_NORM` | `vct_normalize_contacts` | off (presence) | DatasetRunner dock_config scoring |
| `FLEXAIDDS_VCT_ENTROPY_WEIGHT` | `vct_entropy_weight` | `0.0` | DatasetRunner dock_config scoring |
| `FLEXAIDDS_SHARING_ALPHA` | `sharing_alpha` (optional) | pop-scaled `4.0` | DatasetRunner GA section |
| `FLEXAIDDS_BOOM_FRAC` | `boom_frac` (optional) | `1.0` | DatasetRunner GA section |
| `FLEXAIDDS_N_ELITE` | `n_elite` | `1` | DatasetRunner GA section |
| `FLEXAIDDS_USE_SHANNON` | `use_shannon` | off (presence) | DatasetRunner dock_config |
| `FLEXAIDDS_THERMO` | `thermo_enabled` | off (presence) | DatasetRunner dock_config |
| `FLEXAIDDS_T_EFF` | `t_eff` | `0.596` | DatasetRunner thermo_engine |
| `FLEXAIDDS_TENCOM_SCALE` | `tencom_scale` | `1.0` | DatasetRunner thermo_engine |
| `FLEXAIDDS_DATA_DIR` | `data_dir` | empty | DatasetRunner matrix + `--data-dir`; `top.cpp` auto-detect |
| `FLEXAIDDS_CF_WINDOW_SELECTOR` | `cf_window_selector` | `false` | DatasetRunner ctor |
| `FLEXAIDDS_CLUSTER_MEMBER_EMIT` | `cluster_member_emit` | `false` | DatasetRunner ctor |
| `FLEXAIDDS_HBOND_WEIGHT` | `hbond_weight` | `-2.5` | `top.cpp` FA init |

**Estimate:** ~17 unique vars / ~22 call sites migrated; residual in the audit
trio is roughly **~42 getenv call sites** (including `HOME` / `OMP_*` and
double-reads).

## Residual getenv inventory (TODO — later chunks)

### `LIB/DatasetRunner.cpp` (owner: benchmark / DatasetRunner)

| Env var | Purpose | Default / notes | Priority |
|---------|---------|-----------------|----------|
| `FLEXAIDDS_SEED_ELITISM` | Crystal seed elitism gate | mode-dependent | P1 |
| `FLEXAIDDS_FREQSEL` | Frequency-gated pose selector | off | P1 |
| `FLEXAIDDS_FREQSEL_ALPHA` | Freqsel soft weight | `12.0` | P1 |
| `FLEXAIDDS_FREQSEL_RMSD` | Freqsel RMSD bin | `1.5` | P1 |
| `FLEXAIDDS_SEED_ELITISM_DELTA_CF` | Seed elitism CF window | `10.0` | P1 |
| `FLEXAIDDS_HTTP_RETRIES` | RCSB download retries | `3` | P2 |
| `FLEXAIDDS_ORACLE_SITE_DIR` | Oracle sphere directory | unset | P1 |
| `FLEXAIDDS_USE_DP` | Force DensityPeak clustering | off | P2 |
| `FLEXAIDDS_BINARY` / `_BUILD` / `_REPO` | FlexAIDdS binary resolution | path search | P2 |
| `FLEXAIDDS_PRIORITY_TARGETS` | Target priority list | unset | P2 |
| `FLEXAIDDS_IGNORE_CACHE` | Skip completed-cache | off | P2 |
| `FLEXAIDDS_RING_FLEX` | Ring pucker DoF | off | P1 |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | DoF budget scaling mode | `1` | P1 |
| `FLEXAIDDS_BUDGET_SCALE` | Extra budget scale | on | P1 |
| `FLEXAIDDS_FINE_GRID` | Finer angle/dihedral steps | off (presence) | P2 |
| `OMP_NUM_THREADS` | Worker OMP sizing | hardware | leave (OpenMP contract) |
| `HOME` | `~` expansion | `/tmp` | leave (platform) |
| `FLEXAIDDS_COGNATE_SITE` | Cognate redock site inject | off | P1 |
| `FLEXAIDDS_SCORE_NATIVE` | Native CF diagnostic | off | P1 |
| `FLEXAIDDS_MULTI_CLEFT` | Parallel cleft regions | off | P2 |
| `FLEXAIDDS_CONSENSUS_SCORER` | Cross-restart consensus | off | P1 |
| `FLEXAIDDS_HVIB` | tENCoM/H(ω) validator | on unless `=0` | P1 |
| `FLEXAIDDS_THERMO_CSV` | Extra thermo CSV columns | off | P2 |

### `LIB/gaboom.cpp` (owner: GA engine)

| Env var | Purpose | Default / notes | Priority |
|---------|---------|-----------------|----------|
| `FLEXAIDDS_CHAIN_NORM` | Multi-chain receptor norm | off | P2 |
| `OMP_NUM_THREADS` | Default OMP=2 if unset | platform | leave |
| `FLEXAIDDS_USE_SHANNON` | Enable Shannon in fitness | off | P1 (already on ProtocolConfig; wire gaboom) |
| `FLEXAIDDS_INSTREAM_INTERVAL` | In-stream clustering interval | compile default | P2 |
| `FLEXAIDDS_NO_SEC` | Disable SEC early exits | off | P1 |
| `FLEXAIDDS_T_HOT` | Hot anneal temperature | `0.0` | P1 |
| `FLEXAIDDS_BENCHMARK` | Full-gen benchmark mode | off | P1 |
| `FLEXAIDDS_SMFREE_REQUIRE_T` | Require T for SMFREE | off | P2 |

### `LIB/top.cpp` (owner: engine entry)

| Env var | Purpose | Default / notes | Priority |
|---------|---------|-----------------|----------|
| `FLEXAIDDS_GRID_CACHE_DIR` | Grid cache directory | unset | P2 |
| `FLEXAIDDS_LIGAND_BATCH` | Batch ligand directory | unset | P2 |
| `FLEXAIDS_VH_DEBUG` | Vector H-bond debug dump | off | P3 |
| `FLEXAIDDS_CLEFT_SPHERE_FILE` | Explicit cleft spheres | unset | P1 |
| `FLEXAIDDS_ORACLE_SITE` | Oracle site path | unset | P1 |
| `FLEXAIDDS_SCORE_NATIVE` / `FLEXAIDDS_NATIVE_ONLY` | Native scoring path | off | P1 |

### Also still env-backed (out of chunk-1 audit trio but related)

`LIB/config_parser.cpp` still applies `FLEXAIDDS_VCT_ENTROPY_WEIGHT`,
`FLEXAIDDS_N_ELITE`, `FLEXAIDDS_SHARING_ALPHA`, `FLEXAIDDS_BOOM_FRAC`,
`FLEXAIDDS_ENTROPY_WEIGHT`, `FLEXAIDDS_DIVERSITY_MONITORING`, ranking flags.
**Chunk 2 recommendation:** have `apply_config()` accept an optional
`const ProtocolConfig*` so engine-side overrides share the same typed source
as DatasetRunner (env adapter still works for standalone CLI).

## Migration rules (for later PRs)

1. **Do not rewrite DatasetRunner in one PR.** One cluster of related knobs per chunk.
2. Add fields to `ProtocolConfig` first, extend `from_env()` / JSON, then flip call sites.
3. Keep default values byte-identical when env is unset (tests + golden dock_config when practical).
4. Prefer reading `protocol_cfg_` (or a passed-in snapshot) over re-calling `getenv` in loops.
5. Leave platform env (`HOME`, `PATH`, `OMP_*`) as raw getenv unless a strong reason exists.
6. Document each migrated var in the table above and drop it from Residual.

## Usage sketch

```cpp
#include "ProtocolConfig.h"

// Compatibility (benchmarks / CLI):
auto cfg = flexaids::ProtocolConfig::from_env();

// Future: skill / JSON protocol file
// auto cfg = flexaids::ProtocolConfig::from_json(text);

// Engine or runner:
int n = cfg.restarts;
double alpha = cfg.effective_sharing_alpha(pop_base, pop_scaled);
```

## Tests

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_protocol_config -j
ctest --test-dir build -R ProtocolConfigTests --output-on-failure
```
