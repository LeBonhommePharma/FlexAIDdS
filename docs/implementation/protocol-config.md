# ProtocolConfig — typed protocol knobs (getenv consolidation)

**Status:** Chunk 3 (config_parser science knobs + RUN_RECEIPT)  
**Owner path:** `LIB/ProtocolConfig.{h,cpp}`, `LIB/RunReceipt.{h,cpp}`  
**Related audit:** DatasetRunner / gaboom / top / config_parser getenv consolidation.

## Goal

Replace ad-hoc `FLEXAIDDS_*` environment reads with **one typed, serializable
`flexaids::ProtocolConfig`**. Environment variables remain a **compatibility
adapter only** (`ProtocolConfig::from_env()`). Campaigns also write
`RUN_RECEIPT.json` (C0 schema-aligned) with a ProtocolConfig snapshot so mid-run
`setenv` cannot create dual-protocol provenance.

## Snapshot policy

| When | What |
|------|------|
| `DatasetRunner` constructor | `protocol_cfg_ = ProtocolConfig::from_env()` |
| `DatasetRunner::run()` entry | **Re-snapshot** `protocol_cfg_` so long-lived runners pick up env changes after construction |
| `apply_config()` | Single snapshot at entry (or caller-supplied `const ProtocolConfig*`) — **no mid-apply getenv** |
| `GA()` (gaboom) entry | Local `ProtocolConfig::from_env()` for engine-side knobs |
| `top.cpp` site/grid paths | Local snapshot at use site; apply_config receives explicit snapshot |

Chunk 1 only snapshotted the constructor; chunk 2 re-snapshots at `run()`;
chunk 3 eliminates dual getenv inside `apply_config` and emits RUN_RECEIPT.

## RUN_RECEIPT.json schema (C0-aligned)

Written by `DatasetRunner::run()` via `flexaids::write_run_receipt()` to
`<output_dir>/RUN_RECEIPT.json` (also keeps legacy `provenance.json`).

| Key | Type | Notes |
|-----|------|-------|
| `schema_version` | int | `1` |
| `run_id` | string | Dataset name or output basename |
| `started_utc` | string | ISO-8601 UTC |
| `output` | string | Output directory |
| `dataset` | string | Dataset label |
| `mode` | string | `oracle-ceiling` / `defined-cleft-redock` / `autonomous` / `unset` |
| `temperature_K` | number | DockingConfig temperature |
| `pop` | int | GA population |
| `gen` | int | GA generations |
| `restarts` | int | From ProtocolConfig |
| `seed_base` | uint | From ProtocolConfig |
| `seed_elitism` | 0/1 | Mode-aware (oracle ON; defined-cleft/autonomous OFF) |
| `matrix_path` / `matrix_md5` / `matrix_sha256` | string | Scoring matrix |
| `binary_path` / `binary_sha256` | string | FlexAIDdS docking binary |
| `runner_path` / `runner_sha256` | string | Host process when resolvable |
| `git_commit` | string | `git rev-parse HEAD` if available |
| `oracle_site_dir` / `oracle_site_dir_set` | string/bool | Oracle site |
| `protocol_config` | object | Full `ProtocolConfig::to_json()` snapshot |

C0 launch scripts under the iCloud queue may still write a thinner receipt before
start; the C++ path is a superset (adds `protocol_config`, hashes, git, schema).

## Migrated env vars

### Chunk 1

| Env var | Field | Default |
|---------|-------|---------|
| `FLEXAIDDS_SEED_BASE` | `seed_base` | `0` |
| `FLEXAIDDS_RESTARTS` | `restarts` | `5` |
| `FLEXAIDDS_PARALLEL_RESTARTS` | `parallel_restarts` | `restarts > 1` |
| `FLEXAIDDS_VCT_R0` | `vct_r0` | `7.0` |
| `FLEXAIDDS_VCT_NORM` | `vct_normalize_contacts` | off |
| `FLEXAIDDS_VCT_ENTROPY_WEIGHT` | `vct_entropy_weight` | `0.0` |
| `FLEXAIDDS_SHARING_ALPHA` | `sharing_alpha` | pop-scaled `4.0` |
| `FLEXAIDDS_BOOM_FRAC` | `boom_frac` | `1.0` |
| `FLEXAIDDS_N_ELITE` | `n_elite` | `1` |
| `FLEXAIDDS_USE_SHANNON` | `use_shannon` | off |
| `FLEXAIDDS_THERMO` | `thermo_enabled` | off |
| `FLEXAIDDS_T_EFF` | `t_eff` | `0.596` |
| `FLEXAIDDS_TENCOM_SCALE` | `tencom_scale` | `1.0` |
| `FLEXAIDDS_DATA_DIR` | `data_dir` | empty |
| `FLEXAIDDS_CF_WINDOW_SELECTOR` | `cf_window_selector` | `false` |
| `FLEXAIDDS_CLUSTER_MEMBER_EMIT` | `cluster_member_emit` | `false` |
| `FLEXAIDDS_HBOND_WEIGHT` | `hbond_weight` | `-2.5` |

### Chunk 2

| Env var | Field | Default | Call sites |
|---------|-------|---------|------------|
| `FLEXAIDDS_SEED_ELITISM` | `seed_elitism` | `true` | pose selector |
| `FLEXAIDDS_SEED_ELITISM_DELTA_CF` | `seed_elitism_delta_cf` | `10.0` | pose selector |
| `FLEXAIDDS_FREQSEL` | `freqsel` | `false` | pose selector |
| `FLEXAIDDS_FREQSEL_ALPHA` | `freqsel_alpha` | `12.0` | pose selector |
| `FLEXAIDDS_FREQSEL_RMSD` | `freqsel_rmsd` | `1.5` | pose selector |
| `FLEXAIDDS_CONSENSUS_SCORER` | `consensus_scorer` | `false` | DatasetRunner election |
| `FLEXAIDDS_HVIB` | `hvib_enabled` | ON unless `=0` | DatasetRunner |
| `FLEXAIDDS_RING_FLEX` | `ring_flex` | `false` | budget scaling |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | `eval_scale_dihedral` | `1` (`off`→`-1`) | budget scaling |
| `FLEXAIDDS_BUDGET_SCALE` | `budget_scale` | `true` | budget scaling |
| `FLEXAIDDS_FINE_GRID` | `fine_grid` | off (presence) | dock_config |
| `FLEXAIDDS_MULTI_CLEFT` | `multi_cleft` | `0` | dock cmd |
| `FLEXAIDDS_COGNATE_SITE` | `cognate_site` | off (presence) | dock cmd |
| `FLEXAIDDS_SCORE_NATIVE` | `score_native` | off | DatasetRunner + top |
| `FLEXAIDDS_NATIVE_ONLY` | `native_only` | off | top |
| `FLEXAIDDS_USE_DP` | `use_dp` | off (`=1`) | clustering |
| `FLEXAIDDS_IGNORE_CACHE` | `ignore_cache` | off | DatasetRunner |
| `FLEXAIDDS_THERMO_CSV` | `thermo_csv` | off (presence) | CSV export |
| `FLEXAIDDS_ORACLE_SITE_DIR` | `oracle_site_dir` | empty | Astex prep / provenance |
| `FLEXAIDDS_ORACLE_SITE` | `oracle_site` | empty | top grid/site |
| `FLEXAIDDS_CLEFT_SPHERE_FILE` | `cleft_sphere_file` | empty | top |
| `FLEXAIDDS_NO_SEC` | `no_sec` | off (presence) | gaboom |
| `FLEXAIDDS_BENCHMARK` | `benchmark_mode` | off (presence) | gaboom |
| `FLEXAIDDS_T_HOT` | `t_hot` | `0.0` | gaboom anneal |
| `FLEXAIDDS_INSTREAM_INTERVAL` | `instream_interval` | `0` (= compile default) | gaboom |
| `FLEXAIDDS_CHAIN_NORM` | `chain_norm` | off | gaboom |
| `FLEXAIDDS_SMFREE_REQUIRE_T` | `smfree_require_t` | off | gaboom fitness |

### Chunk 3 (config_parser science knobs)

| Env var | Field | Default when unset | Call sites |
|---------|-------|--------------------|------------|
| `FLEXAIDDS_VCT_ENTROPY_WEIGHT` | `vct_entropy_weight` + `vct_entropy_weight_set` | no override | apply_config |
| `FLEXAIDDS_N_ELITE` | `n_elite` + `n_elite_set` | no override | apply_config |
| `FLEXAIDDS_SHARING_ALPHA` | `sharing_alpha` optional | no override | apply_config |
| `FLEXAIDDS_BOOM_FRAC` | `boom_frac` optional | no override | apply_config |
| `FLEXAIDDS_FORCE_CF_RANK_EMISSION` | `force_cf_rank_emission` optional | no override | apply_config |
| `FLEXAIDDS_CLASSIC_ENTROPY_RANKING` | `classic_entropy_ranking` optional | no override | apply_config |
| `FLEXAIDDS_ENTROPY_WEIGHT` | `entropy_weight` optional | no override | apply_config |
| `FLEXAIDDS_DIVERSITY_MONITORING` | `diversity_monitoring` optional | no override | apply_config |

**Ranking algorithm is unchanged** — only the source of the override moves from
raw `getenv` to a single ProtocolConfig snapshot.

## Residual getenv inventory (after chunk 3)

**Audit files residual science-knob getenv: 0 in config_parser.**  
**Audit trio (DatasetRunner / gaboom / top) residual: ~11 platform/infra paths.**

### `LIB/DatasetRunner.cpp` (~7)

| Env var | Purpose | Priority |
|---------|---------|----------|
| `HOME` | `~` expansion | leave (platform) |
| `OMP_NUM_THREADS` | Worker OMP sizing | leave (OpenMP contract) |
| `FLEXAIDDS_HTTP_RETRIES` | RCSB download retries | P2 |
| `FLEXAIDDS_BINARY` / `_BUILD` / `_REPO` | Binary resolution | P2 (infra paths) |
| `FLEXAIDDS_PRIORITY_TARGETS` | Target priority list | P2 |

### `LIB/gaboom.cpp` (~1)

| Env var | Purpose | Priority |
|---------|---------|----------|
| `OMP_NUM_THREADS` | Default OMP=2 if unset | leave |

### `LIB/top.cpp` (~3)

| Env var | Purpose | Priority |
|---------|---------|----------|
| `FLEXAIDDS_GRID_CACHE_DIR` | Grid cache directory | P2 (infra) |
| `FLEXAIDDS_LIGAND_BATCH` | Batch ligand directory | P2 |
| `FLEXAIDS_VH_DEBUG` | Vector H-bond debug dump | P3 |

### `LIB/config_parser.cpp` (~0)

All former FLEXAIDDS_* overrides go through ProtocolConfig.

## Migration rules

1. **Do not rewrite DatasetRunner in one PR.** One cluster of related knobs per chunk.
2. Add fields to `ProtocolConfig` first, extend `from_env()` / JSON, then flip call sites.
3. Keep default values byte-identical when env is unset (tests).
4. Prefer reading `protocol_cfg_` (or a passed-in snapshot) over re-calling `getenv` in loops.
5. Leave platform env (`HOME`, `PATH`, `OMP_*`) as raw getenv unless a strong reason exists.
6. Document each migrated var and drop it from Residual.
7. Do not change pose ranking / clustering order without an explicit request + feature flag.

## Tests

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test_protocol_config -j
ctest --test-dir build -R ProtocolConfigTests --output-on-failure
```
