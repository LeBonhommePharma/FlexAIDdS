# EXPECTED_API — CMA-ES mock unit tests (chunk 6)

Tests target the **chunk1** public API as landed in
`.swarm/cmaes/chunk1_adapter/artifacts/LIB/cmaes_search.{h,cpp}`
(and mirrored at `LIB/cmaes_search.*` after integration).

## Symbols used by the tests

| Symbol | Header | Role in tests |
|--------|--------|----------------|
| `flexaids_cmaes::CmaesConfig` | `cmaes_search.h` | seed, population, max_evals, enable_entropy_trace, sigma0, archive_size |
| `flexaids_cmaes::CmaesResult` | `cmaes_search.h` | best_cf, best_genes, n_evals, status, archive_* |
| `flexaids_cmaes::EntropyTraceSample` | `cmaes_search.h` | gen, H_search, H_energy, F, best_cf, n_evals |
| `cmaes_run_mock` | free fn | pure mock optimizer, no engine |
| `cmaes_fill_chromosomes` | free fn | snapshot archive → `chromosome[]` + contiguous `gene[]` |
| `cmaes_mock_objective` | free fn | not called directly; exercised via `cmaes_run_mock` |

Using-declarations in the header also promote the types to the global namespace:

```cpp
using flexaids_cmaes::CmaesConfig;
using flexaids_cmaes::CmaesResult;
using flexaids_cmaes::EntropyTraceSample;
```

## Not used (dock path)

| Symbol | Note |
|--------|------|
| `cmaes_run_dock` | Live engine path — out of scope for mock unit tests |
| `set_gene_lim`, `set_bins`, `eval_chromosome`, `get_cf_evalue`, `get_apparent_cf_evalue` | Engine seams — stubbed only for link |

## Acceptance expectations

| Test | Expectation |
|------|-------------|
| `MockSeed12345Converges` | seed=12345, dim=8, max_evals≥5000 → `best_cf < 1e-2` (assert `< 1e-4`) |
| `EntropyTraceNonEmpty` | `enable_entropy_trace=true` → non-empty samples; `H_energy` and `F` finite |
| `SnapshotDims` | `cmaes_fill_chromosomes(..., num_genes=dim, ...)` writes `status='n'`, contiguous genes of length `num_genes`, matching archive |

## FLEXAIDS_CMAES_MOCK_ONLY

Chunk1 as landed does **not** yet guard `cmaes_run_dock` with `FLEXAIDS_CMAES_MOCK_ONLY`.
Chunk6 therefore links `tests/cmaes_mock_seams_stub.cpp` so the mock executable does not
need `flexaid_core`. When chunk1 adds:

```cpp
#ifndef FLEXAIDS_CMAES_MOCK_ONLY
int cmaes_run_dock(...) { ... }
#endif
```

set CMake option `CMAES_USE_MOCK_ONLY=ON` (see `CMakeLists_cmaes.fragment`) and drop the stub TU.

## API name drift

None observed relative to the chunk1 prompt (`cmaes_run_mock`, `cmaes_fill_chromosomes`,
`CmaesConfig`, `CmaesResult`, `EntropyTraceSample`). If a later chunk renames symbols,
update this file and the three `TEST(CmaesSearch, …)` cases together.
