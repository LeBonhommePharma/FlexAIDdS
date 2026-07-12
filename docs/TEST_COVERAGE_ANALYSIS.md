# Test Coverage Analysis Report

**Last actualized**: 2026-07-12  
**Inventory method**: live scan of `tests/`, `python/tests/`, and `CMakeLists.txt` in this session (not conversation memory).

## Summary

| Metric | Value (2026-07-12) |
|--------|--------------------|
| C++ GoogleTest source files (`tests/test_*.cpp`) | **69** |
| C++ `TEST` / `TEST_F` / `TEST_P` cases | **~1,570+** |
| CTest targets registered in `CMakeLists.txt` | **66+ executables** (was 63 before orphan registration) |
| Python pytest modules (`python/tests/` + root `tests/*.py`) | **50+** |
| Python `def test_*` functions | **~1,100+** |
| Combined automated cases | **~2,700+** |
| Code coverage CI | [`.github/workflows/coverage.yml`](../.github/workflows/coverage.yml) (lcov on `LIB/*`) |

**Status intent**: all registered C++ and pure-Python suites must pass before push (`AGENTS.md`).

---

## How coverage is measured

1. **Unit / integration counts** — GoogleTest macros + pytest `test_*` functions (this document).
2. **Line coverage (C++)** — `FLEXAIDS_ENABLE_COVERAGE=ON` + lcov extract of `LIB/*` in the Coverage workflow.
3. **Component mapping** — whether a module has a *dedicated* test file (not whether every line is hit).

Historical claim “4.7% component coverage / 52 C++ tests” is **obsolete** and was removed from decision-making. Prefer the tables below.

---

## Test infrastructure

### C++ (GoogleTest)

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
ctest --test-dir build --output-on-failure
```

- Framework: GoogleTest via CMake `FetchContent`
- Heavy targets link engine sources (`gaboom`, `BindingMode`, NATURaL, Shannon, …) + `tests/stubs.cpp`
- Header-only / light targets use `flexaids_add_unit_test` where possible

### Python (pytest)

```bash
# pure-Python surface (no C++ _core required for most tests)
PYTHONPATH=python python3.11 -m pytest python/tests/ -q --tb=line

# skill / root tests
PYTHONPATH=python python3.11 -m pytest tests/test_flexaid_skill.py tests/test_astex_entropy_sdf_utils.py -q
```

- Marker `@requires_core` skips gracefully when `_core` is absent
- Diagnostics modules may load via `importlib` to avoid hard numpy dependency on import

### Coverage CI

See `.github/workflows/coverage.yml`: Clang Debug + `FLEXAIDS_ENABLE_COVERAGE`, ctest, lcov on `*/LIB/*` (excludes tests, python, benchmarks, `_deps`).

---

## Well-tested components

### C++ suites with large case counts (approximate)

| Suite | Focus |
|-------|--------|
| `test_statmech` | Partition function, log-sum-exp, WHAM/TI edges |
| `test_hardware_dispatch` | ShannonThermoStack backends |
| `test_process_ligand` | Ligand processing pipeline |
| `test_dataset_runner` | C++ DatasetRunner helpers |
| `test_knowledge_pool` | Knowledge pool / transfer |
| `test_gaboom` / `test_ga_*` | GA core, operators, validation, population |
| `test_binding_mode_*` | Clustering, statmech, vibrational, I/O, advanced |
| `test_natural` / nascent chain | NATURaL co-translational paths |
| `test_vcontacts*` / soft wall / matrices | Scoring geometry and matrices |
| `test_mol2_sdf_reader` | MOL2 / SDF I/O |
| `test_cavity_detect` / `test_cleft_cavity` | Site detection |
| `test_buffer_safety` | Bounded config string copies |

### Python packages with dedicated tests

| Module | Test file(s) |
|--------|----------------|
| `thermodynamics`, Kirchhoff helpers | `test_thermodynamics*`, `test_kirchhoff`, `test_py_statmech` |
| `schemas.thermo_audit` | `test_thermo_schema` |
| `diagnostics.ranking_bias` | `test_ranking_bias`, `test_ranking_bias_audit_script` |
| `dataset_runner` | `test_dataset_runner*` |
| `dataset_adapters` | `test_dataset_adapters` |
| `truncate_chain` | `test_truncate_chain` |
| `_fallback_types` | `test_fallback_types` |
| `models` / `results` / `io` | `test_models*`, `test_results*`, `test_io`, `test_pdb_io` |
| `docking`, `encom`, `tencm` | `test_docking`, `test_encom`, `test_tencm*` |
| `ml_rescore`, `optimize`, `figures` | matching `test_*` |
| CLI / import | `test_cli`, `test_import_fallback`, `test_version` |

---

## Strengths

1. **High absolute case count** across C++ and Python.
2. **Core scientific paths** (StatMech, BindingMode thermo, GA, scoring geometry) have multi-file suites.
3. **Security-sensitive parsing** has regression tests (`test_buffer_safety`).
4. **Benchmark orchestration** has dry-run / path / packaging tests.
5. **Audit schema** (entropy.help A1.1) has identity and validation tests.
6. **CI breadth**: multi-platform ctest, pure-Python jobs, coverage, sanitizers, license scan.

---

## Gaps and priorities

### Recently closed (this actualization)

| Gap | Action |
|-----|--------|
| Orphan C++ sources not in CMake | Registered `test_binding_mode_io`, `test_ga_population`, `test_production_blockers` |
| No tests for `dataset_adapters` | Added `python/tests/test_dataset_adapters.py` |
| No tests for `truncate_chain` | Added `python/tests/test_truncate_chain.py` |
| No tests for `_fallback_types` | Added `python/tests/test_fallback_types.py` |
| Thin thermo ledger / buffer / schema / ranking_bias | Expanded cases |

### Remaining high priority

1. **Line coverage trends** — track lcov % on `LIB/` over time; fail CI only after a baseline is published.
2. **PyMOL plugin** — no dedicated automated suite (manual / optional dependency).
3. **GPU backends** — CUDA/Metal shaders need hardware-gated tests beyond CPU dispatch.
4. **Legacy C parsers** — continue buffer-safety expansion beyond GIST keys.
5. **End-to-end docking smoke** — keep lightweight; full Astex stays in benchmark skill, not unit CI.

### Medium priority

- Property-based tests beyond the example module (`hypothesis` optional).
- MPI / distributed backends under CI with stubs.
- Documentation snippet tests (doctest) for public Python API.

---

## Component → test map (core engine)

| Component | Primary tests |
|-----------|----------------|
| `LIB/statmech.*` | `test_statmech`, `test_thermo_ledger` |
| `LIB/BindingMode.*` | `test_binding_mode_*` |
| `LIB/gaboom.*` / GA | `test_gaboom`, `test_ga_*`, `test_production_blockers`, `test_shannon_ga`, `test_entropy_ga` |
| `LIB/Vcontacts.*` | `test_vcontacts`, `test_vcontacts_geometry` |
| `LIB/encom.*` / tENCoM | `test_encom`, `test_tencom_*` |
| ShannonThermoStack | `test_hardware_dispatch`, `test_unified_dispatch` |
| Cavity / cleft | `test_cavity_detect`, `test_cleft_cavity` |
| MOL2/SDF | `test_mol2_sdf_reader` |
| Ring / sugar / chiral | `test_ligand_ring_flex`, `test_sugar_pucker`, `test_ring_conformer_library`, `test_chiral_center` |
| Grand partition | `test_grand_partition`, `test_multi_site_gpf` |
| DiFT | `test_dift` (+ Python `test_dift`) |
| NATURaL | `test_natural`, `test_nascent_chain_scheduler` |

---

## Maintenance rules

1. New `.cpp` under `LIB/` → add/extend a GoogleTest target in `CMakeLists.txt` in the same PR.
2. New Python module under `python/flexaidds/` → add `python/tests/test_<module>.py` (or extend an existing suite).
3. Never leave `tests/test_*.cpp` unregistered in CMake.
4. Before claiming “coverage improved”, re-run inventory counts and at least the affected suites.
5. Keep scientific terminology exact: CF scoring proxy vs thermodynamic ledger (`AGENTS.md`).

---

## Related docs

- [`docs/TESTING.md`](TESTING.md) — how to run suites and interpret markers
- [`docs/VALIDATED_CAPABILITIES.md`](VALIDATED_CAPABILITIES.md) — what Core 1.0 claims require automation
- [`docs/dev/thermo_invariants.md`](dev/thermo_invariants.md) — thermodynamic identities under test
- [`AGENTS.md`](../AGENTS.md) — verification before “done”

## Overall assessment

**Core engine and Python analysis surface: strong.**  
**Release-grade line coverage tracking and plugin/GPU paths: still incomplete.**  
Prefer continuous actualization of this file after large test PRs rather than one-off snapshots.
