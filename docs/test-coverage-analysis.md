# Test Coverage Analysis — FlexAIDdS

**Date:** 2026-03-11
**Scope:** C++ (GoogleTest) and Python (pytest) test suites

## Overall Numbers

| Area | Test Files | Test Cases | Estimated Coverage |
|------|-----------|------------|-------------------|
| **C++ (GoogleTest)** | 5 | ~145 | ~35% of LIB/ by lines |
| **Python (pytest)** | 16 | ~321 | ~70% of flexaidds/ |

---

## What's Well-Tested

The **thermodynamics stack** has excellent coverage across both languages:

- **StatMechEngine** (C++ & Python) — partition functions, free energy, entropy, heat capacity, Boltzmann weights, WHAM, thermodynamic integration — ~95% coverage with edge cases and numerical stability tests
- **BindingMode integration** — lazy caching, Boltzmann weighting, delta-G — ~75%
- **ShannonThermoStack** — continuous/discrete entropy, hardware dispatch, torsional entropy — ~85%
- **tENCoM differential engine** — PDB reading, mode comparison, nucleic acid support — ~85%
- **Python I/O & models** — PDB REMARK parsing, data classes, result loading — excellent coverage

---

## Critical Gaps (Priority 1)

### 1. Genetic Algorithm — `gaboom.cpp` (1,700+ lines, NO TESTS)

The core evolutionary search engine has zero direct tests. Selection, crossover, mutation, population initialization, convergence behavior, and adaptive operators are all untested. Only scoring output properties are validated indirectly in `test_ga_validation.cpp`.

**Recommendation:** Add `test_gaboom.cpp` testing:
- Chromosome initialization and gene encoding
- Selection, crossover, mutation operator correctness
- Deterministic convergence on a simple fitness landscape
- Population diversity metrics
- Adaptive operator adjustment

### 2. Voronoi Contact Scoring — `Vcontacts.cpp` (523 lines, minimal indirect testing)

The scoring function that drives all fitness evaluation has no direct unit tests. Voronoi polyhedron computation, contact area calculation, and solvent-accessible surface calculations are untested.

**Recommendation:** Add `test_vcontacts.cpp` with:
- Known-geometry test cases (e.g., two spheres at known distance)
- Contact area validation against analytical values
- Engulfing atom correction verification
- SAS calculation benchmarks

### 3. Cavity & Cleft Detection — `CavityDetect/`, `CleftDetector.cpp` (400+ lines, NO TESTS)

Binding site identification is completely untested. Incorrect cavity detection means the ligand searches the wrong region entirely.

**Recommendation:** Add tests with a known protein structure where the binding site is well-characterized.

---

## Major Gaps (Priority 2)

### 4. File I/O & Input Parsing (~2,000+ lines C++, NO TESTS)

- `Mol2Reader.cpp`, `SdfReader.cpp` — ligand file parsing untested
- `fileio.cpp`, `read_*.cpp`, `write_*.cpp` — grid, constraint, and parameter I/O untested
- `assign_*.cpp` — atom typing and radii assignment untested (errors here cascade everywhere)

**Recommendation:** Add round-trip I/O tests and malformed-input tests for each parser.

### 5. Ring Flexibility — `LigandRingFlex/` (400+ lines, NO TESTS)

Non-aromatic ring conformer sampling (critical for sugars and saturated rings) has no tests.

**Recommendation:** Add tests for ring conformer enumeration and sugar pucker classification.

### 6. Chirality Discrimination — `ChiralCenter/` (150+ lines, NO TESTS)

R/S stereocenter handling is untested. Wrong chirality = wrong molecule.

**Recommendation:** Add tests with known R and S enantiomers verifying correct classification.

### 7. Python `visualization.py` (90 lines, 0 tests)

PyMOL integration module has zero test coverage. No validation of graceful degradation when PyMOL is unavailable.

**Recommendation:** Add tests mocking PyMOL imports to validate load functions and fallback behavior.

### 8. ENCoM Full Pipeline (~200 lines C++, 60% coverage)

`load_modes()` file parsing, eigenvalue cutoff filtering, and error handling for malformed mode files are untested. This is the foundation for Phase 3.

**Recommendation:** Add mode file parsing tests with valid and malformed inputs before Phase 3 begins.

---

## Moderate Gaps (Priority 3)

| Area | Issue |
|------|-------|
| **tENCoM Hessian construction** (`tencm.cpp`) | Hessian matrix building and Jacobi eigendecomposition untested (~65% coverage) |
| **NATURaL module** (`NATURaL/`, 600+ lines) | Core planning, scheduler, DualAssembly synthetic backend, Sim C gating, and fail-closed `--real-ga` paths covered by C++ tests; real GA backend wiring and calibration remain unvalidated |
| **OpenMP concurrency** | Race conditions and thread-safety of global RNG state not tested |
| **GPU fallback paths** | CUDA/Metal unavailable → scalar fallback not validated |
| **CLI --csv and --top flags** | Python CLI output formats partially untested |
| **C++ binding integration** | `docking.py` BindingMode C++ paths barely tested |
| **CI breadth** | Only `test_statmech_smoke.py` runs with C++ bindings in CI; full Python+C++ test suite not in CI |

---

## C++ Test File Summary

| Test File | Tests | Focus | Coverage |
|-----------|-------|-------|----------|
| `test_statmech.cpp` | ~60 | StatMechEngine: partition functions, thermodynamics, WHAM, TI | 95% |
| `test_binding_mode_statmech.cpp` | ~15 | BindingMode ↔ StatMechEngine integration | 75% |
| `test_tencom_diff.cpp` | ~20 | PDB reader, tENCoM differential engine, nucleic acids | 85% |
| `test_hardware_dispatch.cpp` | ~40 | Shannon entropy, ShannonThermoStack, hardware dispatch | 85% |
| `test_ga_validation.cpp` | ~5 | GA scoring validation, batch/serial agreement | 50% |

## C++ Source Files Without Any Tests

| File | Lines | Purpose |
|------|-------|---------|
| `gaboom.cpp` | 1,700+ | Genetic algorithm core |
| `Vcontacts.cpp` | 523 | Voronoi contact scoring |
| `CleftDetector.cpp` | ~400 | Binding site detection |
| `CavityDetect/` | ~200 | Cavity detection (Metal GPU) |
| `ChiralCenter/` | ~150 | Stereocenter discrimination |
| `LigandRingFlex/` | ~400 | Ring conformation sampling |
| `NATURaL/` | ~600 | Co-translational assembly; DualAssembly C++ coverage now exists, real GA backend still uncovered |
| `Mol2Reader.cpp` | ~150 | MOL2 file parsing |
| `SdfReader.cpp` | ~150 | SDF file parsing |
| `fileio.cpp` | ~200 | File I/O operations |
| `assign_*.cpp` | ~600 | Atom typing, radii, constraints |
| `read_*.cpp` | ~1,000 | Input parsing |
| `write_*.cpp` | ~700 | Output formatting |

---

## Python Test File Summary

| Test File | Tests | Focus | C++ Required? |
|-----------|-------|-------|---------------|
| `test_docking.py` | 46 | Docking config, Pose/BindingMode/Population | Partial |
| `test_io.py` | 44 | PDB REMARK parsing, key normalization | No |
| `test_models.py` | 43 | PoseResult, BindingModeResult, DockingResult | No |
| `test_py_statmech.py` | 35 | Pure-Python StatMechEngine | No |
| `test_statmech.py` | 28 | C++ StatMechEngine correctness | Yes |
| `test_encom.py` | 27 | ENCoM vibrational entropy | Partial |
| `test_thermodynamics.py` | 21 | Thermodynamics dataclass + integration | Partial |
| `test_cli.py` | 18 | CLI argument parsing, output formats | No |
| `test_results.py` | 17 | Result loading, grouping, metadata | No |
| `test_thermodynamics_dataclass.py` | 8 | Thermodynamics serialization | No |
| `test_version.py` | 8 | Package metadata | No |
| `test_results_loader_models.py` | 3 | Nested directory loading | No |
| `test_results_io.py` | 2 | Result file parsing integration | No |
| `test_import_fallback.py` | 2 | Graceful degradation without C++ | No |
| `test_statmech_smoke.py` | 1 | CI smoke test | Yes |

## Python Modules Without Tests

| Module | Lines | Purpose |
|--------|-------|---------|
| `visualization.py` | 90 | PyMOL integration (0 tests) |

---

## Structural Recommendations

1. **Add CI job for full Python+C++ test suite** — Currently only a 1-test smoke test runs with bindings. The ~68 C++-dependent Python tests should run in CI.

2. **Add integration/regression tests** — All current tests use synthetic data. Add at least one end-to-end test with a real receptor-ligand pair validating the full pipeline (load → dock → score → cluster → thermodynamics).

3. **Add error-path and edge-case tests broadly** — Malformed inputs, NaN energies, empty collections, and negative values are undertested across both C++ and Python.

4. **Add determinism tests** — Cross-platform reproducibility tests with fixed RNG seeds to catch silent numerical drift.

5. **Expand `test_ga_validation.cpp`** — The existing file is a natural home for GA algorithm tests (selection, crossover, mutation) rather than creating a new file.
