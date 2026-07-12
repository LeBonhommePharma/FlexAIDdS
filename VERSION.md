# FlexAID∆S Version History

**Current Version**: 2.0.1
**Release Date**: 2026-07-12
**Python Package**: 2.0.1 (`python/flexaidds/__version__.py`)
**Repository**: [github.com/LeBonhommePharma/FlexAIDdS](https://github.com/LeBonhommePharma/FlexAIDdS)
**License**: Apache-2.0

---

## Weekly Highlights (2026-05-27 to 2026-06-03)

Selected changes that landed on `master` this week.

### Highlights

- **Thermodynamic reporting stack** — landed the `ThermodynamicBreakdown` ledger and component-wise averages, added `ThermodynamicOutput` / `TotalSampledPartitionFunction`, exposed `temperature_scan`, `fit_delta_Cp`, and affinity converters to Python, and closed the loop with targeted `ctest` + `pytest` verification.
- **Structural interpretation and reporting** — added cleft annotation and flexible-residue selection, safe standard-state affinity calibration utilities, and JSON-driven PyMOL / Markdown / plotting reporting.
- **Benchmark automation on Apple Silicon** — hardened the M3 Pro / iCloud campaign path with iCloud-only support, Metal pre-flight guidance, a master launcher + monitor flow, explicit temperature forwarding, stricter best-mode extraction / validation, and `EntryTaskManager.run()` local thread-pool dispatch.
- **Documentation and web surface** — added the entropy.help foundation, tightened the documentation/validation boundary and source-validator determinism, and published `/flexaid` + `/periodic` with homepage teaser cards and navigation links.
- **Validation datasets** — expanded Astex non-native cross-docking inputs from pair definitions to a 20-pair pilot and the full 85-pair set, plus a null-guard fix in `select_flexible_residues`.

### Key PRs

- **Thermodynamics** — [#212](https://github.com/LeBonhommePharma/FlexAIDdS/pull/212), [#222](https://github.com/LeBonhommePharma/FlexAIDdS/pull/222), [#223](https://github.com/LeBonhommePharma/FlexAIDdS/pull/223), [#224](https://github.com/LeBonhommePharma/FlexAIDdS/pull/224), [#225](https://github.com/LeBonhommePharma/FlexAIDdS/pull/225), [#228](https://github.com/LeBonhommePharma/FlexAIDdS/pull/228)
- **Reporting and preprocessing** — [#215](https://github.com/LeBonhommePharma/FlexAIDdS/pull/215), [#217](https://github.com/LeBonhommePharma/FlexAIDdS/pull/217), [#218](https://github.com/LeBonhommePharma/FlexAIDdS/pull/218)
- **Docs and validation** — [#220](https://github.com/LeBonhommePharma/FlexAIDdS/pull/220), [#221](https://github.com/LeBonhommePharma/FlexAIDdS/pull/221), [#229](https://github.com/LeBonhommePharma/FlexAIDdS/pull/229)
- **Benchmark operations** — [#230](https://github.com/LeBonhommePharma/FlexAIDdS/pull/230)
- **Tooling and site surface** — [#231](https://github.com/LeBonhommePharma/FlexAIDdS/pull/231), [#232](https://github.com/LeBonhommePharma/FlexAIDdS/pull/232), [#233](https://github.com/LeBonhommePharma/FlexAIDdS/pull/233)

---

## v2.0.1 (2026-07-12) — Homebrew + packaging hotfix

Production docking matrix and installer fixes for out-of-the-box Homebrew / PyPI installs.

### Critical
- **Ship production VCT matrix** `MC_st0r5.2_6.dat` (md5 `9dc93717…`) at repo root **and** `WRK/` so release tarballs and Homebrew stable builds type atoms correctly (v2.0.0 only had an outdated matrix under `WRK/`, causing mass "unknown atom type" failures).
- **Homebrew formula**: libexec layout + PATH wrappers setting `FLEXAIDDS_DATA_DIR`; Metal off by default; data co-located with binaries.
- **Engine path resolution** (`LIB/top.cpp`): `realpath(argv[0])`, honor `FLEXAIDDS_DATA_DIR`, fall back to `share/flexaidds` / `WRK`.
- **PyPI metadata**: static PEP 621 version so clean builds no longer emit `UNKNOWN-0.0.0` wheels.

### Python
- Package version **2.0.1** (`python/flexaidds/__version__.py`, `python/pyproject.toml`).

## v2.0.0 (2026-04-04) — Stable Release

First stable release of FlexAID∆S, the entropy-driven molecular docking engine. This is a ground-up rewrite of FlexAID combining genetic algorithms with statistical mechanics thermodynamics for accurate binding free energy prediction. 655 commits ahead of the v1.5 legacy tag. All development phases complete.

### Post-release fixes (2026-04-14)

- **GrandPartitionFunction audit** — 17 code quality fixes: `free_energy()` renamed to `F_bound()`, `[[nodiscard]]` on all queries, `=delete` on copy/move, explicit `log_c` storage, `scoped_lock` for thread safety, concentration guard (> 1000 M rejected), `log_intrinsic_selectivity()` and `all_log_zZ()` API additions
- **Build fixes** — resolved merge conflict in `gaboom.cpp`, restored 8 missing GA constants in `ga_constants.h`, renamed `OVERFLOW` → `ERR_OVERFLOW` (macOS `math.h` conflict), fixed `Backend`/`HardwareBackend` type mismatch, added Eigen3 linkage for `test_encom` and `test_ion_handling`
- **Test fixes** — corrected `ExtremeEnergySpreadLogsumexpStable` threshold, fixed `VeryLargeEigenvaluesStiff` physics assumption, corrected `DeltaGRelativeToAnotherModeConsistent` sign convention, rewrote `hash_genes` tests
- **CI hardening** — all 11 workflow files pinned to commit SHAs with version comments
- **Test suite**: 48 C++ tests, all passing (was 46 with 3 pre-existing failures)

### Highlights

- **Shannon entropy (∆S) scoring** — the core innovation: configurational entropy as a first-class scoring term via StatMechEngine
- **Full cross-platform support** — Linux (GCC/Clang), macOS (Clang + Metal), Windows (MSVC)
- **Python bindings** — complete `flexaidds` package with pybind11 C++ bridge and pure-Python fallback
- **GPU acceleration** — CUDA and Metal compute shaders for batch evaluation, entropy histograms, and cavity detection
- **MPI distributed docking** — grid domain decomposition with parallel transport
- **GIST water-displacement scoring** — grid-based explicit solvation via `GISTEvaluator`
- **Directional H-bond scoring** — geometry-aware hydrogen bond potential via `HBondEvaluator`
- **DatasetRunner benchmarking** — automated distributed benchmarking system for docking campaigns
- **80 test targets** — 48 C++ (GoogleTest) + 32 Python (pytest)
- **4 CI workflows** — build matrix, license scanning, performance regression, sanitizers

### Core Thermodynamics

- **StatMechEngine** — canonical ensemble partition function, Helmholtz free energy, Shannon entropy, heat capacity
- **BindingMode** — pose clustering with intra-mode thermodynamics and Boltzmann weights
- **BindingPopulation** — multi-mode free energy aggregation via log-sum-exp
- **Log-sum-exp stability** — numerically stable for energy differences > 100 kT
- **WHAM** — single-window weighted histogram free energy profiles
- **Thermodynamic integration** — lambda-path free energy perturbation
- **Replica exchange** — Metropolis swap criterion for parallel tempering

### Python Bindings

- **pybind11 core** — `StatMechEngine`, `ENCoMEngine`, `Thermodynamics`, `VibrationalEntropy` exposed to Python
- **Data models** — `PoseResult`, `BindingModeResult`, `DockingResult` with JSON/dict round-trip serialization
- **Pure-Python fallback** — full `flexaidds` package works without C++ compilation
- **CLI inspector** — `python -m flexaidds <dir> [--json|--csv|--top N]`
- **GIL-release bindings** — thread-safe parallel C++ calls from Python
- **22 Python modules** including benchmark, energy matrix training, visualization, dataset adapters

### PyMOL/NRGSuite GUI

- **PyMOL plugin** — 8 modules: entropy heatmap, mode animation, ITC comparison, interactive docking
- **14 PyMOL commands** — load, visualize, color, animate binding modes
- **FlexMolView** — experimental multi-language molecule viewer (Python, Swift, TypeScript)

### Hardware Acceleration

- **Unified dispatch** — CUDA > ROCm/HIP > Metal > AVX-512 > AVX2 > OpenMP > scalar
- **CUDA kernels** — batch CF evaluation, Shannon histograms, tENCoM Hessian, FastOPTICS k-NN
- **Metal shaders** — Shannon entropy, cavity detection, CF evaluation (macOS/Apple Silicon)
- **SIMD** — AVX-512 and AVX2 vectorised geometric primitives
- **OpenMP + Eigen3** — thread parallelism and vectorised linear algebra
- **LTO binaries** — link-time optimized `FlexAIDdS` and `tENCoM` executables

### Scoring & Solvation

- **GIST water-displacement scoring** — grid-based explicit solvation term (`GISTEvaluator`, `GISTGrid`)
- **Directional H-bond scoring** — geometry-aware hydrogen bond potential (`HBondEvaluator`, `hbond_potential.h`)
- **256x256 energy matrix** — extended atom typing with continuous training pipeline
- **Metal ion scoring** — crystallographic VdW radii for 20+ ion types in Voronoi CF
- **Structural water** — ordered waters (B < 20 A^2) participate in contact scoring

### Genetic Algorithm Engine

- **GA diversity** — population entropy monitoring and adaptive diversity pressure (`ga_diversity.h`)
- **GAContext** — structured GA run context for reproducible experiments
- **SMFREE fitness** — Shannon entropy as direct GA fitness component with adaptive pressure
- **Molecular interaction fields** — MIF-guided docking, reference ligand seeding, Spectrophore descriptors
- **Memory safety fixes** — roulette OOB, crossover/mutate UB, adapt_prob division-by-zero guards
- **GA re-entrancy** — thread-safe GA for parallel pipeline routing

### Parallel & Distributed

- **MPI transport** — grid domain decomposition with parallel docking
- **DistributedBackend / ThreadBackend / GPUContextPool** — unified parallel execution layer
- **AtomSoA / VoronoiCFBatch_SoA** — Structure-of-Arrays layouts for SIMD-friendly batching
- **DatasetRunner** — automated distributed benchmarking system (`dataset_runner.py`)
- **Bonhomme Fleet** — distributed docking across Apple devices via iCloud

### Structural Modules

- **tENCoM** — torsional elastic network model for backbone vibrational entropy
- **ShannonThermoStack** — combined configurational + vibrational entropy with HW acceleration
- **LigandRingFlex** — non-aromatic ring conformer sampling (chair/boat/twist, sugar pucker)
- **ChiralCenter** — explicit R/S stereocenter discrimination with energy penalty
- **CavityDetect** — SURFNET gap-sphere cavity detection with Metal GPU support
- **NATURaL** — co-translational/co-transcriptional assembly with NucleationDetector
- **CleftDetector** — binding site identification with ion/water awareness
- **ENCoM** — elastic network contact model for vibrational entropy (∆S_vib)

### Input & Configuration

- **Multi-format input** — PDB, CIF/mmCIF, MOL2, SDF, SMILES (auto 3D build)
- **JSON config** — single config file with sensible defaults for all parameters

### CI/CD & Quality

- **4 CI workflows** — `ci.yml` (multi-platform build matrix), `license-scan.yml` (GPL/AGPL detection), `perf.yml` (benchmark regression), `sanitizers.yml` (ASan/UBSan)
- **Multi-platform CI** — Linux GCC, Linux Clang, macOS Clang, Windows MSVC
- **Benchmark harnesses** — CASF-2016, CrossDock, LIT-PCBA (`tests/benchmarks/`)
- **Smoke validation bundle** — manifest and run script (`benchmarks/smoke/`)

### Experimental

- **Swift package** — macOS/iOS actors wrapping StatMechEngine and ENCoM
- **TypeScript SDK** — PWA dashboard with Mol* 3D viewer
- **TurboQuant** — vector quantization for contact and energy vector compression
- **ml_rescore.py** — machine-learning rescoring module
- **optimize.py** — docking parameter optimization

### Codebase Statistics

- 655 commits since v1.5 tag
- 40K+ lines C++, 19K+ lines Python, 10K lines Swift, 4K lines TypeScript, 12K+ lines tests
- 48 C++ test files (GoogleTest), 32 Python test files (pytest) — 80 total
- CI: Linux GCC, Linux Clang, macOS Clang, Windows MSVC
- 4 CI workflows: `ci.yml`, `license-scan.yml`, `perf.yml`, `sanitizers.yml`

---

## v1.9.0-beta (2026-04-04) — Pre-release

Feature-complete pre-release. See [v2.0.0](#v200-2026-04-04--stable-release) for the full changelog — v1.9.0-beta was the final beta before this stable release.

---

## v1.5 (2026-03-06) — Legacy tag

Original FlexAID release with Voronoi contact function scoring. Tagged as `v1.5` in git history.

- Genetic algorithm docking engine
- Voronoi contact function (VCT) and sphere approximation (SPH)
- Dead-end elimination for torsion pruning
- Legacy INP/GA config file format
- No thermodynamics, no entropy, no hardware acceleration

---

## Versioning Policy

FlexAID∆S follows [Semantic Versioning](https://semver.org/):

- **MAJOR**: Breaking API changes
- **MINOR**: New features, backward-compatible
- **PATCH**: Bug fixes, backward-compatible

Pre-release tags: `-alpha` (unstable), `-beta` (feature-complete, testing), `-rc` (release candidate)
