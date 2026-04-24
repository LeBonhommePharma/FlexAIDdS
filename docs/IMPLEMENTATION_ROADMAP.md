# FlexAIDΔS Implementation Roadmap

**Complete 5-Phase Development Plan**  
**Status as of March 23, 2026**: Phases 1–3 ✅ Complete | Phase 4 🚧 Active | Phase 5 ✅ Complete  
**Repository**: [lmorency/FlexAIDdS](https://github.com/LeBonhommePharma/FlexAIDdS) · [NRGlab/FlexAIDdS](https://github.com/NRGlab/FlexAIDdS)  
**Branch**: `master`  
**Codebase**: 508 commits · 40K lines C++ · 19K lines Python · 10K lines Swift · 4K lines TypeScript · 12K lines tests

---

## Executive Summary

FlexAIDΔS bridges the **30-year entropy gap** in molecular docking by computing true thermodynamic free energies (ΔG = ΔH − TΔS) via Shannon information theory and statistical mechanics. This roadmap documents the complete implementation strategy across 5 development phases:

1. **Phase 1**: Core Thermodynamics ✅ **COMPLETE** (March 2026)
2. **Phase 2**: Python Bindings ✅ **COMPLETE** (March 2026)
3. **Phase 3**: NRGSuite/PyMOL GUI ✅ **COMPLETE** (March 2026)
4. **Phase 4**: Grid Optimization & Docking Intelligence 🚧 **ACTIVE** (March–May 2026)
5. **Phase 5**: Hardware Acceleration ✅ **COMPLETE** (March 2026)

**Shannon's Energy Collapse Metric**: Every optimization targets minimal information entropy waste in the computational pipeline.

---

## Phase 1: Core Thermodynamics ✅ COMPLETE

**Timeline**: January–March 2026  
**Status**: ✅ **100% Complete** (March 9, 2026)  
**Commits**: 47 commits, 8,432 lines added

### Objectives

1. **Statistical Mechanics Engine** → Canonical partition function, Helmholtz free energy, Shannon entropy
2. **Binding Mode Abstraction** → Ensemble clustering, intra-mode thermodynamics
3. **Global Population** → Multi-mode free energy aggregation
4. **Unit Test Suite** → Comprehensive validation (GoogleTest)

### Deliverables

#### 1.1 Statistical Mechanics Core (`LIB/statmech.{h,cpp}`)

```cpp
class StatMechEngine {
public:
    Thermodynamics get_thermodynamics() const;
    double get_partition_function() const;  // Z
    double get_free_energy() const;         // F = -kT ln(Z)
    double get_mean_energy() const;         // ⟨E⟩
    double get_entropy() const;             // S = -k Σ p_i ln(p_i)
    double get_heat_capacity() const;       // C_v
    
    double relative_free_energy(const StatMechEngine& other) const;
    std::vector<double> wham_profile(const std::vector<double>& coords);
    double thermodynamic_integration(const std::vector<double>& lambda_path);
    bool replica_swap_accept(const StatMechEngine& other, double T_self, double T_other);
};
```

**Key Features**:
- Log-sum-exp stability — prevents overflow for large energy differences
- Lazy evaluation — partition function computed once, cached
- WHAM integration — single-window free energy profiles
- TI support — thermodynamic integration over λ reaction coordinate
- OpenMP/Eigen acceleration for matrix operations

**Validation**: ✅ Exact agreement with analytical results (harmonic oscillator), numerical stability for ΔE > 100 kT, entropy bounds verified

#### 1.2 Binding Mode Abstraction (`LIB/BindingMode.{h,cpp}`)

```cpp
class BindingMode {
public:
    Thermodynamics get_thermodynamics();
    double get_free_energy();
    std::vector<double> get_boltzmann_weights();
    double get_shannon_entropy();
    std::vector<std::vector<double>> get_deltaG_matrix();
    StatMechEngine& get_stat_mech();
};
```

**Key Features**:
- Dual API — legacy (`compute_energy/enthalpy/entropy`) + modern (`get_thermodynamics`)
- Shannon entropy per binding mode
- ΔG interaction matrices for inter-residue analysis
- FastOPTICS/DBSCAN clustering integration
- VibrationalEntropy via tENCoM normal modes

#### 1.3 Global Population (`LIB/BindingMode.h` — BindingPopulation)

```cpp
class BindingPopulation {
public:
    Thermodynamics get_global_thermodynamics();
    double get_selectivity();
    std::vector<BindingMode*> get_modes();
    double get_shannon_entropy();
    std::vector<std::vector<double>> get_deltaG_matrix();
};
```

- Multi-mode free energy via log-sum-exp aggregation
- Selectivity metric — dominant mode probability vs entropy
- Population-level Shannon entropy and ΔG matrices

#### 1.4 Metal Cavity Detection (`LIB/CavityDetect/`)

- `CavityDetect.h` — grid-based solvent-accessible cavity identification
- `SpatialGrid.h` — fast neighbour lookup for cavity analysis
- `CavityDetectMetalBridge.h` — Metal GPU bridge (macOS)

#### 1.5 Unit Test Suite

32 C++ test files, 30 Python test files — **27/27 C++ tests passing** (verified March 23, 2026):

| Test Target | Coverage |
|:------------|:---------|
| `test_statmech` | StatMechEngine, partition function, entropy bounds |
| `test_binding_mode_statmech` | BindingMode thermodynamics, Boltzmann weights |
| `test_binding_mode_advanced` | Multi-mode population aggregation |
| `test_binding_mode_vibrational` | tENCoM vibrational entropy integration |
| `test_gaboom` | GA core: QuickSort, fitness_stats, selection |
| `test_ga_core` | Genetic algorithm operators |
| `test_ga_validation` | GA convergence, BatchResult validation |
| `test_entropy_ga` | SMFREE fitness model, Boltzmann GA blending |
| `test_vcontacts` | Voronoi contact surface computation |
| `test_encom` | Elastic Network Contact Model |
| `test_tencom_diff` | tENCoM differential entropy |
| `test_tencom_entropy_diff` | Vibrational entropy differentials |
| `test_fast_optics` | FastOPTICS clustering |
| `test_statmech` | Statistical mechanics core |
| `test_json_config` | JSON parser, config defaults |
| `test_hardware_dispatch` | ShannonThermoStack dispatch |
| `test_hardware_detect_dispatch` | Hardware detection + dispatch |
| `test_unified_dispatch` | AVX-512 geometric primitives |
| `test_ptm_attachment` | Post-translational modification |
| `test_cleft_cavity` | Cleft/cavity detection |
| `test_cavity_detect` | Grid-based cavity analysis |
| `test_ion_handling` | Ion typing, radii assignment |
| `test_ring_conformer_library` | Ligand ring flexibility |
| `test_sugar_pucker` | Sugar ring conformations |
| `test_chiral_center` | Chiral center detection |
| `test_cube_grid` | Cube grid decomposition |
| `test_parallel_dock` | Distributed docking infrastructure |
| `test_binding_residues` | MIF-based binding site ID |
| `test_mif_grid` | Molecular interaction field |
| `test_reflig_spectrophore` | Reference ligand + Spectrophore |
| `test_mol2_sdf_reader` | Mol2/SDF file parsing |
| `test_soft_contact_matrix` | 256×256 energy matrix I/O |

### Phase 1 Lessons Learned

- Lazy evaluation + log-sum-exp = no numerical surprises at any temperature
- Dual API (legacy + modern) prevented breaking existing GA pipeline during migration
- Shannon entropy as optimization metric — not just output — collapses fitness landscape dimensionality

---

## Phase 2: Python Bindings ✅ COMPLETE

**Timeline**: February–March 2026  
**Status**: ✅ **100% Complete** (March 23, 2026)  
**Module count**: 22 Python modules in `python/flexaidds/`

### Objectives

1. **pybind11 Core Bindings** → Full C++ API exposed to Python
2. **Data Model Round-Trip** → JSON/dict serialization for all result types
3. **Analysis Tools** → Benchmarking, visualization, energy matrix training
4. **CLI Interface** → Command-line tools for batch operations

### Deliverables

#### 2.1 Core Bindings (`python/bindings/core_bindings.cpp`)

pybind11 bindings exposing:
- `StatMechEngine`, `Thermodynamics`, `VibrationalEntropy`
- `BindingMode`, `BindingPopulation`
- `PoseResult`, `BindingModeResult`, `DockingResult`
- Shannon entropy, partition functions, Boltzmann weights

#### 2.2 Python Data Models (`python/flexaidds/models.py`)

```python
@dataclass
class DockingResult:
    modes: List[BindingModeResult]
    thermodynamics: Thermodynamics
    
    def to_dict(self) -> dict: ...
    @classmethod
    def from_dict(cls, d: dict) -> 'DockingResult': ...
    @classmethod
    def from_json(cls, path: str) -> 'DockingResult': ...
    def __repr__(self) -> str: ...  # Concise REPL representation
```

All data model classes (`PoseResult`, `BindingModeResult`, `DockingResult`, `Thermodynamics`) support:
- `__repr__` for debugging/REPL experience
- `from_dict()` / `to_dict()` round-trip serialization
- `from_json()` for file-based deserialization

#### 2.3 Analysis & Benchmarking

| Module | Purpose |
|:-------|:--------|
| `benchmark.py` | Comparative benchmark: FlexAIDΔS vs Boltz-2 |
| `boltz2.py` | Boltz-2 NIM client for structure/affinity prediction |
| `energy_matrix.py` | 256×256 soft contact energy matrix operations |
| `energy_matrix_cli.py` | CLI for energy matrix training/evaluation |
| `continuous_training.py` | Continuous training pipeline for energy matrices |
| `dataset_adapters.py` | Dataset adapters (PDBbind, CASF, custom) |
| `train_256x256.py` | 256×256 matrix training driver |
| `visualization.py` | Binding mode visualization, Mol* viewer integration |
| `supercluster.py` | SuperCluster hierarchical clustering |
| `tencom_results.py` | FlexModeResult, FlexPopulationResult parsing |
| `tencm.py` | tENCoM Python interface |
| `encom.py` | ENCoM elastic network interface |

#### 2.4 Package Structure

```
python/flexaidds/
├── __init__.py              # Unified exports, fallback types
├── __main__.py              # CLI entry point
├── __version__.py           # Version tracking
├── _fallback_types.py       # Pure-Python fallbacks (no C++)
├── models.py                # Data model classes
├── docking.py               # Docking orchestration
├── io.py                    # File I/O utilities
├── results.py               # Result parsing
├── thermodynamics.py        # Thermodynamic calculations
├── visualization.py         # Mol* viewer, plots
├── benchmark.py             # Comparative benchmarking
├── boltz2.py                # Boltz-2 NIM client
├── energy_matrix.py         # Energy matrix operations
├── energy_matrix_cli.py     # Energy matrix CLI
├── continuous_training.py   # Continuous training pipeline
├── dataset_adapters.py      # Dataset adapters
├── train_256x256.py         # 256×256 training driver
├── supercluster.py          # SuperCluster
├── tencom_results.py        # tENCoM result parsing
├── tencm.py                 # tENCoM interface
├── encom.py                 # ENCoM interface
└── updater.py               # Self-update utility
```

---

## Phase 3: NRGSuite/PyMOL GUI ✅ COMPLETE

**Timeline**: February–March 2026  
**Status**: ✅ **100% Complete** (March 23, 2026)

### Objectives

1. **PyMOL Plugin Hardening** → Bug fixes, performance, async operations
2. **Entropy Visualization** → Heatmaps, animation, ITC-style plots
3. **Interactive Docking** → Real-time docking from PyMOL

### Deliverables

#### 3.1 PyMOL Plugin (`pymol_plugin/`)

| Module | Features |
|:-------|:---------|
| `__init__.py` | Plugin registration, initialization |
| `gui.py` | Main GUI panel with mode/population controls |
| `entropy_heatmap.py` | Shannon entropy projected onto structure surface |
| `mode_animation.py` | Binding mode transition animation |
| `itc_comparison.py` | ITC-style thermogram comparison plots |
| `interactive_docking.py` | Real-time docking with pose streaming |
| `results_adapter.py` | Bridge between C++ results and PyMOL objects |
| `visualization.py` | CGO renderer, Kabsch alignment, NumPy acceleration |

**Hardening** (merged `check-pymol-nrgsuite`):
- 3 critical bug fixes, 3 high-priority, 5 medium
- NumPy bulk conversion for atom coordinates
- Async docking to prevent GUI freezing
- Multi-ligand ITC comparison support

---

## Phase 4: Grid Optimization & Docking Intelligence 🚧 ACTIVE

**Timeline**: March–May 2026  
**Status**: 🚧 **~75% Complete** (March 23, 2026)

### Objectives

1. **Parallel Grid-Decomposed Docking** → Domain decomposition for distributed docking
2. **Molecular Interaction Fields** → MIF-guided docking optimization
3. **Entropy-Driven GA** → Shannon entropy in genetic algorithm fitness
4. **Voronoi Hydration** → Explicit solvation entropy (planned)

### Deliverables

#### 4.1 Parallel Docking Infrastructure ✅ COMPLETE

| Component | File | Purpose |
|:----------|:-----|:--------|
| `GridDecomposer` | `LIB/GridDecomposer.{h,cpp}` | Domain decomposition of docking grid into sub-volumes |
| `ParallelDock` | `LIB/ParallelDock.{h,cpp}` | Distributed docking coordinator with partition function aggregation |
| `SharedPosePool` | `LIB/SharedPosePool.{h,cpp}` | Lock-free shared pose pool across parallel workers |
| `MPITransport` | `LIB/MPITransport.{h,cpp}` | MPI communication layer for distributed execution |

#### 4.2 Molecular Interaction Fields ✅ COMPLETE

| Component | File | Purpose |
|:----------|:-----|:--------|
| `MIFGrid` | `LIB/MIFGrid.h` | Molecular interaction field grid computation |
| `RefLigSeed` | `LIB/RefLigSeed.h` | Reference ligand pose seeding from co-crystal data |
| `Spectrophore` | `LIB/Spectrophore.h` | Spectrophore molecular descriptors (pybind11 exposed) |
| `BindingResidues` | `LIB/BindingResidues.h` | Auto-identification of key binding residues from MIF scores |

Wired into GA pipeline: MIF/RefLig/GridPrio integrated into `gaboom.cpp`, auto-flex key binding residues by default.

#### 4.3 SMFREE Entropy-Driven GA Fitness ✅ COMPLETE

```cpp
// SMFREE: Statistical Mechanics FREE energy fitness model
// Boltzmann weight blending in GA selection — Shannon entropy as 
// direct fitness component, not just post-hoc analysis
struct SMFREEConfig {
    double entropy_weight;      // Blend factor for S_Shannon in fitness
    double boltzmann_beta;      // Inverse temperature for selection pressure
    bool   adaptive_pressure;   // Auto-tune β during GA evolution
};
```

- New chromosome fields for entropy tracking
- Boltzmann weight blending in tournament selection
- Config parser integration (`config_defaults.h`, `config_parser.cpp`)

#### 4.4 GA Bug Fixes ✅ COMPLETE

- QuickSort descending order fix (was producing incorrect rankings)
- `fitness_stats` loop boundary fix
- Swift FXGA.mm cleanup (better variable names, removed unnecessary `const_cast`)

#### 4.5 Voronoi Hydration 🔜 PLANNED

- Voronoi tessellation of solvent shell (`LIB/Vcontacts.{h,cpp}` — foundation exists)
- Interface water detection and classification
- Empirical calibration against ITC data
- Integration into BindingMode entropy calculation

#### 4.6 256×256 Energy Matrix System ✅ COMPLETE

| Component | File | Purpose |
|:----------|:-----|:--------|
| `soft_contact_matrix.h` | `LIB/soft_contact_matrix.h` | 256-type soft contact matrix storage/lookup |
| `atom_typing_256.h` | `LIB/atom_typing_256.h` | Extended atom typing scheme (256 types) |
| `shannon_matrix_scorer.h` | `LIB/shannon_matrix_scorer.h` | Shannon-informed matrix scoring |
| Continuous training | `python/flexaidds/continuous_training.py` | Online training pipeline |
| Dataset adapters | `python/flexaidds/dataset_adapters.py` | PDBbind/CASF/custom dataset loading |

---

## Phase 5: Hardware Acceleration ✅ COMPLETE

**Timeline**: February–March 2026  
**Status**: ✅ **100% Complete** (March 23, 2026)

### Objectives

1. **Multi-backend GPU** → CUDA + Metal acceleration
2. **SIMD Vectorization** → AVX2/AVX-512 for CPU-bound kernels
3. **Thread Parallelism** → OpenMP + Eigen3 integration
4. **Unified Dispatch** → Runtime backend selection
5. **Cross-platform Build** → Linux/macOS/Windows

### Deliverables

#### 5.1 Hardware Detection (`LIB/hardware_detect.{h,cpp}`)

Runtime detection of available acceleration backends:
- CPU feature flags (AVX2, AVX-512, SSE4.2)
- CUDA device enumeration and capability query
- Metal device detection (macOS)
- OpenMP thread count auto-tuning

#### 5.2 CUDA Kernels ✅ COMPLETE

| Kernel | File | Purpose |
|:-------|:-----|:--------|
| CF batch evaluation | `LIB/cuda_eval.cu` | Contact function batch GPU computation |
| Shannon histograms | `LIB/ShannonThermoStack/shannon_cuda.cu` | Entropy histogram accumulation |
| TENCoM Hessian | `LIB/tENCoM/tencm_cuda.cu` | Contact discovery + Hessian assembly |
| FastOPTICS k-NN | `LIB/gpu_fast_optics.cu` | GPU-accelerated k-nearest-neighbour for clustering |

`gpu_fast_optics.cu` architecture:
- Grid: N threadblocks (one per query point)
- Block: 256 threads (cooperative scan)
- Shared memory query loading + warp-level priority queue
- Host wrapper: `gpu_foptics_knn()` — upload, launch, download

#### 5.3 Metal Shaders ✅ COMPLETE

| Component | File | Purpose |
|:----------|:-----|:--------|
| CF kernel | `LIB/metal_eval.mm` | Metal GPU CF evaluation (macOS) |
| Shannon bridge | `LIB/ShannonThermoStack/ShannonMetalBridge.{h,mm}` | Metal bridge for Shannon entropy |
| Shader | `LIB/ShannonThermoStack/shannon_metal.metal` | Metal shader language kernels |

#### 5.4 SIMD Vectorization ✅ COMPLETE

| Component | File | Features |
|:----------|:-----|:---------|
| Distance kernels | `LIB/simd_distance.h` | AVX-512/AVX2 pairwise distance computation |
| Dispatch macros | `LIB/hardware_dispatch.{h,cpp}` | Runtime SIMD backend selection |

#### 5.5 OpenMP + Eigen3 ✅ COMPLETE

- OpenMP parallelism in StatMechEngine, tENCoM, Voronoi contacts
- Eigen3 vectorized linear algebra for Hessian assembly, normal modes
- Thread-safe pose pool with OpenMP critical sections

#### 5.6 Unified Dispatch Layer ✅ COMPLETE

```cpp
// LIB/HardwareDispatch.{h,cpp}
class HardwareDispatch {
    BackendType select_backend();      // CUDA > Metal > AVX-512 > AVX2 > scalar
    void dispatch_cf_batch(...);       // Route to best available backend
    void dispatch_entropy_histogram(...);
    void dispatch_distance_matrix(...);
    BenchmarkResult benchmark_all();   // Measure all backends, report throughput
};
```

- `implement-priority-todos-ckxCq`: AVX-512 geometric primitives + benchmark suite
- `implement-todo-item-NUVCa`: OpenMP/Eigen acceleration in statmech + tENCoM
- Runtime backend selection with automatic fallback chain

#### 5.7 Cross-Platform Build ✅ COMPLETE

**CMake Build Options**:
```cmake
option(FLEXAIDS_USE_CUDA    "Enable CUDA GPU evaluation"          OFF)
option(FLEXAIDS_USE_METAL   "Enable Metal GPU acceleration"       OFF)
option(FLEXAIDS_USE_AVX2    "Enable AVX2 SIMD acceleration"       ON)
option(FLEXAIDS_USE_AVX512  "Enable AVX-512 SIMD acceleration"    OFF)
option(FLEXAIDS_USE_OPENMP  "Enable OpenMP thread parallelism"    ON)
option(FLEXAIDS_USE_EIGEN   "Enable Eigen3 vectorised algebra"    ON)
option(FLEXAIDS_USE_256_MATRIX "256×256 soft contact matrix"      ON)
option(FLEXAIDS_USE_MPI     "MPI distributed parallel docking"    OFF)
option(BUILD_PYTHON_BINDINGS "Python bindings via pybind11"       OFF)
option(BUILD_FLEXAIDDS_FAST "Ultra-fast docking exe (LTO+native)" ON)
option(BUILD_TESTING        "Unit tests (requires GoogleTest)"    OFF)
option(BUILD_SWIFT_BRIDGE   "Swift bridge (macOS only)"           OFF)
```

**CMakePresets.json**: Pre-configured presets for Linux GCC, Linux Clang, macOS Clang, Windows MSVC (Debug/Release).

**Windows MSVC Support**:
- MSVC ≥ 19.40 (Visual Studio 2022 17.10+) for C++26 (via `/std:c++latest`)
- `_CRT_SECURE_NO_WARNINGS`, `_USE_MATH_DEFINES`, `NOMINMAX` across all targets
- `flexaids_configure_msvc_test()` helper for consistent test target configuration
- CI: Windows job in GitHub Actions

**Security**: Buffer overflow audit and fixes across 18 source files (`snprintf`, `safe_remark_cat` replacing unsafe `sprintf`/`strcat`).

---

## Cross-Cutting: Apple/Swift Integration

**Status**: ✅ **COMPLETE** (March 23, 2026) — 35 Swift source files, 9,730 lines

### Swift Modules (`swift/Sources/`)

| Module | Purpose |
|:-------|:--------|
| **FleetScheduler** | Distributed docking job scheduling across Apple devices |
| `FleetScheduler.swift` | Job distribution, device capability matching |
| `FleetAggregator.swift` | Merge distributed ChunkResult populations |
| `WorkChunk.swift` | Work unit definition for fleet distribution |
| `DeviceCapability.swift` | Device hardware profiling |
| **FlexAIDCore** | C++ bridge for Swift |
| `FXGA.mm` | Objective-C++ bridge to GA engine |
| **FlexAIDdS** | Main Swift application |
| `DockingRunner.swift` | High-level docking orchestration |
| **Intelligence** | AI-powered docking analysis (15 modules) |
| `IntelligencePipeline.swift` | Unified analysis pipeline |
| `IntelligenceOracle.swift` | Central intelligence coordinator |
| `BindingModeNarrator.swift` | Natural language binding mode descriptions |
| `CampaignJournalist.swift` | Docking campaign summarization |
| `CleftAssessor.swift` | Binding cleft quality assessment |
| `ConvergenceCoach.swift` | GA convergence guidance |
| `FleetExplainer.swift` | Fleet operation explanations |
| `HealthEntropyInsight.swift` | Health-entropy correlation insights |
| `LigandFitCritic.swift` | Ligand fit quality assessment |
| `SelectivityAnalyst.swift` | Binding selectivity analysis |
| `ThermoReferee.swift` | Thermodynamic validation |
| `VibrationalInterpreter.swift` | Vibrational mode interpretation |
| **HealthIntegration** | Apple Health data correlation |
| `BindingEntropyScore.swift` | Entropy-health score mapping |
| **MediaIntegration** | Apple ecosystem integration |
| `FitnessRecommender.swift` | Fitness recommendations from docking state |
| `MusicKitManager.swift` | Ambient audio for long docking runs |

### TypeScript Modules (`typescript/`)

| Component | Purpose |
|:----------|:--------|
| `IntelligenceEngine.ts` | Client-side intelligence analysis engine |
| `IntelligencePanel.tsx` | React panel for intelligence insights |
| `RefereePanel.tsx` | Thermodynamic referee verdict display |
| `FleetDashboard.tsx` | Fleet scheduling and monitoring dashboard |
| `MolstarViewer.tsx` | Mol* 3D structure viewer integration |
| **Analyzers** (`intelligence/`) | |
| `BindingModeAnalyzer.ts` | Binding mode quality analysis |
| `CleftAnalyzer.ts` | Cleft feature analysis |
| `ConvergenceAnalyzer.ts` | GA convergence monitoring |
| `PoseQualityAnalyzer.ts` | Pose quality scoring |
| `SelectivityAnalyzer.ts` | Selectivity evaluation |
| **Shared types** (`packages/shared/src/`) | TypeScript type definitions for all intelligence/fleet data structures |

---

## Cross-Cutting: Testing Strategy

### C++ Tests (GoogleTest)

- 32 test files, 27 test targets — **all passing** (March 23, 2026)
- Coverage: StatMech, BindingMode, GA, tENCoM, Voronoi, FOPTICS, hardware dispatch, cavity detection, grid decomposition, parallel docking, MIF, ion handling, ring conformers, chiral centers, file readers

### Python Tests

- 30 test files in `python/tests/`
- CI pipeline: `test_io.py`, `test_cli.py`, `test_docking.py`, `test_results.py`, `test_energy_matrix.py`, `test_train_256x256.py`, `test_tencm.py`, `test_tencom_results.py`, `test_visualization.py`

### Swift Tests

- `FleetAggregatorTests.swift` — distributed result merging
- `IntelligenceFeatureTests.swift` — intelligence module integration
- `IntelligencePipelineTests.swift` — pipeline end-to-end
- `ThermoRefereeTests.swift` — thermodynamic validation

### TypeScript Tests

- `IntelligenceAnalyzers.test.ts` — analyzer unit tests
- `IntelligenceEngineReferee.test.ts` — referee integration
- `IntelligenceFeatures.test.ts` — feature parity with Swift

### CI/CD

- GitHub Actions: Linux GCC, Linux Clang, macOS Clang, Windows MSVC
- Python test matrix across CI
- GitHub Pages: automatic site stats update + Actions deployment

---

## Cross-Cutting: Additional Modules

### NATURaL (`LIB/NATURaL/`)

- `NATURaLDualAssembly.h/.cpp` — Co-translational folding simulation
- `RibosomeElongation.h/.cpp` — Ribosome elongation modeling
- `TransloconInsertion.h/.cpp` — Membrane protein insertion

### tENCoM (`LIB/tENCoM/`)

- `tencm.{h,cpp}` — Thermal ENtropy Contact Model
- `tencom_diff.h` — Differential entropy computation
- `tencom_output.h` — Result formatting
- `pdb_calpha.h` — Cα extraction
- `tencm_metal.h` — Metal GPU bridge
- `tencm_cuda.cu` — CUDA GPU acceleration

### Ligand Flexibility (`LIB/LigandRingFlex/`)

- `LigandRingFlex.h` — Ligand ring flexibility engine
- `RingConformerLibrary.h` — Ring conformer database
- `SugarPucker.h` — Sugar ring conformation handling

### PTM Attachment (`LIB/PTMAttachment/`)

- `PTMAttachment.h` — Post-translational modification modeling

### Cleft Detection (`LIB/CleftDetector.{h,cpp}`)

- Geometric cleft identification for docking site selection

### File Readers

- `Mol2Reader.{h,cpp}` — Tripos Mol2 format
- `SdfReader.{h,cpp}` — MDL SDF format

---

## Licensing & Compliance

| | |
|:--|:--|
| **License** | Apache-2.0 |
| **Accepted contributions** | Apache-2.0, BSD, MIT, MPL-2.0 |
| **Not accepted** | GPL / AGPL — see [clean-room policy](../docs/licensing/clean-room-policy.md) |
| **CLA** | Required for all contributions |

---

## Risk Management

### Technical Risks

| Risk | Mitigation | Status |
|:-----|:-----------|:-------|
| CUDA/Metal portability | Unified dispatch with automatic fallback | ✅ Resolved |
| Windows build complexity | CMakePresets + MSVC CI job | ✅ Resolved |
| Buffer overflows in C legacy code | Security audit + `snprintf`/`safe_remark_cat` | ✅ Resolved |
| Energy matrix overfitting | Continuous training + dataset adapters | ✅ Mitigated |
| Voronoi CGAL licensing | Use in-house Vcontacts (Apache-2.0) | ✅ Resolved |

### Scientific Risks

| Risk | Mitigation |
|:-----|:-----------|
| Entropy approximation accuracy | Validate against ITC experimental data |
| Hydration entropy parameterization | Planned: empirical calibration against ITC-187 |
| Scoring function transferability | Benchmark against CASF-2016, DUD-E |

---

## Timeline Summary

```
Jan 2026  ──────── Phase 1: Core Thermodynamics ─────────── ✅ COMPLETE
Feb 2026  ──────── Phase 5: Hardware Acceleration ────────── ✅ COMPLETE
          ──────── Phase 2: Python Bindings ──────────────── ✅ COMPLETE
          ──────── Phase 3: NRGSuite/PyMOL GUI ──────────── ✅ COMPLETE
Mar 2026  ──────── Branch consolidation (20+ branches) ──── ✅ March 23
          ──────── Phase 4: Grid/Docking Intelligence ───── 🚧 75%
Apr 2026  ──────── Phase 4: Voronoi hydration ───────────── 🔜 PLANNED
May 2026  ──────── Phase 4: ITC calibration ─────────────── 🔜 PLANNED
          ──────── Manuscript: FlexAIDΔS paper ──────────── 🔜 In preparation
```

---

## Appendix: Key Equations

### Shannon Entropy

$$S = -k_B \sum_i p_i \ln p_i$$

where $p_i = e^{-\beta E_i} / Z$ and $Z = \sum_i e^{-\beta E_i}$

### Helmholtz Free Energy

$$F = -k_B T \ln Z = \langle E \rangle - TS$$

### Total Binding Free Energy

$$\Delta G_{bind} = -k_B T \ln \left( \sum_m e^{-\beta F_m} \right)$$

summing over all binding modes $m$ with intra-mode free energies $F_m$

### Selectivity

$$\text{Selectivity} = \frac{p_{dominant}}{1 - p_{dominant}} \cdot \frac{1}{S_{population}}$$

### SMFREE Fitness

$$f_{SMFREE}(i) = (1 - w_S) \cdot CF_i + w_S \cdot (-T \cdot S_i)$$

where $w_S$ is the entropy weight and $S_i$ is the Shannon entropy contribution of pose $i$
