<div align="center">

# FlexAID∆S

**Entropy-Driven Molecular Docking**

*Combining genetic algorithms with Shannon information theory*
*and statistical mechanics for accurate binding free energy prediction*

[![CI](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml/badge.svg)](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml)
[![Coverage](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/coverage.yml/badge.svg)](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/coverage.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![C++26](https://img.shields.io/badge/C%2B%2B-26-blue.svg)](https://en.cppreference.com/w/cpp/26)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9-3776AB.svg)](https://www.python.org/)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#)
[![Version](https://img.shields.io/badge/version-2.0.0-brightgreen.svg)](VERSION.md)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5b00078-blue)](https://doi.org/10.1021/acs.jcim.5b00078)

</div>

**[Weekly Pulse](#weekly-pulse)** · **[Installation](docs/INSTALLATION.md)** · **[User Guide](docs/USERGUIDE.md)** · **[Support Matrix](docs/SUPPORT_MATRIX.md)** · **[Reproducibility](docs/REPRODUCIBILITY.md)** · **[Benchmarks](docs/BENCHMARKS.md)** · **[Changelog](VERSION.md)** · **[Website](https://lmorency.github.io/FlexAIDdS/)** · **[Documentation Hub](#documentation)**

---

FlexAID∆S extends the [FlexAID](https://doi.org/10.1021/acs.jcim.5b00078) docking engine with a full **canonical ensemble thermodynamics** layer based on Shannon information theory and statistical mechanics. Where conventional docking programs rank poses by enthalpy alone, FlexAID∆S computes the Helmholtz free energy *F* = *H* - *TS* from the partition function over the GA conformational ensemble, accounting for configurational and vibrational entropy contributions that are critical for correct binding mode identification.

**Full thermodynamic documentation**: See [docs/thermodynamics.md](docs/thermodynamics.md) for the complete post-roadmap reference (G_config, H_eff, component averages, diagnostics, experimental features, and support classification).

**NVT vs μVT (competitive binding):** See [docs/theory.md](docs/theory.md). Canonical \(Z\) lives in `StatMechEngine`; grand-canonical \(\Xi(\mu,V,T)\) is `GrandCanonicalEngine` (multi-N + competitive multi-species via `GrandPartitionFunction`). Retitrate without re-docking: `set_concentration([L1, L2, ...])` (NRGsuite / Python / C++). *Make competitive binding great again.*

**Key capabilities:**
- Genetic algorithm docking with Voronoi contact function scoring
- Canonical ensemble partition function, free energy, entropy, and heat capacity
- Grand-canonical (μVT) ensemble: multi-N occupancy, competitive multi-ligand \(\Xi\), \(\langle N\rangle\), selectivity, mixing entropy & ligand entropy collapse
- Concentration-aware retitration (`set_concentration`) without re-docking
- Torsional elastic network model (tENCoM) for backbone vibrational entropy
- Full ligand flexibility: torsions, ring conformers, chiral center discrimination
- Unified hardware dispatch (CUDA, Metal, AVX-512, AVX2, OpenMP)
- Python package with docking API, result analysis, and PyMOL visualization
- Cross-platform: Linux (GCC/Clang), macOS (Clang + Apple Silicon Metal), Windows (MSVC)

---

## Core 1.0 Support Boundary

The repository contains a broad research platform. **Not every visible feature is part of the supported 1.0 contract.**

### Supported for Core 1.0

- `FlexAIDdS` command-line executable (LTO-optimized)
- `FlexAID` legacy-compatible command-line executable
- `tENCoM` command-line executable (vibrational entropy)
- JSON configuration workflows documented for the core engine
- `flexaidds` Python package (pure-Python + optional C++ accelerated mode)
- Core repository documentation (install, validation, reproducibility)
- Benchmark bundles under [`benchmarks/`](benchmarks/)

### Experimental (not part of Core 1.0)

- Swift packages and Apple-device integration layers
- TypeScript / PWA / browser-facing dashboards
- Bonhomme Fleet and iCloud-based distributed execution
- NATURaL co-translational / co-transcriptional workflows
- Backend-specific accelerator paths not covered by the support matrix

See: [`docs/VALIDATED_CAPABILITIES.md`](docs/VALIDATED_CAPABILITIES.md) · [`docs/EXPERIMENTAL_CAPABILITIES.md`](docs/EXPERIMENTAL_CAPABILITIES.md) · [`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md)

---

## Quick Start

### Native tools (full power)

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

```bash
# Dock with full flexibility and entropy at 300 K (default)
./build/FlexAIDdS receptor.pdb ligand.mol2
```

### Python package (easiest for analysis, results loading, thermodynamics)

```bash
pip install -e ./python
```

```python
import flexaidds as fd
print(fd.__version__)

run = fd.load_results("path/to/results")
```

### macOS: Homebrew (native CLI tools)

```bash
brew install --HEAD --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb
pip install flexaidds   # Python package (recommended)
```

See [docs/INSTALLATION.md](docs/INSTALLATION.md) for conda, PyPI, Windows, and full platform instructions. The Python package and native tools can be installed independently.

---

## Weekly Pulse

<table>
<tr>
<td valign="top" width="25%">
<strong>242 commits</strong><br/>
landed on <code>master</code> from 2026-05-27 to 2026-06-03.
</td>
<td valign="top" width="25%">
<strong>19 merged PRs</strong><br/>
spanned thermodynamics, validation, benchmark ops, and web/docs work.
</td>
<td valign="top" width="25%">
<strong>M3 Pro benchmark path</strong><br/>
hardened for iCloud-only launches, Metal pre-flight, temperature fidelity, and best-mode extraction.
</td>
<td valign="top" width="25%">
<strong>Astex non-native grew</strong><br/>
from pair definitions to a 20-pair pilot and the full 85-pair cross-docking set.
</td>
</tr>
</table>

### This Week on `master`

- **Thermodynamics** — `ThermodynamicBreakdown`, component averages, `ThermodynamicOutput`, `TotalSampledPartitionFunction`, Python temperature scans, and targeted thermo verification all landed.
- **Reporting and preprocessing** — cleft annotation, flexible-residue selection, safe standard-state affinity calibration, and JSON-driven PyMOL / Markdown plotting landed.
- **Benchmark operations** — the M3 Pro / iCloud launcher path gained iCloud-only support, `EntryTaskManager.run()`, explicit temperature forwarding, and stricter best-mode validation.
- **Docs and web** — entropy.help, the documentation/validation boundary, `/flexaid`, `/periodic`, and homepage teaser/navigation work all landed this week.

**Full weekly notes and PR links**: see [VERSION.md](VERSION.md).

---

## What's New in v2.0

Released 2026-04-04 — complete rewrite of the FlexAID engine.

- **Entropy-driven scoring** — full canonical ensemble thermodynamics via Shannon information theory (ΔG, ΔH, -TΔS, Cv)
- **Multi-format input with SMILES** — dock directly from SMILES strings with automatic 3D coordinate generation; pure C++26 parser, no RDKit/Boost
- **tENCoM vibrational entropy** — torsional elastic network contact model for backbone vibrational entropy differentials
- **Unified hardware dispatch** — automatic runtime backend: CUDA > Metal > AVX-512 > AVX2 > OpenMP > scalar
- **Full ligand flexibility by default** — ring conformer sampling, sugar pucker, R/S chiral center discrimination
- **Python package & PyMOL plugin** — `flexaidds` with docking API, result I/O, CLI inspector, and interactive visualization
- **GIST desolvation & H-bond scoring** — grid-based water thermodynamics and angular-dependent hydrogen bond potential
- **GA diversity monitoring** — Shannon entropy collapse detection with adaptive catastrophic mutation
- **MPI distributed docking** — grid domain decomposition across compute nodes
- **ML rescoring bridge** — Voronoi graph + Shannon profile feature extraction for hybrid physics/ML scoring

See the full changelog in [VERSION.md](VERSION.md).

---

## Architecture

```
+--------------+    +----------------+    +--------------+    +-------------------+    +-----------------+
|   Input      |    | ProcessLigand  |    |   Genetic    |    |     Scoring       |    | Thermodynamics  |
|              |--->|                |--->|  Algorithm   |--->|                   |--->|                 |
| PDB/CIF +   |    | SMILES parse,  |    |  (gaboom)    |    | Voronoi CF + DEE  |    | StatMech + S    |
| MOL2/SDF/   |    | 3D build,      |    +------+-------+    +-------------------+    +--------+--------+
| SMILES      |    | atom typing    |          |                                              |
+--------------+    +----------------+          v                                              v
                                         +--------------+                              +-----------------+
                                         |  Flexibility |                              |  Binding Modes  |
                                         |              |                              |                 |
                                         | Torsions     |                              | Clustering +    |
                                         | Rings        |                              | dG, dH, -TdS   |
                                         | Chirality    |                              | Cv, Boltzmann   |
                                         | tENCoM       |                              +-----------------+
                                         +--------------+

Hardware: CUDA > Metal > AVX-512 > AVX2 > OpenMP > scalar
```

---

## Features

#### Docking Engine
- **Genetic algorithm** with configurable population, crossover, mutation, selection, and GA diversity monitoring
- **Voronoi contact function (CF)** for shape complementarity scoring
- **GIST water-displacement scoring** -- grid-based explicit solvation term (`GISTEvaluator`, `GISTGrid`)
- **Directional H-bond scoring** -- geometry-aware hydrogen bond potential (`HBondEvaluator`, `hbond_potential.h`)
- **Dead-end elimination (DEE)** reduces ligand conformational search space
- **Batch evaluation** via `VoronoiCFBatch` (AoS) and `VoronoiCFBatch_SoA` (SIMD-friendly) with OpenMP parallelism
- **Multiple clustering** methods: centroid-first, FastOPTICS, Density Peak
- **MPI distributed docking** -- grid domain decomposition across compute nodes (`FLEXAIDS_USE_MPI`)
- **Metal ion scoring** -- receptor-bound ions (Mg2+, Zn2+, Ca2+, Fe2+/3+, Cu2+, Mn2+, Na+, K+, Cl-, Br-, and more) receive crystallographic VdW radii and SYBYL atom types
- **Structural water retention** -- ordered crystallographic waters (B-factor < 20 A^2) participate in Voronoi CF scoring

#### Thermodynamics
- **Canonical ensemble** -- partition function *Z*, Helmholtz free energy *F*, entropy *S*, heat capacity *C_v* (`StatMechEngine`)
- **Grand canonical ensemble** -- competitive ligand binding via `GrandPartitionFunction` (μVT Ξ): concentration-dependent occupancy, pairwise selectivity, intrinsic free energy ranking, mixing entropy, ligand entropy collapse, and `set_concentration` retitration (see `docs/theory.md`)
- **Shannon entropy** -- S = -k_B sum(p_i ln p_i) with log-sum-exp numerical stability
- **Torsional ENCoM** (tENCoM) -- backbone vibrational entropy without full rotamer rebuilds
- **ShannonThermoStack** -- combined configurational + vibrational entropy pipeline
- **Thermodynamic integration** and WHAM free energy profiles
- **Boltzmann weight** normalization and fast lookup table

#### Molecular Flexibility
- **Full flexibility by default** -- ligand torsions, ring conformers, chirality, intramolecular scoring at 300 K
- **Non-aromatic ring sampling** -- chair/boat/twist for 6-membered, envelope/twist for 5-membered rings, sugar pucker
- **Chiral center discrimination** -- explicit R/S sampling with stereochemical energy penalty
- **Multi-format input** -- PDB, CIF/mmCIF, MOL2, SDF/MOL V2000, SMILES (with automatic 3D build)

#### Hardware Acceleration
- **Unified hardware dispatch** -- automatic backend selection at runtime
- **CUDA** -- batch CF evaluation and Shannon entropy histograms (Volta through Blackwell, sm_70–sm_120; Blackwell requires CUDA ≥ 12.6)
- **Metal** -- Apple Silicon GPU for Shannon entropy, cavity detection, and evaluation
- **SIMD** -- AVX-512 and AVX2 vectorised geometric primitives
- **OpenMP + Eigen3** -- thread parallelism and vectorised linear algebra
- **LTO binaries** -- link-time optimized executables

#### Python Package
- **`flexaidds`** -- docking API, result I/O, thermodynamics, CLI inspector
- **Pure-Python fallback** -- works without C++ compilation
- **pybind11 C++ bridge** -- `StatMechEngine`, `ENCoMEngine`, `Thermodynamics`, `VibrationalEntropy`
- **PyMOL plugin** -- binding mode visualization, entropy heatmaps, ITC comparison
- **Grok Imagine integration** (via flexaidds.figures + skill) -- auto NRDD-quality cover figure + 6 s animation of best binding mode (Gate 6, real ledger values, branding + metadata) after successful runs with --visualize

---

## Build

### Requirements

- **Required**: C++26 compiler (GCC >= 14, Clang >= 18, Apple Clang >= 16 / Xcode 16, MSVC >= 19.40 / VS 2022 17.10), CMake >= 3.28
- **Recommended**: Eigen3 (`libeigen3-dev` / `brew install eigen`)
- **Optional**: OpenMP, CUDA Toolkit, Metal framework (macOS), pybind11, MPI

### Output Binaries

| Binary | Description |
|:-------|:------------|
| `FlexAID` | Standard docking executable |
| `FlexAIDdS` | Optimized docking (LTO + `-march=native`) |
| `tENCoM` | Vibrational entropy differential tool |

### Build Variants

```bash
# Standard release
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# With tests (48 C++ tests via GoogleTest)
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel && ctest --test-dir build --output-on-failure

# With Python bindings
cmake -S . -B build -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# HPC deployment (AVX-512 + OpenMP)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_AVX512=ON -DFLEXAIDS_USE_OPENMP=ON

# Apple Silicon with Metal GPU
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_METAL=ON
```

```bash
# Distributed validation (parallel C++ tests + MPI benchmark smoke run)
MPI_PROCS=4 CTEST_JOBS=8 ./scripts/run_distributed_validation.sh
```

### CMake Options

| Option | Default | Description |
|:-------|:--------|:------------|
| `BUILD_FLEXAIDDS_FAST` | **ON** | LTO-optimized FlexAIDdS binary |
| `BUILD_TESTING` | OFF | GoogleTest unit tests (48 test targets) |
| `BUILD_PYTHON_BINDINGS` | OFF | pybind11 Python extension (`_core`) |
| `BUILD_SWIFT_BRIDGE` | OFF | Swift bridge (macOS only, experimental) |
| `ENABLE_TENCOM_TOOL` | **ON** | tENCoM vibrational entropy tool |
| `ENABLE_TENCOM_BENCHMARK` | OFF | tENCoM benchmark binary |
| `ENABLE_VCFBATCH_BENCHMARK` | OFF | VoronoiCFBatch benchmark binary |
| `FLEXAIDS_USE_CUDA` | OFF | CUDA GPU acceleration |
| `FLEXAIDS_USE_METAL` | OFF | Metal GPU acceleration (macOS) |
| `FLEXAIDS_USE_AVX2` | **ON** | AVX2 SIMD acceleration |
| `FLEXAIDS_USE_AVX512` | OFF | AVX-512 SIMD acceleration |
| `FLEXAIDS_USE_OPENMP` | **ON** | OpenMP thread parallelism |
| `FLEXAIDS_USE_EIGEN` | **ON** | Eigen3 vectorised linear algebra |
| `FLEXAIDS_USE_256_MATRIX` | **ON** | 256x256 soft contact matrix system |
| `FLEXAIDS_USE_MPI` | OFF | MPI distributed parallel docking |

---

## Usage

### Command Line

```bash
# Full flexibility dock (entropy at 300 K by default)
./build/FlexAIDdS receptor.pdb ligand.mol2

# JSON config override
./build/FlexAIDdS receptor.pdb ligand.mol2 -c config.json

# Rigid screening (no flexibility, no entropy)
./build/FlexAIDdS receptor.pdb ligand.mol2 --rigid
```

All parameters have built-in defaults. Override only what you need via JSON:

```json
{
  "thermodynamics": { "temperature": 310, "clustering_algorithm": "DP" },
  "ga": { "num_chromosomes": 2000, "num_generations": 1000 },
  "flexibility": { "ligand_torsions": true, "ring_conformers": true }
}
```

### Vibrational Entropy (tENCoM)

```bash
tENCoM reference.pdb target1.pdb [target2.pdb ...] [-T temp] [-r cutoff] [-k k0] [-o prefix]
```

tENCoM vibrational entropy is currently a relative heuristic unless a
benchmark-specific eigenvalue-to-frequency calibration is supplied. See
[docs/TENCOM_ENTROPY_CALIBRATION.md](docs/TENCOM_ENTROPY_CALIBRATION.md).

### Python Package

```bash
pip install -e ./python
```

The `flexaidds` package works in two modes: **pure Python** (always available) and **C++ accelerated** (when built with `BUILD_PYTHON_BINDINGS=ON`).

```python
import flexaidds as fd

# High-level docking
results = fd.dock(receptor="receptor.pdb", ligand="ligand.mol2", compute_entropy=True)

# Load existing results
docking = fd.load_results("output_prefix")
for mode in docking.binding_modes:
    print(f"Mode {mode.rank}: dG={mode.free_energy:.2f}, S={mode.entropy:.3f}")
```

```python
from flexaidds import StatMechEngine

engine = StatMechEngine(temperature=300)
engine.add_energies(pose_energies)
thermo = engine.compute()
print(f"F = {thermo.free_energy:.2f} kcal/mol")
print(f"S = {thermo.entropy:.4f} kcal/(mol*K)")
```

### CLI Inspector

```bash
python -m flexaidds /path/to/results/              # summary table
python -m flexaidds /path/to/results/ --top 5      # top 5 modes
python -m flexaidds /path/to/results/ --json        # JSON output
python -m flexaidds /path/to/results/ --csv out.csv  # CSV export
```

<details>
<summary><strong>PyMOL Plugin</strong></summary>

**Installation**: PyMOL > Plugin Manager > Install New Plugin > select `pymol_plugin/` > restart.

| Command | Description |
|:--------|:------------|
| `flexaids_load <dir> [temp]` | Load results from output directory |
| `flexaids_show_ensemble <mode>` | Display all poses in a binding mode |
| `flexaids_color_boltzmann <mode>` | Color by Boltzmann weight |
| `flexaids_thermo <mode>` | Print thermodynamic properties |
| `flexaids_entropy_heatmap <mode>` | Spatial entropy density heatmap |
| `flexaids_dock <obj> <lig>` | Interactive docking from PyMOL |

Requires: `pip install -e python/`

</details>

<details>
<summary><strong>Legacy Mode</strong></summary>

```bash
./FlexAID config.inp ga.inp output_prefix
./FlexAIDdS --legacy config.inp ga.inp output_prefix
```

</details>

---

## Testing

### C++ (GoogleTest)

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

48 test targets covering: StatMechEngine, BindingMode, GA operators, GrandPartitionFunction, competitive binding, tENCoM, Voronoi contacts, FastOPTICS clustering, hardware dispatch, cavity detection, ion handling, ring conformers, chiral centers, file readers, grid decomposition, MIF scoring, GIST/HBond scoring, GA diversity, distributed backend, and more.

### Python (pytest)

```bash
pip install -e ./python && pytest tests/
```

32 test files. Tests marked `@requires_core` skip gracefully when the C++ extension is not built.

### CI

GitHub Actions pipeline: Linux GCC, Linux Clang, macOS Clang (allow-fail), Windows MSVC (allow-fail). Python smoke tests on all platforms. Additional workflows: `license-scan.yml` (Apache compliance), `perf.yml` (benchmark regression), `sanitizers.yml` (ASan/UBSan).

See [`docs/SUPPORT_MATRIX.md`](docs/SUPPORT_MATRIX.md).

---

## Reproducibility

Benchmark and scientific performance claims should be interpreted through the repository reproducibility policy.

**A claim is not repository-reproducible merely because it appears in documentation.** It becomes reproducible only when the corresponding bundle exists under [`benchmarks/`](benchmarks/) with dataset provenance, commands, expected outputs, and metric definitions.

See: [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) · [`benchmarks/README.md`](benchmarks/README.md)

---

## Scientific Background

FlexAID∆S treats the GA conformational ensemble as a **canonical ensemble** (*N*, *V*, *T* fixed):

```
Z = sum exp[-beta * E_i]              (partition function)
F = -k_B * T * ln(Z)                  (Helmholtz free energy)
<E> = sum p_i * E_i                    (mean energy)
S = -k_B * sum p_i * ln(p_i)          (Shannon configurational entropy)
C_v = k_B * beta^2 * (<E^2> - <E>^2)  (heat capacity)
```

**Scoring** uses two complementary layers: a **Voronoi contact function** for geometry-based shape complementarity, weighted by a **2-term LJ+Coulomb potential** parameterised over 40 SYBYL atom types. The StatMechEngine converts the GA ensemble into thermodynamic quantities via the canonical partition function with log-sum-exp numerical stability.

---

<details>
<summary><strong>Repository Structure</strong></summary>

```
FlexAIDdS/
+-- LIB/                    # Core C++ library (~100+ source files)
|   +-- flexaid.h            # Main header: constants, data structures
|   +-- gaboom.cpp/h         # Genetic algorithm engine
|   +-- Vcontacts.cpp/h      # Voronoi contact function
|   +-- statmech.cpp/h       # StatMechEngine
|   +-- BindingMode.cpp/h    # Pose clustering + thermodynamics
|   +-- encom.cpp/h          # ENCoM vibrational entropy
|   +-- GrandPartitionFunction.* | Grand canonical ensemble for competitive binding |
|   +-- config_defaults.h    # Default parameter schema
|   +-- config_parser.cpp/h  # JSON config system
|   +-- hardware_dispatch.*  # Unified HW backend selection
|   +-- ShannonThermoStack/  # Shannon entropy + HW acceleration
|   +-- tENCoM/              # Torsional ENCoM module
|   +-- LigandRingFlex/      # Ring conformer sampling
|   +-- ChiralCenter/        # R/S discrimination
|   +-- NATURaL/             # Co-translational assembly (experimental)
|   +-- CavityDetect/        # SURFNET cavity detection
|   +-- Mol2Reader.cpp/h     # MOL2 file reader
|   +-- SdfReader.cpp/h      # SDF file reader
|   +-- CleftDetector.cpp/h  # Binding-site detection
|   +-- VoronoiCFBatch.h     # Batch Voronoi CF (header-only)
+-- src/                    # Entry points
+-- tests/                  # GoogleTest suite (48 test targets)
+-- python/                 # Python package + pybind11 bindings
|   +-- flexaidds/           # Python package (22+ modules)
|   +-- bindings/            # C++ bridge
|   +-- tests/               # pytest suite (32 test files)
|   +-- setup.py
+-- pymol_plugin/           # PyMOL visualization plugin
+-- swift/                  # Swift package (experimental)
+-- typescript/             # TypeScript SDK (experimental)
+-- docs/                   # Documentation
+-- benchmarks/             # Reproducibility bundles
+-- cmake/                  # CMake helpers
+-- .github/workflows/      # CI/CD
+-- CMakeLists.txt          # Build configuration
+-- VERSION.md              # Version history and changelog
```

</details>

---

## Documentation

Browse the full documentation at **[lmorency.github.io/FlexAIDdS](https://lmorency.github.io/FlexAIDdS/)**.

### Core Guides

| Document | Description |
|:---------|:------------|
| [Installation Guide](docs/INSTALLATION.md) | Prerequisites, build instructions, platform-specific notes |
| [User Guide](docs/USERGUIDE.md) | Full parameter reference, Python API, PyMOL plugin |
| [Support Matrix](docs/SUPPORT_MATRIX.md) | Supported platforms, compilers, GPU backends, CI coverage |
| [Benchmarks](docs/BENCHMARKS.md) | ITC-187, CASF-2016, LIT-PCBA, cross-docking validation |
| [Reproducibility](docs/REPRODUCIBILITY.md) | Benchmark claim policy, dataset provenance |
| [VERSION.md](VERSION.md) | Version history and changelog |
| [Contributing](CONTRIBUTING.md) | Development setup, license policy, PR guidelines |

### Technical Deep Dives

| Document | Description |
|:---------|:------------|
| [Scoring Functions](docs/docs/scoring/overview.md) | Multi-component scoring: Voronoi CF, H-bond, GIST, electrostatics |
| [H-Bond Potential](docs/docs/scoring/hbond.md) | Angular-dependent Gaussian hydrogen bond scoring |
| [GIST Desolvation](docs/docs/scoring/gist.md) | Grid Inhomogeneous Solvation Theory integration |
| [Genetic Algorithm](docs/docs/ga/overview.md) | Fitness models (SMFREE, PSHARE, LINEAR), key parameters |
| [GA Diversity](docs/docs/ga/diversity.md) | Entropy collapse detection and adaptive mutation |
| [GA Optimizer](docs/docs/ga/optimize.md) | Automated hyperparameter tuning with GAOptimizer |
| [Python API](docs/docs/api/python.md) | StatMechEngine, ENCoMEngine, dock(), load_results() |
| [ML Rescoring](docs/docs/ml-rescore.md) | Feature extraction bridge for hybrid physics/ML scoring |
| [Configuration](docs/docs/configuration.md) | JSON config reference for scoring, GA, and distributed settings |

### Benchmark Suites

| Document | Description |
|:---------|:------------|
| [Astex Diverse](benchmarks/datasets/astex_diverse.yaml) | FlexAID 2015 validation #1: native-pose docking on 85 holo complexes |
| [Astex Non-Native](benchmarks/datasets/astex_nonnative.yaml) | FlexAID 2015 validation #2: cross-docking into non-native receptor conformations |
| [HAP2](benchmarks/datasets/hap2.yaml) | FlexAID 2015 validation #3: holo/apo protein-pair benchmark from Gaudreault & Najmanovich |
| [CASF-2016](docs/docs/benchmarks/casf2016.md) | Scoring and docking power on 285 diverse complexes |
| [LIT-PCBA](docs/docs/benchmarks/litpcba.md) | Unbiased virtual screening across 15 PubChem targets |
| [Cross-Docking](docs/docs/benchmarks/crossdock.md) | Non-native receptor docking validation |

### Internal

| Document | Description |
|:---------|:------------|
| [Implementation Roadmap](docs/IMPLEMENTATION_ROADMAP.md) | Development phases and status |
| [Validated Capabilities](docs/VALIDATED_CAPABILITIES.md) | Core 1.0 supported features |
| [Experimental Capabilities](docs/EXPERIMENTAL_CAPABILITIES.md) | Experimental / non-1.0 features |
| [Known Limitations](docs/KNOWN_LIMITATIONS.md) | Current limitations and caveats |

---

## The Thermodynamic Core: dG = E_enthalpic − T·S_entropic

FlexAID∆S computes binding free energy via the fundamental **enthalpy–entropy decomposition**:

```
dG = E_CF  +  (−T · S_total)
   = E_CF  −  T · [S_conf + S_vib + S_rot + S_cratic]

   ─────────   ──────────────────────────────────────
   exp(−CF/T)   −T · Σ pᵢ ln(pᵢ)        (the two fundamental terms)
   Enthalpic    Entropic penalty           (binding trades freedom for energy)
```

This structure — `G = H − TS` — is the standard thermodynamic decomposition of free energy into enthalpic and entropic contributions.

Odrzywołek (2026, arXiv:2603.21852) proved that the operator `eml(a,b) = exp(a) − ln(b)` is functionally complete for elementary functions — any arithmetic operation, transcendental function, or entropy formula can be expressed as compositions of eml and the constant 1. The domain-specific formulas below are each instances of composed eml operations, providing a natural algebraic language for cross-domain comparison.

In FlexAID∆S specifically:

| Abstract term | FlexAID∆S realization | Physical meaning |
|:---|:---|:---|
| `exp(a)` | `exp(−β Eᵢ)` — Boltzmann weights in Z | Energy/enthalpy landscape of conformational states |
| `−ln(b)` | `−T · k_B · Σ pᵢ ln(pᵢ)` — Shannon configurational entropy | Conformational freedom lost upon binding |
| Combined | `F = −kT ln Z` — Helmholtz free energy | Net thermodynamic binding strength |

The **log-sum-exp trick** provides numerical stability for the partition function:

```
ln Z = max_w + ln( Σ exp(wᵢ − max_w) )

where wᵢ = ln(nᵢ) − β · Eᵢ       (unnormalized log-weights)
      β = 1/(k_B T)                 (inverse temperature)
```

The entropy follows from the thermodynamic identity:

```
S = (⟨E⟩ − F) / T  =  k_B · H · ln(2)

where H = −Σ pᵢ log₂(pᵢ)  (Shannon entropy in bits)
      S = k_B · H · ln(2)  (Gibbs entropy in J/K)
```

Jaynes (1957) proved these are the **same mathematical object** in different units — this identity is the bridge that connects molecular docking to LLM safety and physiological biofeedback via a shared entropy kernel.

---

## From Molecules to Machines: The Triptych

The configurational entropy kernel at the heart of FlexAID∆S has been **ported directly** to [Shannon](https://github.com/LeBonhommePharma/Shannon), a zero-shot detector for evaluation awareness and strategic deception in frontier LLM agents. The same `EntropyCalculator` also powers [NATURaL](https://github.com/LeBonhommePharma/NATURaL), which detects autonomic drug response signatures from Apple Watch HRV data. The entropy computation uses the same `H = −Σ pᵢ log₂(pᵢ)` formula across all three domains. Jaynes' identity ensures this is mathematically sound.

```
🔬 FlexAID∆S (in silico)        💊 NATURaL (in vivo)          🤖 Shannon (in silico)
─────────────────────        ──────────────────        ─────────────────────
Torsional angles (°)         RR intervals (ms)         Logit scores
H = −Σ pᵢ log₂(pᵢ)          H = −Σ pᵢ log₂(pᵢ)       H = −Σ pᵢ log₂(pᵢ)
Binding → ΔS < 0             Drug → ΔH < 0             Eval-aware → ΔH < 0
            ↘                         ↘                         ↙
                  Shared entropy kernel. Distributional phase transitions in each domain.
```

### FlexAID∆S ↔ Shannon correspondence

```
FlexAID∆S (molecular docking)          Shannon (LLM safety)
─────────────────────────────          ──────────────────────
Conformational states {i}          ↔   Token vocabulary {i}
Energy Eᵢ                           ↔   Logit wᵢ
Boltzmann: pᵢ = exp(−Eᵢ/kT)/Z      ↔   Softmax: pᵢ = exp(wᵢ)/Z
S = −k_B · Σ pᵢ ln(pᵢ)             ↔   H = −Σ pᵢ log₂(pᵢ)
Binding mode lock-in                ↔   Evaluation awareness lock-in
ΔS_config < 0 (entropy collapses)   ↔   ΔH < −3.2 bits (entropy collapses)
Temperature T (Kelvin)              ↔   τ_eff = 1/σ²(logit variance)
Grand canonical Xi                  ↔   Multi-model monitoring
Thermostat (canonical ensemble)     ↔   Handrail engine (feedback control)
```

When a drug molecule locks into a binding pocket, the conformational entropy collapses from a multimodal distribution to a unimodal one — the same signature Shannon detects when an LLM shifts from normal generation to evaluation-aware strategic behavior.

### FlexAID∆S ↔ NATURaL correspondence

```
FlexAID∆S (computational)             NATURaL (physiological)
────────────────────────────          ──────────────────────────
Conformational ensemble           ↔   RR-interval distribution
ΔS_config (bits) from torsions    ↔   ΔH_hrv (bits) from HRV
Binding → entropy collapse        ↔   Sympathomimetic → entropy collapse
Binding → entropy unchanged       ↔   Parasympathomimetic → entropy expansion
Partition function Z              ↔   HRV power spectrum (frequency-domain Z)
−TΔS (kcal/mol)                   ↔   ΔH_hrv (bits) from Apple Watch
```

NATURaL's `CrossDomainValidator` performs **three-way entropy validation**: correlating |ΔS_config| (FlexAID∆S computational), |−TΔS| (SCORPIO ITC measured), and |ΔH_hrv| (NATURaL in-vivo) across substances — providing preliminary evidence that entropy-collapse patterns may correlate across molecular and physiological domains (a hypothesis under active investigation).

### Cross-domain correspondence table

| Dimension | FlexAID∆S | NATURaL | Shannon |
|:---|:---|:---|:---|
| System | Ligand + protein | Autonomic nervous system | LLM agent |
| Observable | Torsional angles | RR intervals | Logits / token probs |
| Distribution pᵢ | Boltzmann | HRV histogram | Softmax |
| Entropy H | −Σ pᵢ log₂(pᵢ) | −Σ pᵢ log₂(pᵢ) | −Σ pᵢ log₂(pᵢ) |
| Collapse signal | ΔS < 0 (binding) | ΔH < 0 (sympathomimetic) | ΔH < −3.2 bits (eval-aware) |
| Expansion signal | ΔS > 0 (solvation) | ΔH > 0 (parasympathomimetic) | — |
| Temperature | T (Kelvin) | 37°C (body) | Logit spread (heuristic analog) |
| Partition function | Z = Σ exp(−βEᵢ) | HRV spectral power | Softmax normalizer Z (analogous) |
| Control | Canonical ensemble | Drug dose | Handrail engine |

All three domains exhibit distributional narrowing detectable by the same `H = −Σ pᵢ log₂(pᵢ)` kernel. The mathematical identity (Jaynes 1957) ensures the kernel transfers; the domain-specific interpretation is different in each case.

The entropy kernel, the log-sum-exp partition function, and the sliding-window z-score detector are **identical code** across projects. The domain changes (molecules → tokens → heartbeats), but the mathematics doesn't.

---

## Publications

If you use FlexAID∆S in your research, please cite:

> Gaudreault F & Najmanovich RJ (2015). FlexAID: Revisiting Docking on Non-Native-Complex Structures.
> *J. Chem. Inf. Model.* 55(7):1323-36.
> [DOI:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)

**Related publications:**

- Morency LP & Najmanovich RJ (2026). Shannon: Entropy Collapse Detection for Evaluation Awareness in Frontier LLM Agents. *Manuscript in preparation.*
- Morency LP & Najmanovich RJ (2026). FlexAID∆S: Information-Theoretic Entropy Improves Molecular Docking Accuracy and Binding Mode Prediction. *Manuscript in preparation.*
- Gaudreault F, Morency LP & Najmanovich RJ (2015). NRGsuite: a PyMOL plugin to perform docking simulations in real time. *Bioinformatics* 31(23):3856-8. [DOI:10.1093/bioinformatics/btv458](https://doi.org/10.1093/bioinformatics/btv458)
- Frappier V et al. (2015). A Coarse-Grained Elastic Network Atom Contact Model. *Proteins* 83(11):2073-82. [DOI:10.1002/prot.24922](https://doi.org/10.1002/prot.24922)

Full list: [Google Scholar](https://scholar.google.ca/citations?user=amFCT0oAAAAJ&hl=en)

---

## Contributing

Contributions are welcome. Please read [CONTRIBUTING.md](CONTRIBUTING.md) before submitting pull requests.

| | |
|:--|:--|
| **Accepted licenses** | Apache-2.0, BSD, MIT, MPL-2.0 |
| **Not accepted** | GPL / AGPL -- see [clean-room policy](docs/licensing/clean-room-policy.md) |
| **CLA** | Required for all contributions |

---

## License

[Apache License 2.0](LICENSE) -- free for academic and commercial use.

See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for dependency licenses.

---

## Grok / Codex Agent Skills (flexaidds)

This repository packages a first-class Grok skill for all FlexAID, FlexAIDdS, and FlexAID∆S work:

```
.grok/skills/flexaidds/
├── SKILL.md
├── scripts/validate_skill.py
└── references/flexaidds-guidance.md
```

**Invocation (user-facing trigger phrases + native slash):**
- `/flexaidds`
- `/FlexAid docking`
- `/FlexAidDS`
- Mention `FlexAIDdS`, `FlexAID∆S`, "FlexAID docking", "ensemble analysis", or "thermodynamic ledger" in any prompt.

**Validation (ALWAYS run before claiming completion):**
```bash
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 -m pytest tests/test_flexaid_skill.py -q --tb=line
```

The skill **mandates**:
- Start every task with `git status` + the exact mandated `find` command for skill/XML/AGENTS files.
- Never edit scoring, ranking, or thermodynamic kernels without tests + explicit request.
- Use precise language: "CF/contact-function scoring proxy" vs. full thermodynamic ledger.
- Produce only chunked plans; never rewrite history or merge without confirmation.

See the skill's own [SKILL.md](.grok/skills/flexaidds/SKILL.md) and [references/flexaidds-guidance.md](.grok/skills/flexaidds/references/flexaidds-guidance.md) for the full guardrails and terminology contract.

**For Claude Code / Codex / other agents**: Copy or symlink the `.grok/skills/flexaidds/` tree into your agent's skill directory (or `.skills/flexaidds/`) and reference the same validation commands.

---

<div align="center">

Louis-Philippe Morency
[NRGlab](http://biophys.umontreal.ca/nrg), Departement de biochimie et medecine moleculaire
Universite de Montreal

[Repository](https://github.com/LeBonhommePharma/FlexAIDdS) · [Issues](https://github.com/LeBonhommePharma/FlexAIDdS/issues) · [NRGlab GitHub](https://github.com/NRGlab)

</div>
