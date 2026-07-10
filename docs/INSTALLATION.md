# Installation Guide

Complete build and installation instructions for FlexAID∆S on all supported platforms.

---

## Prerequisites

### Required

| Dependency | Minimum Version | Notes |
|:-----------|:----------------|:------|
| C++ compiler | GCC ≥ 14, Clang ≥ 18, Apple Clang ≥ 16 (Xcode 16), MSVC ≥ 19.40 (VS 2022 17.10) | C++26 support required |
| CMake | ≥ 3.28 | Build system (needed for `CXX_STANDARD 26`) |
| Python | ≥ 3.9 | For the `flexaidds` Python package |

### Optional

| Dependency | Purpose | Install |
|:-----------|:--------|:--------|
| Eigen3 | Vectorised linear algebra | `apt install libeigen3-dev` / `brew install eigen` |
| OpenMP | Thread parallelism | `apt install libomp-dev` / `brew install libomp` |
| CUDA Toolkit | NVIDIA GPU acceleration | [developer.nvidia.com](https://developer.nvidia.com/cuda-toolkit) |
| ROCm/HIP | AMD GPU acceleration (MI100/MI200/MI300) | [rocm.docs.amd.com](https://rocm.docs.amd.com/) |
| Metal framework | Apple GPU acceleration | Included with Xcode on macOS |
| pybind11 | Python ↔ C++ bindings | `pip install pybind11[global]` |
| Ninja | Faster builds | `apt install ninja-build` / `brew install ninja` |

> **Note**: No RDKit or Boost dependency is required. The ProcessLigand module (SMILES parsing, ring perception, aromaticity detection, 3D coordinate building) is implemented in pure C++26 + Eigen.

---

## Quick Install

### Python package via pip (easiest for analysis, results, thermodynamics)

> **Status (2026-07):** `flexaidds` is **not yet published on the public PyPI index**.
> Use the GitHub install below until the first PyPI release ships (GitHub Actions
> workflow `.github/workflows/pypi-release.yml` + trusted publishing). After that,
> `pip install flexaidds` will work as usual.

**Recommended today — install from GitHub (no clone required):**
```bash
pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
# later upgrade:
pip install --upgrade --force-reinstall \
  "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
# or: python -m flexaidds --self-update   # uses GitHub Releases + git fallback
```

From a local checkout (development):
```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
pip install -e ./python
```

Once published on PyPI:
```bash
pip install flexaidds
pip install --upgrade flexaidds
```

- Installs the `flexaidds` package (and the `flexaidds` + `flexaidds-benchmark` CLIs).
- The native C++ extension (`_core`) is **optional**. If build tools/Eigen are missing **or** compilation fails, the package still installs and works in pure-Python mode (with fallbacks). Set `FLEXAIDDS_SKIP_CORE=1` to force pure-Python.
- Works from sdist, git, local checkout, and (after first publish) PyPI.

Verify:
```bash
python -c "import flexaidds as fd; print(fd.__version__, 'HAS_CORE=', getattr(fd, 'HAS_CORE_BINDINGS', False))"
python -m flexaidds --help || echo "CLI works via python -m"
python -m flexaidds --check-update
```

A GitHub Actions workflow (`.github/workflows/pypi-release.yml`) handles building wheels + sdist and publishing via trusted publishing on release (or manual `workflow_dispatch` to TestPyPI/PyPI).

### Full native tools + Python (CMake)

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j $(nproc)
```

This produces the main executables (`FlexAIDdS`, `FlexAID`, `tENCoM`, ...).

### Conda

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
conda env create -f environment.yml
conda activate flexaidds
```

Inside the env you get the Python package (editable). You can additionally build the full native tools with the usual cmake commands (conda compilers/eigen are already available).

To build a conda package:
```bash
conda-build conda   # produces a flexaidds conda package
```

See the table below for per-platform dependency commands.

---

## Platform-Specific Instructions

### Linux (Ubuntu/Debian)

```bash
# Install all dependencies
sudo apt-get update
sudo apt-get install -y cmake ninja-build libeigen3-dev libomp-dev g++ python3-dev

# Build
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

### macOS (Apple Silicon & Intel)

#### Easy install via Homebrew (recommended for native CLI tools)

```bash
# Stable tagged release (v2.0.0+)
brew install --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb

# Or latest development tip
brew install --HEAD --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb

# Update later
brew reinstall --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb
# or, for HEAD builds:
brew reinstall --HEAD --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb
```

This installs `FlexAIDdS`, `tENCoM`, `FlexAID` + required data files.

Then (recommended — GitHub until PyPI publish):
```bash
pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
# after first PyPI release: pip install flexaidds
```

#### Build from source (full control)

```bash
brew install cmake ninja libomp eigen
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

> **Apple Silicon note**: AVX2/AVX-512 flags are automatically disabled on ARM64. Metal GPU acceleration is available with `-DFLEXAIDS_USE_METAL=ON`.

### Windows (Visual Studio 2022)

```cmd
REM Open "x64 Native Tools Command Prompt for VS 2022"
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git && cd FlexAIDdS
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_OPENMP=OFF -DFLEXAIDS_USE_EIGEN=OFF
cmake --build build --parallel
```

> **Windows note**: Install Eigen via `choco install eigen` or manually set `CMAKE_PREFIX_PATH`. OpenMP support on Windows requires additional configuration.

---

## Python Package Only (pip / conda) — no full C++ build required

Most users who want to **analyze results**, load docking runs, use thermodynamic models, or run benchmarks can use just the Python package:

```bash
# pip (from repo root)
pip install -e ./python

# conda (recommended for complex dependency environments)
conda env create -f environment.yml
conda activate flexaidds
```

The package provides:
- `flexaidds.load_results`, data models, CLI (`python -m flexaidds`)
- Pure-Python + optional accelerated `StatMechEngine`, thermodynamics
- Dataset runners, visualization helpers, etc.

The heavy compiled `_core` extension (for maximum speed on StatMech/ENCoM) is built automatically during `pip install -e ./python` **if** a compiler + Eigen + pybind11 are present. Otherwise it silently falls back (see `HAS_CORE_BINDINGS`).

See also `python/README.md`.

## Homebrew (macOS native tools)

The formula provides the high-performance native executables with data files staged correctly:

```bash
# Stable release
brew install --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb

# Development / latest
brew install --HEAD --formula https://raw.githubusercontent.com/LeBonhommePharma/FlexAIDdS/master/Formula/flexaidds.rb

# After install, add the Python package (GitHub until PyPI publish):
pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
# after first PyPI release: pip install flexaidds
```

The formula lives at `Formula/flexaidds.rb`. It supports both a stable `url`/`sha256` (tagged release) and `head` for tip-of-tree. Bottles may be added later.

When cutting a new stable release, update the formula’s `url`, `sha256`, and version together:
```bash
# Example maintainer steps
TAG=v2.0.1
curl -sL "https://github.com/LeBonhommePharma/FlexAIDdS/archive/refs/tags/${TAG}.tar.gz" | shasum -a 256
# paste sha256 into Formula/flexaidds.rb, commit, push
```

---

## Releasing & Distribution

- **Python package (PyPI)**: The workflow `.github/workflows/pypi-release.yml` builds sdist + wheels (cibuildwheel) and publishes on GitHub Release using trusted publishing. Run the workflow manually (`workflow_dispatch`) for TestPyPI. The optional `_core` extension is built when possible; otherwise pure-Python wheels/sdists still publish.
  - Required once per environment: configure PyPI / TestPyPI **trusted publishing** for the `pypi` and `testpypi` GitHub Environments.
  - Smoke checks: sdist install + `import flexaidds` run in CI before publish.

- **Homebrew**: Update `Formula/flexaidds.rb` (`url` + `sha256` for stable; `head` tracks `master`) when cutting releases. Users install via the raw formula URL (or a personal tap).

- **Native binaries**: Existing release workflow attaches platform archives on tag.

- **In-package updater**: `python -m flexaidds --check-update` / `--self-update` uses the GitHub Releases API and `pip install --upgrade flexaidds`.

---

## Build Options

All CMake options with defaults:

| Option | Default | Description |
|:-------|:--------|:------------|
| `BUILD_FLEXAIDDS_FAST` | **ON** | LTO-optimized FlexAIDdS binary |
| `ENABLE_TENCOM_TOOL` | **ON** | tENCoM vibrational entropy tool |
| `FLEXAIDS_USE_CUDA` | OFF | NVIDIA GPU acceleration (Volta → Blackwell) |
| `FLEXAIDS_USE_ROCM` | OFF | AMD GPU acceleration (MI100/MI200/MI300) |
| `FLEXAIDS_USE_METAL` | OFF | Apple GPU acceleration (macOS only) |
| `FLEXAIDS_USE_AVX2` | **ON** | AVX2 SIMD (auto-disabled on ARM) |
| `FLEXAIDS_USE_AVX512` | OFF | AVX-512 SIMD acceleration |
| `FLEXAIDS_USE_OPENMP` | **ON** | OpenMP thread parallelism |
| `FLEXAIDS_USE_EIGEN` | **ON** | Eigen3 vectorised linear algebra |
| `FLEXAIDS_USE_256_MATRIX` | **ON** | 256×256 soft contact matrix system |
| `BUILD_PYTHON_BINDINGS` | OFF | pybind11 Python extension (`_core`) |
| `BUILD_TESTING` | OFF | GoogleTest unit tests |
| `ENABLE_TENCOM_BENCHMARK` | OFF | Standalone tENCoM benchmark binary |
| `ENABLE_VCFBATCH_BENCHMARK` | OFF | VoronoiCFBatch benchmark binary |

---

## Build Variants

### Standard Release (default)

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j $(nproc)
```

### With Python Bindings

```bash
cmake .. -DBUILD_PYTHON_BINDINGS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . -j $(nproc)
```

### With Tests

```bash
cmake .. -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . -j $(nproc) && ctest --test-dir .
```

### CUDA GPU Acceleration

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_CUDA=ON
cmake --build . -j $(nproc)
```

Requires: CUDA Toolkit installed, NVIDIA GPU with compute capability ≥ 7.0 (Volta through Blackwell; Blackwell requires CUDA ≥ 12.6).

### Metal GPU Acceleration (macOS)

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_METAL=ON
cmake --build . -j $(nproc)
```

Requires: macOS with Xcode installed.

### ROCm/HIP GPU Acceleration (AMD)

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_ROCM=ON
cmake --build . -j $(nproc)
```

Requires: ROCm toolkit installed. Supports AMD Instinct MI100 (gfx908), MI200 (gfx90a), and MI300 (gfx942). The runtime dispatch priority is: CUDA > ROCm > Metal > AVX-512 > AVX-2 > OpenMP > scalar.

### HPC Deployment (maximum performance)

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DFLEXAIDS_USE_AVX512=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_CUDA=ON
cmake --build . -j $(nproc)
```

Build once on the target architecture for best `-march=native` optimization.

---

## Verifying the Installation

### C++ Binaries

```bash
./build/FlexAIDdS --version
./build/tENCoM --help
```

### Python Package

```bash
python -c "import flexaidds; print(flexaidds.__version__)"
```

### Test Suite

```bash
# C++ tests
ctest --test-dir build --output-on-failure

# Python tests
cd python && pytest tests/ -q
```

---

## Troubleshooting

### CMake cannot find Eigen3

```bash
# Ubuntu/Debian
sudo apt-get install libeigen3-dev

# macOS
brew install eigen

# Or set the path manually
cmake .. -DCMAKE_PREFIX_PATH=/path/to/eigen
```

### OpenMP not found on macOS

```bash
brew install libomp
cmake .. -DCMAKE_PREFIX_PATH="$(brew --prefix)"
```

### CUDA not detected

Ensure `nvcc` is in your PATH:
```bash
export PATH=/usr/local/cuda/bin:$PATH
cmake .. -DFLEXAIDS_USE_CUDA=ON
```

### ROCm not detected

Ensure the ROCm toolkit is installed and `hipcc` is in your PATH:
```bash
export PATH=/opt/rocm/bin:$PATH
cmake .. -DFLEXAIDS_USE_ROCM=ON
```

Verify your GPU is supported: `rocminfo | grep gfx`. Supported targets: gfx908, gfx90a, gfx942.

### Python bindings import error

If `import flexaidds` works but C++ functions are unavailable:
```bash
# Rebuild with bindings and set the output directory
cmake .. -DBUILD_PYTHON_BINDINGS=ON \
    -DCMAKE_LIBRARY_OUTPUT_DIRECTORY=$(pwd)/../python/flexaidds
cmake --build . --target _core
```

### Windows: LNK1104 / linker errors

Use the "x64 Native Tools Command Prompt for VS 2022" — not the default terminal. This ensures `cl.exe` and link paths are properly configured.

---

## Next Steps

- [User Guide](USERGUIDE.md) — full parameter reference and usage examples
- [Benchmarks](BENCHMARKS.md) — performance and accuracy data
- [Contributing](../CONTRIBUTING.md) — development setup and guidelines
