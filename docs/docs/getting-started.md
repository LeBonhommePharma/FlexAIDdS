# Getting Started

**Last actualized**: 2026-07-12

## Prerequisites

- C++26 compiler (GCC ≥ 14, Clang ≥ 18, Apple Clang ≥ 16 / Xcode 16, or MSVC ≥ 19.40)
- CMake ≥ 3.28
- Python ≥ 3.9 (3.11 recommended for the analysis package)
- Optional: Eigen3, OpenMP, CUDA toolkit, Metal (macOS), pybind11

Full platform notes: [`docs/INSTALLATION.md`](../INSTALLATION.md).

## Build the native engines

From the repository root:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
```

With tests:

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### Useful CMake options

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_TESTING` | OFF | GoogleTest suite |
| `BUILD_PYTHON_BINDINGS` | OFF | pybind11 `_core` extension |
| `FLEXAIDS_USE_AVX2` | ON | AVX2 SIMD |
| `FLEXAIDS_USE_AVX512` | OFF | AVX-512 SIMD |
| `FLEXAIDS_USE_CUDA` | OFF | CUDA evaluation |
| `FLEXAIDS_USE_METAL` | OFF | Metal (macOS) |
| `FLEXAIDS_USE_OPENMP` | ON | OpenMP |
| `FLEXAIDS_USE_EIGEN` | ON | Eigen3 linear algebra |
| `ENABLE_TENCOM_TOOL` | OFF | Standalone tENCoM CLI extras |

Primary executables after a successful build:

- `FlexAIDdS` — modern docking + entropy path
- `FlexAID` — legacy-compatible docking entry
- `tENCoM` — vibrational entropy tool (when enabled)

## Install the Python package

From the repository root (not from inside `python/` unless you adjust paths):

```bash
pip install -e ./python
```

Verify:

```bash
python -c "import flexaidds as fd; print(fd.__version__, fd.HAS_CORE_BINDINGS)"
python -m flexaidds --help
```

Run pure-Python tests:

```bash
PYTHONPATH=python python -m pytest python/tests/ -q
```

See [`docs/TESTING.md`](../TESTING.md) for the full test matrix.

## First docking runs

```bash
# Zero-config docking (defaults include flexibility + entropy temperature)
./build/FlexAIDdS receptor.pdb ligand.mol2

# Inspect an existing results directory
python -m flexaidds /path/to/results/ --top 5
```

Python high-level dock helper (requires a working engine binary on `PATH` or via `binary=`):

```python
import flexaidds as fd

population = fd.dock("receptor.pdb", "ligand.mol2", temperature=298.15)
```

Loading completed runs (no engine required):

```python
import flexaidds as fd

result = fd.load_results("path/to/results")
print(result)  # DockingResult with modes / poses / REMARK ledger fields
```

## Runtime data (matrices + definition files)

Production docking needs interaction matrices (`MC_*.dat`) and definition files (`AMINO.def`, …) next to the binary. The skill helper ensures them:

```bash
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
```

## Next steps

- [Configuration](configuration.md)
- [Python API](api/python.md)
- [Thermodynamics](thermodynamics.md)
- [User Guide](../USERGUIDE.md) (full CLI parameter reference)
