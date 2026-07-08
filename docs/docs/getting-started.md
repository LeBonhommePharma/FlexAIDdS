# Getting Started

## Prerequisites

- C++26 compiler (GCC ≥ 14, Clang ≥ 18, Apple Clang ≥ 16 / Xcode 16)
- CMake ≥ 3.28
- Python ≥ 3.9

## Build

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . -j $(nproc)
```

### Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `BUILD_TESTING` | OFF | Enable GoogleTest suite |
| `FLEXAIDS_USE_AVX512` | OFF | AVX-512 SIMD acceleration |
| `FLEXAIDS_USE_CUDA` | OFF | CUDA GPU evaluation |
| `FLEXAIDS_USE_METAL` | OFF | Metal GPU (macOS) |
| `FLEXAIDS_USE_MPI` | OFF | MPI distributed docking |

## Python Package

```bash
cd python
pip install -e ./python
pytest tests/
```

## Usage

```bash
# CLI
./FlexAIDdS receptor.pdb ligand.mol2

# Python
import flexaidds
result = flexaidds.dock("receptor.pdb", "ligand.mol2")
```
