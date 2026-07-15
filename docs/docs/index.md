# FlexAIDΔS Documentation

**FlexAIDΔS** (FlexAID with ΔS Entropy) is an entropy-driven molecular docking engine combining genetic algorithms with statistical mechanics thermodynamics.

## Key Features

- **Entropy-driven scoring**: Computes full thermodynamic properties (ΔG, ΔH, -TΔS, Cv) via partition functions
- **Voronoi contact function**: Soft scoring based on surface complementarity
- **Shannon configurational entropy**: Hardware-accelerated (CUDA/Metal/AVX-512)
- **Vibrational entropy**: ENCoM elastic network model for receptor flexibility
- **Angular H-bond potential**: Gaussian bell function for directional hydrogen bonds
- **GIST desolvation**: Grid-based water thermodynamics from MD trajectories
- **GA diversity monitoring**: Entropy collapse detection with adaptive mutation
- **ML rescoring bridge**: Feature extraction for hybrid physics/ML scoring

## Quick Start

```bash
# Build
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Python
pip install -e ./python
python -m flexaidds /path/to/results/
```

## Architecture

See the [Getting Started](getting-started.md) guide for installation and usage.
