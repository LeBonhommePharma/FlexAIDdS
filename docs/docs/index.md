# FlexAIDΔS Documentation

**FlexAIDΔS** (FlexAID with ΔS entropy) is an entropy-driven molecular docking engine that combines genetic algorithms with statistical mechanics.

> **Source of truth for agent workflow**: repository root [`AGENTS.md`](../../AGENTS.md).  
> **Last documentation actualization**: 2026-07-12.

## What this site covers

This MkDocs tree (`docs/docs/`) is the navigable technical reference for:

- Getting started (build + Python package)
- Configuration and scoring
- Genetic algorithm behavior
- Python API surface
- Benchmark suite docs
- Thermodynamics terminology (scoring proxy vs ledger)
- Testing and support boundaries

Repository-root Markdown under `docs/` (User Guide, Installation, VALIDATED_CAPABILITIES, etc.) remains the full product documentation set. Prefer both when writing papers or validation packages.

## Key capabilities

- **GA docking** with Voronoi contact-function (**CF**) scoring as the **search ranking proxy**
- **Canonical ensemble thermodynamics** via `StatMechEngine` (logZ, G/F, H_eff, S, Cv)
- **Shannon configurational entropy** with hardware dispatch (CUDA / Metal / AVX / OpenMP / scalar)
- **Vibrational entropy** via ENCoM / tENCoM
- **Ligand flexibility**: torsions, ring conformers, chiral centers
- **Python package** `flexaidds` for results I/O, thermo analysis, DatasetRunner, ML rescoring bridge
- **Cross-platform** Linux / macOS / Windows (see Support Matrix)

## Scientific guardrail (read first)

| Concept | Role |
|---------|------|
| CF / Voronoi contact function | **Scoring proxy** used during GA search and many pose ranks |
| Thermodynamic ledger (F, H, −TS, Cv) | Computed on the **ensemble** after search (StatMech / BindingMode) |
| Experimental ΔG / Kd | Only with calibration + validated protocol — never equate raw CF to affinity |

## Quick start

```bash
# Native engines
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/FlexAIDdS receptor.pdb ligand.mol2

# Python analysis package
pip install -e ./python
python -m flexaidds /path/to/results/
```

## Documentation map

| Page | Contents |
|------|----------|
| [Getting Started](getting-started.md) | Prerequisites, build options, first runs |
| [Configuration](configuration.md) | JSON config reference |
| [Scoring Overview](scoring/overview.md) | CF, H-bond, GIST, components |
| [GA Overview](ga/overview.md) | Fitness models and GA parameters |
| [Python API](api/python.md) | Package modules and examples |
| [Thermodynamics](thermodynamics.md) | Ensemble quantities and support class |
| [Architecture](architecture.md) | Pipeline layers and directories |
| [Testing](testing.md) | How to run and interpret tests |
| [Support Boundary](support-boundary.md) | Core 1.0 vs experimental |

## External product docs (repo root `docs/`)

- [Installation](../INSTALLATION.md)
- [User Guide](../USERGUIDE.md)
- [Support Matrix](../SUPPORT_MATRIX.md)
- [Thermodynamics (full)](../thermodynamics.md)
- [Testing Guide](../TESTING.md)
- [Test Coverage Analysis](../TEST_COVERAGE_ANALYSIS.md)
- [Validated Capabilities](../VALIDATED_CAPABILITIES.md)
- [Known Limitations](../KNOWN_LIMITATIONS.md)
