# flexaidds Python Package

Python bindings and analysis helpers for the FlexAID∆S thermodynamic core.

The current Python package now covers two complementary layers:

- `Thermodynamics`: free energy, mean energy, entropy, heat capacity, and energy dispersion from the canonical ensemble engine
- `StatMechEngine`: canonical ensemble builder and analyzer exposed through the compiled extension
- `load_results(...)`: a read-only loader for existing docking result directories, grouping PDB outputs into binding-mode ensembles from `REMARK` metadata and filename heuristics
- `PoseResult`, `BindingModeResult`, `DockingResult`: immutable Python-side data models for downstream analysis

Higher-level live docking orchestration is still intentionally staged behind the ongoing C++ integration work, but read-only ensemble inspection is now available immediately.

## Installation

### With pip (recommended for most users)

```bash
# Recommended today (package not yet on public PyPI):
pip install "git+https://github.com/LeBonhommePharma/FlexAIDdS.git#subdirectory=python"
python -m flexaidds --check-update       # or --self-update
# Console scripts after install: flexaidds, flexaidds-benchmark

# From the repo root (development)
pip install -e ./python

# After the first PyPI release:
pip install flexaidds
pip install --upgrade flexaidds
```

This installs the **Python analysis package** only (not the native `FlexAIDdS` docking binary — use Homebrew or a CMake build for that). The compiled acceleration (`_core`) is built if the required build dependencies (pybind11, Eigen headers, a C++ compiler) are present. Otherwise the package installs successfully in **pure-Python fallback mode** (full result loading, models, pure-Python thermodynamics, etc. still work). Set `FLEXAIDDS_SKIP_CORE=1` to force pure-Python.

**Native docking tools (separate):**
```bash
brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS
brew install lebonhommepharma/flexaidds/flexaidds
```

You can also do a non-editable install after building a wheel:
```bash
cd python
python -m pip wheel . -w /tmp/wheels
pip install /tmp/wheels/flexaidds-*.whl
```

### With conda

```bash
conda env create -f environment.yml
conda activate flexaidds
# (the env file does an editable pip install of the python package)
```

For a full conda package build (advanced):
```bash
conda-build conda
```

See the top-level `docs/INSTALLATION.md` for complete platform instructions (including when you also want the native `FlexAIDdS` / `tENCoM` executables).

## Quick example: thermodynamics

```python
import flexaidds as fd

engine = fd.StatMechEngine(temperature=300.0)
engine.add_sample(-7.5)
engine.add_sample(-6.0)
engine.add_sample(-5.5, multiplicity=2.0)
thermo = engine.compute()

print("F =", thermo.free_energy)
print("E =", thermo.mean_energy)
print("S =", thermo.entropy)
print("Cv =", thermo.heat_capacity)
print("sigma_E =", thermo.std_energy)
```

## Quick example: load docking results

```python
import flexaidds as fd

run = fd.load_results("path/to/output_dir")
print(run.n_modes)
print(run.binding_modes[0].best_cf)
print(run.binding_modes[0].free_energy)
print(run.to_records())
```

The loader scans PDB-like files recursively, parses `REMARK` lines such as `binding_mode`, `pose_rank`, `CF`, `free_energy`, `enthalpy`, `entropy`, and falls back to filename heuristics like `binding_mode_1_pose_2.pdb` when metadata are sparse.

## Command-line inspection

You can inspect a result directory directly without writing Python code:

```bash
# Summary table
python -m flexaidds path/to/output_dir

# JSON output
python -m flexaidds path/to/output_dir --json

# Export to CSV
python -m flexaidds path/to/output_dir --csv results.csv

# Show only top N binding modes
python -m flexaidds path/to/output_dir --top 5
```

This is useful for smoke-checking ensemble structure before pushing results into notebooks, pandas, or PyMOL-side tooling.

## Scope

This package is designed to stay additive and low-conflict while the deeper C++ refactor lands:

- no edits to `LIB/` sources for the read-only results layer
- no edits to root `CMakeLists.txt`
- no interference with ongoing `BindingMode` / `StatMechEngine` engine-side integration work

That makes it safe to evolve the Python analysis surface in parallel with the core thermodynamic plumbing.
