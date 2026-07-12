# Python API Reference

**Package**: `flexaidds`  
**Last actualized**: 2026-07-12  

Install: `pip install -e ./python` from the repository root.

The package works in **pure-Python mode** without the optional C++ `_core` extension. When `_core` is built (`BUILD_PYTHON_BINDINGS=ON` or a wheel that includes it), `HAS_CORE_BINDINGS` is `True` and accelerated engines replace pure-Python fallbacks for matching symbols.

```python
import flexaidds as fd
print(fd.__version__, fd.HAS_CORE_BINDINGS)
```

## High-level docking and results

### `flexaidds.dock(receptor, ligand, **kwargs)`

Runs the native engine and returns a `BindingPopulation`.

```python
pop = flexaidds.dock(
    "receptor.pdb",
    "ligand.mol2",
    binding_site="auto",
    compute_entropy=True,
    temperature=300.0,
    timeout=3600,
)
```

### `flexaidds.load_results(path)`

Parses a results directory into `DockingResult` (poses, modes, REMARK ledger fields).

```python
from flexaidds import load_results
run = load_results("results/my_campaign/target_state")
```

### Core data models

| Type | Module | Role |
|------|--------|------|
| `PoseResult` | `models` | Single pose + scores / REMARKs |
| `BindingModeResult` | `models` | Cluster / mode aggregate |
| `DockingResult` | `models` | Full run container |
| `Pose` / `BindingMode` / `BindingPopulation` | `docking` | In-memory docking objects |
| `Atom` / `PDBStructure` | `io` | Structure I/O helpers |

## Statistical mechanics

### `StatMechEngine` / `Thermodynamics` / `ThermodynamicBreakdown`

Always importable. Prefer C++ implementations when `HAS_CORE_BINDINGS`.

```python
from flexaidds import StatMechEngine, kB_kcal

engine = StatMechEngine(300.0)
engine.add_sample(-12.0)
engine.add_sample(-10.5, weight=2.0)
thermo = engine.compute()  # free_energy, mean_energy, entropy, ...
```

Utilities:

- `deltaG_standard_to_Kd_M` / `Kd_M_to_deltaG_standard` — standard-state conversion (safe only when calibrated)
- `gibbs_helmholtz_dG`, `kirchhoff_dH`, `kirchhoff_dS` — temperature helpers
- `TemperatureScanPoint`, `DeltaCpFit` — model-derived scan / fit types

### Audit schema (`schemas.thermo_audit`)

Public entropy.help / A1.1 types:

- `ThermodynamicOutput` / `ThermodynamicOutputDC`
- `TotalSampledPartitionFunction` / `TotalSampledPartitionFunctionDC`
- `Provenance` / `ProvenanceDC`
- `make_total_sampled_output(...)` — factory that enforces F = −kT logZ identities

```python
from flexaidds import make_total_sampled_output

out = make_total_sampled_output(
    logZ=-40.0,
    mean_energy=-15.0,
    temperature_K=298.15,
    n_samples=1000,
    git_sha="abc123",
    timestamp="2026-07-12T00:00:00Z",
    gate_results={"gate6_crosscheck": {"passed": True}},
)
out.validate()
```

## Vibrational and torsional entropy

| Symbol | Module | Notes |
|--------|--------|-------|
| `ENCoMEngine`, `NormalMode`, `VibrationalEntropy` | `encom` | Elastic network vibrational entropy |
| `TorsionalENM`, `compute_shannon_entropy`, … | `tencm` | Torsional ENM + Shannon stack helpers |
| `parse_tencom_pdb`, `parse_tencom_json` | `tencom_results` | tENCoM result parsers |
| `DiFTEngine`, `score_torsional`, … | `dift` | Dihedral Fourier thermodynamics |

## Scoring, matrices, and ML bridge

| Symbol | Role |
|--------|------|
| `EnergyMatrix`, `encode_256_type`, `parse_dat_file` | Interaction matrix I/O |
| `VoronoiGraphExtractor`, `ShannonProfileExtractor`, `FeatureBuilder`, `MLRescorer` | ML rescoring feature bridge |
| `GAOptimizer` | GA hyperparameter search (does not alter production ranking by default) |

## Benchmarking and datasets

| Symbol | Role |
|--------|------|
| `DatasetRunner`, `DatasetConfig`, … | Campaign orchestration (`flexaidds.dataset_runner`) |
| `dataset_adapters.create_adapter` | PDBbind / ITC-187 / MOAD / BindingDB / ChEMBL / DUD-E / DEKOIS |
| `normalize_affinity` | Convert Kd/Ki/IC50/pKd/… → ΔG (kcal/mol) for training pipelines |
| `run_benchmark`, `BenchmarkSystem`, … | Method comparison helpers (`benchmark`) |

## Diagnostics and reporting

| Symbol | Role |
|--------|------|
| `diagnostics.ranking_bias` | CF vs RMSD ranking pathology metrics (post-hoc on REMARK CF) |
| `generate_pymol_script`, `generate_markdown_report`, `write_all_reports` | Reporting from `DockingResult` JSON only |
| `prepare_publication_figures`, NRDD cover helpers | Post-hoc figures (Gate 6 gated; never affects scores) |
| `truncate_chain.write_extended_ca_chain` | Synthetic nascent-chain PDBs for co-translational pipelines |

## CLI

```bash
python -m flexaidds <results_dir> [--json|--csv|--top N]
```

Updater helpers: `check_for_updates` (package self-update path).

## Fallback types (`_fallback_types`)

When `_core` is absent, pure-Python stubs keep imports working:

`WHAMBin`, `TIPoint`, `Replica`, `State`, `BoltzmannLUT`, `TemperatureScanPoint`, `DeltaCpFit`.

These are **API-compatible placeholders**, not full C++ numerics replacements for every algorithm.

## Terminology in API docs

- Methods that rank on CF report **contact-function scores**, not free energies.
- `free_energy` / `G_config` fields from ensemble analysis are **ensemble-derived** quantities.
- Do not publish raw docking scores as experimental Kd without the calibration path and validation docs.

## Tests covering this surface

See `python/tests/` and [`docs/TESTING.md`](../TESTING.md). Representative modules: `test_thermodynamics*`, `test_thermo_schema`, `test_results*`, `test_dataset_adapters`, `test_fallback_types`, `test_ranking_bias`.
