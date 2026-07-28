# Thermodynamic Source Map

## Core Canonical Ensemble

- `LIB/statmech.h/.cpp`
  - `statmech::kB_kcal`, `statmech::kB_SI`
  - `statmech::State`
  - `statmech::Thermodynamics` legacy output
  - `statmech::ThermodynamicBreakdown` explicit ledger
  - `statmech::StatMechEngine::compute()` canonical ensemble foundation
  - `statmech::StatMechEngine::compute_breakdown()` ledger adapter
  - `statmech::StatMechEngine::boltzmann_weights()` normalized pose probabilities
  - `statmech::StatMechEngine::boltzmann_pmf()` single-window Boltzmann PMF; not true multi-window WHAM

## Binding Mode Energy And Ranking

- `LIB/BindingMode.h/.cpp`
  - `Pose::CF` contact-function effective scoring energy
  - `Pose::receptor_strain`
  - `Pose::total_energy() = CF + receptor_strain`
  - `BindingMode::rebuild_engine()` populates `StatMechEngine` from `Pose::total_energy()`
  - `BindingMode::compute_energy()` legacy ranking energy with vibrational/NATURaL corrections
  - `BindingMode::get_thermodynamics()` legacy compatibility struct
  - `BindingMode::get_thermodynamic_breakdown()` explicit configurational/correction ledger
  - `BindingPopulation::EnergyComparator` sorts by cached legacy energy; unchanged by the ledger

## Grand Canonical And Selectivity

- `LIB/GrandPartitionFunction.h/.cpp`
  - `GrandPartitionFunction::add_ligand()`
  - `GrandPartitionFunction::log_Xi()`
  - `GrandPartitionFunction::binding_probability()`
  - `GrandPartitionFunction::empty_probability()`
  - `GrandPartitionFunction::selectivity()` apparent concentration-aware selectivity
  - `GrandPartitionFunction::log_intrinsic_selectivity()` concentration-independent selectivity
  - `TargetServer` / sessions carry conc_M (P3); per-receptor Ξ computation from ensemble log_Z
- Emission (P3+): DatasetRunner emits `grand_summary` (dict) → `*_grand_summary.csv` (ligand,log_Z,conc_M,log_Xi,p_bind) + JSON; C++ --conc wires default_conc_M; competition_example.yaml provides per-ligand conc_M
- Python: `flexaidds.grand_canonical.compute_grand_partition`, `_PyGrandPartitionFunction`, `LigandRank`; roundtrips in test_grand_canonical + grand_synthetic fixtures (exact 3L cases verified)
- MultiSiteGPF: product of per-site GPFs + cooperativity (tested in test_multi_site_gpf.cpp); non-breaking, documented for future multi-cleft integration (no auto-wiring from CavityDetect yet)

## Vibrational, NATURaL, GIST, H-bond, And Strain Terms

- `LIB/encom.h`, `LIB/encom.cpp`, `LIB/tENCoM/*`: vibrational entropy and model-scale corrections.
- `LIB/ShannonThermoStack/ShannonThermoStack.h/.cpp`: combined configurational/vibrational stack.
- `LIB/NATURaL/*`, `FA_Global::natural_deltaG` in `LIB/flexaid.h`: NATURaL correction.
- `LIB/GISTEvaluator.h/.cpp`, `LIB/GISTGrid.h`: GIST/desolvation scoring inputs.
- `LIB/HBondEvaluator.*`, `LIB/hbond_potential.h`: directional H-bond scoring terms.
- `Pose::receptor_strain` in `LIB/BindingMode.h`: receptor conformer strain penalty.

## Output Surfaces

- PDB/REMARK output: `LIB/BindingMode.cpp`, `LIB/FOPTICS.cpp`, `LIB/cluster.cpp`, `LIB/BinarySnapshot.cpp`.
- Python dataclasses and JSON/CSV: `python/flexaidds/models.py`.
- Python result parsing: `python/flexaidds/results.py`, `python/flexaidds/io.py`.
- PyMOL plugin/reporting: `pymol_plugin/gui.py`, `pymol_plugin/results_adapter.py`, `pymol_plugin/visualization.py`.

## Python And Bindings

- `python/flexaidds/thermodynamics.py`
  - `Thermodynamics` legacy dataclass
  - `ThermodynamicBreakdown` ledger dataclass
  - `_PyStatMechEngine` pure-Python fallback
  - `StatMechEngine.compute_breakdown()` Python API
- `python/flexaidds/_core.cpp`: standalone pybind module.
- `python/bindings/core_bindings.cpp`: CMake/full-GA pybind module.

## Ambiguities Flagged

- Legacy `Thermodynamics.free_energy` can mean pure configurational free energy in `StatMechEngine`, but corrected total free energy in `BindingMode::get_thermodynamics()`.
- README/docs historically mention WHAM; current implementation is single-window Boltzmann PMF unless a real multi-window WHAM layer is added.
- `CF` and `H_eff` are effective scoring-energy quantities, not physical calorimetric enthalpies without calibration.
- Existing PyMOL labels use `DeltaG/free_energy`; reporting should distinguish legacy totals from `ThermodynamicBreakdown` fields when available.

## Core Versus Experimental Status

- Core: canonical ensemble in `StatMechEngine`, `ThermodynamicBreakdown` after tests pass.
- Core diagnostic: component averages once component completeness metadata exists.
- Diagnostic only: compensation metrics.
- Experimental: joint receptor-ligand ensemble, affinity calibration until benchmarked, temperature scan, model-derived DeltaCp, cleft annotation.
- Visualization only: PyMOL/reporting utilities.
