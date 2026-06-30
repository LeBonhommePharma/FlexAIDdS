# Astex Entropy Benchmark

Clean, focused benchmark for thermodynamic affinity on the Astex Diverse Set (native and non-native modes) using your existing data.

**Focus**: ITC-style thermodynamic decomposition (ΔG, ΔH, −TΔS), entropy/enthalpy index, and re-ranking poses from external tools using your entropy metrics (Shannon energy collapse, tENCoM, thermodynamic scoring). Not for autonomous blind docking power.

Success definition (for reference poses):
- RMSD ≤ 2.0 Å **and** passes PoseBusters validation.

Re-ranking uses your entropy metrics to produce better thermo-aware ranking and decomposition.

## Minimal Dependencies (MacBook Pro)

```bash
# Recommended: use conda/mamba for rdkit + openbabel on Apple Silicon
conda create -n astex_entropy python=3.11
conda activate astex_entropy

conda install -c conda-forge rdkit openbabel pandas pyyaml typer matplotlib seaborn

# PoseBusters (for validation)
pip install posebusters

# Optional for better openbabel python
pip install openbabel-wheel
```

External tools you must have in PATH (or adjust paths in code):
- `vina` (AutoDock Vina)
- `rbdock`, `rbcavity` (rDock)
- `boltz` (Boltz-2, `pip install boltz` or from source; supports protein-ligand)

Your FlexAIDdS / tENCoM entropy tools (for existing metrics):
- Built `tENCoM` (or python/flexaidds with bindings for shannon/thermo stack)
- Add to PATH or set FLEXAIDDS_BIN in env.

Data already in repo:
- `benchmarks/datasets/astex_diverse.yaml`
- `benchmarks/datasets/astex_nonnative.yaml`
- `benchmarks/astex_diverse/...`
- `benchmarks/astex_nonnative/...`

## Quick Start on MacBook Pro

```bash
cd /path/to/your/FlexAIDdS

# 1. Data prep (creates inputs for Vina/rDock/Boltz from the yamls)
python -m benchmarks.astex_entropy.cli data_prep --out-dir benchmarks/astex_entropy/data

# 2. Run external tools (native mode example)
python -m benchmarks.astex_entropy.cli run --mode native --out-dir results/astex_entropy/native

# 3. Rescore poses from a tool with your entropy metrics + ITC/ΔH/ΔS focus + entropy/enthalpy index
python -m benchmarks.astex_entropy.cli rescore --poses_from vina --mode native --out-dir results/astex_entropy/native --tool vina

# Non-native
python -m benchmarks.astex_entropy.cli run --mode non_native --out-dir results/astex_entropy/non_native
python -m benchmarks.astex_entropy.cli rescore --poses_from rdock --mode non_native --out-dir results/astex_entropy/non_native --tool rdock
```

See `cli.py --help` for all options.

## What it does

- Loads targets from your yamls (native holo or non-native pairs).
- Prepares receptor/ligand inputs for Vina (pdbqt + box), rDock (.prm + .mol2), Boltz-2 (yaml input for complex prediction).
- Runs the tools to generate poses.
- For each pose:
  - Compute RMSD to reference (crystal ligand pose).
  - Run PoseBusters.
  - Run your entropy metrics: Shannon collapse, tENCoM (vibrational), thermodynamic scoring.
  - Compute ΔH / ΔS decomposition (or proxies via your stack).
  - Compute entropy/enthalpy index (e.g. sign(TΔS) vs ΔH, or |TΔS| / |ΔG| fraction).
- Re-ranks the poses using a thermo score (e.g. predicted ΔG or entropy-weighted).
- Reports per-tool, per-mode:
  - "PB-success" rate (RMSD≤2 + PB pass) for top-1 after re-rank vs original tool ranking.
  - Average predicted ΔH, −TΔS, entropy/enthalpy index on successful vs all poses.
  - Plots: decomposition, index vs RMSD, success before/after re-rank.

This lets you see how your entropy/thermo re-ranking affects quality and thermodynamic character of poses from Vina/rDock/Boltz on Astex.

## Directory Layout (after running)

```
benchmarks/astex_entropy/
  cli.py
  data_prep.py
  runners.py
  entropy.py          # wrappers + ΔH/ΔS + index
  evaluate.py
  README.md
  requirements.txt

results/astex_entropy/
  native/
    <pdb>/
      vina/poses.sdf
      rdock/...
      boltz/...
      rescore_vina/  # re-ranked + metrics json
  non_native/...
```

## Notes for MacBook Pro (M-series)

- Use conda for rdkit/openbabel to get arm64 builds.
- Vina: `brew install vina` or download binary.
- rDock: compile or use conda `conda install -c bioconda rdock`.
- Boltz: `pip install boltz` (may need torch).
- For entropy: build the project once (`cmake -B build ...`), then `export PATH=$PWD/build:$PATH` for tENCoM etc. Or use python -m flexaidds if bindings installed.
- Run with `caffeinate` for long runs: `caffeinate python -m ...`
- Results are deterministic enough for comparison; use --seed if available in tools.

## Extending for ITC

If you have ITC experimental ΔH/ΔS/ΔG for some targets (via your fetch_itc_data or calibrate), the evaluate can load them and compute correlations of predicted decomposition vs experimental for the re-ranked poses. See comments in entropy.py / evaluate.py.

This is ready to drop in. Run data_prep first.

Good luck with the thermo-focused benchmarking this morning!
