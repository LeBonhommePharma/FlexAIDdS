# DatasetRunner: Thermo Focus & Reproducible Bulletproof Implementation

## Critical Direction (from author Le Bonhomme Pharma)
- **Forget autonomous blind docking**. It is NOT the competition or target.
- Autonomous blind docking (full blind, no oracle, chasing docking power on blind sets like TIER-3) has **never** been the goal.
- **Target**: Thermodynamic affinity prediction on ITC data.
  - Accurate ΔG (scoring_power_pearson_r, rmse vs experimental).
  - Full decomposition: deltaH, deltaS (from ITC ground truth).
  - Entropy contribution: entropy_rescue_rate (Shannon Energy Collapse), how -TΔS term rescues correct modes.
  - Entropy/enthalpy index for characterization and re-ranking.
  - Ablation studies on entropy components (Shannon, tENCoM, vib, etc.).
- Use ITC-187, BindingDB_ITC_subset, SCORPIO_ITC as primary.
- For Astex: use native/nonnative for pose quality in thermo context, but evaluate on thermo metrics + entropy index, not docking power.
- Use external tools (Vina, rDock, Boltz-2) + re-rank with **existing entropy metrics** (Shannon energy collapse, tENCoM, thermodynamic scoring).
- PoseBusters: **only for validity filter** — a pose is "good for thermo" if RMSD ≤ 2.0 Å **AND** passes all PB checks. Not as autonomous success rate.
- Reproducible: always record git_sha, binary_sha, host, full provenance, temperature (298K for ITC), exact commits (e.g. legacy 94.1% self-dock at 8196829f35a2bf065919ccd1508f62f00059895d; recent cross-dock at aef5c4b77c9167d9693c4d59f14d2e717eb82cdf or 525cf811...).
- Bulletproof: error isolation, resume, per-target artifacts, CI-ready, type hints, logging, validation. No batching commits without push. Verify with actual runs.

## New Supporting Module: benchmarks/astex_entropy/
Clean, MacBook-Pro native benchmark (drop-in ready).

**Minimal deps only**: rdkit, pandas, pyyaml, typer, matplotlib, seaborn, openbabel (+ posebusters for PB).

**Features**:
- Native + non_native using existing astex_diverse.yaml + astex_nonnative.yaml.
- Runs Vina, rDock, Boltz-2 to generate poses.
- For **every pose**: run entropy metrics → re-rank.
- Thermo focus: deltaH/deltaS decomp, entropy/enthalpy index.
- Success filter: RMSD≤2.0 + PoseBusters (all_passed).
- Typer CLI: data_prep, run --mode native|non_native, rescore --poses_from vina|rdock|boltz.
- Production but clean/simple (no gates/abstractions).

**HOW-TO (MacBook Pro)**:
```bash
# setup
conda create -n astex_entropy python=3.11
conda activate ...
conda install -c conda-forge rdkit openbabel pandas pyyaml matplotlib seaborn
pip install "typer[all]" posebusters

# (tools in PATH: vina, rdock, boltz; entropy bins or flexaidds python)

python -m benchmarks.astex_entropy.cli data_prep --out-dir benchmarks/astex_entropy/data
python -m benchmarks.astex_entropy.cli run --mode native --out-dir results/astex_entropy/native
python -m benchmarks.astex_entropy.cli rescore --poses_from vina --mode native --out-dir results/astex_entropy/native
# similar for non_native, other tools
```

See full README in benchmarks/astex_entropy/README.md for details, thermo reporting, ITC hooks.

## Integration with DatasetRunner (this module)
- Use `python/flexaidds/dataset_runner/` for **internal FlexAIDdS** runs on thermo datasets.
- Use `benchmarks/astex_entropy/` for **external tools + entropy re-rank + PB filter**.
- To make unified automated reproducible:
  - Discover from benchmarks/datasets/*.yaml (already done; thermo yamls updated to scoring_power + entropy_rescue only, no docking_power).
  - For "other" / hybrid: support loading external poses and rescoring (extend run_dataset or add rescore path).
  - Always: git_sha, temperature=298.0 for ITC, full provenance in reports.
  - Metrics: prioritize scoring_power_*, entropy_rescue_rate, add thermo decomp if exp data available.
  - Bulletproof extensions: per-entry resume (already), error continue (added), validate configs.
  - Repro: pin to commits, record binary, support --dry-run, --resume.
  - External integration: call into astex_entropy runners/entropy for Vina etc. poses, then feed to PoseScore.

## Key Provenance (for reproducible runs)
- Legacy 94.1% Astex self-dock (old protocol): 8196829f35a2bf065919ccd1508f62f00059895d
- Recent cross-dock provenance: aef5c4b77c9167d9693c4d59f14d2e717eb82cdf (or similar 525cf811...)
- Current workspace for thermo work: 938b2a3c (on feat/thermoaffinity-suite-v21)
- v48 etc. are internal June 2026 campaign labels (not software versions); success rates lower on modern protocols (e.g. v43 ~81% native).
- Recent external pairs: ~0-1% on non-native (as expected).

## Updated Thermo Datasets (in benchmarks/datasets/ and python copy)
- itc187.yaml, bindingdb_itc.yaml, scorpio.yaml: focus scoring_power_* + entropy_rescue_rate. Baselines adjusted for realistic thermo.
- Astex yamls: keep for native/nonnative pose sources, but use in thermo context.

## Making DatasetRunner Bulletproof & Automated
- Error isolation: wrap runs, continue, mark failed (implemented).
- Reproducible artifacts: per-target json, full report with sha.
- Support other datasets: add YAML + implement loader/runner if external.
- Integrate astex_entropy: option to "rescore_external" using the module's logic.
- Thermo ledger: temperature param, exact for ITC.
- Validation: check metrics, baselines, data presence.
- CI: use dry-run, tier1.
- No autonomous: configs for thermo sets exclude docking_power; docs emphasize ITC focus.

See also:
- benchmarks/astex_entropy/ (new clean external thermo benchmark)
- previous ThermoAffinitySuite updates (PR #239 direction)
- BENCHMARK_STANDARD.md (tiers, but de-emphasize blind for thermo)
- python/flexaidds/ (entropy impls: shannon, tENCoM, thermo)

For the worker: extend to orchestrate both internal + external rescore in one reproducible pipeline, always with provenance, focus thermo metrics.

Run with actual execution to verify before claims.
