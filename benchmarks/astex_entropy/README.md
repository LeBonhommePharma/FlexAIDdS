# Astex Entropy Benchmark

Clean Astex Diverse Set benchmark for LP's MacBook Pro 14" M3 Pro / 18 GB RAM:

- Modes: `native` and `non_native`
- Pose generators: FlexAIDdS, AutoDock Vina, rDock, Boltz-2
- Rescoring: Shannon energy collapse, tENCoM, and thermodynamic `G_bind`
- Success definition: `RMSD <= 2.0 A` and PoseBusters passes
- Required validators: PoseBusters and tENCoM/Eigen. No PoseBusters, no
  tENCoM, no benchmark claim.

This module is intentionally small. It writes manifests, tool inputs, poses,
rescored CSVs, and plots under:

```text
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy
```

## Install

From the FlexAIDdS checkout:

```bash
cd /Users/lp.more/Projects/FlexAIDdS

/opt/homebrew/bin/python3.12 -m venv .venv-astex-entropy
source .venv-astex-entropy/bin/activate
python -m pip install --upgrade pip

# Python dependencies only: rdkit, pandas, pyyaml, typer, matplotlib, seaborn, openbabel.
python -m pip install -r benchmarks/astex_entropy/requirements.txt
```

The default `config.yaml` is pinned to local tool paths under the checkout:

```text
/Users/lp.more/Projects/FlexAIDdS/build_lto/benchmark_datasets
/Users/lp.more/Projects/FlexAIDdS/.tools/bin/vina
/Users/lp.more/Projects/FlexAIDdS/.tools/rdock/bin/rbcavity
/Users/lp.more/Projects/FlexAIDdS/.tools/rdock/bin/rbdock
/Users/lp.more/Projects/FlexAIDdS/build_lto/cavity_detect_cli
/Users/lp.more/Projects/Get_Cleft/Get_Cleft
/Users/lp.more/Projects/FlexAIDdS/.venv-posebusters/bin/bust
/Users/lp.more/Projects/FlexAIDdS/build_lto/tencom_entropy_diff
/Users/lp.more/Projects/FlexAIDdS/.venv-boltz/bin/boltz
/opt/homebrew/bin/obabel
```

Install the external command-line tools:

```bash
brew install open-babel

cd /Users/lp.more/Projects/FlexAIDdS

mkdir -p .tools/bin
curl -L --fail \
  --output .tools/bin/vina \
  https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_mac_aarch64
chmod +x .tools/bin/vina

# rDock should expose these two binaries:
test -x .tools/rdock/bin/rbcavity
test -x .tools/rdock/bin/rbdock

/opt/homebrew/bin/python3.12 -m venv .venv-posebusters
.venv-posebusters/bin/python -m pip install --upgrade pip
.venv-posebusters/bin/python -m pip install posebusters

/opt/homebrew/bin/python3.12 -m venv .venv-boltz
.venv-boltz/bin/python -m pip install --upgrade pip
.venv-boltz/bin/python -m pip install boltz

/opt/homebrew/bin/obabel -V
/Users/lp.more/Projects/FlexAIDdS/build_lto/benchmark_datasets --help
/Users/lp.more/Projects/FlexAIDdS/build_lto/cavity_detect_cli --help
/Users/lp.more/Projects/FlexAIDdS/build_lto/tencom_entropy_diff --help || true
/Users/lp.more/Projects/FlexAIDdS/.tools/bin/vina --help
/Users/lp.more/Projects/FlexAIDdS/.venv-posebusters/bin/bust --help
NUMBA_CACHE_DIR="/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/cache/numba" \
  /Users/lp.more/Projects/FlexAIDdS/.venv-boltz/bin/python \
  /Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_entropy/boltz_cli.py predict --help
```

Boltz-2 is configured for the 18 GB MacBook Pro constraint: CPU accelerator,
one device, one dataloader worker, one parallel sample, and one diffusion
sample by default. The affinity head is also capped to one diffusion sample so
`properties: affinity` does not silently fall back to Boltz's heavier defaults.
It uses the existing `~/.boltz` model/CCD cache. The local `boltz_cli.py`
wrapper disables Numba disk caching inside the Boltz process to avoid the
cache-locator crash in this venv package path.

## One-command Benchmark

Use the orchestrator for normal benchmarking. It prepares the shared manifest,
runs the selected pose generators, rescoring each generated pose set with
Shannon collapse, tENCoM, thermodynamic `G_bind`, RMSD, and PoseBusters.
The orchestrator preflights PoseBusters and `tencom_entropy_diff` before doing
anything else; both are mandatory.

## Methodology

The benchmark protocol from now on is:

1. Prepare native targets from `astex_diverse.yaml` and non-native targets from
   `astex_nonnative.yaml`.
2. Define the FlexAIDdS binding site with FlexAIDdS's native `CavityDetector`,
   the in-tree SURFNET/GetCleft-equivalent cavity detector. The primary protocol
   selects the detected cleft nearest the original target ligand centroid and
   writes that sphere file as `cleft_sphere_file`.
3. Keep FlexAIDdS pose-blinded and no-native-seed: `include_oracle_site: false`,
   one `benchmark_datasets` worker, one FlexAIDdS child, `restarts: 1`, and
   `parallel_restarts: 0`.
4. Use CPU/OpenMP for FlexAIDdS GA chromosome fitness, even from a Metal-enabled
   build. The current Metal evaluator is translation-only over raw genes and is
   not parity-safe with the CPU `ic2cf`/`vcfunction` path.
5. Run the selected pose generators: FlexAIDdS, Vina, rDock, and/or Boltz-2.
6. Rescore every emitted pose with Shannon energy collapse, tENCoM/Eigen, and
   thermodynamic `G_bind`.
7. Validate every pose with PoseBusters.
8. Report success only from `success_pb`, which means `RMSD <= 2.0 A` and
   PoseBusters passed. RMSD-only success is not a benchmark result.

The primary result is the main occupied cavity because it matches the ligand
site definition used by the original Astex-style redocking benchmark. A
multi-cavity run is useful as a blind-search ablation, not as the primary
headline number. To run that diagnostic, copy `config.yaml` and change:

```yaml
tools:
  flexaidds:
    cavity:
      detector: native
      protocol: multi_cavity
      top_cavities: 3
```

Multi-cavity child targets are named like `1G9V__clf1`, `1G9V__clf2`, and
`1G9V__clf3`. Rescoring groups them back under the original target through
`analysis_target_id` so success rates remain target-level, not cavity-level.

Native one-target smoke:

```bash
python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds --max-targets 1 --skip-rescore
```

Native full head-to-head:

```bash
python -m benchmarks.astex_entropy orchestrate --mode native --tools flexaidds,vina,rdock,boltz
```

Non-native smoke:

```bash
python -m benchmarks.astex_entropy orchestrate --mode non_native --tools flexaidds --max-targets 3 --download-missing --skip-rescore
```

Run native and non-native in one pass:

```bash
python -m benchmarks.astex_entropy orchestrate --mode all --tools flexaidds,vina,rdock,boltz --download-missing
```

The standalone script form is equivalent:

```bash
python -m benchmarks.astex_entropy.orchestrate --mode native --tools flexaidds,vina
```

Each orchestrator run writes:

```text
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/orchestrator_runs/<run_id>/orchestrator_summary.json
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/orchestrator_runs/<run_id>/orchestrator_summary.md
```

## Prepare Data

Native Astex uses the existing prepared checkout layout:

`/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_diverse/astex_diverse/<PDB>/<PDB>_apo.pdb`

Run a fast native smoke:

```bash
python -m benchmarks.astex_entropy data_prep --mode native --max-targets 5
```

Prepare all native targets found in `astex_diverse_set.csv`:

```bash
python -m benchmarks.astex_entropy data_prep --mode native
```

Prepare non-native pairs. Add `--download-missing` to fetch missing PDB files
from RCSB into the iCloud work dir under `data/non_native/`.

```bash
python -m benchmarks.astex_entropy data_prep --mode non_native --max-targets 3 --download-missing
```

For a one-target wiring smoke:

```bash
python -m benchmarks.astex_entropy data_prep --mode native --max-targets 1 --force
python -m benchmarks.astex_entropy run --mode native --tools flexaidds
python -m benchmarks.astex_entropy rescore --mode native --poses_from flexaidds
python -m benchmarks.astex_entropy run --mode native --tools vina
python -m benchmarks.astex_entropy rescore --mode native --poses_from vina
python -m benchmarks.astex_entropy run --mode native --tools rdock
python -m benchmarks.astex_entropy rescore --mode native --poses_from rdock
```

## Run Pose Generators

Native, all four tools:

```bash
python -m benchmarks.astex_entropy run --mode native
```

FlexAIDdS only, with M3 Pro defaults from `config.yaml`:

```bash
python -m benchmarks.astex_entropy run --mode native --tools flexaidds
```

The generated command uses:

```text
FLEXAIDDS_RESTARTS=1
FLEXAIDDS_PARALLEL_RESTARTS=0
--threads 1 --omp-threads 6 --mode autonomous --ga-generations 500 --ga-population 1000
```

Boltz CPU inference can be much slower than Vina/rDock and may download a large
model cache under `~/.boltz`. For an interactive smoke, run it separately:

```bash
python -m benchmarks.astex_entropy run --mode native --tools boltz
```

Native, only Vina:

```bash
python -m benchmarks.astex_entropy run --mode native --tools vina
```

Non-native tier-1 smoke:

```bash
python -m benchmarks.astex_entropy run --mode non_native --tools flexaidds
python -m benchmarks.astex_entropy run --mode non_native --tools vina,rdock,boltz
```

Outputs:

```text
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/poses/native_flexaidds_poses.csv
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/poses/native_vina_poses.csv
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/poses/native_rdock_poses.csv
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/poses/native_boltz_poses.csv
```

Each `run` invocation rewrites the current mode/tool pose CSVs. Use a separate
`work_dir` in `config.yaml` if you want to keep smoke and full outputs side by
side.

Use `--dry-run` to verify manifests, generated FlexAIDdS commands, tool inputs,
and configured executables without launching FlexAIDdS, Vina, rDock, or Boltz.
Dry runs do not overwrite existing pose CSVs.

## Rescore

Rescore FlexAIDdS poses in native mode:

```bash
python -m benchmarks.astex_entropy rescore --mode native --poses_from flexaidds
```

Rescore Vina poses in native mode:

```bash
python -m benchmarks.astex_entropy rescore --mode native --poses_from vina
```

Rescore rDock or Boltz poses:

```bash
python -m benchmarks.astex_entropy rescore --mode native --poses_from rdock
python -m benchmarks.astex_entropy rescore --mode native --poses_from boltz
```

Non-native:

```bash
python -m benchmarks.astex_entropy rescore --mode non_native --poses_from vina
python -m benchmarks.astex_entropy rescore --mode non_native --poses_from rdock
python -m benchmarks.astex_entropy rescore --mode non_native --poses_from boltz
python -m benchmarks.astex_entropy rescore --mode non_native --poses_from flexaidds
```

Rescore outputs:

```text
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/rescored/native/vina/rescored_poses.csv
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/rescored/native/vina/report.md
/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/astex_entropy/rescored/native/vina/gbind_vs_rmsd.png
```

## Scoring Columns

`rescored_poses.csv` includes:

- `analysis_target_id`: original target ID used for target-level ranking and
  multi-cavity aggregation
- `rmsd_A`: RDKit best RMSD to the reference ligand
- `posebusters_all_passed`: PoseBusters validation result
- `success_pb`: `rmsd_A <= 2.0` and PoseBusters passed
- `H_vct_proxy`: external tool score converted to lower-is-better energy
- `ensemble_shannon_nats`: Shannon entropy over the tool pose ensemble
- `shannon_energy_collapse`: pose surprisal, `-ln(p_i)`
- `TdS_shannon`: `kB * T * (-ln(p_i))`
- `TdS_vib`: tENCoM vibrational correction parsed from `tencom_entropy_diff`
- `G_bind`: `H_vct_proxy + TdS_shannon - TdS_vib`
- `rank_entropy`: per-target rank after entropy thermodynamic rescoring

## Notes

- PoseBusters and tENCoM are hard requirements for `rescore`; if either command
  is missing, rescoring fails instead of emitting fake successes.
- `tencom_entropy_diff` is treated as the Eigen-backed tENCoM validator. If it
  is missing, rebuild FlexAIDdS with Eigen/tENCoM support before benchmarking.
- Metal-enabled builds are allowed, but `UnifiedHardwareDispatch` must not
  select Metal for `FITNESS_EVAL` until the Metal path reconstructs the same
  internal-coordinate geometry and CF terms as the CPU scorer. Metal can remain
  available for parity-safe non-GA kernels.
- FlexAIDdS cavity files are generated by `build_lto/cavity_detect_cli` by
  default. Set `tools.flexaidds.cavity.detector: getcleft` only for an explicit
  external NRGlab Get_Cleft comparison.
- FlexAIDdS defaults to `main_occupied_cavity`. The manifest writes
  `cleft_sphere_file` for every target, and `benchmark_datasets` passes it to
  FlexAIDdS as `FLEXAIDDS_CLEFT_SPHERE_FILE`.
- FlexAIDdS runs through the configured `benchmark_datasets` binary, writes its
  raw output under the iCloud work dir, then harvests emitted pose PDBs into the
  same pose CSV format as Vina/rDock/Boltz.
- The module defaults to `/Users/lp.more/Projects/FlexAIDdS` for source data
  because this benchmark workspace is not the main Git checkout.
- Edit `config.yaml` if your binaries or source checkout live elsewhere.
