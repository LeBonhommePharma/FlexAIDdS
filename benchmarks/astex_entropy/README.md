# Astex Entropy Benchmark

Clean Astex Diverse Set benchmark for a MacBook Pro:

- Modes: `native` and `non_native`
- Pose generators: AutoDock Vina, rDock, Boltz-2
- Rescoring: Shannon energy collapse, tENCoM, and thermodynamic `G_bind`
- Success definition: `RMSD <= 2.0 A` and PoseBusters passes

This module is intentionally small. It writes manifests, tool inputs, poses,
rescored CSVs, and plots under `results/astex_entropy/`.

## Install

From the real FlexAIDdS checkout:

```bash
cd /Users/lp.more/Projects/FlexAIDdS

/opt/homebrew/bin/python3.12 -m venv .venv-astex-entropy
source .venv-astex-entropy/bin/activate
python -m pip install --upgrade pip

# Python dependencies only.
python -m pip install -r benchmarks/astex_entropy/requirements.txt
```

The default `config.yaml` is pinned to local tool paths under the checkout:

```text
.tools/bin/vina
.tools/rdock/bin/rbcavity
.tools/rdock/bin/rbdock
.venv-posebusters/bin/bust
.venv-boltz/bin/boltz
/opt/homebrew/bin/obabel
```

Install the external command-line tools:

```bash
brew install open-babel

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
.tools/bin/vina --help
.venv-posebusters/bin/bust --help
NUMBA_CACHE_DIR=results/astex_entropy/cache/numba .venv-boltz/bin/boltz predict --help
```

Boltz-2 is configured for the 18 GB MacBook Pro constraint: CPU accelerator,
one device, one dataloader worker, one parallel sample, and one diffusion
sample by default. It also uses `NUMBA_CACHE_DIR` under `results/astex_entropy`
to avoid Numba cache-locator crashes in the venv package path.

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
from RCSB into `results/astex_entropy/data/non_native/`.

```bash
python -m benchmarks.astex_entropy data_prep --mode non_native --max-targets 3 --download-missing
```

For a one-target wiring smoke:

```bash
python -m benchmarks.astex_entropy data_prep --mode native --max-targets 1 --force
python -m benchmarks.astex_entropy run --mode native --tools vina
python -m benchmarks.astex_entropy rescore --mode native --poses_from vina
python -m benchmarks.astex_entropy run --mode native --tools rdock
python -m benchmarks.astex_entropy rescore --mode native --poses_from rdock
```

## Run Pose Generators

Native, all three tools:

```bash
python -m benchmarks.astex_entropy run --mode native
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
python -m benchmarks.astex_entropy run --mode non_native --tools vina,rdock,boltz
```

Outputs:

```text
results/astex_entropy/poses/native_vina_poses.csv
results/astex_entropy/poses/native_rdock_poses.csv
results/astex_entropy/poses/native_boltz_poses.csv
```

Each `run` invocation rewrites the current mode/tool pose CSVs. Use a separate
`work_dir` in `config.yaml` if you want to keep smoke and full outputs side by
side.

Use `--dry-run` to verify manifests, tool inputs, and configured executables
without launching Vina, rDock, or Boltz. Dry runs do not overwrite existing pose
CSVs.

## Rescore

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
```

Rescore outputs:

```text
results/astex_entropy/rescored/native/vina/rescored_poses.csv
results/astex_entropy/rescored/native/vina/report.md
results/astex_entropy/rescored/native/vina/gbind_vs_rmsd.png
```

## Scoring Columns

`rescored_poses.csv` includes:

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
- The module defaults to `/Users/lp.more/Projects/FlexAIDdS` for source data
  because some agent sessions start in the separate benchmark workspace under
  `Documents/PhD/Programs/FlexAIDdS`.
- Edit `config.yaml` if your binaries or source checkout live elsewhere.
