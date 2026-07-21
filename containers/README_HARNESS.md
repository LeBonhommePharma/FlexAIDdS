# Locked-arch Apptainer harness (VALIDATION G — local half)

Copyright 2026 Le Bonhomme Pharma · SPDX-License-Identifier: Apache-2.0

This directory ships the **recipe only**. Building the `.sif` requires
Apptainer/Singularity on an **x86_64 Linux** host (Alliance/Narval login node,
CI Linux runner, or a local Linux box). macOS hosts land the recipes and print
the intended `sbatch` / smoke commands; they do **not** build the image.

| Artifact | Role |
|----------|------|
| `flexaidds_locked_x86_64.def` | Apptainer definition (Ubuntu 22.04, GCC≥14, CMake≥3.28, FAST binary) |
| `flexaidds_locked_x86_64.sif` | Built image (not committed; produce locally / on cluster) |
| `../scripts/cmaes_ab_manifest.json` | Smoke complex **1G9V** + GA vs CMA-ES **2e6** eval budget |
| `../scripts/validate_cmaes_manifest.py` | Manifest schema gate (exit 0/1) |
| `../scripts/narval_cmaes_array.sh` | Slurm array (`#SBATCH`, `${CC_ACCOUNT}`, `apptainer exec`) |

Labels baked into the image: `arch=x86_64`, `toolchain=locked`.

---

## 0. Preconditions

```bash
# On the Linux build host
command -v apptainer || command -v singularity
uname -m   # must report x86_64 (or amd64)

# Point at a FlexAIDdS checkout (NO absolute path baked into the recipe)
export FLEXAIDDS_SRC="${FLEXAIDDS_SRC:-$PWD}"   # repo root containing CMakeLists.txt
test -f "${FLEXAIDDS_SRC}/CMakeLists.txt"

# Optional: Alliance account for Narval submit (remote half)
# export CC_ACCOUNT=rrg-yourlab   # required by narval_cmaes_array.sh
```

**Source policy:** the definition does **not** `%files` copy or `%post` clone.
Sources must be **bind-mounted** at `/opt/flexaidds/src` during `apptainer build`.

---

## 1. Build the SIF

From this `containers/` directory (or pass absolute/relative paths to the `.def`):

```bash
# Recommended: bind the live checkout so %post can compile BUILD_FLEXAIDDS_FAST=ON
apptainer build \
  --bind "${FLEXAIDDS_SRC}:/opt/flexaidds/src" \
  flexaidds_locked_x86_64.sif \
  flexaidds_locked_x86_64.def
```

Singularity-compatible spelling:

```bash
singularity build \
  --bind "${FLEXAIDDS_SRC}:/opt/flexaidds/src" \
  flexaidds_locked_x86_64.sif \
  flexaidds_locked_x86_64.def
```

Notes:

- Build on an x86_64 host in the **same CPU family** as the run nodes when
  possible (CMake FAST enables `-march=native` on Linux).
- Rootless builds may need `apptainer build --fakeroot …` depending on site
  policy; on Alliance clusters prefer the site module (`module load apptainer`).
- After a successful bake, record provenance:
  ```bash
  sha256sum flexaidds_locked_x86_64.sif | tee flexaidds_locked_x86_64.sif.sha256
  apptainer exec flexaidds_locked_x86_64.sif cat /opt/flexaidds/share/TOOLCHAIN.txt
  apptainer exec flexaidds_locked_x86_64.sif cat /opt/flexaidds/share/BINARY.txt
  apptainer test flexaidds_locked_x86_64.sif
  ```

Without the bind mount, the image still installs the **locked toolchain** but
ships only a stub binary (exit 127). Always bind sources for a usable SIF.

---

## 2. Smoke dock (one complex, bind mounts)

Eval budget for a **container smoke** is intentionally tiny (prove the binary
and mounts). Full eval-matched **2e6** A/B is driven by the manifest / array
script (items E / G full run), not this smoke.

Prepare a short JSON config (example fields only; paths are **inside** the
container after bind):

```json
{
  "ga": { "num_chromosomes": 1000, "num_generations": 10 },
  "thermodynamics": { "temperature": 300 },
  "reference_ligand": {
    "file": "/work/benchmarks/astex_diverse/astex_diverse/1G9V/1G9V_ligand.sdf",
    "seed_fraction": 0,
    "pose_seed_enabled": false
  }
}
```

Smoke command (repo bound at `/work`; SIF path is relative to your working dir):

```bash
export FLEXAIDDS_SRC="${FLEXAIDDS_SRC:-$PWD}"
export SIF="${SIF:-flexaidds_locked_x86_64.sif}"
mkdir -p "${FLEXAIDDS_SRC}/out/harness_smoke_1G9V"

apptainer exec \
  --bind "${FLEXAIDDS_SRC}:/work" \
  "${SIF}" \
  env FLEXAIDDS_DATA_DIR=/opt/flexaidds/data \
      FLEXAID_SEED=12345 \
      FLEXAIDDS_NO_SEC=1 \
      FLEXAIDDS_RESTARTS=1 \
      FLEXAIDDS_SEARCH="${FLEXAIDDS_SEARCH:-ga}" \
  FlexAIDdS \
    /work/benchmarks/astex_diverse/astex_diverse/1G9V/1G9V_apo.pdb \
    /work/benchmarks/astex_repro/smoke/1G9V/1G9V_dockin.sdf \
    -c /work/out/harness_smoke_1G9V/dock_config_smoke.json \
    -o /work/out/harness_smoke_1G9V/smoke
```

CMA-ES arm (once the engine is wired and `FLEXAIDDS_SEARCH=cmaes` is honoured):

```bash
FLEXAIDDS_SEARCH=cmaes apptainer exec \
  --bind "${FLEXAIDDS_SRC}:/work" \
  "${SIF}" \
  env FLEXAIDDS_DATA_DIR=/opt/flexaidds/data \
      FLEXAID_SEED=12345 \
      FLEXAIDDS_NO_SEC=1 \
      FLEXAIDDS_RESTARTS=1 \
      FLEXAIDDS_SEARCH=cmaes \
  FlexAIDdS \
    /work/benchmarks/astex_diverse/astex_diverse/1G9V/1G9V_apo.pdb \
    /work/benchmarks/astex_repro/smoke/1G9V/1G9V_dockin.sdf \
    -c /work/out/harness_smoke_1G9V/dock_config_smoke.json \
    -o /work/out/harness_smoke_1G9V/smoke_cmaes
```

`FLEXAIDDS_DATA_DIR` defaults to `/opt/flexaidds/data` inside the image
(`%environment`); the energy matrix `MC_st0r5.2_6.dat` and `AMINO.def` /
`NUCLEOTIDES.def` are staged there at build time.

---

## 3. Manifest + validator (local, no Apptainer required)

```bash
# From artifacts/ (or any CWD — pass the manifest path)
python3 scripts/validate_cmaes_manifest.py scripts/cmaes_ab_manifest.json
# exit 0 = schema OK; exit 1 = schema / budget / path errors
```

The manifest encodes **1G9V** and the eval-matched budget
`population=1000 × generations=2000 = 2_000_000` for both **ga** and **cmaes**
arms. Repo-relative input paths only (no machine-specific absolutes).

---

## 4. Narval / Alliance array (remote half — print or submit)

On a **Compute Canada / Alliance login node** only:

```bash
export CC_ACCOUNT=rrg-yourlab          # required
export FLEXAIDDS_SRC=$SCRATCH/FlexAIDdS
export SIF=$FLEXAIDDS_SRC/.swarm/cmaes/chunk5_harness/artifacts/containers/flexaidds_locked_x86_64.sif
export MANIFEST=$FLEXAIDDS_SRC/.swarm/cmaes/chunk5_harness/artifacts/scripts/cmaes_ab_manifest.json

module load apptainer
cd "$FLEXAIDDS_SRC"
sbatch .swarm/cmaes/chunk5_harness/artifacts/scripts/narval_cmaes_array.sh
```

If this host is **not** a login node: **do not** `sbatch`. Print the command
and leave G remote half OPEN (see `validation_evidence/harness/`).

Array tasks:

| `SLURM_ARRAY_TASK_ID` | Arm | Env |
|----------------------|-----|-----|
| 0 | GA | `FLEXAIDDS_SEARCH=ga` |
| 1 | CMA-ES | `FLEXAIDDS_SEARCH=cmaes` |

Both arms use the same 2e6 eval budget from the manifest.

---

## 5. What closes VALIDATION G (local half)

1. `.sif` built on Linux → record `sha256` (e.g. `G_sif.sha256`).
2. In-container smoke dock → real RMSD/CF log on disk.
3. Collapse fingerprint on the smoke entropy/trace (chunk4
   `analysis/collapse_fingerprint.py`) → invariant on disk.
4. Manifest validator exit 0; array script present with `#SBATCH` +
   `${CC_ACCOUNT}`.

Remote half (Narval `sbatch` + cross-arch fingerprint) is separate and must not
be faked from a macOS laptop.
