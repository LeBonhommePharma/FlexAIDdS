#!/usr/bin/env bash
# =============================================================================
# narval_cmaes_array.sh — Alliance / Narval Slurm array for GA vs CMA-ES A/B
#
# Copyright 2026 Le Bonhomme Pharma
# SPDX-License-Identifier: Apache-2.0
#
# Submit (on a Compute Canada / Alliance login node ONLY):
#   export CC_ACCOUNT=rrg-yourlab
#   export FLEXAIDDS_SRC=$SCRATCH/FlexAIDdS          # repo root on cluster FS
#   export SIF=$FLEXAIDDS_SRC/.swarm/cmaes/chunk5_harness/artifacts/containers/flexaidds_locked_x86_64.sif
#   export MANIFEST=$FLEXAIDDS_SRC/.swarm/cmaes/chunk5_harness/artifacts/scripts/cmaes_ab_manifest.json
#   module load apptainer
#   sbatch path/to/narval_cmaes_array.sh
#
# Array map (matches cmaes_ab_manifest.json arms[].array_task_id):
#   0 → GA      (FLEXAIDDS_SEARCH=ga)
#   1 → CMA-ES  (FLEXAIDDS_SEARCH=cmaes)
#
# Eval budget: population=1000 × generations=2000 = 2e6 per arm (manifest).
# No machine-specific absolute paths are hard-coded; resolve via env / $SLURM_SUBMIT_DIR.
# =============================================================================

#SBATCH --job-name=flexaidds-cmaes-ab
#SBATCH --account=${CC_ACCOUNT}
#SBATCH --array=0-1
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=out/cmaes_ab/slurm_%A_%a.out
#SBATCH --error=out/cmaes_ab/slurm_%A_%a.err

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve paths (env overrides; never hard-code /Users/... or other host abs)
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HARNESS_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"          # .../artifacts
# Prefer explicit FLEXAIDDS_SRC; else walk up from harness; else submit dir.
if [[ -z "${FLEXAIDDS_SRC:-}" ]]; then
  if [[ -f "${HARNESS_ROOT}/../../../CMakeLists.txt" ]]; then
    FLEXAIDDS_SRC="$(cd "${HARNESS_ROOT}/../../.." && pwd)"
  else
    FLEXAIDDS_SRC="${SLURM_SUBMIT_DIR:-${PWD}}"
  fi
fi
export FLEXAIDDS_SRC

MANIFEST="${MANIFEST:-${SCRIPT_DIR}/cmaes_ab_manifest.json}"
SIF="${SIF:-${HARNESS_ROOT}/containers/flexaidds_locked_x86_64.sif}"
CONFIG_DIR="${CONFIG_DIR:-${FLEXAIDDS_SRC}/out/cmaes_ab/config}"
OUT_ROOT="${OUT_ROOT:-${FLEXAIDDS_SRC}/out/cmaes_ab}"

mkdir -p "${OUT_ROOT}" "${CONFIG_DIR}" \
  "${FLEXAIDDS_SRC}/out/cmaes_ab/1G9V"

# ---------------------------------------------------------------------------
# Account + tooling checks
# ---------------------------------------------------------------------------
if [[ -z "${CC_ACCOUNT:-}" ]]; then
  echo "ERROR: CC_ACCOUNT is unset. Export your Alliance account, e.g.:" >&2
  echo "  export CC_ACCOUNT=rrg-yourlab" >&2
  exit 2
fi

if ! command -v apptainer >/dev/null 2>&1 && ! command -v singularity >/dev/null 2>&1; then
  # Alliance modules
  if command -v module >/dev/null 2>&1; then
    module load apptainer 2>/dev/null || module load singularity 2>/dev/null || true
  fi
fi

if command -v apptainer >/dev/null 2>&1; then
  RT=apptainer
elif command -v singularity >/dev/null 2>&1; then
  RT=singularity
else
  echo "ERROR: neither apptainer nor singularity found on PATH" >&2
  exit 3
fi

if [[ ! -f "${SIF}" ]]; then
  echo "ERROR: SIF not found: ${SIF}" >&2
  echo "Build on x86_64 Linux (see containers/README_HARNESS.md):" >&2
  echo "  apptainer build --bind \"\${FLEXAIDDS_SRC}:/opt/flexaidds/src\" \\" >&2
  echo "      flexaidds_locked_x86_64.sif flexaidds_locked_x86_64.def" >&2
  exit 4
fi

if [[ ! -f "${MANIFEST}" ]]; then
  echo "ERROR: manifest not found: ${MANIFEST}" >&2
  exit 5
fi

# Optional schema gate (pure Python; no third-party deps)
if command -v python3 >/dev/null 2>&1; then
  python3 "${SCRIPT_DIR}/validate_cmaes_manifest.py" "${MANIFEST}"
fi

# ---------------------------------------------------------------------------
# Parse arm from array id (stdlib python — portable on cluster nodes)
# ---------------------------------------------------------------------------
TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
export TASK_ID MANIFEST

eval "$(python3 - <<'PY'
import json, os, sys
mid = int(os.environ["TASK_ID"])
man = json.load(open(os.environ["MANIFEST"], encoding="utf-8"))
arm = None
for a in man["arms"]:
    if int(a.get("array_task_id", -1)) == mid:
        arm = a
        break
if arm is None:
    # fallback: index into arms list
    arms = man["arms"]
    if mid < 0 or mid >= len(arms):
        print(f"echo ERROR: array task {mid} out of range; exit 6", file=sys.stderr)
        sys.exit(6)
    arm = arms[mid]
eb = man["eval_budget"]
cx = man["complex"]
env = man.get("environment", {})
search = arm["env"].get("FLEXAIDDS_SEARCH", arm.get("search", "ga"))
print(f"ARM_ID={arm['id']!r}")
print(f"SEARCH={search!r}")
print(f"OUT_PREFIX={arm['output_prefix']!r}")
print(f"POP={int(eb['population'])}")
print(f"GEN={int(eb['generations'])}")
print(f"TARGET_EVALS={int(eb['target_evals'])}")
print(f"COMPLEX={cx['code']!r}")
print(f"RECEPTOR={cx['receptor']!r}")
print(f"LIGAND={cx['ligand']!r}")
print(f"REFLIG={cx['reference_ligand']!r}")
for k, v in env.items():
    # export safe KEY=value pairs for bash
    print(f"export {k}={str(v)!r}")
print(f"export FLEXAIDDS_SEARCH={search!r}")
PY
)"

echo "=== Narval CMA-ES A/B task ==="
echo "job=${SLURM_JOB_ID:-local} array=${TASK_ID} arm=${ARM_ID} search=${SEARCH}"
echo "src=${FLEXAIDDS_SRC}"
echo "sif=${SIF}"
echo "budget pop=${POP} gen=${GEN} target_evals=${TARGET_EVALS}"
echo "complex=${COMPLEX}"
echo "runtime=${RT}"

# ---------------------------------------------------------------------------
# Materialize a dock config with the 2e6 budget (repo-relative paths → /work)
# ---------------------------------------------------------------------------
CONFIG_JSON="${CONFIG_DIR}/dock_${COMPLEX}_${ARM_ID}_2e6.json"
python3 - <<PY
import json
from pathlib import Path
cfg = {
    "flexibility": {
        "force_rigid": False,
        "intramolecular": True,
        "permeability": 0.9,
        "soft_wall_cutoff": 0.4,
        "receptor_rotamer_prep": True,
    },
    "optimization": {"grid_spacing": 0.375},
    "scoring": {
        "normalize_area": True,
        "vct_dist_weight_r0": 7,
        "vct_normalize_contacts": False,
        "hbond_enabled": True,
        "hbond_search_enabled": True,
        "hbond_rank_enabled": False,
        "metal_coord_enabled": True,
        "sas_weight": 1.0,
        "tencom_weight": 0.0,
        "vct_entropy_weight": 0,
    },
    "seeding": {"mif_enabled": True},
    "reference_ligand": {
        "file": "/work/${REFLIG}",
        "seed_fraction": 0,
        "pose_seed_enabled": False,
        "k_nearest": 10,
    },
    "coarse_init": {
        "enabled": True,
        "grid_step": 3.0,
        "n_seeds": 25,
        "n_orientations": 16,
    },
    "thermodynamics": {
        "temperature": 300,
        "clustering_algorithm": "CF",
        "cluster_rmsd": 2.0,
    },
    "ga": {
        "num_chromosomes": int("${POP}"),
        "num_generations": int("${GEN}"),
        "crossover_rate": 0.8,
        "mutation_rate": 0.03,
        "diversity_monitoring": True,
        "adaptive": True,
        "adaptive_k": [0.95, 0.1, 1.0, 0.05],
        "sharing_alpha": 4,
        "boom_inject_interval": 100,
        "boom_inject_fraction": 0,
        "n_elite": 1,
        "fitness_model": "SMFREE",
    },
}
Path("${CONFIG_JSON}").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
print("wrote", "${CONFIG_JSON}")
PY

# Ensure parent of output prefix exists (host path; mirrored under /work)
mkdir -p "${FLEXAIDDS_SRC}/$(dirname "${OUT_PREFIX}")"

# ---------------------------------------------------------------------------
# apptainer exec — bind repo at /work; data dir from image %environment
# ---------------------------------------------------------------------------
set -x
"${RT}" exec \
  --bind "${FLEXAIDDS_SRC}:/work" \
  "${SIF}" \
  env \
    FLEXAIDDS_DATA_DIR="${FLEXAIDDS_DATA_DIR:-/opt/flexaidds/data}" \
    FLEXAID_SEED="${FLEXAID_SEED:-12345}" \
    FLEXAIDDS_NO_SEC="${FLEXAIDDS_NO_SEC:-1}" \
    FLEXAIDDS_RESTARTS="${FLEXAIDDS_RESTARTS:-1}" \
    FLEXAIDDS_SEARCH="${FLEXAIDDS_SEARCH}" \
    OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}" \
  /opt/flexaidds/bin/FlexAIDdS \
    "/work/${RECEPTOR}" \
    "/work/${LIGAND}" \
    -c "/work/out/cmaes_ab/config/dock_${COMPLEX}_${ARM_ID}_2e6.json" \
    -o "/work/${OUT_PREFIX}"
set +x

echo "=== task ${TASK_ID} arm=${ARM_ID} finished ==="
echo "outputs under: ${FLEXAIDDS_SRC}/${OUT_PREFIX}*"
