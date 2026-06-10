#!/usr/bin/env bash
# =============================================================================
# run_rdock_astex.sh — rDock self-docking over the Astex Diverse 85 set
# =============================================================================
# Head-to-head companion to the FlexAIDdS v27 Astex benchmark.
#
# For each of the 85 complexes this script:
#   1. Prepares the receptor   <CODE>_apo.pdb  -> receptor.mol2  (protonated)
#   2. Prepares the ref ligand <CODE>_ligand.sdf -> ref.sd       (protonated)
#   3. Writes a cavity system .prm using rDock's reference-ligand site mapper
#   4. Runs `rbcavity` to carve the binding-site cavity (.as)
#   5. Runs `rbdock -n 10` to emit 10 docked poses scored by RbtInterIdxSF
#
# Success criterion (applied later by parse_rdock_results.py): sub-2 A
# Hungarian RMSD of the TOP-1 emitted pose vs the crystal ligand.
#
# This script does NOT compute RMSD — it only produces the raw rDock SD output.
# Run parse_rdock_results.py afterwards to score and build the comparison CSV.
#
# Requirements: rDock (rbcavity, rbdock) on PATH with $RBT_ROOT set; Open Babel.
# If rDock is missing the script prints install instructions and exits 2.
# =============================================================================
set -u

# ---- Paths -------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASTEX_DIR="${REPO_ROOT}/benchmarks/astex_diverse/astex_diverse"
OUT_ROOT="${RDOCK_OUT:-${HOME}/flexaidds_benchmark_results/rdock_astex}"
N_POSES="${RDOCK_N_POSES:-10}"
PH="${RDOCK_PH:-7.4}"          # protonation pH for obabel
CAVITY_RADIUS="${RDOCK_RADIUS:-6.0}"   # ref-ligand site mapper radius (A)

# ---- Tool checks -------------------------------------------------------------
need_install=0
if ! command -v rbdock >/dev/null 2>&1 || ! command -v rbcavity >/dev/null 2>&1; then
    need_install=1
fi
if [ "${need_install}" -eq 1 ]; then
    cat <<'EOF'
============================================================================
 rDock is NOT installed (rbdock / rbcavity not found on PATH).
----------------------------------------------------------------------------
 Verified package landscape (checked on this machine, Apple Silicon / arm64):
   * conda-forge : NO rdock package.
   * bioconda    : linux-64 ONLY (rdock=24.04.204_legacy). No osx-64, no
                   osx-arm64 build -> cannot `conda install` it natively here.
   * Docker/podman: not installed on this machine.

 So on this macOS arm64 host there are two real paths:

 (A) Native source build (no Docker needed) — VERIFIED on this machine:
       # Full step-by-step (incl. 3 required libc++ source patches) is in
       #   scripts/RDOCK_SETUP.md
       brew install gsl cppunit popt
       git clone https://github.com/CBDD/rDock.git ~/Software/rDock
       # ...apply the 3 patches from RDOCK_SETUP.md, then:
       cd ~/Software/rDock && make build -j4 CXX=clang++ \
         CXX_EXTRA_FLAGS="-D_LIBCPP_ENABLE_CXX17_REMOVED_FEATURES \
           -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION \
           -D_LIBCPP_ENABLE_CXX17_REMOVED_BINDERS"
       export RBT_ROOT=~/Software/rDock
       export PATH="$RBT_ROOT/bin:$PATH"
       export DYLD_LIBRARY_PATH="$RBT_ROOT/lib:$DYLD_LIBRARY_PATH"

 (B) Run the prebuilt linux-64 bioconda package in a Linux container/VM:
       # install Docker Desktop or colima first, then:
       docker run --platform linux/amd64 -it \
         -v "$PWD":/work continuumio/miniconda3 bash
       #   (inside) conda install -n base -c bioconda -c conda-forge rdock
       # Bind-mount this repo so the Astex inputs are reachable, run this
       # script inside the container, then parse results back on the host.

 After install, set RBT_ROOT and re-run. The script auto-detects rbdock.
============================================================================
EOF
    exit 2
fi

if ! command -v obabel >/dev/null 2>&1; then
    echo "ERROR: Open Babel (obabel) is required for receptor/ligand prep." >&2
    exit 3
fi

# $RBT_ROOT supplies the standard scoring fn + docking protocol files
: "${RBT_ROOT:?RBT_ROOT must be set (points at the rDock install root)}"
DOCK_PRM="${RBT_DOCK_PRM:-${RBT_ROOT}/data/scripts/dock.prm}"
if [ ! -f "${DOCK_PRM}" ]; then
    echo "ERROR: rDock docking protocol not found: ${DOCK_PRM}" >&2
    echo "       Set RBT_DOCK_PRM to the correct dock.prm path." >&2
    exit 3
fi

mkdir -p "${OUT_ROOT}"
echo "rDock      : $(command -v rbdock)"
echo "RBT_ROOT   : ${RBT_ROOT}"
echo "dock.prm   : ${DOCK_PRM}"
echo "Astex dir  : ${ASTEX_DIR}"
echo "Output dir : ${OUT_ROOT}"
echo "Poses/lig  : ${N_POSES}   pH: ${PH}   cavity radius: ${CAVITY_RADIUS} A"
echo

CODES=$(ls -d "${ASTEX_DIR}"/*/ 2>/dev/null | xargs -n1 basename)
total=$(echo "${CODES}" | wc -w | tr -d ' ')
echo "Found ${total} complexes."
echo

ok=0; skipped=0; failed=0
for code in ${CODES}; do
    src="${ASTEX_DIR}/${code}"
    rec_pdb="${src}/${code}_apo.pdb"
    lig_sdf="${src}/${code}_ligand.sdf"
    work="${OUT_ROOT}/${code}"
    mkdir -p "${work}"

    if [ ! -f "${rec_pdb}" ]; then
        echo "[${code}] SKIP — no apo receptor PDB"; skipped=$((skipped+1)); continue
    fi
    if [ ! -f "${lig_sdf}" ]; then
        # 1TW6 = Smac AVPI peptide ligand, no extracted SDF (see project memory)
        echo "[${code}] SKIP — no ligand SDF (peptide/unsupported ligand)"
        skipped=$((skipped+1)); continue
    fi

    # 1) Receptor PDB -> protonated mol2.
    #    No Gasteiger charges: rDock assigns its own atom types and computes
    #    the InterIdxSF terms internally, and Gasteiger fails to converge on
    #    metalloproteins (e.g. HEM in 1G9V) -> 0 molecules. -p adds H at pH.
    rec_mol2="${work}/receptor.mol2"
    obabel "${rec_pdb}" -O "${rec_mol2}" -p "${PH}" \
        >/dev/null 2>"${work}/obabel_receptor.log" || {
        echo "[${code}] FAIL — receptor mol2 conversion"; failed=$((failed+1)); continue; }

    # 2) Reference / docking ligand -> protonated SD (rDock needs explicit H).
    ref_sd="${work}/ref.sd"
    obabel "${lig_sdf}" -O "${ref_sd}" -p "${PH}" \
        >/dev/null 2>"${work}/obabel_ligand.log" || {
        echo "[${code}] FAIL — ligand SD conversion"; failed=$((failed+1)); continue; }

    # 3) Cavity system .prm (reference-ligand site mapper around the crystal pose)
    prm="${work}/system.prm"
    cat > "${prm}" <<EOF
RBT_PARAMETER_FILE_V1.00
TITLE astex_${code}

RECEPTOR_FILE receptor.mol2
RECEPTOR_FLEX 3.0

SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL ref.sd
    RADIUS ${CAVITY_RADIUS}
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
    GRIDSTEP 0.5
END_SECTION

SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION
EOF

    # 4) Carve cavity (.as). Run from work dir so relative paths in .prm resolve.
    #    NB: this rDock revision uses the cxxopts CLI — write=-W, dump=-d
    #    (the legacy single-dash `-was` is rejected by the new parser).
    ( cd "${work}" && rbcavity -r "system.prm" -W -d \
        >"${work}/rbcavity.log" 2>&1 ) || {
        echo "[${code}] FAIL — rbcavity"; failed=$((failed+1)); continue; }

    # 5) Dock N poses. Output -> docked.sd (rbdock appends .sd).
    ( cd "${work}" && rbdock -i "ref.sd" -o "docked" -r "system.prm" \
        -p "${DOCK_PRM}" -n "${N_POSES}" \
        >"${work}/rbdock.log" 2>&1 ) || {
        echo "[${code}] FAIL — rbdock"; failed=$((failed+1)); continue; }

    if [ -f "${work}/docked.sd" ]; then
        echo "[${code}] OK"; ok=$((ok+1))
    else
        echo "[${code}] FAIL — no docked.sd produced"; failed=$((failed+1))
    fi
done

echo
echo "=========================================================="
echo "rDock Astex run complete:  OK=${ok}  SKIP=${skipped}  FAIL=${failed}  (of ${total})"
echo "Raw poses in: ${OUT_ROOT}/<CODE>/docked.sd"
echo "Next: python scripts/parse_rdock_results.py --rdock-dir ${OUT_ROOT}"
echo "=========================================================="
