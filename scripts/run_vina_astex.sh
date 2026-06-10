#!/usr/bin/env bash
# =============================================================================
# run_vina_astex.sh — AutoDock Vina self-docking over the Astex Diverse 85 set
# =============================================================================
# Third engine in the FlexAIDdS / rDock / Vina head-to-head (mirrors
# run_rdock_astex.sh). For each complex this script:
#   1. Receptor <CODE>_apo.pdb -> receptor.pdbqt   (rigid, protonated)
#   2. Ligand   <CODE>_ligand.sdf -> ligand.pdbqt  (protonated)
#   3. Docking box: centred on the crystal-ligand centroid, sized to the
#      ligand's max extent + 10 A padding per axis
#   4. `vina --exhaustiveness 8` -> out.pdbqt (ranked poses, MODEL 1 = best)
#
# Success criterion (applied by parse_vina_results.py): sub-2 A Hungarian RMSD
# of the TOP-1 pose (MODEL 1) vs the crystal ligand — identical to rDock/FlexAIDdS.
#
# Requirements: AutoDock Vina CLI (`vina`) + Open Babel. If Vina is missing the
# script prints install instructions and exits 2.
# =============================================================================
set -u

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ASTEX_DIR="${REPO_ROOT}/benchmarks/astex_diverse/astex_diverse"
OUT_ROOT="${VINA_OUT:-${HOME}/flexaidds_benchmark_results/vina_astex}"
EXHAUST="${VINA_EXHAUSTIVENESS:-8}"
N_MODES="${VINA_NUM_MODES:-10}"
PADDING="${VINA_PADDING:-10.0}"   # A added to ligand extent on each axis
PH="${VINA_PH:-7.4}"

# ---- Tool checks -------------------------------------------------------------
if ! command -v vina >/dev/null 2>&1; then
    cat <<'EOF'
============================================================================
 AutoDock Vina is NOT installed (`vina` CLI not found on PATH).
----------------------------------------------------------------------------
 Install (any one):
   conda install -c conda-forge vina      # provides the `vina` CLI binary
   # or, Python bindings only (no CLI):
   pip install vina                       # then drive via the Python API

 conda-forge ships an osx-arm64 `vina` build, so this is the easy path on
 Apple Silicon. After install, re-run this script.
============================================================================
EOF
    exit 2
fi
if ! command -v obabel >/dev/null 2>&1; then
    echo "ERROR: Open Babel (obabel) is required for receptor/ligand prep." >&2
    exit 3
fi

mkdir -p "${OUT_ROOT}"
echo "vina       : $(command -v vina)  ($(vina --version 2>&1 | head -1))"
echo "Astex dir  : ${ASTEX_DIR}"
echo "Output dir : ${OUT_ROOT}"
echo "exhaust=${EXHAUST}  num_modes=${N_MODES}  padding=${PADDING} A  pH=${PH}"
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
        echo "[${code}] SKIP — no ligand SDF (peptide/unsupported ligand)"
        skipped=$((skipped+1)); continue
    fi

    # 1) Receptor -> rigid pdbqt (-xr keeps it rigid; -p adds H at pH)
    rec_pdbqt="${work}/receptor.pdbqt"
    obabel -ipdb "${rec_pdb}" -opdbqt -O "${rec_pdbqt}" -p "${PH}" -xr \
        >/dev/null 2>"${work}/obabel_receptor.log" || {
        echo "[${code}] FAIL — receptor pdbqt conversion"; failed=$((failed+1)); continue; }

    # 2) Ligand -> pdbqt (-h adds H; obabel sets rotatable bonds / ROOT)
    lig_pdbqt="${work}/ligand.pdbqt"
    obabel -isdf "${lig_sdf}" -opdbqt -O "${lig_pdbqt}" -h \
        >/dev/null 2>"${work}/obabel_ligand.log" || {
        echo "[${code}] FAIL — ligand pdbqt conversion"; failed=$((failed+1)); continue; }

    # 3) Box from crystal-ligand geometry (centroid centre, extent+padding size)
    box=$(python3 - "${lig_sdf}" "${PADDING}" <<'PY'
import sys
path, pad = sys.argv[1], float(sys.argv[2])
xs=ys=zs=None
xmn=ymn=zmn= 1e9; xmx=ymx=zmx=-1e9
with open(path, errors="replace") as fh:
    lines = fh.read().splitlines()
try:
    natoms = int(lines[3][0:3])
except Exception:
    natoms = 0
for a in range(natoms):
    ln = lines[4+a]
    try:
        x=float(ln[0:10]); y=float(ln[10:20]); z=float(ln[20:30])
    except Exception:
        p=ln.split(); x,y,z=float(p[0]),float(p[1]),float(p[2])
    xmn=min(xmn,x); xmx=max(xmx,x)
    ymn=min(ymn,y); ymx=max(ymx,y)
    zmn=min(zmn,z); zmx=max(zmx,z)
cx=(xmn+xmx)/2; cy=(ymn+ymx)/2; cz=(zmn+zmx)/2
sx=(xmx-xmn)+pad; sy=(ymx-ymn)+pad; sz=(zmx-zmn)+pad
print(f"{cx:.3f} {cy:.3f} {cz:.3f} {sx:.3f} {sy:.3f} {sz:.3f}")
PY
)
    if [ -z "${box}" ]; then
        echo "[${code}] FAIL — box computation"; failed=$((failed+1)); continue; fi
    read -r CX CY CZ SX SY SZ <<<"${box}"

    # 4) Dock
    out_pdbqt="${work}/out.pdbqt"
    vina --receptor "${rec_pdbqt}" --ligand "${lig_pdbqt}" \
         --center_x "${CX}" --center_y "${CY}" --center_z "${CZ}" \
         --size_x "${SX}" --size_y "${SY}" --size_z "${SZ}" \
         --exhaustiveness "${EXHAUST}" --num_modes "${N_MODES}" \
         --out "${out_pdbqt}" \
         >"${work}/vina.log" 2>&1 || {
        echo "[${code}] FAIL — vina"; failed=$((failed+1)); continue; }

    if [ -f "${out_pdbqt}" ]; then
        echo "[${code}] OK  (box ${SX}x${SY}x${SZ} A @ ${CX},${CY},${CZ})"; ok=$((ok+1))
    else
        echo "[${code}] FAIL — no out.pdbqt produced"; failed=$((failed+1))
    fi
done

echo
echo "=========================================================="
echo "Vina Astex run complete:  OK=${ok}  SKIP=${skipped}  FAIL=${failed}  (of ${total})"
echo "Raw poses in: ${OUT_ROOT}/<CODE>/out.pdbqt"
echo "Next: python scripts/parse_vina_results.py --vina-dir ${OUT_ROOT}"
echo "=========================================================="
