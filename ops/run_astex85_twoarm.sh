#!/usr/bin/env bash
# =============================================================================
# run_astex85_twoarm.sh — Astex-85 three-arm campaign
#
# Two independent CF.com blow-up mechanisms are under test here (elected pose
# CF.com reaches ~-3200 against a healthy ~-130):
#
#   Type I  (water-driven)   : a water-rich pocket lets the GA bury the ligand
#                              in the crystallographic solvent shell, harvesting
#                              unbounded com from sub-Ångström ligand-O ⋯ HOH-O
#                              contacts (e.g. 1JD0, 156 retained waters).
#                              → smart water retention removes the attack surface.
#   Type II (protein-O crowd): an exposed sp3-O pocket (Ser/Thr/Tyr OH, backbone
#                              carbonyls) does the same with protein oxygens, so
#                              stripping every water does NOT stop it (1SG0: com
#                              -3254, 95 % from a single ligand-C × O.3 pair).
#                              → needs intensive (per-contact) com + a lower clamp.
#
#   Arm A (water-only)  : smart water retention ON, RAW extensive com
#                         (no VCT_NORM, no COM_FLOOR, no thermo).
#                         Tests the water-prep-only hypothesis in isolation.
#   Arm B (PRODUCTION)  : smart water retention ON
#                         + FLEXAIDDS_VCT_NORM=1 + FLEXAIDDS_COM_FLOOR=500
#                         + thermodynamics (THERMO_SCORE=1, T_EFF=0.596,
#                           SOFTBETA_ELECTION=1). The combined fix.
#   Arm C (com-only)    : ALL waters stripped (keep_structural_waters:false)
#                         + FLEXAIDDS_VCT_NORM=1 + FLEXAIDDS_COM_FLOOR=500,
#                         no thermo. Tests com-taming without any water fix.
#
# NOTE ON ARM B ≈ ARM C + SMART-WATER — the thermo selection gate is INERT:
# gaboom.cpp:1341 is printf-only, and the ranking that picks the elected pose is
# the QuickSort at gaboom.cpp:705, finalized BEFORE the thermodynamic compute
# runs. So THERMO_SCORE/T_EFF/SOFTBETA_ELECTION do not move which pose is
# elected. Any Arm B − Arm C difference is therefore attributable to WATER
# RETENTION (B keeps bridging waters, C strips all), NOT to thermodynamics.
# Read the B-vs-C delta as the water-prep effect on top of the com fix.
#
# VCT_NORM + COM_FLOOR are deliberately NOT applied to every arm: Arm A must run
# with raw com to test the Codex water-only hypothesis, and Arm C must run with
# no water fix to test the com-taming hypothesis. Applying the fixes everywhere
# would collapse the ablation.
#
# Arms run SEQUENTIALLY on purpose: benchmark_datasets instances share one cache
# directory and concurrent runs corrupt it.
#
# Reproducibility: OMP_NUM_THREADS=1 and serial restarts
# (FLEXAIDDS_PARALLEL_RESTARTS=0) make each worker bit-deterministic; the 6
# workers are independent processes on different targets, so cross-target
# parallelism does not affect any single target's result.
# =============================================================================
set -uo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
ENGINE="${REPO}/build/FlexAIDdS"
RUNNER="${REPO}/build/benchmark_datasets"
STAMP="$(date +%Y%m%d_%H%M%S)"
ROOT="${HOME}/flexaidds_benchmark_results/astex85_threearm_${STAMP}"
CACHE="${HOME}/.flexaidds/benchmarks"
LOG="${ROOT}/campaign.log"

# Survive session teardown (v7 multi-worker SIGTERM kill). setsid puts us in our
# own process group; only HUP is ignored so the campaign stays deliberately killable.
trap '' HUP

mkdir -p "${ROOT}"
xattr -w com.apple.fileprovider.ignore#P 1 "${ROOT}" 2>/dev/null || true
touch "${ROOT}/.metadata_never_index" 2>/dev/null || true

exec >>"${LOG}" 2>&1

echo "=== Astex-85 three-arm campaign ${STAMP} ==="
echo "engine : ${ENGINE}"
echo "sha256 : $(shasum -a 256 "${ENGINE}" | cut -d' ' -f1)"
echo "commit : $(cd "${REPO}" && git rev-parse HEAD)"
echo "root   : ${ROOT}"

# Single-instance guard.
if pgrep -f "benchmark_datasets --benchmark" >/dev/null 2>&1; then
    echo "[ABORT] a benchmark_datasets run is already active — refusing to share the cache"
    exit 1
fi

# ── Reproducibility / protocol env common to every arm ────────────────────────
export OMP_NUM_THREADS=1
export FLEXAIDDS_PARALLEL_RESTARTS=0     # serial restarts → deterministic
export FLEXAID_SEED=42
export FLEXAIDDS_SEED_BASE=42
export FLEXAIDDS_BINARY="${ENGINE}"
export FLEXAIDDS_ORACLE_SITE_DIR="${REPO}/benchmarks/astex_diverse/astex_diverse"

# ── Mechanism logger ──────────────────────────────────────────────────────────
# After an arm finishes, walk its per-target elected poses and pull the CF.com
# channel and the single largest per-pair contact contribution out of the pose
# REMARKs. Emitting this alongside RMSD turns the campaign into a direct
# mechanism test: we can see "did com blow up, and on which contact type?"
# rather than inferring it from RMSD. Healthy com ≈ -130; blow-up ≳ -1000.
com_summary() {
    local out="$1"
    local csv="${out}/com_mechanism.csv"
    echo "pdb_id,cf_com,cf_sas,cf_wal,cf_total,dominant_contact_type,dominant_contact_value" > "${csv}"
    local pose
    for pose in "${out}"/*/elected_pose.pdb; do
        [ -f "${pose}" ] || continue
        local pdb; pdb="$(basename "$(dirname "${pose}")")"
        local com sas wal tot
        com="$(awk '/^REMARK CF.com=/{print $2}' "${pose}" | sed 's/.*=//' | head -1)"
        sas="$(awk '/^REMARK CF.sas=/{print $2}' "${pose}" | sed 's/.*=//' | head -1)"
        wal="$(awk '/^REMARK CF.wal=/{print $2}' "${pose}" | sed 's/.*=//' | head -1)"
        tot="$(awk -F= '/^REMARK CF=/{print $2}' "${pose}" | head -1)"
        # Dominant contact-type pair, if the engine emitted per-pair REMARKs.
        # Format: "REMARK contact <t1>-<t2> ... <value>". Pick the row whose
        # last field has the largest magnitude.
        local dom
        dom="$(awk '/^REMARK contact /{v=$NF; a=v<0?-v:v; if(a>best){best=a; line=$3" "$NF}} END{if(line!="")print line}' "${pose}")"
        local dtype dval
        dtype="$(echo "${dom}" | awk '{print $1}')"
        dval="$(echo "${dom}" | awk '{print $2}')"
        echo "${pdb},${com:-NA},${sas:-NA},${wal:-NA},${tot:-NA},${dtype:-NA},${dval:-NA}" >> "${csv}"
    done
    echo "  [MECHANISM] wrote $(wc -l < "${csv}") rows → ${csv}"
}

run_arm() {
    local name="$1"; shift
    local out="${ROOT}/${name}"
    mkdir -p "${out}"
    echo ""
    echo "=== ARM ${name} starting $(date -u +%FT%TZ) ==="
    env | grep -E '^FLEXAIDDS_(THERMO_SCORE|T_EFF|SOFTBETA_ELECTION|SMART_WATER|STRIP_ALL_WATERS|VCT_NORM|COM_FLOOR)=' | sort || true
    caffeinate -i "${RUNNER}" \
        --benchmark astex \
        --output "${out}" \
        --cache  "${CACHE}" \
        --threads 6 \
        --omp-threads 1 \
        --job-timeout-seconds 3600 \
        --force
    echo "=== ARM ${name} finished rc=$? $(date -u +%FT%TZ) ==="
    com_summary "${out}"
}

# ── Arm A: smart water ON, raw extensive com (Codex water-only hypothesis) ─────
(
  unset FLEXAIDDS_THERMO_SCORE FLEXAIDDS_T_EFF FLEXAIDDS_SOFTBETA_ELECTION
  unset FLEXAIDDS_VCT_NORM FLEXAIDDS_COM_FLOOR FLEXAIDDS_STRIP_ALL_WATERS
  export FLEXAIDDS_SMART_WATER=1
  run_arm "armA_smartwater_rawcom"
)

# ── Arm B: PRODUCTION — smart water + com fix + thermodynamics ─────────────────
(
  unset FLEXAIDDS_STRIP_ALL_WATERS
  export FLEXAIDDS_SMART_WATER=1
  export FLEXAIDDS_VCT_NORM=1
  export FLEXAIDDS_COM_FLOOR=500
  export FLEXAIDDS_THERMO_SCORE=1
  export FLEXAIDDS_T_EFF=0.596
  export FLEXAIDDS_SOFTBETA_ELECTION=1
  run_arm "armB_production"
)

# ── Arm C: strip ALL waters + com fix, no thermo (com-taming hypothesis) ───────
(
  unset FLEXAIDDS_THERMO_SCORE FLEXAIDDS_T_EFF FLEXAIDDS_SOFTBETA_ELECTION FLEXAIDDS_SMART_WATER
  export FLEXAIDDS_STRIP_ALL_WATERS=1
  export FLEXAIDDS_VCT_NORM=1
  export FLEXAIDDS_COM_FLOOR=500
  run_arm "armC_stripwater_comfix"
)

echo ""
echo "=== CAMPAIGN COMPLETE $(date -u +%FT%TZ) ==="
touch "${ROOT}/.CAMPAIGN_DONE"
