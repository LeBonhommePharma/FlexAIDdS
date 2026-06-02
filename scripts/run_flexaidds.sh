#!/usr/bin/env bash
# =============================================================================
# run_flexaidds.sh — Reproducibility wrapper for FlexAIDdS single docking runs
#
# Part of the flexaidds-docking skill (CREATE_BUNDLE support added in the
# 1-reproducibility-bundle-implementation workstream).
#
# Purpose:
#   - Provide a stable, auditable entry point for docking invocations from
#     the skill, agents, or manual use.
#   - Always emit a machine-readable reproducibility JSON sidecar
#     (git commit, binary SHA256, input hashes, provenance, thermo + confidence).
#   - Capture the complete raw stdout/stderr of the FlexAIDdS engine (.raw.txt).
#   - (Optional, opt-in) When CREATE_BUNDLE=1 and the run succeeds, emit a
#     single .tar.gz archive containing exactly the artifacts needed for
#     third-party replay and audit.
#
# Usage:
#   CREATE_BUNDLE=1 SKIP_REBUILD=1 \
#   FLEXAID_BINARY=/path/to/FlexAID \
#   bash scripts/run_flexaidds.sh receptor.pdb ligand.mol2 \
#        [--outdir DIR] [--temperature 298.15] [--seed 42] [--visualize]
#
#   # Or let it auto-create a timestamped results dir under $HOME/flexaidds_results
#   CREATE_BUNDLE=1 bash scripts/run_flexaidds.sh 1stp.pdb biotin.mol2
#
#   # With Imagine figure + animation (best mode, Gate 6, NRDD aesthetic):
#   FLEXAIDDS_SOURCE=/path/to/FlexAIDdS SKIP_REBUILD=1 \
#   bash scripts/run_flexaidds.sh 1stp biotin.mol2 --temperature 298.15 -o results/test_run --visualize
#
# Environment variables (all optional, non-breaking defaults):
#   CREATE_BUNDLE=1          → after success, also write run_YYYYMMDD_HHMMSS_XXXX.tar.gz
#   SKIP_REBUILD=1           → passed through; recorded in provenance
#   FLEXAID_BINARY=...       → explicit path to the C++ binary (preferred)
#   FLEXAIDDS_SOURCE=...     → repo root (used by run_metadata.py for git info)
#   RESULTS_DIR=...          → force a specific output directory (advanced)
#   RUN_ID=...               → optional explicit short identifier for bundle name
#   VISUALIZE=1 or --visualize → after success+Gate 6, auto-prep results/figures/ prompts for Grok Imagine (NRDD cover + 6s anim of best mode)
#
# Output layout (inside RESULTS_DIR):
#   reproducibility.json     ← full record from run_metadata.create_run_record
#   dock.raw.txt             ← complete unfiltered engine output
#   (optional) inputs/       ← copies of local receptor + ligand when used
#   run_*.tar.gz             ← only when CREATE_BUNDLE=1 and success
#
# Bundle contract (when CREATE_BUNDLE=1):
#   - Created ONLY on success (non-zero exit or missing artifacts → no bundle)
#   - Filename: run_YYYYMMDD_HHMMSS_<short-run-id>.tar.gz
#   - Contains: reproducibility.json, dock.raw.txt, inputs/* (local files only)
#   - Logged with the exact phrases required by the spec:
#       “Creating reproducibility bundle…”
#       “Bundle created at …”
#
# Design principles:
#   - Minimal: the bundle logic is <60 LOC of pure bash + one python helper call.
#   - Robust: never aborts the main run; best-effort on hashes, git, etc.
#   - No behavior change whatsoever when CREATE_BUNDLE is unset or 0.
#   - Works on macOS (shasum) and Linux (sha256sum); falls back to python.
#   - Zero new Python dependencies beyond what the skill already ships.
#
# Apache-2.0 · FlexAIDdS reproducibility tooling
# =============================================================================

set -euo pipefail

# ─── Colours & logging (consistent with other repo scripts) ───────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'

info()  { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { printf "${GREEN}[ OK ]${NC}  %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*" >&2; }
fail()  { printf "${RED}[FAIL]${NC}  %s\n" "$*" >&2; }

# ─── Resolve script & repo roots ──────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# ─── Defaults & argument parsing ──────────────────────────────────────────────
RECEPTOR=""
LIGAND=""
OUTDIR=""
TEMPERATURE="298.15"
SEED="42"
VISUALIZE="${VISUALIZE:-0}"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir|-o)    OUTDIR="$2"; shift ;;
        --temperature)  TEMPERATURE="$2"; shift ;;
        --seed)         SEED="$2"; shift ;;
        --visualize|-v) VISUALIZE=1 ;;
        -h|--help)
            sed -n '2,80p' "$0" | sed 's/^# \?//'
            exit 0
            ;;
        -*)
            EXTRA_ARGS+=("$1")
            ;;
        *)
            if [[ -z "${RECEPTOR}" ]]; then
                RECEPTOR="$1"
            elif [[ -z "${LIGAND}" ]]; then
                LIGAND="$1"
            else
                EXTRA_ARGS+=("$1")
            fi
            ;;
    esac
    shift
done

if [[ -z "${RECEPTOR}" || -z "${LIGAND}" ]]; then
    fail "Usage: $0 <receptor.pdb> <ligand.mol2> [options]"
    exit 2
fi

# ─── Determine RESULTS_DIR (timestamped by default) ───────────────────────────
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [[ -n "${RESULTS_DIR:-}" ]]; then
    # Caller forced a location
    :
elif [[ -n "${OUTDIR}" ]]; then
    RESULTS_DIR="${OUTDIR}"
else
    RESULTS_DIR="${HOME}/flexaidds_results/run_${TIMESTAMP}"
fi
mkdir -p "${RESULTS_DIR}"

RAW_PATH="${RESULTS_DIR}/dock.raw.txt"
JSON_PATH="${RESULTS_DIR}/reproducibility.json"

# Short identifier used both for provenance and for bundle filename
SHORT_RUN_ID="${RUN_ID:-$(echo "${TIMESTAMP}" | tr -d '_')}"   # fallback; will be refined below

# ─── Locate Python and the skill helpers ──────────────────────────────────────
if command -v python3 >/dev/null 2>&1; then
    PYTHON=python3
elif command -v python >/dev/null 2>&1; then
    PYTHON=python
else
    fail "python3/python not found in PATH"
    exit 1
fi

# Make the local python package importable if we are inside the repo
export PYTHONPATH="${REPO_ROOT}/python:${PYTHONPATH:-}"

# ─── Helper: robust SHA256 (macOS + Linux + python fallback) ──────────────────
compute_sha256() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "file-not-found"
        return
    fi
    if command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$file" | awk '{print $1}'
    elif command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$file" | awk '{print $1}'
    else
        "$PYTHON" -c '
import hashlib, sys
print(hashlib.sha256(open(sys.argv[1], "rb").read()).hexdigest())
' "$file"
    fi
}

# ─── Helper: detect whether a path is a real local file we should archive ─────
is_local_regular_file() {
    local p="$1"
    [[ -f "$p" && ! "$p" =~ ^- && ! "$p" =~ ^(https?|ftp|s3|gs):// && "$p" != "/dev/"* && "$p" != "<"* ]]
}

# ─── Helper: generate an 8-char short run id (pure, no network) ───────────────
generate_short_run_id() {
    if command -v openssl >/dev/null 2>&1; then
        openssl rand -hex 4
    else
        # Portable fallback — still 8 hex chars
        "$PYTHON" -c '
import secrets
print(secrets.token_hex(4))
' 2>/dev/null || echo "${TIMESTAMP//_/}" | tail -c 9
    fi
}

# ─── NEW: Bundle creation (only called on success path when CREATE_BUNDLE=1) ──
# Fully self-contained, documented, and minimal.
create_reproducibility_bundle() {
    local results_dir="$1"
    local json_file="$2"
    local raw_file="$3"
    local receptor_path="$4"
    local ligand_path="$5"

    local ts
    ts="$(date +%Y%m%d_%H%M%S)"
    local short_id
    short_id="$(generate_short_run_id)"
    # Update global for consistency if caller wants it
    SHORT_RUN_ID="${short_id}"

    local bundle_name="run_${ts}_${short_id}.tar.gz"
    local bundle_path="${results_dir}/${bundle_name}"

    printf "Creating reproducibility bundle...\n"

    # Use a private staging directory so the tar contains a clean, minimal tree
    local stage
    stage="$(mktemp -d "${results_dir}/.bundle.XXXXXX")"
    # We intentionally do NOT trap-clean here; we want the bundle even on weird exits

    mkdir -p "${stage}/inputs"

    # 1. The canonical reproducibility JSON (already contains everything required)
    if [[ -f "${json_file}" ]]; then
        cp -p "${json_file}" "${stage}/reproducibility.json"
    else
        warn "reproducibility.json missing — bundle will be incomplete"
    fi

    # 2. The complete raw engine transcript
    if [[ -f "${raw_file}" ]]; then
        cp -p "${raw_file}" "${stage}/dock.raw.txt"
    else
        warn "dock.raw.txt missing — bundle will be incomplete"
    fi

    # 3. Local input files only (never URLs, never pipes, never missing files)
    if is_local_regular_file "${receptor_path}"; then
        cp -p "${receptor_path}" "${stage}/inputs/$(basename "${receptor_path}")"
    fi
    if is_local_regular_file "${ligand_path}"; then
        cp -p "${ligand_path}" "${stage}/inputs/$(basename "${ligand_path}")"
    fi

    # Create a deterministic, portable tarball (no owner/mode variance that matters)
    tar -czf "${bundle_path}" -C "${stage}" .

    # Best-effort cleanup of the transient stage dir
    rm -rf "${stage}" 2>/dev/null || true

    printf "Bundle created at %s\n" "${bundle_path}"
    printf "  (sha256: %s  size: %s)\n" \
        "$(compute_sha256 "${bundle_path}")" \
        "$(du -h "${bundle_path}" 2>/dev/null | awk '{print $1}')"
}

# ─── Locate the FlexAIDdS binary (or fall back to python entry point) ─────────
FLEXAID_BIN="${FLEXAID_BINARY:-}"
if [[ -z "${FLEXAID_BIN}" ]]; then
    # Common build locations inside the repo
    for cand in \
        "${REPO_ROOT}/BIN/FlexAID" \
        "${REPO_ROOT}/build/FlexAID" \
        "${REPO_ROOT}/build-test/FlexAID" \
        "$(command -v FlexAID 2>/dev/null || true)"
    do
        if [[ -n "$cand" && -x "$cand" ]]; then
            FLEXAID_BIN="$cand"
            break
        fi
    done
fi

# ─── Execute the actual docking run (capture EVERYTHING to .raw.txt) ──────────
# We deliberately do not use "set -e" around the run so we can inspect the exit
# code and still produce the JSON + optional bundle on partial success.
info "Receptor : ${RECEPTOR}"
info "Ligand   : ${LIGAND}"
info "Results  : ${RESULTS_DIR}"
info "Binary   : ${FLEXAID_BIN:-<python-flexaidds-dock>}"
info "CREATE_BUNDLE=${CREATE_BUNDLE:-0}"

RUN_EXIT=0
{
    echo "=== FlexAIDdS run started $(date -Iseconds) ==="
    echo "RECEPTOR=${RECEPTOR}"
    echo "LIGAND=${LIGAND}"
    echo "TEMPERATURE=${TEMPERATURE}"
    echo "SEED=${SEED}"
    echo "FLEXAID_BINARY=${FLEXAID_BIN:-python}"
    echo "CREATE_BUNDLE=${CREATE_BUNDLE:-0}"
    echo "=== engine output begins ==="
} > "${RAW_PATH}"

if [[ -n "${FLEXAID_BIN}" && -x "${FLEXAID_BIN}" ]]; then
    # Legacy binary path (captured)
    set +e
    "${FLEXAID_BIN}" \
        --receptor "${RECEPTOR}" \
        --ligand "${LIGAND}" \
        --temperature "${TEMPERATURE}" \
        --seed "${SEED}" \
        "${EXTRA_ARGS[@]}" \
        >> "${RAW_PATH}" 2>&1
    RUN_EXIT=$?
    set -e
else
    # Preferred modern path: use the Python high-level API (produces structured + raw)
    # We still capture a .raw.txt for the bundle contract.
    set +e
    "$PYTHON" -c '
import sys, os, datetime
sys.path.insert(0, os.environ.get("PYTHONPATH", "").split(":")[0] or ".")
try:
    import flexaidds as fd
except Exception as e:
    print("ERROR: could not import flexaidds:", e, file=sys.stderr)
    sys.exit(2)

rec = sys.argv[1]
lig = sys.argv[2]
temp = float(sys.argv[3])
seed = int(sys.argv[4])

print(f"[python] FlexAIDdS dock start {datetime.datetime.utcnow().isoformat()}Z", flush=True)
pop = fd.dock(receptor=rec, ligand=lig, temperature=temp, seed=seed, compute_entropy=True)
print(f"[python] produced {len(pop.binding_modes)} binding mode(s)", flush=True)

# Also emit a tiny raw-style summary so .raw.txt is never empty
for i, m in enumerate(pop.rank_by_free_energy()[:3]):
    print(f"  mode {i}: ΔG={m.free_energy:.3f} kcal/mol  n_poses={m.n_poses}")
' "${RECEPTOR}" "${LIGAND}" "${TEMPERATURE}" "${SEED}" \
    >> "${RAW_PATH}" 2>&1
    RUN_EXIT=$?
    set -e
fi

echo "=== run exit code: ${RUN_EXIT} ===" >> "${RAW_PATH}"

# ─── Build the reproducibility JSON (always attempted on any completion) ──────
# Uses the skill's own run_metadata module so the schema stays in sync.
info "Building reproducibility record..."

# We export a few vars so the python snippet can see them without fragile argv tricks
export RECEPTOR_ARG="${RECEPTOR}"
export LIGAND_ARG="${LIGAND}"
export RUN_EXIT
export RESULTS_DIR

"$PYTHON" - <<'PYEOF' > "${JSON_PATH}" 2>/dev/null || true
import os, sys, json, datetime, hashlib
from pathlib import Path

def _sha(p):
    try:
        return hashlib.sha256(Path(p).read_bytes()).hexdigest() if p and Path(p).exists() else None
    except Exception:
        return None

# The script already put REPO_ROOT/python on PYTHONPATH
try:
    from flexaidds.run_metadata import create_run_record
except Exception:
    # Ultra-minimal fallback so the bundle contract (JSON always present) is still honored
    def create_run_record(docking_results=None, **_):
        prov = {
            "timestamp": datetime.datetime.utcnow().isoformat() + "Z",
            "hostname": os.uname().nodename if hasattr(os, "uname") else "unknown",
            "kernel": "unknown",
            "python_version": sys.version.split()[0],
            "cmake_version": None,
            "flexaidds_git": {"commit": None, "dirty": None, "branch": None, "error": "run_metadata unavailable"},
            "binary": {"path": os.environ.get("FLEXAID_BINARY"), "sha256": None, "error": "unavailable"},
            "used_skip_rebuild": os.environ.get("SKIP_REBUILD", "0").lower() in ("1", "true", "yes"),
            "skill": "flexaidds-docking/run_flexaidds.sh"
        }
        rec = {
            "_summary": f"FlexAIDdS run (fallback) {prov['timestamp']}",
            "overall_run_confidence": {"score": 0.0, "level": "low", "values_extracted": 0, "reason": "run_metadata import failed"},
            "provenance": prov,
            "input_hashes": {
                "receptor": _sha(os.environ.get("RECEPTOR_ARG")),
                "ligand": _sha(os.environ.get("LIGAND_ARG"))
            }
        }
        if docking_results:
            rec.update(docking_results)
        return rec

docking_blob = {
    "receptor": os.environ.get("RECEPTOR_ARG", ""),
    "ligand": os.environ.get("LIGAND_ARG", ""),
    "exit_code": int(os.environ.get("RUN_EXIT", "0") or 0),
    "results_dir": os.environ.get("RESULTS_DIR", ""),
    "raw_output": "dock.raw.txt"
}

rec = create_run_record(
    docking_results=docking_blob,
    binary_path=Path(os.environ.get("FLEXAID_BINARY")) if os.environ.get("FLEXAID_BINARY") else None,
    repo_root=Path(os.environ.get("FLEXAIDDS_SOURCE", "")) or None
)

# Guarantee the input_hashes block exists (required by bundle spec)
if "input_hashes" not in rec:
    rec["input_hashes"] = {
        "receptor": _sha(os.environ.get("RECEPTOR_ARG")),
        "ligand": _sha(os.environ.get("LIGAND_ARG"))
    }

print(json.dumps(rec, indent=2, sort_keys=False))
PYEOF

# If the python step produced literally nothing, write a last-resort valid record
if [[ ! -s "${JSON_PATH}" ]]; then
    cat > "${JSON_PATH}" <<EOF
{
  "_summary": "FlexAIDdS run ${TIMESTAMP} (emergency record)",
  "overall_run_confidence": {"score": 0.0, "level": "low", "values_extracted": 0, "reason": "JSON writer produced no output"},
  "provenance": {
    "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
    "hostname": "$(hostname)",
    "kernel": "unknown",
    "python_version": "$($PYTHON --version 2>/dev/null | awk '{print $2}')",
    "flexaidds_git": {"commit": null, "dirty": null, "error": "emergency fallback"},
    "binary": {"path": "${FLEXAID_BIN:-}", "sha256": null},
    "used_skip_rebuild": ${SKIP_REBUILD:-false},
    "skill": "flexaidds-docking/run_flexaidds.sh"
  },
  "receptor": "${RECEPTOR}",
  "ligand": "${LIGAND}",
  "exit_code": ${RUN_EXIT},
  "input_hashes": {
    "receptor": "$(compute_sha256 "${RECEPTOR}")",
    "ligand": "$(compute_sha256 "${LIGAND}")"
  }
}
EOF
fi


# ─── Success gate + optional bundle (THE ONLY PLACE THAT TOUCHES CREATE_BUNDLE) ─
SUCCESS=false
if [[ "${RUN_EXIT}" -eq 0 && -s "${JSON_PATH}" ]]; then
    SUCCESS=true
    ok "Docking run completed successfully."
else
    warn "Docking run reported non-zero exit or produced no JSON (exit=${RUN_EXIT})."
fi

# === CREATE_BUNDLE feature (strictly opt-in, strictly success-only) ============
if [[ "${CREATE_BUNDLE:-0}" == "1" ]]; then
    if [[ "${SUCCESS}" == "true" && -f "${JSON_PATH}" && -f "${RAW_PATH}" ]]; then
        create_reproducibility_bundle \
            "${RESULTS_DIR}" \
            "${JSON_PATH}" \
            "${RAW_PATH}" \
            "${RECEPTOR}" \
            "${LIGAND}"
    else
        warn "CREATE_BUNDLE=1 was set but run did not succeed or artifacts are missing — no bundle written (non-breaking)."
    fi
fi
# ===============================================================================

# === VISUALIZE / Imagine figure gate (P1 of grok-imagine-figure-rendering) =====
# After success + (when requested) after the bundle gate, invoke the pure-Python
# figures prep. This writes results/figures/{prompt_cover.txt, prompt_animation.txt,
# figure_metadata.json} with real ΔG/ΔH/-TΔS from the ledger + Gate 6 status.
# The skill agent (or caller) is then expected to feed the prompts to image_gen /
# video_gen / image_edit and place the materialized assets in the same dir.
# Only runs on SUCCESS path. Non-breaking when VISUALIZE=0 (default).
if [[ "${SUCCESS}" == "true" && "${VISUALIZE:-0}" == "1" ]]; then
    info "VISUALIZE=1 set — preparing publication cover + 6s animation prompts (Gate 6 aware)..."
    # Reuse the same python that was located earlier; RESULTS_DIR is already absolute
    "$PYTHON" -c '
import os, sys
from pathlib import Path
sys.path.insert(0, os.environ.get("PYTHONPATH", "").split(":")[0] or ".")
try:
    from flexaidds.figures import prepare_publication_figures
except Exception as e:
    print("[figures] import failed:", e, file=sys.stderr)
    sys.exit(0)  # never break the run

rd = Path(os.environ.get("RESULTS_DIR", "."))
res = prepare_publication_figures(
    rd,
    visualize=True,
    require_gate6=True,
    use_pymol_base=False,
    force=False,
)
if res.get("proceeded"):
    print("[figures] prepared:", res.get("figures_dir"))
    print("[figures] cover prompt:", res.get("cover_prompt_path"))
    print("[figures] animation prompt:", res.get("animation_prompt_path"))
    print("[figures] metadata:", res.get("metadata_path"))
    print("[figures] gate6_passed:", res.get("gate6_passed"))
    print("Agent: now call image_gen / video_gen with the prompts above (aspect 3:2 for cover, duration=6 for anim). Save outputs inside the figures/ dir as cover_best_mode.png + animation_6s.mp4. Use image_edit for banner/equation polish if needed.")
else:
    print("[figures] skipped (", res.get("skipped", "unknown"), ")")
' 2>&1 || warn "figures preparation encountered an error (non-fatal; run continues)"
fi
# ===============================================================================

# Final status for callers / CI
if [[ "${SUCCESS}" == "true" ]]; then
    ok "Results + reproducibility JSON written to: ${RESULTS_DIR}"
    echo "${RESULTS_DIR}"   # machine-readable last line for scripts that want to consume it
    exit 0
else
    fail "Run failed (see ${RAW_PATH} and ${JSON_PATH})"
    exit "${RUN_EXIT}"
fi
