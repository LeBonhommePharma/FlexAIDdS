#!/bin/bash
#
# launch_full_benchmark.sh — Canonical launcher for full FlexAIDdS benchmark runs
# on M3 Pro (iCloud-only results storage).
#
# This is the ONLY supported way to start production full first runs
# (especially Non-Native or at custom temperatures like 298 K / 310 K).
#
# Usage:
#   python3 .grok/skills/flexaidds/scripts/launch_full_benchmark.sh \
#       <dataset> <temperature> <results_subdir_name>
#
# The script enforces:
#   - Full skill ritual (validate + ensure, no --quick for production)
#   - Correct PATH + benchmark_datasets availability
#   - Strict pre-flight (especially for heavy Non-Native runs)
#   - Early run_status.json with PID, temperature, command, timestamps
#   - Proper detached execution + reliable logging
#   - All artifacts land only under $FLEXAIDDS_RESULTS (iCloud)
#
# After launch it prints simple monitoring commands.
#
# This script lives inside the flexaidds skill so it is always available
# and versioned with the skill.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"   # repo root

# --- 1. Arguments (check before env so --help-style usage works) ------------
if [ $# -ne 3 ]; then
    echo "Usage: $0 <dataset> <temperature> <results_subdir_name>"
    echo "Example: $0 astex_nonnative 310 astex_nonnative_310K"
    echo ""
    echo "Requires benchmark_datasets on PATH. Set FLEXAIDDS_BUILD or source ~/.flexaidds_env"
    exit 1
fi

# --- 2. Resolve build (reject stale env, pin SHA256) -------------------------
DATASET="$1"
TEMPERATURE="$2"
SUBDIR_NAME="$3"

if [ -f "${HOME}/.flexaidds_env" ]; then
    # shellcheck disable=SC1091
    source "${HOME}/.flexaidds_env"
fi

echo "=== [Skill] Resolving active FlexAIDdS build (SHA-pinned) ==="
if ! eval "$(python3 "$SCRIPT_DIR/resolve_build.py" --export-shell --repo-root "$REPO_ROOT")"; then
    echo "FATAL: resolve_build.py failed — no valid build tree."
    echo "       Fix: python3 $SCRIPT_DIR/resolve_build.py --sync-env"
    exit 1
fi

python3 "$SCRIPT_DIR/resolve_build.py" --write-pin --repo-root "$REPO_ROOT" >/dev/null

if ! command -v benchmark_datasets >/dev/null 2>&1; then
    echo "FATAL: benchmark_datasets not found after resolve_build."
    echo "       Expected at: ${FLEXAIDDS_RUNNER:-$FLEXAIDDS_BUILD/benchmark_datasets}"
    exit 1
fi

FLEXAIDDS_RUNNER="${FLEXAIDDS_RUNNER:-$FLEXAIDDS_BUILD/benchmark_datasets}"
FLEXAIDDS_BINARY="${FLEXAIDDS_BINARY:-$FLEXAIDDS_BUILD/FlexAIDdS}"

RUNNING_PIDS="$(pgrep -f "benchmark_datasets" 2>/dev/null || true)"
if [ -n "$RUNNING_PIDS" ]; then
    echo "FATAL: benchmark_datasets already running (PIDs: $RUNNING_PIDS)"
    echo "       Monitor or stop the existing run before launching a duplicate."
    exit 1
fi

echo "  build_dir:       $FLEXAIDDS_BUILD"
echo "  engine_sha256:   $FLEXAIDDS_ENGINE_SHA256"
echo "  runner_sha256:   $FLEXAIDDS_RUNNER_SHA256"

FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$HOME/flexaidds_results}"

# --- 3. Full Skill Ritual (mandatory for production) ------------------------
echo "=== [Skill] Full ritual for production full benchmark run ==="
python3 "$SCRIPT_DIR/validate_skill.py"
python3 "$SCRIPT_DIR/ensure_docking_data.py"   # no --quick

# --- 4. Strict pre-flight ---------------------------------------------------
echo "=== Pre-flight checks ==="

if [[ "$TEMPERATURE" != "298" && "$TEMPERATURE" != "310" ]]; then
    echo "WARNING: Temperature $TEMPERATURE is not one of the two standard values (298/310)."
fi

# For Non-Native / heavy datasets, require the full data set
if [[ "$DATASET" == *"nonnative"* || "$DATASET" == *"non_native"* ]]; then
    echo "Heavy dataset ($DATASET) detected — requiring complete support files..."
    MISSING=()
    for f in Lovell_LIB.dat rotobs.lst SYBYL_emat.dat; do
        if [ ! -f "$FLEXAIDDS_BUILD/$f" ] && [ ! -f "$FLEXAIDDS_ICLOUD/build/$f" ]; then
            MISSING+=("$f")
        fi
    done
    if [ ${#MISSING[@]} -gt 0 ]; then
        echo "FATAL: Missing required files for full Non-Native run: ${MISSING[*]}"
        echo "       Re-run ensure_docking_data.py without --quick or use a full data checkout."
        exit 1
    fi
fi

# --- Metal / Hardware Acceleration pre-flight (M3 Pro & Apple Silicon) -------
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Apple platform (M3 Pro etc.) — Metal hardware acceleration pre-flight..."
    METAL_COUNT=0
    if [ -d "$FLEXAIDDS_BUILD" ]; then
        METAL_COUNT=$(ls "$FLEXAIDDS_BUILD"/*.metallib 2>/dev/null | wc -l | tr -d ' ')
    fi
    if [ "$METAL_COUNT" -gt 0 ]; then
        echo "  ✓ $METAL_COUNT .metallib shaders present (ShannonEntropy, CavityDetect, MetalRMSD, TurboQuant, etc.)"
        echo "  → UnifiedHardwareDispatch / ShannonThermoStack / tENCoM / CavityDetect will use Metal when kernels execute."
    else
        echo "  ⚠ No .metallib shaders found — Metal acceleration will fall back to CPU paths."
    fi
    if command -v system_profiler >/dev/null 2>&1; then
        METAL_INFO=$(system_profiler SPDisplaysDataType 2>/dev/null | grep -E 'Metal Support|Chipset Model' | head -2 | xargs)
        echo "  Runtime: $METAL_INFO"
    fi
    echo "  Build was configured with FLEXAIDS_USE_METAL=ON (see CMakeCache.txt in build dir)."

    # Binary-level linkage check (stronger confirmation)
    if [ -x "$FLEXAIDDS_BINARY" ]; then
        if otool -L "$FLEXAIDDS_BINARY" 2>/dev/null | grep -qi metal; then
            echo "  ✓ $FLEXAIDDS_BINARY links against Metal.framework — Metal acceleration linked in."
        else
            echo "  ⚠ $FLEXAIDDS_BINARY does not appear to link Metal.framework directly."
        fi
    fi

    # Ready-to-use verification command for this run's logs (for the user to run later)
    echo "  To verify Metal usage in this run's logs after the C++ kernels execute:"
    echo "    grep -iE 'metal|backend|dispatch|shannon|using metal' \"\$LOG_FILE\" | tail -20"
fi

# --- 5. Prepare output directory on iCloud ----------------------------------
TS=$(date +%s)
OUT_DIR="$FLEXAIDDS_RESULTS/full-${TEMPERATURE}K-${SUBDIR_NAME}-${TS}"
mkdir -p "$OUT_DIR"

CONFIG_FILE="$OUT_DIR/config_${TEMPERATURE}.json"
cat > "$CONFIG_FILE" << EOF
{
  "temperature": $TEMPERATURE,
  "launched_via": "flexaidds skill launcher",
  "timestamp": "$(date -Iseconds)"
}
EOF

# --- 6. Early structured status (the key bulletproofing item) ---------------
STATUS_FILE="$OUT_DIR/run_status.json"
python3 - <<PYEOF
import json, os
from datetime import datetime, timezone
payload = {
    "status": "launched",
    "wrapper_pid": $$,
    "dataset": "$DATASET",
    "temperature": $TEMPERATURE,
    "start_time": datetime.now(timezone.utc).isoformat(),
    "command": "$0 $*",
    "output_dir": "$OUT_DIR",
    "build_dir": os.environ.get("FLEXAIDDS_BUILD", ""),
    "engine_path": os.environ.get("FLEXAIDDS_BINARY", ""),
    "runner_path": os.environ.get("FLEXAIDDS_RUNNER", ""),
    "engine_sha256": os.environ.get("FLEXAIDDS_ENGINE_SHA256", ""),
    "runner_sha256": os.environ.get("FLEXAIDDS_RUNNER_SHA256", ""),
}
with open("$STATUS_FILE", "w") as f:
    json.dump(payload, f, indent=2)
PYEOF

echo "Early status written: $STATUS_FILE"

# --- 7. Launch (properly detached + logged) ---------------------------------
LOG_FILE="$OUT_DIR/binary.log"
ERR_FILE="$OUT_DIR/stderr.log"

echo "=== Launching real benchmark (this can take a long time) ==="
echo "Dataset:     $DATASET"
echo "Temperature: $TEMPERATURE K"
echo "Output:      $OUT_DIR"
echo "Log:         $LOG_FILE"

# Portable detachment (macOS + Linux):
# - nohup for survival after logout
# - no 'setsid' (Linux-only; it was causing the command/redirection to fail silently on macOS M3 Pro)
# - </dev/null for clean stdin
# - disown to fully detach from job control (best effort)
CACHE_DIR="${FLEXAIDDS_CACHE:-${FLEXAIDDS_CACHE_ROOT:-$HOME/.flexaidds/benchmarks}}"
mkdir -p "$CACHE_DIR"

nohup "$FLEXAIDDS_RUNNER" \
    --benchmark "$DATASET" \
    --temperature "$TEMPERATURE" \
    --output "$OUT_DIR/${DATASET}_${TEMPERATURE}" \
    --threads 1 \
    --omp-threads "${OMP_NUM_THREADS:-6}" \
    --cache "$CACHE_DIR" \
    --job-timeout-seconds 7200 \
    >> "$LOG_FILE" 2>> "$ERR_FILE" </dev/null &

CHILD_PID=$!
disown $CHILD_PID 2>/dev/null || true

# Record the real child PID
python3 - <<PYEOF
import json
with open("$STATUS_FILE") as f:
    data = json.load(f)
data["benchmark_runner_pid"] = $CHILD_PID
with open("$STATUS_FILE", "w") as f:
    json.dump(data, f, indent=2)
PYEOF

echo ""
echo "=== LAUNCH COMPLETE (using flexaidds skill launcher) ==="
echo ""
echo "Monitor with:"
echo "  tail -f $LOG_FILE"
echo "  tail -f $ERR_FILE"
echo ""
echo "Health / PID:"
echo "  cat $STATUS_FILE"
echo ""
echo "When finished, run the usual post-run verification:"
echo "  python3 .grok/skills/flexaidds/scripts/validate_skill.py"
echo "  python benchmarks/re-dock/icloud_fs_check.py --path $OUT_DIR"
echo ""
echo "You can safely log out. The real work is detached."
echo ""

# Mark as running
python3 - <<PYEOF
import json, time
with open("$STATUS_FILE") as f:
    data = json.load(f)
data["status"] = "running"
data["launcher_finished"] = time.time()
with open("$STATUS_FILE", "w") as f:
    json.dump(data, f, indent=2)
PYEOF

# --- P1: Enhanced post-run guidance for "best BindingMode" exact answer + trap ---
echo ""
echo "=== [P1] Post-run commands for best BindingMode validity + extract (after child done): ==="
echo "  python3 .grok/skills/flexaidds/scripts/summarize_campaign.py \"$OUT_DIR\" --verbose"
echo "  # (when extended) --extract-best-mode to get top free_energy BindingMode + full thermo ledger (the exact requested answer)"
echo "  python3 .grok/skills/flexaidds/scripts/validate_skill.py"
echo "  python benchmarks/re-dock/icloud_fs_check.py --path \"$OUT_DIR\""
echo "  cat \"$OUT_DIR/run_status.json\""
echo "  grep -iE 'metal|backend|dispatch|shannon|using metal|success rate|RMSD|Binding Mode' \"$LOG_FILE\" | tail -30 || true"
echo ""

# Trap: ensure status file reflects launcher exit (non-fatal, for cases where launcher is monitored). Robust to spaces in iCloud paths.
trap 'python3 -c "
import json, time, os, sys
sf = os.environ.get(\"STATUS_FILE\", \"\")
if sf and os.path.exists(sf):
  try:
    with open(sf) as f: data = json.load(f)
    if data.get(\"status\") in (\"launched\", \"running\"):
      data[\"status\"] = \"launcher_exited\"
      data[\"launcher_end_time\"] = time.time()
      with open(sf, \"w\") as f: json.dump(data, f, indent=2)
  except Exception:
    pass
" ' EXIT INT TERM

# Note: full --extract-best and strict validity gate implemented in summarize (P1.2) and dedicated helper.

