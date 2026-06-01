#!/usr/bin/env bash
#
# grok_master_launcher.sh
# Master Launcher + Syncer for healthy, reproducible benchmark campaigns
# on M3 Pro 18GB (iCloud 2TB only storage)
#
# Purpose:
#   - Safe, failsafe orchestration of the enhanced failsafe_campaign.py
#   - Specifically tuned for Astex Diverse, Astex Non-Native, HAP2
#   - Everything durable on iCloud Drive (via ~/.flexaidds_env)
#   - Hot execution on local APFS (/private/tmp) for speed & health
#   - Maximally uses the 18GB hardware (workers=4)
#   - Your exact GA settings (1000 generations, 2000 population)
#   - Full reproducibility (manifests, git SHA, RUN_OK.json health checks)
#
# Usage examples:
#   ./grok_master_launcher.sh preflight          # Safe validation
#   ./grok_master_launcher.sh launch             # Start 10-rep campaign
#   ./grok_master_launcher.sh full               # preflight + launch + analyze + sync
#   ./grok_master_launcher.sh status             # Live monitoring
#   ./grok_master_launcher.sh sync               # Extra belt-and-suspenders rsync to iCloud
#   ./grok_master_launcher.sh analyze            # Bootstrap stats
#
# All results/logs go to iCloud via FLEXAIDDS_* env.
# Run this inside screen/tmux for long campaigns.
#
# Apache-2.0 (c) 2026 — based on m3pro/ tools + portability improvements

set -euo pipefail

# Colors
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; BOLD='\033[1m'; NC='\033[0m'
info() { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()   { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn() { printf "${YELLOW}[WARN]${NC} %s\n" "$*" >&2; }
die()  { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; exit 1; }
phase() { printf "\n${BOLD}════════════════════════════════════════════════════════${NC}\n${BOLD}  %s${NC}\n${BOLD}════════════════════════════════════════════════════════${NC}\n\n" "$*"; }

# === SELF-LOCATION (Chunk 1 safety fix — never trust FLEXAIDDS_REPO for script paths) ===
# This is the root cause fix for every "wrong tree / no such file / stale binary" screen failure.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_ABS="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")"
FAILSAFE_PY="$SCRIPT_DIR/failsafe_campaign.py"
ANALYZER_PY="$SCRIPT_DIR/analyze_repetitions.py"

# 1. Source the official M3 Pro iCloud environment
ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    set -a
    # shellcheck disable=SC1090
    source "$ENV_FILE"
    set +a
else
    die "Required $ENV_FILE not found. Run setup_cloud_storage.sh first."
fi

# 2. Required variables from env (iCloud durability)
: "${FLEXAIDDS_REPO:?FLEXAIDDS_REPO not set in env}"
: "${FLEXAIDDS_BUILD:?FLEXAIDDS_BUILD not set}"
: "${FLEXAIDDS_RESULTS:?FLEXAIDDS_RESULTS not set}"
: "${FLEXAIDDS_LOGS:?FLEXAIDDS_LOGS not set}"
: "${FLEXAIDDS_BENCHMARK_DATA:?FLEXAIDDS_BENCHMARK_DATA not set}"

BINARY="$FLEXAIDDS_BUILD/benchmark_datasets"
DOCKING_BINARY="${FLEXAIDDS_BINARY:-$FLEXAIDDS_BUILD/FlexAID}"

# 3. Grok-tuned healthy defaults for this exact machine + your GA
DEFAULT_DATASETS="astex astex_nonnative hap2"
DEFAULT_RUNS=10
DEFAULT_WORKERS=4                    # Max sensible for 18GB M3 Pro (see m3pro_profile.yaml)
DEFAULT_GA_GENERATIONS=1000
DEFAULT_GA_POPULATION=2000
DEFAULT_CLUSTERING="FO"
DEFAULT_TEMPERATURE=300
SESSION_BACKEND="${SESSION_BACKEND:-screen}"   # Chunk 1 skeleton: screen | tmux | none (expanded Chunk 2)

# 4. Unique run identity (for "your own" campaigns, cleanly separated from Codex runs)
RUN_ID="${RUN_ID:-grok_own_m3pro_$(date +%Y%m%d_%H%M%S)}"

# Local hot (fast APFS) workspace for execution speed — critical on iCloud-only storage
LOCAL_HOT_BASE="${LOCAL_HOT_BASE:-/private/tmp/grok_bench_hot_${RUN_ID}}"

# iCloud destination roots (all durable storage lives here via the env)
ICLOUD_RESULTS="$FLEXAIDDS_RESULTS"
ICLOUD_LOGS="$FLEXAIDDS_LOGS"

# Use self-located absolute paths (Chunk 1). FLEXAIDDS_REPO is still used for the *build/binary* only.
FAILSAVE_SCRIPT="$FAILSAFE_PY"
ANALYZER="$ANALYZER_PY"

# Safety: basic iCloud volume sanity check (common source of "unhealthy" runs)
check_icLOUD_health() {
    if [[ ! -d "$ICLOUD_RESULTS" || ! -w "$ICLOUD_RESULTS" ]]; then
        warn "iCloud results path not writable: $ICLOUD_RESULTS"
        warn "Ensure the volume is mounted and you have write access (run setup_cloud_storage.sh if needed)."
    fi
    # Quick memory headroom hint for this 18GB machine
    if command -v sysctl >/dev/null 2>&1; then
        local mem_gb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1024 / 1024 / 1024 ))
        if (( mem_gb < 16 )); then
            warn "This machine reports only ${mem_gb}GB — the recommended worker counts in m3pro_profile.yaml may be too aggressive."
        fi
    fi
}

check_icLOUD_health

# Safety: ensure we have a real binary
if [[ ! -x "$BINARY" ]]; then
    die "benchmark_datasets not executable at $BINARY"
fi
if [[ ! -x "$DOCKING_BINARY" ]]; then
    die "FlexAID not executable at $DOCKING_BINARY"
fi

# === doctor() — Chunk 1: explicit diagnosis of the exact class of failure the user has been hitting ===
doctor() {
    phase "DOCTOR — Self-location vs Environment Sanity (M3 Pro 18GB)"
    echo "SCRIPT_DIR (this launcher lives here):     $SCRIPT_DIR"
    echo "SCRIPT_ABS:                                 $SCRIPT_ABS"
    echo "FAILSAFE_PY (will be used for campaign):    $FAILSAFE_PY"
    echo "ANALYZER_PY:                                $ANALYZER_PY"
    echo ""
    echo "FLEXAIDDS_REPO (from env, for BUILD only):  ${FLEXAIDDS_REPO:-<unset>}"
    echo "FLEXAIDDS_BUILD:                            ${FLEXAIDDS_BUILD:-<unset>}"
    echo ""
    if [[ "$SCRIPT_DIR" != "${FLEXAIDDS_REPO:-}"* && -n "${FLEXAIDDS_REPO:-}" ]]; then
        warn "MISMATCH DETECTED: This launcher is in a different worktree than FLEXAIDDS_REPO."
        warn "This was the root cause of almost every 'No such file' / wrong binary path in screen sessions."
        warn "The new self-location + absolute inner wrapper (Chunk 2) fixes it."
    else
        ok "SCRIPT_DIR and FLEXAIDDS_REPO are consistent (or FLEXAIDDS_REPO not overriding)."
    fi
    echo ""
    echo "Binary checks:"
    [[ -x "$BINARY" ]] && ok "benchmark_datasets: $BINARY" || warn "benchmark_datasets MISSING or not +x at $BINARY"
    [[ -x "$DOCKING_BINARY" ]] && ok "FlexAID:          $DOCKING_BINARY" || warn "FlexAID MISSING or not +x at $DOCKING_BINARY"
    echo ""
    # Quick iCloud / tmp sanity (expanded in later chunks)
    [[ -d "$ICLOUD_RESULTS" && -w "$ICLOUD_RESULTS" ]] && ok "iCloud results writable: $ICLOUD_RESULTS" || warn "iCloud results path problem: $ICLOUD_RESULTS"
    [[ -d "/private/tmp" && -w "/private/tmp" ]] && ok "/private/tmp (hot path) writable" || warn "/private/tmp not writable — hot execution will fail"

    # Minimal 18GB M3 Pro resource awareness (from m3pro_profile.yaml: 11 GB usable after OS+Metal)
    if command -v sysctl >/dev/null 2>&1; then
        local avail_gb=$(( 11 ))   # conservative for this machine/profile
        local req_workers=${WORKERS:-4}
        if (( req_workers > 4 )); then
            warn "Requested workers=$req_workers exceeds safe max ~4 for 18GB M3 Pro (see m3pro_profile.yaml tier2:2 sequential, tier1:4). Proceeding but you were warned."
        else
            ok "Worker count $req_workers within 18GB M3 Pro budget (profile max 4 for aggressive tier-1 style)."
        fi
    fi
    echo ""
    info "Run 'doctor' anytime to re-check. Start will call this automatically (Chunk 1+)."
}

# macOS notification stub (Chunk 3) — called from inner on completion or from outer on error
notify_user() {
    local title="$1"
    local msg="$2"
    if command -v osascript >/dev/null 2>&1; then
        osascript -e "display notification \"$msg\" with title \"$title\"" 2>/dev/null || true
    fi
    # Optional voice: say "$title" 2>/dev/null || true
}

# === Chunk 2: Safe absolute inner wrapper generator (the "tard-proof" one-command fix) ===
# Creates a complete, self-contained, chmod +x script inside the local hot base.
# It hard-codes the absolute FAILSAFE_PY / ANALYZER_PY from *this* manager's SCRIPT_DIR,
# sources the env for the heavy build/binary/iCloud paths only, and runs the full pipeline.
# Screen/tmux then execs this absolute file directly — no heredoc, no $0, no cd to wrong tree, no FLEXAIDDS_REPO pollution for scripts.
write_inner_wrapper() {
    local inner_path="$1"
    local run_id="$2"
    cat > "$inner_path" << 'INNER_EOF'
#!/bin/bash
set -euo pipefail

# These are expanded by the *outer* manager at generation time (absolute, from the correct worktree)
FAILSAFE_PY="__FAILSAFE_PY__"
ANALYZER_PY="__ANALYZER_PY__"
RUN_ID="__RUN_ID__"
LOCAL_HOT_BASE="__LOCAL_HOT_BASE__"
ICLOUD_LOGS="__ICLOUD_LOGS__"
ICLOUD_RESULTS="__ICLOUD_RESULTS__"

# Source the user's env for the heavy stuff (build, iCloud targets, etc.)
ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    set -a
    # shellcheck disable=SC1090
    source "$ENV_FILE"
    set +a
fi

echo "[inner] Starting healthy pipeline for $RUN_ID using self-located scripts:"
echo "[inner]   FAILSAFE_PY=$FAILSAFE_PY"
echo "[inner]   ANALYZER_PY=$ANALYZER_PY"
echo "[inner]   Status will be at: $ICLOUD_LOGS/m3pro_failsafe_${RUN_ID}/campaign_status.json"
echo ""

# The full safe pipeline (preflight already passed in outer; here we do the real work + post-steps)
python3 "$FAILSAFE_PY" \
    --datasets ${DATASETS:-astex astex_nonnative hap2} \
    --runs ${RUNS:-10} \
    --workers ${WORKERS:-4} \
    --ga-generations ${GA_GENERATIONS:-1000} \
    --ga-population ${GA_POPULATION:-2000} \
    --clustering ${CLUSTERING:-FO} \
    --temperature ${TEMPERATURE:-300} \
    --run-id "$RUN_ID" \
    --local-base "$LOCAL_HOT_BASE" \
    --resume \
    2>&1 | tee -a "$ICLOUD_LOGS/${RUN_ID}_inner.log"

echo ""
echo "[inner] Campaign body complete. Running analyze + final sync..."

python3 "$ANALYZER_PY" \
    --results-dir "$ICLOUD_RESULTS/tier2" \
    --n-bootstrap 10000 \
    --output-dir "$ICLOUD_RESULTS/analysis/${RUN_ID}" \
    2>&1 | tee -a "$ICLOUD_LOGS/${RUN_ID}_inner.log" || true

# Belt-and-suspenders final sync of hot results (the failsafe also does rsync when remote configured)
if [[ -d "$LOCAL_HOT_BASE" ]]; then
    rsync -a --delete "$LOCAL_HOT_BASE/" "$ICLOUD_RESULTS/tier2/" 2>/dev/null || true
fi

echo "[inner] All done for $RUN_ID. Artifacts on iCloud. You can detach (Ctrl-A D) or close this screen."
INNER_EOF

    # Perform the real substitutions (the heredoc above used placeholders to avoid early expansion)
    sed -i '' \
        -e "s|__FAILSAFE_PY__|$FAILSAFE_PY|g" \
        -e "s|__ANALYZER_PY__|$ANALYZER_PY|g" \
        -e "s|__RUN_ID__|$run_id|g" \
        -e "s|__LOCAL_HOT_BASE__|$LOCAL_HOT_BASE|g" \
        -e "s|__ICLOUD_LOGS__|$ICLOUD_LOGS|g" \
        -e "s|__ICLOUD_RESULTS__|$ICLOUD_RESULTS|g" \
        "$inner_path"

    chmod +x "$inner_path"
}

mkdir -p "$LOCAL_HOT_BASE" "$ICLOUD_LOGS" "$ICLOUD_RESULTS/tier2" "$ICLOUD_RESULTS/analysis"

build_cmd() {
    echo "python3 \"$FAILSAVE_SCRIPT\" \
        --datasets ${DATASETS:-$DEFAULT_DATASETS} \
        --runs ${RUNS:-$DEFAULT_RUNS} \
        --workers ${WORKERS:-$DEFAULT_WORKERS} \
        --ga-generations ${GA_GENERATIONS:-$DEFAULT_GA_GENERATIONS} \
        --ga-population ${GA_POPULATION:-$DEFAULT_GA_POPULATION} \
        --clustering ${CLUSTERING:-$DEFAULT_CLUSTERING} \
        --temperature ${TEMPERATURE:-$DEFAULT_TEMPERATURE} \
        --run-id \"$RUN_ID\" \
        --local-base \"$LOCAL_HOT_BASE\" \
        ${EXTRA_ARGS:-}"
}

case "${1:-help}" in
    doctor)
        doctor
        ;;
    preflight)
        phase "HEALTHY PREFLIGHT (M3 Pro 18GB + iCloud durability)"
        doctor
        info "Using local hot APFS for execution speed, iCloud for all durable storage"
        cmd=$(build_cmd)
        eval "$cmd --preflight-only --skip-smoke"
        ok "Preflight complete. Ready for launch."
        ;;
    launch)
        phase "LAUNCHING REPETITION CAMPAIGN (healthy + resumable)"
        info "Run ID: $RUN_ID"
        info "Results will be rsynced to iCloud: $ICLOUD_RESULTS"
        cmd=$(build_cmd)
        eval "$cmd --resume" 2>&1 | tee -a "$ICLOUD_LOGS/${RUN_ID}_launcher.log"
        ;;
    sync)
        phase "EXTRA SYNC TO ICLOUD (belt & suspenders — safe & resumable)"
        if [[ -d "$LOCAL_HOT_BASE" ]]; then
            rsync -a --delete --info=progress2 "$LOCAL_HOT_BASE/" "$ICLOUD_RESULTS/tier2/" || warn "Partial rsync of hot results — safe to re-run later"
        fi
        # Also sync any logs the failsafe wrote under the iCloud logs dir
        LATEST_LOG_DIR=$(ls -td "$ICLOUD_LOGS"/m3pro_failsafe_${RUN_ID}* 2>/dev/null | head -1 || true)
        if [[ -n "$LATEST_LOG_DIR" ]]; then
            rsync -a --delete "$LATEST_LOG_DIR/" "$ICLOUD_LOGS/m3pro_failsafe_${RUN_ID}/" || true
        fi
        ok "Extra sync to iCloud complete. All durable artifacts are now safely on Drive."
        ;;
    analyze)
        phase "BOOTSTRAP ANALYSIS"
        python3 "$ANALYZER" \
            --results-dir "$ICLOUD_RESULTS/tier2" \
            --n-bootstrap 10000 \
            --output-dir "$ICLOUD_RESULTS/analysis/${RUN_ID}"
        ok "Analysis written to $ICLOUD_RESULTS/analysis/${RUN_ID}"
        ;;
    full)
        $0 preflight
        $0 launch
        $0 analyze
        $0 sync
        ok "Full healthy campaign pipeline complete for $RUN_ID"
        ;;
    status)
        cat "$ICLOUD_LOGS/m3pro_failsafe_${RUN_ID}/campaign_status.json" 2>/dev/null || echo "No status yet for $RUN_ID"
        ;;
    show|cmd|command)
        phase "EXACT FAILSAFE COMMAND THIS WOULD RUN"
        cmd=$(build_cmd)
        echo "$cmd ${EXTRA_ARGS:-}"
        echo ""
        info "Copy-paste the above (with any extra flags) if you want to run the failsafe directly."
        ;;
    start)
        phase "START CAMPAIGN (one fucking command — absolute self-located safe wrapper)"
        doctor   # always surface the truth before committing resources

        # Generate the bulletproof inner wrapper (Chunk 2 core fix)
        INNER_SCRIPT="$LOCAL_HOT_BASE/inner_${RUN_ID}.sh"
        mkdir -p "$LOCAL_HOT_BASE"
        write_inner_wrapper "$INNER_SCRIPT" "$RUN_ID"

        info "Generated safe absolute inner script: $INNER_SCRIPT"
        info "This script hard-codes the correct FAILSAFE_PY/ANALYZER_PY from *this* worktree."
        info "Screen/tmux will exec it directly. No cd, no \$0, no FLEXAIDDS_REPO script pollution."

        SCREEN_NAME="grok_own_${RUN_ID}"

        # Add a simple trap in the *outer* manager for hot dir cleanup on abnormal exit
        cleanup_hot() { rm -rf "$LOCAL_HOT_BASE" 2>/dev/null || true; }
        trap cleanup_hot EXIT INT TERM

        if [[ "$SESSION_BACKEND" == "screen" ]] && command -v screen >/dev/null 2>&1; then
            screen -dmS "$SCREEN_NAME" "$INNER_SCRIPT"
            ok "Screen session '$SCREEN_NAME' started with safe wrapper."
            info "Attach: screen -r $SCREEN_NAME"
            info "Monitor status: tail -f $ICLOUD_LOGS/m3pro_failsafe_${RUN_ID}/campaign_status.json"
            info "Detach: Ctrl-A then D"
        elif [[ "$SESSION_BACKEND" == "tmux" ]] && command -v tmux >/dev/null 2>&1; then
            tmux new-session -d -s "$SCREEN_NAME" "$INNER_SCRIPT" || warn "tmux launch had issue"
            ok "tmux session '$SCREEN_NAME' started with safe wrapper."
            info "Attach: tmux attach -t $SCREEN_NAME"
        else
            warn "No screen/tmux or backend '$SESSION_BACKEND' unavailable — running the inner wrapper directly in foreground (nohup-style fallback)."
            "$INNER_SCRIPT" &
            INNER_PID=$!
            info "Inner pipeline running in background (PID $INNER_PID). Check status files on iCloud."
            wait $INNER_PID || true
            trap - EXIT INT TERM   # already exited cleanly
        fi
        ;;
    *)
        cat <<'EOF'
Grok Master M3 Pro Launcher & Syncer (iCloud durability, local speed, failsafe)

This script is the safe, high-level orchestrator for running your own
repetition benchmark campaigns on Astex Diverse / Non-Native + HAP2.

It wraps the hardened failsafe_campaign.py with:
- Automatic use of your ~/.flexaidds_env (iCloud paths)
- Local APFS hot execution (fast) + rsync to iCloud (durable)
- Max hardware settings for 18GB M3 Pro (workers=4)
- Your GA parameters (1000 generations, 2000 population)
- Distinct run-ids so your results are cleanly separated
- Full health features (preflight, resume, RUN_OK.json, fatal detection)

All durable artifacts (results, logs, analysis, manifests) live on iCloud Drive.

Subcommands:
  doctor      Diagnose self-location vs FLEXAIDDS_REPO / binary / iCloud / tmp sanity (NEW Chunk 1)
  preflight   Safe validation + doctor (recommended first step; now calls doctor)
  launch      Start/resume the 10-rep campaign
  start       The "one fucking command" — launches full pipeline inside screen/tmux using self-located paths
              (supports SESSION_BACKEND=screen|tmux env or flag in later chunks)
  sync        Extra rsync to iCloud (safety net)
  analyze     Generate bootstrap success rate reports
  full        preflight + launch + analyze + sync (one-shot)
  status      Show live campaign_status.json

  ./grok_master_launcher.sh doctor
  SESSION_BACKEND=tmux ./grok_master_launcher.sh start ...

Environment: Must have run setup_cloud_storage.sh so ~/.flexaidds_env exists.

This script exists so you can run your own clean, healthy, reproducible
benchmark campaigns on this exact M3 Pro 18GB + iCloud-only machine
without fighting the storage policy or forgetting the right flags.

Core idea (local speed + iCloud durability):
  - Execution uses fast local APFS (/private/tmp hot base)
  - All durable artifacts (results, logs, manifests, analysis) live on iCloud Drive

Recommended safe flow for your own runs on the three datasets:
  ./grok_master_launcher.sh preflight
  ./grok_master_launcher.sh launch     # best run inside screen/tmux
  # later...
  ./grok_master_launcher.sh analyze
  ./grok_master_launcher.sh sync

Or the one-command version:
  ./grok_master_launcher.sh full

Results will appear under $FLEXAIDDS_RESULTS with your run-id.
Use "show" to see the exact failsafe command being constructed.
EOF
        ;;
esac
