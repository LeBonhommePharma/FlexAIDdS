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

# 4. Unique run identity (for "your own" campaigns, cleanly separated from Codex runs)
RUN_ID="${RUN_ID:-grok_own_m3pro_$(date +%Y%m%d_%H%M%S)}"

# Local hot (fast APFS) workspace for execution speed — critical on iCloud-only storage
LOCAL_HOT_BASE="${LOCAL_HOT_BASE:-/private/tmp/grok_bench_hot_${RUN_ID}}"

# iCloud destination roots (all durable storage lives here via the env)
ICLOUD_RESULTS="$FLEXAIDDS_RESULTS"
ICLOUD_LOGS="$FLEXAIDDS_LOGS"

FAILSAVE_SCRIPT="$FLEXAIDDS_REPO/benchmarks/m3pro/failsafe_campaign.py"
ANALYZER="$FLEXAIDDS_REPO/benchmarks/m3pro/analyze_repetitions.py"

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
    preflight)
        phase "HEALTHY PREFLIGHT (M3 Pro 18GB + iCloud durability)"
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
        phase "START CAMPAIGN IN SCREEN (recommended for long runs)"
        SCREEN_NAME="grok_bench_${RUN_ID}"
        if command -v screen >/dev/null 2>&1; then
            info "Launching inside screen session: $SCREEN_NAME"
            screen -dmS "$SCREEN_NAME" bash -c "
                source ~/.flexaidds_env
                cd \"$FLEXAIDDS_REPO\"
                $0 preflight
                $0 launch
                $0 analyze
                $0 sync
                echo 'Campaign complete. Press any key to close this screen.'
                read -n 1
            "
            ok "Screen session '$SCREEN_NAME' started."
            info "Attach with: screen -r $SCREEN_NAME"
            info "Detach with: Ctrl-A then D"
        else
            warn "screen not found. Falling back to regular launch (run this in your own tmux/screen)."
            $0 full
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
  preflight   Safe validation (recommended first step)
  launch      Start/resume the 10-rep campaign
  sync        Extra rsync to iCloud (safety net)
  analyze     Generate bootstrap success rate reports
  full        preflight + launch + analyze + sync (one-shot)
  status      Show live campaign_status.json

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
