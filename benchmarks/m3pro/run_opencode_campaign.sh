#!/usr/bin/env bash
# run_opencode_campaign.sh — opencode (GLM-5.1) benchmark campaign runner
#
# Runs benchmark datasets with run31+ to avoid collision with Claude Code
# PID 96971 which uses run01-30. Uses local SSD for hot I/O, syncs to
# iCloud on completion. Writes fleet status JSON for Bonhomme Dashboard.
#
# Campaign:
#   Phase 1: Astex Diverse — 3 runs (run31-33)
#   Phase 2: Astex Non-Native — 1 run (run31)
#   Phase 3: HAP2 — 3 runs (run31-33)
#
# Usage:
#   ./benchmarks/m3pro/run_opencode_campaign.sh
#   ./benchmarks/m3pro/run_opencode_campaign.sh --resume
#   ./benchmarks/m3pro/run_opencode_campaign.sh --dataset astex
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

info()  { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
die()   { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; exit 1; }
phase() { printf "\n${BOLD}══════════════════════════════════════════════════════════${NC}\n"; printf "${BOLD}  %s${NC}\n" "$*"; printf "${BOLD}══════════════════════════════════════════════════════════${NC}\n\n"; }

ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    source "$ENV_FILE"
fi

REPO="/Users/lp.more/Projects/FlexAIDdS"
BINARY="$REPO/build_opencode/benchmark_datasets"

FLEXAIDDS_FAST=1
export FLEXAIDDS_FAST
export FLEXAIDDS_BUILD="$REPO/build_opencode"
export FLEXAIDDS_REPO="$REPO"
FAST_BASE="${FLEXAIDDS_FAST_BASE:-$HOME/.flexaidds_fast}"
RESULTS="$FAST_BASE/results"
LOGS="${FLEXAIDDS_LOGS:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS/logs}"
CACHE="${FLEXAIDDS_BENCHMARK_DATA:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS/benchmark_data}"
# Prefer FLEXAIDDS_ICLOUD everywhere; avoid baking full iCloud container paths.
ICLOUD_BASE="${FLEXAIDDS_ICLOUD:-${HOME}/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS}"
ICLOUD_RESULTS="${ICLOUD_BASE}/results"
FLEET_STATUS_DIR="${ICLOUD_BASE}"

export FLEXAIDDS_OMP_THREADS="${FLEXAIDDS_OMP_THREADS:-2}"

RESUME=false
SINGLE_DATASET=""
WORKERS=4
GA_GENERATIONS=1000
GA_POPULATION=2000
TEMPERATURE=300
CLUSTERING="FO"

DS_NAMES=(astex astex_nonnative hap2)
DS_RUNS=(3 1 3)
RUN_START=31

while [[ $# -gt 0 ]]; do
    case "$1" in
        --resume)        RESUME=true ;;
        --dataset)       shift; SINGLE_DATASET="$1" ;;
        --workers)       shift; WORKERS="$1" ;;
        --ga-generations) shift; GA_GENERATIONS="$1" ;;
        --ga-population) shift; GA_POPULATION="$1" ;;
        --temperature)   shift; TEMPERATURE="$1" ;;
        --clustering)    shift; CLUSTERING="$1" ;;
        --help|-h)
            echo "Usage: $0 [--resume|--dataset NAME|--workers N]"
            echo "       [--ga-generations N|--ga-population N|--temperature K]"
            exit 0
            ;;
        *) die "Unknown argument: $1" ;;
    esac
    shift
done

if [[ -n "$SINGLE_DATASET" ]]; then
    found=0
    for i in "${!DS_NAMES[@]}"; do
        if [[ "${DS_NAMES[$i]}" == "$SINGLE_DATASET" ]]; then found=1; break; fi
    done
    if [[ "$found" -eq 0 ]]; then
        die "Unknown dataset: $SINGLE_DATASET (valid: ${DS_NAMES[*]})"
    fi
fi

if [[ ! -f "$BINARY" ]]; then
    die "benchmark_datasets binary not found: $BINARY"
fi

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
MASTER_LOG="$LOGS/opencode_campaign_${TIMESTAMP}.log"
mkdir -p "$LOGS" "$RESULTS/tier2" "$RESULTS/analysis"

{
    echo "FlexAIDdS opencode Campaign — $TIMESTAMP"
    echo "=============================================="
    echo ""
    echo "  Runner:            opencode (GLM-5.1)"
    echo "  Run numbers:       start=$RUN_START"
    echo "  Workers:           $WORKERS"
    echo "  GA:                pop=$GA_POPULATION gen=$GA_GENERATIONS"
    echo "  Temperature:       $TEMPERATURE K"
    echo "  Clustering:        $CLUSTERING"
    echo "  Resume:            $RESUME"
    echo "  Binary:            $BINARY"
    echo "  Results (hot):     $RESULTS"
    echo "  Results (iCloud):  $ICLOUD_RESULTS"
    echo "  Cache:             $CACHE"
    echo ""
    echo "  Campaign:"
    for i in "${!DS_NAMES[@]}"; do
        n="${DS_RUNS[$i]}"
        echo "    ${DS_NAMES[$i]}: $n runs (run$(printf '%02d' $RUN_START)..run$(printf '%02d' $((RUN_START + n - 1))))"
    done
    echo ""
    echo "  Hardware:"
    echo "    CPU:   $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Apple Silicon')"
    echo "    Cores: $(sysctl -n hw.ncpu)"
    echo "    RAM:   $(( $(sysctl -n hw.memsize) / 1073741824 )) GB"
    echo ""
    echo "  Git: $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo 'unknown') on $(cd "$REPO" && git branch --show-current 2>/dev/null || echo 'unknown')"
    echo ""
} | tee "$MASTER_LOG"

write_fleet_status() {
    local active_ds="$1"
    local completed_runs="$2"
    local total_runs="$3"
    local status_file="$FLEET_STATUS_DIR/fleet_status_opencode.json"

    local astex_done=0 astex_status="pending"
    local astex_nn_done=0 astex_nn_status="pending"
    local hap2_done=0 hap2_status="pending"

    for i in "${!DS_NAMES[@]}"; do
        local n="${DS_NAMES[$i]}"
        case "$n" in
            astex)           astex_done="${DS_DONE[$i]:-0}"; astex_status="${DS_STAT[$i]:-pending}" ;;
            astex_nonnative) astex_nn_done="${DS_DONE[$i]:-0}"; astex_nn_status="${DS_STAT[$i]:-pending}" ;;
            hap2)            hap2_done="${DS_DONE[$i]:-0}"; hap2_status="${DS_STAT[$i]:-pending}" ;;
        esac
    done

    cat > "$status_file" <<EOF
{
  "runner": "opencode-glm5.1",
  "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "activeDataset": "$active_ds",
  "campaign": {
    "astex": { "completed": $astex_done, "total": ${DS_RUNS[0]:-0}, "status": "$astex_status" },
    "astex_nonnative": { "completed": $astex_nn_done, "total": ${DS_RUNS[1]:-0}, "status": "$astex_nn_status" },
    "hap2": { "completed": $hap2_done, "total": ${DS_RUNS[2]:-0}, "status": "$hap2_status" }
  },
  "metrics": {
    "jobID": "opencode_${TIMESTAMP}",
    "totalChunks": $total_runs,
    "completedChunks": $completed_runs,
    "failedChunks": ${CAMP_FAILED:-0},
    "orphanedChunks": 0,
    "activeDevices": 1,
    "totalTFLOPS": 5.1,
    "estimatedRemainingSeconds": null,
    "timestamp": "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  },
  "devices": [
    {
      "deviceID": "m3pro-opencode",
      "model": "MacBook Pro M3 Pro",
      "estimatedTFLOPS": 5.1,
      "availableMemoryGB": 9,
      "computeWeight": 0.5,
      "thermalState": "nominal"
    }
  ]
}
EOF
    info "Fleet status updated"
}

sync_to_icloud() {
    local ds="$1"
    local src="$RESULTS/tier2/$ds"
    local dst="$ICLOUD_RESULTS/tier2/$ds"

    if [[ ! -d "$src" ]]; then
        warn "No results to sync for $ds"
        return
    fi

    info "Syncing $ds results to iCloud..."
    for run_dir in "$src"/run*; do
        if [[ -d "$run_dir" ]]; then
            run_name=$(basename "$run_dir")
            mkdir -p "$dst/$run_name"
            rsync -a --ignore-existing "$run_dir/" "$dst/$run_name/"
        fi
    done
    ok "Synced $ds to iCloud: $dst"
}

GLOBAL_START=$(date +%s)

DS_DONE=(0 0 0)
DS_STAT=("pending" "pending" "pending")
TOTAL_RUNS=0
TOTAL_DONE=0
CAMP_FAILED=0

for n in "${DS_RUNS[@]}"; do
    TOTAL_RUNS=$(( TOTAL_RUNS + n ))
done

for idx in "${!DS_NAMES[@]}"; do
    ds="${DS_NAMES[$idx]}"
    n_runs="${DS_RUNS[$idx]}"

    if [[ -n "$SINGLE_DATASET" && "$ds" != "$SINGLE_DATASET" ]]; then
        continue
    fi

    phase "Dataset: $ds ($n_runs runs, run$(printf '%02d' $RUN_START)..run$(printf '%02d' $((RUN_START + n_runs - 1))))"

    DS_STAT[$idx]="running"
    write_fleet_status "$ds" "$TOTAL_DONE" "$TOTAL_RUNS"

    DS_START=$(date +%s)
    DS_DIR="$RESULTS/tier2/$ds"
    mkdir -p "$DS_DIR"

    for (( i=0; i<n_runs; i++ )); do
        run_num=$(( RUN_START + i ))
        RUN_TAG=$(printf "%s_run%02d_%s" "$ds" "$run_num" "$TIMESTAMP")
        RUN_DIR="$DS_DIR/run$(printf '%02d' "$run_num")"
        mkdir -p "$RUN_DIR"

        RESULT_CSV="$RUN_DIR/${RUN_TAG}_results.csv"

        if [[ "$RESUME" == true && -f "$RESULT_CSV" && -s "$RESULT_CSV" ]]; then
            ok "[$ds] Run $run_num already exists — skipping"
            DS_DONE[$idx]=$(( ${DS_DONE[$idx]} + 1 ))
            TOTAL_DONE=$(( TOTAL_DONE + 1 ))
            write_fleet_status "$ds" "$TOTAL_DONE" "$TOTAL_RUNS"
            continue
        fi

        RUN_START_TS=$(date +%s)
        info "[$ds] Run $run_num starting..."

        if "$BINARY" \
            --benchmark "$ds" \
            --threads "$WORKERS" \
            --cache "$CACHE" \
            --output "$RUN_DIR" \
            --ga-generations "$GA_GENERATIONS" \
            --ga-population "$GA_POPULATION" \
            --temperature "$TEMPERATURE" \
            --clustering "$CLUSTERING" \
            2>&1 | tee -a "$MASTER_LOG"; then
            RUN_END=$(date +%s)
            RUN_DURATION=$((RUN_END - RUN_START_TS))
            ok "[$ds] Run $run_num complete in ${RUN_DURATION}s"
            DS_DONE[$idx]=$(( ${DS_DONE[$idx]} + 1 ))
            TOTAL_DONE=$(( TOTAL_DONE + 1 ))
            echo "$ds run$(printf '%02d' "$run_num") ${RUN_DURATION}s" >> "$MASTER_LOG"
        else
            RUN_END=$(date +%s)
            RUN_DURATION=$((RUN_END - RUN_START_TS))
            warn "[$ds] Run $run_num FAILED after ${RUN_DURATION}s"
            CAMP_FAILED=$(( CAMP_FAILED + 1 ))
            TOTAL_DONE=$(( TOTAL_DONE + 1 ))
            echo "$ds run$(printf '%02d' "$run_num") FAILED ${RUN_DURATION}s" >> "$MASTER_LOG"
        fi

        write_fleet_status "$ds" "$TOTAL_DONE" "$TOTAL_RUNS"
    done

    DS_END=$(date +%s)
    DS_DURATION=$((DS_END - DS_START))
    DS_STAT[$idx]="complete"
    ok "[$ds] All $n_runs runs complete in ${DS_DURATION}s ($((DS_DURATION / 60))m $((DS_DURATION % 60))s)"
    echo "" >> "$MASTER_LOG"
    echo "$ds: total=${DS_DURATION}s ($n_runs runs)" >> "$MASTER_LOG"

    sync_to_icloud "$ds"
    write_fleet_status "idle" "$TOTAL_DONE" "$TOTAL_RUNS"
done

phase "Campaign Complete"

GLOBAL_END=$(date +%s)
GLOBAL_DURATION=$((GLOBAL_END - GLOBAL_START))

{
    echo ""
    echo "=============================================="
    echo "  OPENCODE CAMPAIGN COMPLETE"
    echo "=============================================="
    echo ""
    echo "  Runner:           opencode (GLM-5.1)"
    echo "  Total wall-clock: ${GLOBAL_DURATION}s ($(( GLOBAL_DURATION / 3600 ))h $(( (GLOBAL_DURATION % 3600) / 60 ))m)"
    echo "  Datasets:         ${#DS_NAMES[@]}"
    echo "  Total runs:       $TOTAL_RUNS"
    echo "  Failed:           $CAMP_FAILED"
    echo "  Results (hot):    $RESULTS/tier2/"
    echo "  Results (iCloud): $ICLOUD_RESULTS/tier2/"
    echo "  Log:              $MASTER_LOG"
    echo "  Fleet status:     $FLEET_STATUS_DIR/fleet_status_opencode.json"
    echo ""
    if [[ $CAMP_FAILED -eq 0 ]]; then
        echo "  All runs successful. Ready for bootstrap analysis."
    else
        echo "  WARNING: $CAMP_FAILED runs failed. Check logs."
    fi
    echo ""
} | tee -a "$MASTER_LOG"

write_fleet_status "complete" "$TOTAL_DONE" "$TOTAL_RUNS"
