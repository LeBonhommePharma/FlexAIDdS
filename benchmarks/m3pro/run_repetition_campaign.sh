#!/usr/bin/env bash
# run_repetition_campaign.sh — 30-run benchmark repetition campaign for FlexAIDdS
#
# Runs each dataset N times (default 30) with independent GA seeds,
# saving per-run JSON results for post-hoc bootstrap analysis.
#
# Usage:
#   ./benchmarks/m3pro/run_repetition_campaign.sh                    # all datasets, 30 runs
#   ./benchmarks/m3pro/run_repetition_campaign.sh --pilot            # 1 run of all datasets
#   ./benchmarks/m3pro/run_repetition_campaign.sh --dataset astex    # single dataset, 30 runs
#   ./benchmarks/m3pro/run_repetition_campaign.sh --runs 50          # custom run count
#   ./benchmarks/m3pro/run_repetition_campaign.sh --resume           # skip already-completed runs
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
phase() { printf "\n${BOLD}════════════════════════════════════════════════════════${NC}\n"; printf "${BOLD}  %s${NC}\n" "$*"; printf "${BOLD}════════════════════════════════════════════════════════${NC}\n\n"; }

ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    source "$ENV_FILE"
else
    die "Environment file not found: $ENV_FILE"
fi

REPO="${FLEXAIDDS_REPO:?not set}"
BINARY="${FLEXAIDDS_BUILD:?not set}/benchmark_datasets"
RESULTS="${FLEXAIDDS_RESULTS:?not set}"
LOGS="${FLEXAIDDS_LOGS:?not set}"
CACHE="${FLEXAIDDS_BENCHMARK_DATA:?not set}"

N_RUNS=30
PILOT=false
RESUME=false
SINGLE_DATASET=""
WORKERS=8
GA_GENERATIONS=2000
GA_POPULATION=1000
TEMPERATURE=300
CLUSTERING="FO"

DATASETS=(
    astex
    astex_nonnative
    hap2
    casf2016
    itc187
    dude37
    muv
    lsd_docking
    erds_specificity
    psychopharm23
    crossdock_modern
)

while [[ $# -gt 0 ]]; do
    case "$1" in
        --pilot)         PILOT=true; N_RUNS=1 ;;
        --resume)        RESUME=true ;;
        --runs)          shift; N_RUNS="$1" ;;
        --dataset)       shift; SINGLE_DATASET="$1" ;;
        --workers)       shift; WORKERS="$1" ;;
        --ga-generations) shift; GA_GENERATIONS="$1" ;;
        --ga-population) shift; GA_POPULATION="$1" ;;
        --temperature)   shift; TEMPERATURE="$1" ;;
        --clustering)   shift; CLUSTERING="$1" ;;
        --help|-h)
            echo "Usage: $0 [--pilot|--runs N|--dataset NAME|--resume|--workers N]"
            echo "       [--ga-generations N|--ga-population N|--temperature K]"
            exit 0
            ;;
        *) die "Unknown argument: $1" ;;
    esac
    shift
done

if [[ -n "$SINGLE_DATASET" ]]; then
    DATASETS=("$SINGLE_DATASET")
fi

if [[ ! -f "$BINARY" ]]; then
    die "benchmark_datasets binary not found: $BINARY"
fi

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
MASTER_LOG="$LOGS/campaign_${TIMESTAMP}.log"
mkdir -p "$LOGS" "$RESULTS/tier2" "$RESULTS/analysis"

{
    echo "FlexAIDdS Repetition Campaign — $TIMESTAMP"
    echo "=============================================="
    echo ""
    echo "  Runs per dataset:   $N_RUNS"
    echo "  Datasets:           ${DATASETS[*]}"
    echo "  Workers:            $WORKERS"
    echo "  GA:                 pop=$GA_POPULATION gen=$GA_GENERATIONS"
    echo "  Temperature:        $TEMPERATURE K"
    echo "  Clustering:         $CLUSTERING"
    echo "  Pilot mode:         $PILOT"
    echo "  Resume:             $RESUME"
    echo "  Results:            $RESULTS"
    echo "  Cache:              $CACHE"
    echo ""
    echo "  Hardware:"
    echo "    CPU:   $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Apple Silicon')"
    echo "    Cores: $(sysctl -n hw.ncpu)"
    echo "    RAM:   $(( $(sysctl -n hw.memsize) / 1073741824 )) GB"
    echo ""
    echo "  Git: $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo 'unknown') on $(cd "$REPO" && git branch --show-current 2>/dev/null || echo 'unknown')"
    echo ""
} | tee "$MASTER_LOG"

GLOBAL_START=$(date +%s)

for ds in "${DATASETS[@]}"; do
    phase "Dataset: $ds ($N_RUNS runs)"

    DS_START=$(date +%s)
    DS_DIR="$RESULTS/tier2/$ds"
    mkdir -p "$DS_DIR"

    for run_num in $(seq 1 "$N_RUNS"); do
        RUN_TAG=$(printf "%s_run%02d_%s" "$ds" "$run_num" "$TIMESTAMP")
        RUN_DIR="$DS_DIR/run$(printf '%02d' "$run_num")"
        mkdir -p "$RUN_DIR"

        RESULT_JSON="$RUN_DIR/${RUN_TAG}.json"
        RESULT_CSV="$RUN_DIR/${RUN_TAG}_results.csv"
        RESULT_MD="$RUN_DIR/${RUN_TAG}_report.md"

        if [[ "$RESUME" == true && -f "$RESULT_CSV" && -s "$RESULT_CSV" ]]; then
            ok "[$ds] Run $run_num/$N_RUNS already exists — skipping"
            continue
        fi

        RUN_START=$(date +%s)
        info "[$ds] Run $run_num/$N_RUNS starting..."

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
            RUN_DURATION=$((RUN_END - RUN_START))
            ok "[$ds] Run $run_num/$N_RUNS complete in ${RUN_DURATION}s"

            echo "$ds run$(printf '%02d' "$run_num") ${RUN_DURATION}s" >> "$MASTER_LOG"
        else
            RUN_END=$(date +%s)
            RUN_DURATION=$((RUN_END - RUN_START))
            warn "[$ds] Run $run_num/$N_RUNS FAILED after ${RUN_DURATION}s"
            echo "$ds run$(printf '%02d' "$run_num") FAILED ${RUN_DURATION}s" >> "$MASTER_LOG"
        fi
    done

    DS_END=$(date +%s)
    DS_DURATION=$((DS_END - DS_START))
    ok "[$ds] All $N_RUNS runs complete in ${DS_DURATION}s ($((DS_DURATION / 60))m $((DS_DURATION % 60))s)"
    echo "" >> "$MASTER_LOG"
    echo "$ds: total=${DS_DURATION}s ($N_RUNS runs)" >> "$MASTER_LOG"
done

phase "Campaign Complete"

GLOBAL_END=$(date +%s)
GLOBAL_DURATION=$((GLOBAL_END - GLOBAL_START))

{
    echo ""
    echo "=============================================="
    echo "  CAMPAIGN COMPLETE"
    echo "=============================================="
    echo ""
    echo "  Total wall-clock: ${GLOBAL_DURATION}s ($(( GLOBAL_DURATION / 3600 ))h $(( (GLOBAL_DURATION % 3600) / 60 ))m)"
    echo "  Datasets:         ${#DATASETS[@]}"
    echo "  Runs per dataset: $N_RUNS"
    echo "  Total runs:       $(( ${#DATASETS[@]} * N_RUNS ))"
    echo "  Results:          $RESULTS/tier2/"
    echo "  Log:              $MASTER_LOG"
    echo ""
    echo "  Next step: run bootstrap analysis"
    echo "    python benchmarks/m3pro/analyze_repetitions.py --results-dir $RESULTS/tier2 --n-bootstrap 10000"
    echo ""
} | tee -a "$MASTER_LOG"
