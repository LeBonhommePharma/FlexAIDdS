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
#   ./benchmarks/m3pro/run_repetition_campaign.sh --check-quality-every 5   # Periodic quality gate (recommended)
#   ./benchmarks/m3pro/run_repetition_campaign.sh --check-quality-every 5 --no-early-exit
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

set -euo pipefail

# Self-location (so we can reliably find analyze_repetitions.py even when called from elsewhere)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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
    set -a
    source "$ENV_FILE"
    set +a
else
    die "Environment file not found: $ENV_FILE"
fi

REPO="${FLEXAIDDS_REPO:?not set}"
BUILD="${FLEXAIDDS_BUILD:?not set}"
BINARY="$BUILD/benchmark_datasets"
DOCKING_BINARY="${FLEXAIDDS_BINARY:-$BUILD/FlexAID}"
RESULTS="${FLEXAIDDS_RESULTS:?not set}"
LOGS="${FLEXAIDDS_LOGS:?not set}"
CACHE="${FLEXAIDDS_BENCHMARK_DATA:?not set}"

export FLEXAIDDS_REPO="$REPO"
export FLEXAIDDS_BUILD="$BUILD"
export FLEXAIDDS_BINARY="$DOCKING_BINARY"
export FLEXAIDDS_RESULTS="$RESULTS"
export FLEXAIDDS_LOGS="$LOGS"
export FLEXAIDDS_BENCHMARK_DATA="$CACHE"

N_RUNS=30
PILOT=false
RESUME=false
SINGLE_DATASET=""
# Workers: 1 = optimal on M3 Pro (6 OMP threads undivided across 6 P-cores).
# Do NOT raise without also halving --omp-threads (see DatasetRunner OMP fix dc5f5ab).
WORKERS=1
# GA spec: benchmark plan §6.1 = 500 gen × 1000 chrom = 510k evals/complex (~7 min/complex).
# Shannon HSC early-stop (H < 1.3863 nats) terminates converged runs before gen 500.
# Override at CLI with --ga-generations 2000 for exploratory hard-landscape runs only.
GA_GENERATIONS=500
GA_POPULATION=1000
TEMPERATURE=300
CLUSTERING="FO"

# Periodic Quality Gate (P1)
CHECK_QUALITY_EVERY=0          # 0 = disabled. Recommended: 5 or 10
EARLY_EXIT_ON_FAIL=true

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
        --check-quality-every)
            shift; CHECK_QUALITY_EVERY="$1" ;;
        --no-early-exit)
            EARLY_EXIT_ON_FAIL=false ;;
        --help|-h)
            cat <<EOF
Usage: $0 [options]

Core options:
  --pilot
  --runs N
  --dataset NAME
  --resume
  --workers N
  --ga-generations N
  --ga-population N
  --temperature K
  --clustering NAME

Quality Gate (P1 - Periodic Early Exit):
  --check-quality-every N     Run analyzer + quality gate every N runs per dataset.
                              If status becomes FAIL, abort remaining runs (saves compute).
                              Recommended values: 5 or 10.
  --no-early-exit             Continue even if the quality gate returns FAIL.

Example (recommended for long campaigns):
  $0 --check-quality-every 5
EOF
            exit 0
            ;;
        *) die "Unknown argument: $1" ;;
    esac
    shift
done

if [[ -n "$SINGLE_DATASET" ]]; then
    DATASETS=("$SINGLE_DATASET")
fi

if [[ ! -x "$BINARY" ]]; then
    die "benchmark_datasets binary not executable: $BINARY"
fi
if [[ ! -x "$DOCKING_BINARY" ]]; then
    die "FlexAID docking binary not executable: $DOCKING_BINARY"
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
    echo "  Runner binary:      $BINARY"
    echo "  Docking binary:     $DOCKING_BINARY"
    if (( CHECK_QUALITY_EVERY > 0 )); then
        echo "  Quality Gate:       Every $CHECK_QUALITY_EVERY runs (early exit on FAIL: $EARLY_EXIT_ON_FAIL)"
    else
        echo "  Quality Gate:       Disabled (use --check-quality-every N to enable)"
    fi
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

        # Periodic quality-gate check (intended implementation):
        # After every N runs (e.g. every 5), we could run a quick partial analysis
        # and call check_campaign_quality_gate (or a lighter version).
        # For now the full check is at the end; we can make it periodic with a flag.
        if (( run_num % 5 == 0 )); then
            # Placeholder for future periodic gate:
            #   ./analyze_repetitions.py --results-dir ... --output-dir ...
            #   check_campaign_quality_gate
            true
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

            # === Periodic Quality Gate Check ===
            if (( CHECK_QUALITY_EVERY > 0 && run_num % CHECK_QUALITY_EVERY == 0 )); then
                info "[$ds] Run $run_num — running periodic quality gate check (every $CHECK_QUALITY_EVERY runs)..."

                # Run analyzer on **only this dataset** with lighter bootstrap for speed
                python3 "$SCRIPT_DIR/analyze_repetitions.py" \
                    --results-dir "$RESULTS/tier2" \
                    --dataset "$ds" \
                    --n-bootstrap 2000 \
                    --output-dir "$RESULTS/analysis" \
                    2>&1 | tee -a "$MASTER_LOG" || true

                check_campaign_quality_gate "$RESULTS/analysis/campaign_status.json"
            fi
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

# === Periodic Quality-Gate (fully implemented) ===
check_campaign_quality_gate() {
    local status_file="${1:-$RESULTS/analysis/campaign_status.json}"

    if [[ ! -f "$status_file" ]]; then
        return 0
    fi

    local status
    status=$(python3 -c '
import json, sys
try:
    data = json.load(open(sys.argv[1]))
    print(data.get("overall_status", "UNKNOWN"))
    print(data.get("actionable", ""))
    rates = {}
    for name, value in (data.get("dataset_success_rates") or {}).items():
        mean = value.get("mean") if isinstance(value, dict) else value
        if mean is not None:
            rates[name] = mean
    if not rates:
        for ds in data.get("details", []):
            sr = ds.get("metrics", {}).get("success_rate", {})
            mean = sr.get("mean")
            if mean is not None:
                rates[ds.get("dataset", "unknown")] = mean
    if rates:
        print("PER_DATASET")
        for name in sorted(rates):
            print(f"{name}={rates[name]:.3f}")
except Exception as e:
    print("ERROR")
    print(str(e))
' "$status_file" 2>/dev/null || echo -e "ERROR\n")

    local overall_status
    local actionable
    local dataset_lines

    overall_status=$(echo "$status" | head -1)
    actionable=$(echo "$status" | sed -n '2p')
    dataset_lines=$(echo "$status" | sed -n '/^PER_DATASET$/,$p' | tail -n +2)

    if [[ "$overall_status" == "FAIL" ]]; then
        warn "════════════════════════════════════════════════════════"
        warn " QUALITY GATE FAILED"
        warn "════════════════════════════════════════════════════════"
        warn "Status file : $status_file"
        warn "Actionable  : $actionable"
        if [[ -n "$dataset_lines" ]]; then
            warn "Per-dataset success rates:"
            while IFS= read -r line; do
                [[ -z "$line" ]] && continue
                warn "  $line"
            done <<< "$dataset_lines"
        fi
        warn ""
        if [[ "$EARLY_EXIT_ON_FAIL" == "true" ]]; then
            warn "Early exit triggered. Remaining runs aborted to save compute."
            exit 2
        else
            warn "Continuing because --no-early-exit was used."
        fi
    elif [[ "$overall_status" == "WARN" ]]; then
        warn "Quality gate: WARN — $actionable"
        if [[ -n "$dataset_lines" ]]; then
            warn "Per-dataset success rates:"
            while IFS= read -r line; do
                [[ -z "$line" ]] && continue
                warn "  $line"
            done <<< "$dataset_lines"
        fi
    fi
}

# Optional: call after the full campaign (or integrate into loop for per-dataset checks)
check_campaign_quality_gate
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
    if (( CHECK_QUALITY_EVERY > 0 )); then
        echo "  (Periodic quality gate was enabled every $CHECK_QUALITY_EVERY runs)"
    fi
    echo ""
} | tee -a "$MASTER_LOG"
