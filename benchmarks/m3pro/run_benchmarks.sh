#!/usr/bin/env bash
# run_benchmarks.sh — Master benchmark runner for MacBook Pro M3 Pro 18GB
#
# Runs all benchmark tiers sequentially with memory-aware worker counts:
#   Phase 1: C++ kernel benchmarks (dispatch, vcfbatch, tencom)
#   Phase 2: Tier-1 dataset benchmark (CASF-2016, 5 targets, ~5 min)
#   Phase 3: Tier-2 dataset benchmarks (all 10 datasets, sequential, hours)
#   Phase 4: Final report consolidation (results already on iCloud)
#
# Usage:
#   ./benchmarks/m3pro/run_benchmarks.sh                # all phases
#   ./benchmarks/m3pro/run_benchmarks.sh --kernels-only  # phase 1 only
#   ./benchmarks/m3pro/run_benchmarks.sh --tier1-only    # phase 2 only
#   ./benchmarks/m3pro/run_benchmarks.sh --tier2-only    # phase 3 only
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

set -euo pipefail

# ─── Colors ──────────────────────────────────────────────────────────────────

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

info()  { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
phase() { printf "\n${BOLD}════════════════════════════════════════════════════════${NC}\n"; printf "${BOLD}  %s${NC}\n" "$*"; printf "${BOLD}════════════════════════════════════════════════════════${NC}\n\n"; }
die()   { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; exit 1; }

# ─── Load environment ────────────────────────────────────────────────────────

ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    source "$ENV_FILE"
else
    die "Environment file not found: $ENV_FILE — run setup_cloud_storage.sh first"
fi

REPO="${FLEXAIDDS_REPO:?not set}"
BUILD="${FLEXAIDDS_BUILD:?not set}"
RESULTS="${FLEXAIDDS_RESULTS:?not set}"
LOGS="${FLEXAIDDS_LOGS:?not set}"
DATA="${FLEXAIDDS_BENCHMARK_DATA:?not set}"
BINARY="${FLEXAIDDS_BINARY:?not set}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ─── Parse arguments ─────────────────────────────────────────────────────────

RUN_KERNELS=true
RUN_TIER1=true
RUN_TIER2=true

for arg in "$@"; do
    case "$arg" in
        --kernels-only) RUN_TIER1=false; RUN_TIER2=false ;;
        --tier1-only)   RUN_KERNELS=false; RUN_TIER2=false ;;
        --tier2-only)   RUN_KERNELS=false; RUN_TIER1=false ;;
        --help|-h)
            echo "Usage: $0 [--kernels-only|--tier1-only|--tier2-only]"
            exit 0
            ;;
        *) die "Unknown argument: $arg" ;;
    esac
done

# ─── Validate build ─────────────────────────────────────────────────────────

if [[ ! -f "$BINARY" ]]; then
    if [[ "$RUN_KERNELS" == true && "$RUN_TIER1" == false && "$RUN_TIER2" == false ]]; then
        warn "FlexAID binary not found — kernels-only will proceed with whatever benchmark_* binaries exist"
    else
        die "FlexAID binary not found: $BINARY — run build_m3pro.sh first (needed for tier-1/2)"
    fi
fi

# ─── Hardware detection ──────────────────────────────────────────────────────

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
MASTER_LOG="$LOGS/benchmark_run_${TIMESTAMP}.log"
mkdir -p "$LOGS"

{
    echo "FlexAIDdS Benchmark Run — $TIMESTAMP"
    echo "========================================"
    echo ""
    echo "Hardware Profile:"
    echo "  CPU:     $(sysctl -n machdep.cpu.brand_string 2>/dev/null || echo 'Apple Silicon')"
    echo "  Cores:   $(sysctl -n hw.ncpu) total ($(sysctl -n hw.perflevel0.logicalcpu 2>/dev/null || echo '?') P + $(sysctl -n hw.perflevel1.logicalcpu 2>/dev/null || echo '?') E)"
    echo "  RAM:     $(( $(sysctl -n hw.memsize) / 1073741824 )) GB"
    echo "  GPU:     $(system_profiler SPDisplaysDataType 2>/dev/null | grep 'Chipset Model' | sed 's/.*: //' || echo 'unknown')"
    echo "  OS:      $(sw_vers -productName 2>/dev/null || echo macOS) $(sw_vers -productVersion 2>/dev/null || echo '')"
    echo "  Storage: iCloud 2TB only (all results written directly here)"
    echo ""
    echo "Git SHA:   $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
    echo "Branch:    $(cd "$REPO" && git branch --show-current 2>/dev/null || echo 'unknown')"
    echo ""
} | tee "$MASTER_LOG"

GLOBAL_START=$(date +%s)

# ─── Phase 1: C++ Kernel Benchmarks ─────────────────────────────────────────

if [[ "$RUN_KERNELS" == true ]]; then
    phase "Phase 1: C++ Kernel Benchmarks"
    mkdir -p "$RESULTS/kernels"

    KERNEL_REPORT="$RESULTS/kernels/report_${TIMESTAMP}.txt"
    : > "$KERNEL_REPORT"

    # --- benchmark_dispatch ---
    if [[ -f "$BUILD/benchmark_dispatch" ]]; then
        info "Running benchmark_dispatch (size=100000, reps=100)..."
        DISPATCH_OUT="$RESULTS/kernels/dispatch_${TIMESTAMP}.txt"
        "$BUILD/benchmark_dispatch" --size 100000 --reps 100 2>&1 | tee "$DISPATCH_OUT" | tee -a "$KERNEL_REPORT"
        ok "Dispatch benchmark complete"
    else
        warn "benchmark_dispatch not found — skipping"
    fi

    echo "" >> "$KERNEL_REPORT"

    # --- benchmark_vcfbatch ---
    if [[ -f "$BUILD/benchmark_vcfbatch" ]]; then
        info "Running benchmark_vcfbatch (pop=200, genes=20)..."
        VCFBATCH_OUT="$RESULTS/kernels/vcfbatch_${TIMESTAMP}.txt"
        "$BUILD/benchmark_vcfbatch" 200 20 2>&1 | tee "$VCFBATCH_OUT" | tee -a "$KERNEL_REPORT"
        ok "VoronoiCFBatch benchmark complete"
    else
        warn "benchmark_vcfbatch not found — skipping"
    fi

    echo "" >> "$KERNEL_REPORT"

    # --- benchmark_tencom ---
    if [[ -f "$BUILD/benchmark_tencom" ]]; then
        info "Running benchmark_tencom..."
        TENCOM_OUT="$RESULTS/kernels/tencom_${TIMESTAMP}.txt"
        "$BUILD/benchmark_tencom" 2>&1 | tee "$TENCOM_OUT" | tee -a "$KERNEL_REPORT"
        ok "tENCoM benchmark complete"
    else
        warn "benchmark_tencom not found — skipping"
    fi

    ok "Phase 1 complete — kernel results in $RESULTS/kernels/"
    echo "" >> "$MASTER_LOG"
    echo "Phase 1 (Kernels): COMPLETE" >> "$MASTER_LOG"
    cat "$KERNEL_REPORT" >> "$MASTER_LOG"

    # (results already written to iCloud)
fi

# ─── Phase 2: Tier-1 Dataset Benchmark ──────────────────────────────────────

if [[ "$RUN_TIER1" == true ]]; then
    phase "Phase 2: Tier-1 Dataset Benchmark (CASF-2016, 5 targets)"
    mkdir -p "$RESULTS/tier1"

    # 4 workers: ~4.5 GB per worker from 18GB total (3GB OS + 4GB Metal + 2.5GB×4 workers ≈ 17GB)
    info "Running CASF-2016 tier-1 with 4 workers..."
    python -m benchmarks.run \
        --dataset casf2016 \
        --tier 1 \
        --workers 4 \
        --results-dir "$RESULTS/tier1" \
        --data-dir "$DATA" \
        --binary "$BINARY" \
        --report-prefix "$RESULTS/tier1/report_${TIMESTAMP}" \
        2>&1 | tee -a "$MASTER_LOG"

    TIER1_EXIT=$?
    if [[ $TIER1_EXIT -eq 0 ]]; then
        ok "Phase 2 complete — no regressions"
    elif [[ $TIER1_EXIT -eq 1 ]]; then
        warn "Phase 2 complete — REGRESSIONS DETECTED"
    else
        warn "Phase 2 failed with exit code $TIER1_EXIT"
    fi

    echo "Phase 2 (Tier-1): exit=$TIER1_EXIT" >> "$MASTER_LOG"
    # (results already written to iCloud)
fi

# ─── Phase 3: Tier-2 Dataset Benchmarks ─────────────────────────────────────

if [[ "$RUN_TIER2" == true ]]; then
    phase "Phase 3: Tier-2 Dataset Benchmarks (10 datasets, sequential)"
    mkdir -p "$RESULTS/tier2"

    # Sequential execution to stay within 18GB RAM.
    # 2 workers per dataset: ~4.5GB per worker (3GB OS + 4GB Metal + 4.5GB×2 ≈ 16GB)
    DATASETS=(astex_diverse astex_nonnative hap2 casf2016 itc187 dude37 muv lsd_docking erds_specificity psychopharm23)
    TIER2_REGRESSIONS=0

    for ds in "${DATASETS[@]}"; do
        info "Running $ds (tier-2, 2 workers, bootstrap=5000)..."
        DS_START=$(date +%s)

        python -m benchmarks.run \
            --dataset "$ds" \
            --tier 2 \
            --workers 2 \
            --bootstrap \
            --n-bootstrap 5000 \
            --results-dir "$RESULTS/tier2" \
            --data-dir "$DATA" \
            --binary "$BINARY" \
            --report-prefix "$RESULTS/tier2/${ds}_${TIMESTAMP}" \
            2>&1 | tee -a "$MASTER_LOG"

        DS_EXIT=$?
        DS_END=$(date +%s)
        DS_DURATION=$((DS_END - DS_START))

        if [[ $DS_EXIT -eq 0 ]]; then
            ok "$ds complete in ${DS_DURATION}s — no regressions"
        elif [[ $DS_EXIT -eq 1 ]]; then
            warn "$ds complete in ${DS_DURATION}s — REGRESSION DETECTED"
            TIER2_REGRESSIONS=$((TIER2_REGRESSIONS + 1))
        else
            warn "$ds failed (exit=$DS_EXIT) after ${DS_DURATION}s"
        fi

        echo "  $ds: exit=$DS_EXIT, duration=${DS_DURATION}s" >> "$MASTER_LOG"

        # Mirror after each dataset for incremental backup
        # (results already written to iCloud)
    done

    echo "" >> "$MASTER_LOG"
    echo "Phase 3 (Tier-2): ${#DATASETS[@]} datasets, $TIER2_REGRESSIONS regressions" >> "$MASTER_LOG"
fi

# ─── Phase 4: Final Report + Blocking Mirror ────────────────────────────────

phase "Phase 4: Final Report + Sync"

GLOBAL_END=$(date +%s)
GLOBAL_DURATION=$((GLOBAL_END - GLOBAL_START))

{
    echo ""
    echo "========================================"
    echo "  BENCHMARK RUN COMPLETE"
    echo "========================================"
    echo ""
    echo "  Total wall-clock: ${GLOBAL_DURATION}s ($(( GLOBAL_DURATION / 60 ))m $(( GLOBAL_DURATION % 60 ))s)"
    echo "  Results:          $RESULTS"
    echo "  Log:              $MASTER_LOG"
    echo ""
} | tee -a "$MASTER_LOG"

echo ""
echo "================================================================"
echo "  All done. Results are on iCloud 2TB only:"
echo "    iCloud:       $RESULTS"
echo "    Master log:   $MASTER_LOG"
echo "================================================================"
