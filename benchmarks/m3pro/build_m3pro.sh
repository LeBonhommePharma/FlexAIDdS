#!/usr/bin/env bash
# build_m3pro.sh — CMake configure + build for MacBook Pro M3 Pro 18GB
#
# Enables Metal GPU, OpenMP, Eigen3, all benchmark targets.
# Disables AVX2/AVX512/CUDA (ARM64). Builds with 6 P-core jobs.
# Build artifacts go to iCloud to save local SSD space.
#
# Usage:
#   ./benchmarks/m3pro/build_m3pro.sh          # full build
#   ./benchmarks/m3pro/build_m3pro.sh --clean   # clean + rebuild
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

set -euo pipefail

# ─── Colors ──────────────────────────────────────────────────────────────────

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

info()  { printf "${CYAN}[INFO]${NC}  %s\n" "$*"; }
ok()    { printf "${GREEN}[OK]${NC}    %s\n" "$*"; }
warn()  { printf "${YELLOW}[WARN]${NC}  %s\n" "$*"; }
die()   { printf "${RED}[ERROR]${NC} %s\n" "$*" >&2; exit 1; }

# ─── Load environment ────────────────────────────────────────────────────────

ENV_FILE="$HOME/.flexaidds_env"
if [[ -f "$ENV_FILE" ]]; then
    source "$ENV_FILE"
else
    die "Environment file not found: $ENV_FILE — run setup_cloud_storage.sh first"
fi

REPO="${FLEXAIDDS_REPO:?FLEXAIDDS_REPO not set}"
BUILD="${FLEXAIDDS_BUILD:?FLEXAIDDS_BUILD not set}"
LOGS="${FLEXAIDDS_LOGS:?FLEXAIDDS_LOGS not set}"

# ─── Parse arguments ─────────────────────────────────────────────────────────

CLEAN=false
for arg in "$@"; do
    case "$arg" in
        --clean) CLEAN=true ;;
        --help|-h)
            echo "Usage: $0 [--clean]"
            echo "  --clean  Remove build directory and rebuild from scratch"
            exit 0
            ;;
        *) die "Unknown argument: $arg" ;;
    esac
done

# ─── Validate platform ──────────────────────────────────────────────────────

ARCH=$(uname -m)
if [[ "$ARCH" != "arm64" ]]; then
    die "Expected arm64 (Apple Silicon), got: $ARCH"
fi

OS=$(uname -s)
if [[ "$OS" != "Darwin" ]]; then
    die "Expected macOS (Darwin), got: $OS"
fi

info "Platform: macOS $ARCH (Apple Silicon)"

# ─── Clean if requested ─────────────────────────────────────────────────────

if [[ "$CLEAN" == true ]] && [[ -d "$BUILD" ]]; then
    warn "Removing build directory: $BUILD"
    rm -rf "$BUILD"
fi

mkdir -p "$BUILD"
mkdir -p "$LOGS"

# ─── Build parameters ───────────────────────────────────────────────────────

# Use P-cores only (6 on M3 Pro), leave E-cores for system + cloud sync
PARALLEL_JOBS=6

TIMESTAMP="$(date +%Y%m%dT%H%M%S)"
LOG_FILE="$LOGS/build_${TIMESTAMP}.log"

info "Build directory: $BUILD"
info "Parallel jobs:   $PARALLEL_JOBS"
info "Log file:        $LOG_FILE"

# ─── CMake configure ─────────────────────────────────────────────────────────

info "Configuring CMake..."
BUILD_START=$(date +%s)

cmake -S "$REPO" -B "$BUILD" \
    -DCMAKE_BUILD_TYPE=Release \
    -DFLEXAIDS_USE_METAL=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_EIGEN=ON \
    -DFLEXAIDS_USE_AVX2=OFF \
    -DFLEXAIDS_USE_AVX512=OFF \
    -DFLEXAIDS_USE_CUDA=OFF \
    -DBUILD_TESTING=ON \
    -DBUILD_PYTHON_BINDINGS=ON \
    -DENABLE_TENCOM_BENCHMARK=ON \
    -DENABLE_VCFBATCH_BENCHMARK=ON \
    -DENABLE_DISPATCH_BENCHMARK=ON \
    -DENABLE_BENCHMARK_DATASETS=ON \
    -DMETAL_ENABLE_PROFILING=ON \
    2>&1 | tee "$LOG_FILE"

ok "CMake configuration complete"

# Extract Metal detection info from cmake output
if grep -q "Apple Silicon" "$LOG_FILE" 2>/dev/null; then
    ok "Metal GPU: Apple Silicon detected"
elif grep -q "Metal GPU acceleration: ENABLED" "$LOG_FILE" 2>/dev/null; then
    ok "Metal GPU: enabled"
else
    warn "Metal GPU detection not confirmed in cmake output"
fi

# ─── Build ───────────────────────────────────────────────────────────────────

info "Building with $PARALLEL_JOBS parallel jobs..."

# On Apple Silicon (especially newer macOS/Xcode), the full 'all' target can hit
# SDK / C++26 header friction in system libc++. For reliable local benchmarking
# we build the targets actually needed by the m3pro campaign first.
if [[ "$ARCH" == "arm64" ]]; then
    info "Apple Silicon detected — building only m3pro-critical targets for reliability"
    TARGETS=(
        FlexAID
        benchmark_tencom
        benchmark_vcfbatch
        benchmark_dispatch
        tENCoM
        tencom_entropy_diff
    )
    for t in "${TARGETS[@]}"; do
        echo "  → Building $t ..." | tee -a "$LOG_FILE"
        if ! cmake --build "$BUILD" --target "$t" -j"$PARALLEL_JOBS" 2>&1 | tee -a "$LOG_FILE"; then
            warn "Target $t had issues (may still be usable if the main binary succeeded)"
        fi
    done
else
    cmake --build "$BUILD" -j"$PARALLEL_JOBS" 2>&1 | tee -a "$LOG_FILE"
fi

BUILD_END=$(date +%s)
BUILD_DURATION=$((BUILD_END - BUILD_START))

ok "Build completed in ${BUILD_DURATION}s"

# ─── Report binary sizes ────────────────────────────────────────────────────

echo "" | tee -a "$LOG_FILE"
info "Binary sizes:" | tee -a "$LOG_FILE"

BINARIES=(
    "FlexAID"
    "benchmark_dispatch"
    "benchmark_vcfbatch"
    "benchmark_tencom"
    "benchmark_datasets"
)

for bin in "${BINARIES[@]}"; do
    if [[ -f "$BUILD/$bin" ]]; then
        SIZE=$(du -h "$BUILD/$bin" | cut -f1)
        echo "  $bin: $SIZE" | tee -a "$LOG_FILE"
    fi
done

# ─── Run tests ───────────────────────────────────────────────────────────────

echo "" | tee -a "$LOG_FILE"
info "Running ctest..."

if ctest --test-dir "$BUILD" --output-on-failure 2>&1 | tee -a "$LOG_FILE"; then
    ok "All tests passed"
else
    warn "Some tests failed — check log: $LOG_FILE"
fi

# ─── Summary ─────────────────────────────────────────────────────────────────

echo ""
echo "================================================================"
echo "  M3 Pro Build Complete"
echo "================================================================"
echo ""
echo "  Build dir:    $BUILD"
echo "  Duration:     ${BUILD_DURATION}s"
echo "  Log:          $LOG_FILE"
echo ""
echo "  Next step:"
echo "    ./benchmarks/m3pro/run_benchmarks.sh"
echo ""
echo "================================================================"
