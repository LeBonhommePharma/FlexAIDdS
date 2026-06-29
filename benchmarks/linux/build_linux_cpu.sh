#!/usr/bin/env bash
# build_linux_cpu.sh — Release AVX2+OpenMP build for perf-swarm microbenchmarks.
#
# Usage:
#   ./benchmarks/linux/build_linux_cpu.sh
#   ./benchmarks/linux/build_linux_cpu.sh --cuda
#   ./benchmarks/linux/build_linux_cpu.sh --build-dir /path/to/build
#
# Apache-2.0 (c) 2026 NRGlab, Universite de Montreal

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BUILD_DIR="${REPO}/build"
USE_CUDA=OFF

while [[ $# -gt 0 ]]; do
    case "$1" in
        --cuda) USE_CUDA=ON; shift ;;
        --build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [--cuda] [--build-dir DIR]"
            exit 0
            ;;
        *) echo "Unknown arg: $1" >&2; exit 1 ;;
    esac
done

if [[ "$(uname -s)" != "Linux" ]]; then
    echo "ERROR: Linux x86-64 required (got $(uname -s))" >&2
    exit 1
fi

if [[ "${USE_CUDA}" == "ON" ]]; then
    if ! command -v nvcc >/dev/null 2>&1; then
        echo "ERROR: --cuda requires nvcc in PATH" >&2
        exit 1
    fi
fi

: "${CC:=gcc-14}"
: "${CXX:=g++-14}"

mkdir -p "${BUILD_DIR}"

cmake -S "${REPO}" -B "${BUILD_DIR}" -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_COMPILER="${CC}" \
    -DCMAKE_CXX_COMPILER="${CXX}" \
    -DFLEXAIDS_USE_CUDA="${USE_CUDA}" \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_AVX2=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DBUILD_TESTING=ON \
    -DENABLE_TENCOM_BENCHMARK=ON \
    -DENABLE_VCFBATCH_BENCHMARK=ON

cmake --build "${BUILD_DIR}" -j "$(nproc)" \
    --target benchmark_tencom benchmark_vcfbatch

echo "Built:"
ls -lh "${BUILD_DIR}/benchmark_tencom" "${BUILD_DIR}/benchmark_vcfbatch"