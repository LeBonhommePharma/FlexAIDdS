#!/usr/bin/env bash
# build_vcontacts_bisect.sh — dual HEAD builds for Vcontacts/SoA bisect (smoke-12)
#
# Produces:
#   build_bisect_soa_off/FlexAIDdS  — HEAD sources, FLEXAIDS_USE_SOA_DISTANCES=OFF
#   build_bisect_soa_on/FlexAIDdS   — HEAD sources, FLEXAIDS_USE_SOA_DISTANCES=ON (default)
#
# Compare against v131_safe worktree (82ad51f4 + sulfo + holo, pre-27e68e51).
#
# Usage:
#   bash scripts/build_vcontacts_bisect.sh [--repo PATH]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${1:-${FLEXAIDDS_REPO:-$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)}}"

build_one() {
    local label="$1"
    local soa_flag="$2"
    local build_dir="${REPO}/build_bisect_${label}"

    echo "=== configure ${label} (SOA=${soa_flag}) → ${build_dir} ==="
    cmake -S "${REPO}" -B "${build_dir}" -G Ninja \
        -DCMAKE_BUILD_TYPE=Release \
        -DFLEXAIDS_USE_SOA_DISTANCES="${soa_flag}" \
        -DFLEXAIDS_USE_OPENMP=ON \
        -DBUILD_TESTING=OFF

    echo "=== build ${label} ==="
    cmake --build "${build_dir}" -j"$(sysctl -n hw.ncpu 2>/dev/null || echo 8)" \
        --target FlexAIDdS benchmark_datasets

    local bin="${build_dir}/FlexAIDdS"
    local sha
    sha="$(shasum -a 256 "${bin}" | awk '{print $1}')"
    echo "  ${label}: ${bin} sha256=${sha}"
}

git -C "${REPO}" rev-parse --short HEAD
build_one "soa_off" "OFF"
build_one "soa_on" "ON"
echo "Vcontacts bisect builds ready."