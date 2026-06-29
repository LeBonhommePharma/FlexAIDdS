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
#   bash scripts/build_vcontacts_bisect.sh [--repo PATH] [--force]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
# Script location wins over stale FLEXAIDDS_REPO from caller env.
REPO="${FLEXAIDDS_GIT_ROOT:-${DEFAULT_REPO}}"
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --force) FORCE=1; shift ;;
        --repo)
            shift
            [[ $# -gt 0 ]] || { echo "ERROR: --repo requires a path" >&2; exit 2; }
            REPO="$1"
            shift
            ;;
        --repo=*) REPO="${1#--repo=}"; shift ;;
        -h|--help)
            sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0
            ;;
        *)
            if [[ -d "$1" ]]; then
                REPO="$1"
                shift
            else
                echo "ERROR: unknown argument: $1" >&2
                exit 2
            fi
            ;;
    esac
done

NPROC="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

build_one() {
    local label="$1"
    local soa_flag="$2"
    local build_dir="${REPO}/build_bisect_${label}"
    local bin="${build_dir}/FlexAIDdS"
    local runner="${build_dir}/benchmark_datasets"

    if [[ "${FORCE}" -eq 1 ]] && [[ -d "${build_dir}" ]]; then
        echo "Removing ${build_dir} (--force)"
        rm -rf "${build_dir}"
    fi

    if [[ -x "${bin}" ]] && [[ -x "${runner}" ]] && [[ -f "${build_dir}/MC_st0r5.2_6.dat" ]] && [[ "${FORCE}" -eq 0 ]]; then
        local sha
        sha="$(shasum -a 256 "${bin}" | awk '{print $1}')"
        echo "  ${label}: artifacts present — skipping rebuild (${bin} sha256=${sha})"
        return 0
    fi

    echo "=== configure ${label} (SOA=${soa_flag}) → ${build_dir} ==="
    cmake -S "${REPO}" -B "${build_dir}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FLEXAIDDS_FAST=ON \
        -DFLEXAIDS_USE_SOA_DISTANCES="${soa_flag}" \
        -DFLEXAIDS_USE_OPENMP=ON \
        -DFLEXAIDS_USE_AVX2=OFF \
        -DFLEXAIDS_USE_METAL=OFF \
        -DFLEXAIDS_USE_CUDA=OFF \
        -DBUILD_TESTING=OFF

    echo "=== build ${label} ==="
    cmake --build "${build_dir}" -j"${NPROC}" \
        --target FlexAIDdS benchmark_datasets

    for data_file in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def; do
        if [[ ! -f "${build_dir}/${data_file}" ]] && [[ -f "${REPO}/${data_file}" ]]; then
            cp "${REPO}/${data_file}" "${build_dir}/${data_file}"
        fi
    done

    local sha
    sha="$(shasum -a 256 "${bin}" | awk '{print $1}')"
    echo "  ${label}: ${bin} sha256=${sha}"
}

echo "Vcontacts bisect repo: ${REPO}"
git -C "${REPO}" rev-parse --short HEAD
build_one "soa_off" "OFF"
build_one "soa_on" "ON"
echo "Vcontacts bisect builds ready."