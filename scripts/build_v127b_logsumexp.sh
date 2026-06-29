#!/usr/bin/env bash
# build_v127b_logsumexp.sh — a4056163 worktree (logsumexp only, pre H-bond/VCT patch)
#
# Binary at a4056163 = Option B logsumexp-stable boltzmann_composite WITHOUT
# ba5364d3 H-bond cosine gate + C.ar stacking discriminator.
#
# Usage: bash scripts/build_v127b_logsumexp.sh [--force]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
WORKTREE="${REPO}/../FlexAIDdS_v127b_repro"
COMMIT="a4056163"
BUILD="${WORKTREE}/build_lto"

FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        -h|--help)
            head -15 "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
    esac
done

if [[ -e "${WORKTREE}/.git" ]]; then
    echo "Worktree exists — HEAD $(git -C "${WORKTREE}" rev-parse --short HEAD)"
else
    echo "Creating worktree at ${COMMIT}..."
    git -C "${REPO}" worktree add "${WORKTREE}" "${COMMIT}"
fi

if [[ "${FORCE}" -eq 1 ]] && [[ -d "${BUILD}" ]]; then
    rm -rf "${BUILD}"
fi

if [[ -x "${BUILD}/FlexAIDdS" ]] && [[ -x "${BUILD}/benchmark_datasets" ]]; then
    echo "Artifacts present — skipping rebuild"
    exit 0
fi

NPROC="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

cmake -B "${BUILD}" -S "${WORKTREE}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_FLEXAIDDS_FAST=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_AVX2=OFF \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_CUDA=OFF

cmake --build "${BUILD}" --target FlexAIDdS benchmark_datasets -j "${NPROC}"
echo "v127b build OK @ ${COMMIT}"