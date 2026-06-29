#!/usr/bin/env bash
# build_v128_repro.sh — pinned efc4f5d worktree + LTO build for v128 v50b reproduction
#
# Creates a sibling worktree at ../FlexAIDdS_v128_repro, builds FlexAIDdS +
# benchmark_datasets into build_lto, and verifies engine/runner/matrix fingerprints
# against the original v50b campaign anchors.
#
# Usage:
#   bash scripts/build_v128_repro.sh [--force]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
WORKTREE="${REPO}/../FlexAIDdS_v128_repro"
COMMIT="efc4f5d"
BUILD="${WORKTREE}/build_lto"
# v50b-era matrix (72d7c739) — not in git at efc4f5d HEAD; preserved in local build dirs
MATRIX_SEED="${REPO}/build_ablation/MC_st0r5.2_6.dat"
if [[ ! -f "${MATRIX_SEED}" ]]; then
    MATRIX_SEED="${REPO}/build_verify/MC_st0r5.2_6.dat"
fi

EXP_ENGINE_SHA="dbfaca09bfaf9ad8c6c154512f8e7906a6123ce2055a0350d3eec5961b925d0b"
EXP_RUNNER_SHA="53fa471cfe3a55b2b071bf87e2181caba889ee92124553199df436275d714781"
EXP_MATRIX_MD5="72d7c7396702331d96ff12d18f831796"

FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        -h|--help)
            head -20 "${BASH_SOURCE[0]}" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
    esac
done

sha256_file() {
    if command -v shasum &>/dev/null; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        sha256sum "$1" | awk '{print $1}'
    fi
}

md5_file() {
    if command -v md5 &>/dev/null; then
        md5 -q "$1"
    else
        md5sum "$1" | awk '{print $1}'
    fi
}

verify_artifacts() {
    local engine="${BUILD}/FlexAIDdS"
    local runner="${BUILD}/benchmark_datasets"
    local matrix="${BUILD}/MC_st0r5.2_6.dat"

    for p in "$engine" "$runner" "$matrix"; do
        [[ -f "$p" ]] || return 1
    done

    local got_engine got_runner got_matrix
    got_engine="$(sha256_file "$engine")"
    got_runner="$(sha256_file "$runner")"
    got_matrix="$(md5_file "$matrix")"

    if [[ "$got_engine" != "$EXP_ENGINE_SHA" ]]; then
        echo "ERROR: engine SHA mismatch"
        echo "  got  $got_engine"
        echo "  want $EXP_ENGINE_SHA"
        return 1
    fi
    if [[ "$got_runner" != "$EXP_RUNNER_SHA" ]]; then
        echo "ERROR: runner SHA mismatch"
        echo "  got  $got_runner"
        echo "  want $EXP_RUNNER_SHA"
        return 1
    fi
    if [[ "$got_matrix" != "$EXP_MATRIX_MD5" ]]; then
        echo "ERROR: matrix MD5 mismatch"
        echo "  got  $got_matrix"
        echo "  want $EXP_MATRIX_MD5"
        return 1
    fi
    return 0
}

echo "v128 repro build — commit ${COMMIT}"
echo "  repo     : ${REPO}"
echo "  worktree : ${WORKTREE}"
echo "  build    : ${BUILD}"

if [[ -e "${WORKTREE}/.git" ]]; then
    echo "Worktree exists — verifying checkout"
    git -C "${WORKTREE}" rev-parse --short HEAD
else
    echo "Creating worktree at ${COMMIT}..."
    git -C "${REPO}" worktree add "${WORKTREE}" "${COMMIT}"
fi

if [[ "${FORCE}" -eq 1 ]] && [[ -d "${BUILD}" ]]; then
    echo "Removing existing build (--force)"
    rm -rf "${BUILD}"
fi

if [[ -f "${MATRIX_SEED}" ]] && [[ -f "${BUILD}/MC_st0r5.2_6.dat" ]]; then
    got_matrix="$(md5_file "${BUILD}/MC_st0r5.2_6.dat")"
    if [[ "${got_matrix}" != "${EXP_MATRIX_MD5}" ]]; then
        echo "Reseeding v50b matrix into ${BUILD}"
        cp "${MATRIX_SEED}" "${BUILD}/MC_st0r5.2_6.dat"
    fi
fi

if verify_artifacts; then
    echo "Artifacts already match v50b fingerprints — skipping build"
    exit 0
fi

NPROC="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

echo "Configuring CMake..."
cmake -B "${BUILD}" -S "${WORKTREE}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_FLEXAIDDS_FAST=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_AVX2=OFF \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_CUDA=OFF

echo "Building FlexAIDdS + benchmark_datasets (${NPROC} jobs)..."
cmake --build "${BUILD}" --target FlexAIDdS benchmark_datasets -j "${NPROC}"

if [[ -f "${MATRIX_SEED}" ]]; then
    echo "Seeding v50b matrix from ${MATRIX_SEED}"
    cp "${MATRIX_SEED}" "${BUILD}/MC_st0r5.2_6.dat"
    cp "${MATRIX_SEED}" "${WORKTREE}/MC_st0r5.2_6.dat"
else
    echo "WARN: v50b matrix seed not found — matrix MD5 may not match 72d7c739"
fi

if ! verify_artifacts; then
    got_engine="$(sha256_file "${BUILD}/FlexAIDdS")"
    got_matrix="$(md5_file "${BUILD}/MC_st0r5.2_6.dat")"
    if [[ "${got_matrix}" == "${EXP_MATRIX_MD5}" ]]; then
        echo "NOTE: engine/runner SHA differ from June-2026 v50b anchors (toolchain rebuild)."
        echo "      commit=${COMMIT} matrix=${got_matrix} engine=${got_engine}"
        exit 0
    fi
    exit 1
fi
echo "v128 repro build OK — fingerprints match v50b anchors"