#!/usr/bin/env bash
# build_vcontacts_commit_bisect.sh — single-commit revert builds on HEAD (post-27e68e51 revert)
#
# Restores pre-commit file versions (git checkout SHA^ -- <files>) to avoid revert conflicts.
#
# Variants:
#   revert_f9c80fe5  — SoA double sqrdist parity
#   revert_d2295cf0  — thread-local atoms for parallel coord cache
#   revert_d4d68592  — PR4 scalar-identical Vcontacts loop
#
# Usage:
#   bash scripts/build_vcontacts_commit_bisect.sh revert_f9c80fe5
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${FLEXAIDDS_GIT_ROOT:-$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)}"
NPROC="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

VARIANT="${1:-}"
SHA=""
FILES=()
EXTRA_FILES=()
EXTRA_SHA=""
case "${VARIANT}" in
    revert_f9c80fe5)
        SHA="f9c80fe5"
        FILES=(CMakeLists.txt LIB/CMakeLists.txt LIB/Vcontacts.cpp)
        ;;
    revert_d2295cf0)
        SHA="d2295cf0"
        # calc_rmsd_chrom signature must match cluster/gaboom.h (66cecbb9 added out_n_atoms).
        FILES=(LIB/DensityPeak_Cluster.cpp LIB/FastOPTICS_cluster.cpp LIB/calc_rmsd_chrom.cpp)
        EXTRA_FILES=(LIB/gaboom.h LIB/cluster.cpp)
        EXTRA_SHA="66cecbb9^"
        ;;
    revert_d4d68592)
        SHA="d4d68592"
        FILES=(LIB/Vcontacts.cpp LIB/gaboom.cpp scripts/run_soa_parity_gate.py)
        ;;
    *)
        echo "Usage: $0 {revert_f9c80fe5|revert_d2295cf0|revert_d4d68592}" >&2
        exit 2
        ;;
esac

WT="${REPO}/../FlexAIDdS_${VARIANT}"
BUILD="${WT}/build_bisect"
BIN="${BUILD}/FlexAIDdS"

if [[ ! -e "${WT}/.git" ]]; then
    echo "Creating worktree ${WT} @ master..."
    git -C "${REPO}" worktree add -B "bisect_${VARIANT}" "${WT}" master
else
    git -C "${WT}" fetch --quiet origin master 2>/dev/null || true
    git -C "${WT}" checkout -f bisect_"${VARIANT}" 2>/dev/null || git -C "${WT}" checkout -f master
    git -C "${WT}" reset --hard origin/master 2>/dev/null || git -C "${WT}" reset --hard master
fi

echo "Restoring pre-${SHA} versions of: ${FILES[*]}"
git -C "${WT}" checkout "${SHA}^" -- "${FILES[@]}"
if [[ -n "${EXTRA_SHA}" && ${#EXTRA_FILES[@]} -gt 0 ]]; then
    echo "Restoring link-compatible ${EXTRA_SHA} versions of: ${EXTRA_FILES[*]}"
    git -C "${WT}" checkout "${EXTRA_SHA}" -- "${EXTRA_FILES[@]}"
fi

echo "Worktree HEAD+patch: $(git -C "${WT}" rev-parse --short HEAD) (file-level undo ${SHA})"

if [[ -d "${BUILD}" ]]; then
    echo "Removing stale ${BUILD} for clean link"
    rm -rf "${BUILD}"
fi

# Avoid getcwd failures if the caller's cwd was inside a removed build tree.
cd "${REPO}"

cmake -S "${WT}" -B "${BUILD}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_FLEXAIDDS_FAST=ON \
    -DFLEXAIDS_USE_SOA_DISTANCES=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_AVX2=OFF \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_CUDA=OFF \
    -DBUILD_TESTING=OFF

# Runner links from main-repo build_lto (shared DatasetRunner); only engine varies.
cmake --build "${BUILD}" -j"${NPROC}" --target FlexAIDdS
MAIN_RUNNER="${REPO}/build_lto/benchmark_datasets"
if [[ ! -x "${MAIN_RUNNER}" ]]; then
    echo "Building shared runner in ${REPO}/build_lto..."
    cmake -S "${REPO}" -B "${REPO}/build_lto" \
        -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_FLEXAIDDS_FAST=ON \
        -DBUILD_TESTING=OFF \
        -DFLEXAIDS_USE_METAL=OFF
    cmake --build "${REPO}/build_lto" -j"${NPROC}" --target benchmark_datasets
fi

for data_file in MC_st0r5.2_6.dat AMINO.def NUCLEOTIDES.def; do
    if [[ ! -f "${BUILD}/${data_file}" ]] && [[ -f "${WT}/${data_file}" ]]; then
        cp "${WT}/${data_file}" "${BUILD}/${data_file}"
    fi
done

echo "${VARIANT} ready: ${BIN}"
shasum -a 256 "${BIN}" | awk '{print "sha256="$1}'