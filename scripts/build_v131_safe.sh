#!/usr/bin/env bash
# build_v131_safe.sh — v127-safe Lane A binary (pre-perf-regression + sulfo + holo data)
#
# Worktree @ 82ad51f4 (v127 launch commit, before 27e68e51 gaboom parallel / SoA churn),
# then cherry-pick:
#   04ff1735  sulfonamide –SO2NH– → N.am/S.O SybylTyper remap
#   bf8cf1d2  1G9V_apo + 1TW6_holo receptor data
#
# Build: Release LTO-style fast build, OpenMP ON, SoA distances OFF.
#
# Usage: bash scripts/build_v131_safe.sh [--force]
#
# Copyright 2026 Le Bonhomme Pharma. Apache-2.0.

set -euo pipefail

REPO="/Users/lp.more/Projects/FlexAIDdS"
WORKTREE="${FLEXAIDDS_V131_WORKTREE:-${REPO}/../FlexAIDdS_v131_safe}"
BASE_COMMIT="82ad51f4"
CHERRY_PICKS=(04ff1735 bf8cf1d2)
BUILD="${WORKTREE}/build_lto"
MANIFEST="${WORKTREE}/.v131_safe_manifest"

FORCE=0
for arg in "$@"; do
    case "$arg" in
        --force) FORCE=1 ;;
        -h|--help)
            sed -n '2,14p' "${BASH_SOURCE[0]}" | sed 's/^# \?//'
            exit 0
            ;;
    esac
done

expected_manifest() {
    local head
    head="$(git -C "${WORKTREE}" rev-parse HEAD)"
    cat > "${MANIFEST}" <<EOF
base=${BASE_COMMIT}
head=${head}
cherry_picks=${CHERRY_PICKS[*]}
soa_distances=OFF
built_at=$(date -u +%Y-%m-%dT%H:%M:%SZ)
EOF
}

apply_cherry_picks() {
    local cp rev
    for cp in "${CHERRY_PICKS[@]}"; do
        if git -C "${WORKTREE}" merge-base --is-ancestor "${cp}" HEAD; then
            echo "  already has ${cp}"
            continue
        fi
        echo "  cherry-pick ${cp}..."
        git -C "${WORKTREE}" cherry-pick "${cp}"
    done
}

if [[ -e "${WORKTREE}/.git" ]]; then
    echo "Worktree exists — HEAD $(git -C "${WORKTREE}" rev-parse --short HEAD)"
    apply_cherry_picks
else
    echo "Creating worktree @ ${BASE_COMMIT}..."
    git -C "${REPO}" worktree add "${WORKTREE}" "${BASE_COMMIT}"
    apply_cherry_picks
fi

expected_manifest
echo "Manifest: $(cat "${MANIFEST}" | tr '\n' ' ')"

if [[ "${FORCE}" -eq 1 ]] && [[ -d "${BUILD}" ]]; then
    echo "Removing ${BUILD} (--force)"
    rm -rf "${BUILD}"
fi

if [[ -x "${BUILD}/FlexAIDdS" ]] && [[ -x "${BUILD}/benchmark_datasets" ]] && [[ "${FORCE}" -eq 0 ]]; then
    echo "Artifacts present — skipping rebuild"
    echo "  FlexAIDdS: ${BUILD}/FlexAIDdS"
    exit 0
fi

NPROC="$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 4)"

cmake -B "${BUILD}" -S "${WORKTREE}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_FLEXAIDDS_FAST=ON \
    -DFLEXAIDS_USE_OPENMP=ON \
    -DFLEXAIDS_USE_AVX2=OFF \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_CUDA=OFF \
    -DFLEXAIDS_USE_SOA_DISTANCES=OFF

cmake --build "${BUILD}" --target FlexAIDdS benchmark_datasets -j "${NPROC}"

echo "v131_safe build OK"
echo "  worktree: ${WORKTREE}"
echo "  HEAD:     $(git -C "${WORKTREE}" rev-parse --short HEAD)"
echo "  binary:   ${BUILD}/FlexAIDdS"
if [[ -x "${BUILD}/FlexAIDdS" ]]; then
    shasum -a 256 "${BUILD}/FlexAIDdS" | awk '{print "  sha256:   "$1}'
fi