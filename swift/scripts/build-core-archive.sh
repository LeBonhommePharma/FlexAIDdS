#!/usr/bin/env bash
#
# build-core-archive.sh — build the REAL FlexAID C++ engine for the Swift package.
#
# Why this exists
# ---------------
# Sources/FlexAIDCore/*.mm is a thin Objective-C++ bridge over the actual
# engine (statmech::StatMechEngine, BindingMode/BindingPopulation, ENCoM,
# tENCoM, the Shannon stack, GA, read_input, ic2cf). Those implementations live
# in <repo>/LIB and are driven by the root CMakeLists (149 translation units).
#
# SwiftPM refuses to compile sources outside the package root, so the engine
# cannot be added to a SwiftPM target directly. This script builds the genuine
# CMake `flexaid_core` objects and archives them so `swift test` links the real
# implementation.
#
# What this script must never do: emit stubs, synthesize symbols, or hand
# `-undefined dynamic_lookup` to the linker. A green test bought that way is a
# fabricated result and worse than a reported link failure.
#
# Usage:
#   swift/scripts/build-core-archive.sh          # writes swift/.build/cxxcore/lib
#   FLEXAIDDS_CORE_JOBS=4 swift/scripts/build-core-archive.sh
#
# The build directory is deliberately under swift/.build so it can never touch
# <repo>/build, which is reserved for running benchmark campaigns.
#
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
package_root="$(cd "${script_dir}/.." && pwd)"
repo_root="$(cd "${package_root}/.." && pwd)"

build_dir="${FLEXAIDDS_CORE_BUILD_DIR:-${package_root}/.build/cxxcore}"
# NOTE: not "${build_dir}/lib" — macOS filesystems are case-insensitive and
# that path collides with CMake's own LIB/ subdirectory.
lib_dir="${build_dir}/swiftlink"
jobs="${FLEXAIDDS_CORE_JOBS:-2}"

case "${build_dir}" in
    "${repo_root}/build"|"${repo_root}/build/"*)
        echo "refusing to build into ${repo_root}/build (campaign output)" >&2
        exit 1
        ;;
esac

generator=(-G Ninja)
command -v ninja >/dev/null 2>&1 || generator=()

cmake -S "${repo_root}" -B "${build_dir}" "${generator[@]}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TESTING=OFF \
    -DBUILD_PYTHON_BINDINGS=OFF \
    -DFLEXAIDS_USE_METAL=OFF \
    -DFLEXAIDS_USE_CUDA=OFF \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON

cmake --build "${build_dir}" --target flexaid_core -j "${jobs}"

object_dir="${build_dir}/LIB/CMakeFiles/flexaid_core.dir"
if [ ! -d "${object_dir}" ]; then
    echo "expected CMake object directory not found: ${object_dir}" >&2
    exit 1
fi

mkdir -p "${lib_dir}"
archive="${lib_dir}/libflexaid_core.a"
rm -f "${archive}"

# flexaid_core is an OBJECT library, so CMake produces loose .o files rather
# than an archive. Collect exactly those objects — nothing is generated here.
object_list="$(mktemp)"
trap 'rm -f "${object_list}"' EXIT
find "${object_dir}" -name '*.o' -print0 | sort -z > "${object_list}"
if [ ! -s "${object_list}" ]; then
    echo "no objects found under ${object_dir}" >&2
    exit 1
fi

if [ "$(uname -s)" = "Darwin" ]; then
    # Apple's libtool tolerates duplicate object basenames across LIB/ subdirs.
    xargs -0 libtool -static -o "${archive}" < "${object_list}"
else
    xargs -0 ar rcs "${archive}" < "${object_list}"
    ranlib "${archive}"
fi

# Runtime dependencies of the archive itself (OpenMP, ...). Package.swift reads
# this file instead of hardcoding a machine-specific toolchain path.
link_flags="${lib_dir}/flexaid_core.link"
: > "${link_flags}"
omp_lib="$(sed -n 's/^OpenMP_[a-zA-Z0-9_]*_LIBRARY:FILEPATH=//p' "${build_dir}/CMakeCache.txt" | head -1)"
if [ -n "${omp_lib}" ] && [ -f "${omp_lib}" ]; then
    printf -- '-L%s\n-lomp\n' "$(dirname "${omp_lib}")" >> "${link_flags}"
fi

echo "archive:    ${archive}"
echo "link flags: $(tr '\n' ' ' < "${link_flags}")"
echo
echo "swift test will now link the real engine. Re-run this script after"
echo "changing anything under LIB/."
