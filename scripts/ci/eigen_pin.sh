#!/usr/bin/env bash
# =============================================================================
# eigen_pin.sh — read the pinned Eigen version out of the single source of
# truth, and optionally fetch that exact release.
#
# Why this exists: every workflow used to install Eigen from a package
# manager (apt libeigen3-dev, brew eigen, choco eigen), which meant CI
# compiled whatever version the distro happened to ship.  Measured, that was
# an Eigen 3.x on Linux against a 5.x on macOS -- two different major
# versions of the self-adjoint eigensolver across one matrix.
#
# No workflow may carry an Eigen version literal.  They call this script, and
# this script reads cmake/FlexAIDEigen.cmake, so the version exists in
# exactly one place in the repository.
#
# Usage:
#   eigen_pin.sh version           -> e.g. 5.0.1
#   eigen_pin.sh sha256            -> pinned tarball SHA-256
#   eigen_pin.sh url               -> pinned tarball URL
#   eigen_pin.sh fetch [destdir]   -> download + verify + extract;
#                                     prints the include dir on stdout
#
# Written for bash 3.2 (macOS ships only 3.2.57): no ${var,,}, no mapfile,
# no wait -n, no associative arrays.
#
# Copyright 2026 Le Bonhomme Pharma. Licensed under Apache-2.0.
# =============================================================================
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
pin_file="$repo_root/cmake/FlexAIDEigen.cmake"

if [ ! -f "$pin_file" ]; then
    echo "eigen_pin.sh: cannot find $pin_file" >&2
    exit 1
fi

# Parse `set(FLEXAIDS_EIGEN_VERSION "x.y.z")` / `set(FLEXAIDS_EIGEN_SHA256 "...")`.
eigen_version="$(sed -n 's/^[[:space:]]*set([[:space:]]*FLEXAIDS_EIGEN_VERSION[[:space:]]*"\([^"]*\)".*/\1/p' "$pin_file" | head -1)"
eigen_sha256="$(sed -n 's/^[[:space:]]*set([[:space:]]*FLEXAIDS_EIGEN_SHA256[[:space:]]*"\([^"]*\)".*/\1/p' "$pin_file" | head -1)"

if [ -z "$eigen_version" ] || [ -z "$eigen_sha256" ]; then
    echo "eigen_pin.sh: failed to parse the Eigen pin out of $pin_file" >&2
    echo "eigen_pin.sh: version='$eigen_version' sha256='$eigen_sha256'" >&2
    exit 1
fi

eigen_url="https://gitlab.com/libeigen/eigen/-/archive/${eigen_version}/eigen-${eigen_version}.tar.gz"

# Portable SHA-256 of a file: coreutils sha256sum on Linux, shasum on macOS.
sha256_of() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

cmd="${1:-version}"

case "$cmd" in
    version) printf '%s\n' "$eigen_version" ;;
    sha256)  printf '%s\n' "$eigen_sha256" ;;
    url)     printf '%s\n' "$eigen_url" ;;
    fetch)
        dest="${2:-${RUNNER_TEMP:-/tmp}/flexaidds-eigen}"
        src="$dest/eigen-${eigen_version}"
        if [ -f "$src/Eigen/Dense" ]; then
            printf '%s\n' "$src"
            exit 0
        fi
        mkdir -p "$dest"
        tarball="$dest/eigen-${eigen_version}.tar.gz"
        # -fL: fail on HTTP error instead of saving an error page, follow
        # redirects. A saved error page would otherwise fail the hash check
        # with a confusing message.
        curl -fLsS -o "$tarball" "$eigen_url"
        got="$(sha256_of "$tarball")"
        if [ "$got" != "$eigen_sha256" ]; then
            echo "eigen_pin.sh: SHA-256 mismatch for $eigen_url" >&2
            echo "  expected: $eigen_sha256" >&2
            echo "  found   : $got" >&2
            exit 1
        fi
        tar xzf "$tarball" -C "$dest"
        if [ ! -f "$src/Eigen/Dense" ]; then
            echo "eigen_pin.sh: extracted tree has no Eigen/Dense at $src" >&2
            exit 1
        fi
        printf '%s\n' "$src"
        ;;
    *)
        echo "eigen_pin.sh: unknown command '$cmd'" >&2
        echo "usage: eigen_pin.sh {version|sha256|url|fetch [destdir]}" >&2
        exit 2
        ;;
esac