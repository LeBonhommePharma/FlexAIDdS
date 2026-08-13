#!/usr/bin/env bash
# Profile-guided optimization helper. DEFAULT is unused: cmake FLEXAID_PGO=off.
# This script never flips a production/default build. Opt-in only:
#   bash scripts/run_flexaid_pgo.sh generate   # instrument + tell you to run a dock
#   bash scripts/run_flexaid_pgo.sh use        # rebuild consuming *.gcda / *.profraw
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="${1:-off}"
BUILD="${ROOT}/build_pgo"
case "$MODE" in
  off)
    echo "FLEXAID_PGO=off (default). Nothing to do."
    exit 0
    ;;
  generate)
    cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_FLEXAIDDS_FAST=OFF -DFLEXAID_PGO=generate
    cmake --build "$BUILD" --target FlexAIDdS -j4
    echo "Instrumented binary: $BUILD/FlexAIDdS"
    echo "Run a representative dock (e.g. 1G9V, FLEXAID_SEED=12345), then:"
    echo "  $0 use"
    ;;
  use)
    cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_FLEXAIDDS_FAST=OFF -DFLEXAID_PGO=use
    cmake --build "$BUILD" --target FlexAIDdS -j4
    echo "PGO-use binary: $BUILD/FlexAIDdS"
    ;;
  *)
    echo "usage: $0 off|generate|use" >&2
    exit 2
    ;;
esac
