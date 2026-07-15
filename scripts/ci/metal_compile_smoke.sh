#!/usr/bin/env bash
# Metal shader compile smoke — no full docking, no OBJCXX link.
# Exit codes:
#   0  shaders compiled (or metallib linked when available)
#   2  Metal toolchain not available (caller may treat as skip)
#   1  compile failure
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${1:-${ROOT}/build-metal-smoke}"
mkdir -p "${OUT_DIR}"

require_cmd() {
  command -v "$1" >/dev/null 2>&1
}

METALC=""
if require_cmd xcrun; then
  if METALC="$(xcrun -f metal 2>/dev/null)" && [[ -n "${METALC}" && -x "${METALC}" ]]; then
    :
  else
    METALC=""
  fi
fi

if [[ -z "${METALC}" ]]; then
  echo "SKIP: Metal shader compiler not found (xcrun -f metal)."
  echo "Full Metal link/runtime gate requires a self-hosted runner labeled self-hosted-m3."
  echo "See docs/CI_RELEASE_GATES.md and .github/workflows/metal-self-hosted.yml"
  exit 2
fi

METALLIB=""
METAL_DIR="$(dirname "${METALC}")"
if [[ -x "${METAL_DIR}/metallib" ]]; then
  METALLIB="${METAL_DIR}/metallib"
elif require_cmd metallib; then
  METALLIB="$(command -v metallib)"
fi

echo "Metal compiler: ${METALC}"
echo "metallib: ${METALLIB:-<not found — .air only>}"
echo "Output: ${OUT_DIR}"

SHADERS=(
  "LIB/CavityDetect/CavityDetect.metal"
  "LIB/ShannonThermoStack/shannon_metal.metal"
  "LIB/TurboQuant.metal"
  "LIB/MetalRMSD.metal"
  "LIB/gpu_fast_optics_metal.metal"
)

AIR_FILES=()
failed=0
for rel in "${SHADERS[@]}"; do
  src="${ROOT}/${rel}"
  if [[ ! -f "${src}" ]]; then
    echo "WARN: missing shader ${rel} — skipping"
    continue
  fi
  base="$(basename "${rel}" .metal)"
  air="${OUT_DIR}/${base}.air"
  echo "=== metal -c ${rel} ==="
  if ! "${METALC}" -c "${src}" -o "${air}"; then
    echo "FAIL: compile ${rel}"
    failed=1
    continue
  fi
  AIR_FILES+=("${air}")
  echo "OK: ${air}"
done

if [[ "${failed}" -ne 0 ]]; then
  echo "Metal compile smoke FAILED"
  exit 1
fi

if [[ ${#AIR_FILES[@]} -eq 0 ]]; then
  echo "FAIL: no Metal shaders found to compile"
  exit 1
fi

if [[ -n "${METALLIB}" ]]; then
  lib="${OUT_DIR}/flexaidds_smoke.metallib"
  echo "=== metallib -> ${lib} ==="
  "${METALLIB}" -o "${lib}" "${AIR_FILES[@]}"
  echo "OK: ${lib}"
else
  echo "NOTE: metallib not available; .air compile-only smoke is sufficient for hosted CI."
fi

echo "Metal compile smoke: SUCCESS (${#AIR_FILES[@]} shaders)"
exit 0
