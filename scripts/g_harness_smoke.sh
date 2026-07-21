#!/usr/bin/env bash
# VALIDATION G local smoke — host-testable half + optional container build.
# SPDX-License-Identifier: Apache-2.0
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="${G_HARNESS_OUT:-$ROOT/validation_evidence/build_ab/G_harness}"
mkdir -p "$OUT"
LOG="$OUT/g_smoke.log"
: >"$LOG"
exec > >(tee -a "$LOG") 2>&1

echo "=== G harness smoke $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "host=$(hostname) uname=$(uname -srm)"
echo "apptainer=$(command -v apptainer || echo MISSING)"
echo "singularity=$(command -v singularity || echo MISSING)"
echo "docker=$(command -v docker || echo MISSING)"
echo "sbatch=$(command -v sbatch || echo MISSING)"

# 1) Manifest
python3 scripts/validate_cmaes_manifest.py scripts/cmaes_ab_manifest.json \
  | tee "$OUT/G_manifest_validate.txt"
echo "PASS:manifest"

# 2) Recipe presence
test -f containers/flexaidds_locked_x86_64.def
test -f containers/Dockerfile.locked
test -f containers/README_HARNESS.md
test -f scripts/narval_cmaes_array.sh
test -f analysis/collapse_fingerprint.py
echo "PASS:recipes"

# 3) Fingerprint INVARIANT (self-compare on real 2e6 CMA-ES trace if present)
TRACE="$ROOT/validation_evidence/build_ab/E_ab_2e6/cmaes/cmaes_2e6_cmaes_entropy.csv"
if [[ ! -f "$TRACE" ]]; then
  TRACE="$ROOT/validation_evidence/build_ab/C_cmaes_smoke/cmaes_smoke_cmaes_entropy.csv"
fi
if [[ ! -f "$TRACE" ]]; then
  # synthetic fallback from mock testdata
  TRACE="$ROOT/analysis/testdata/mock_trace.csv"
fi
test -f "$TRACE"
python3 analysis/collapse_fingerprint.py "$TRACE" --out "$OUT/fp_a.json"
python3 analysis/collapse_fingerprint.py "$TRACE" --out "$OUT/fp_b.json"
CMP_OUT=$(python3 analysis/collapse_fingerprint.py --compare "$OUT/fp_a.json" "$OUT/fp_b.json" --tol 1e-6)
echo "$CMP_OUT" | tee "$OUT/G_fingerprint_invariant.txt"
echo "$CMP_OUT" | grep -q INVARIANT
echo "PASS:fingerprint_invariant trace=$TRACE"

# 4) Optional Apptainer .sif
if command -v apptainer >/dev/null 2>&1; then
  echo "=== building sif ==="
  apptainer build --force \
    --bind "$ROOT:/opt/flexaidds/src" \
    "$OUT/flexaidds_locked_x86_64.sif" \
    containers/flexaidds_locked_x86_64.def
  shasum -a 256 "$OUT/flexaidds_locked_x86_64.sif" | tee "$OUT/G_sif.sha256"
  echo "PASS:sif_build"
else
  echo "SKIP:sif_build (apptainer missing)"
  cat >"$OUT/G_sif_build_command.txt" <<EOF
STATUS=NOT_EXECUTED
REASON=apptainer_missing
# On Linux x86_64:
export FLEXAIDDS_SRC=$ROOT
apptainer build --bind "\${FLEXAIDDS_SRC}:/opt/flexaidds/src" \\
  $OUT/flexaidds_locked_x86_64.sif containers/flexaidds_locked_x86_64.def
EOF
fi

# 5) Optional Docker image
if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
  echo "=== docker build ==="
  docker build -f containers/Dockerfile.locked -t flexaidds:locked-g-smoke "$ROOT"
  echo "PASS:docker_build"
else
  echo "SKIP:docker_build (daemon unavailable)"
fi

# 6) Host smoke dock (software stack) if binary present
BIN="${FLEXAIDDS_BIN:-$ROOT/.swarm/cmaes/orchestrator/build_fast/FlexAIDdS}"
if [[ -x "$BIN" ]]; then
  REC="$ROOT/benchmarks/astex_diverse/astex_diverse/1G9V/1G9V_apo.pdb"
  LIG="$ROOT/benchmarks/astex_repro/smoke/1G9V/1G9V_dockin.sdf"
  CFG="$ROOT/validation_evidence/build_ab/C_cmaes_smoke/dock_config.json"
  if [[ -f "$REC" && -f "$LIG" && -f "$CFG" ]]; then
    export FLEXAIDDS_SEARCH=cmaes
    export FLEXAIDDS_CMAES_MAX_EVALS="${FLEXAIDDS_CMAES_MAX_EVALS:-2000}"
    export FLEXAIDDS_DATA_DIR="$(dirname "$BIN")"
    "$BIN" "$REC" "$LIG" -c "$CFG" -o "$OUT/host_smoke" >"$OUT/host_smoke.log" 2>&1 || true
    if rg -q 'backend=cmaes' "$OUT/host_smoke.log"; then
      echo "PASS:host_cmaes_smoke"
      if [[ -f "$OUT/host_smoke_cmaes_entropy.csv" ]]; then
        python3 analysis/collapse_fingerprint.py "$OUT/host_smoke_cmaes_entropy.csv" \
          --out "$OUT/host_smoke_fp.json"
        echo "PASS:host_smoke_fingerprint"
      fi
    else
      echo "WARN:host_cmaes_smoke (no backend=cmaes line; see host_smoke.log)"
    fi
  else
    echo "SKIP:host_smoke (inputs missing)"
  fi
else
  echo "SKIP:host_smoke (binary missing at $BIN)"
fi

# 7) Narval print-only
cat >"$OUT/G_narval_submit_command.txt" <<'EOF'
STATUS=NOT_EXECUTED
REASON=host_is_not_login_node_unless_sbatch_present
# Intended:
# sbatch --account=${CC_ACCOUNT} scripts/narval_cmaes_array.sh
EOF
if command -v sbatch >/dev/null 2>&1; then
  echo "NOTE: sbatch present — not auto-submitting (safety)"
fi

echo "=== G harness smoke DONE ==="
echo "artifacts_dir=$OUT"
