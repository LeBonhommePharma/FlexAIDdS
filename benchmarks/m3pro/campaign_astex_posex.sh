#!/bin/bash
# Campaign: Astex Diverse (GCE validation) + PoseX restart
# Date: 2026-07-08
# Purpose: (1) Validate grand canonical ensemble on Astex tier-1
#          (2) Restart PoseX with corrected pop_size=900

set -e

# Resolve repo root from this script (never hard-code machine paths).
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
LAUNCHER="${REPO}/benchmarks/m3pro/grok_master_launcher.sh"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${REPO}/results/campaign_logs"

mkdir -p "$LOG_DIR"

echo "=========================================="
echo "CAMPAIGN: Astex Diverse (GCE) + PoseX"
echo "Timestamp: $TIMESTAMP"
echo "=========================================="

# Phase 1: Astex Diverse tier-1 with GCE validation
echo ""
echo "PHASE 1: Running Astex Diverse tier-1"
echo "  - Validates grand canonical ensemble features"
echo "  - Expected runtime: ~30-45 min on M3 Pro with 4 workers"
echo "  - Dataset: astex_diverse (19 complexes, tier-1)"
echo ""

ASTEX_LOG="${LOG_DIR}/astex_diverse_gce_${TIMESTAMP}.log"
bash "$LAUNCHER" \
  --datasets astex_diverse \
  --tier 1 \
  --run-id "astex_diverse_gce_${TIMESTAMP}" \
  2>&1 | tee "$ASTEX_LOG"

ASTEX_STATUS=$?
if [ $ASTEX_STATUS -eq 0 ]; then
  echo "✓ Astex Diverse tier-1 PASSED"
else
  echo "✗ Astex Diverse tier-1 FAILED (exit code $ASTEX_STATUS)"
  echo "Check log: $ASTEX_LOG"
  exit $ASTEX_STATUS
fi

# Phase 2: PoseX cross-docking with corrected pop_size=900
echo ""
echo "PHASE 2: Restarting PoseX cross-docking"
echo "  - Configuration: pop_size=900 (reduced from 2000 to fit M3 Pro Metal)"
echo "  - Dataset: posex_cd (1312 pairs across ~250 targets)"
echo "  - Expected runtime: ~3-5 hours with 4 workers"
echo "  - Success metric: RMSD ≤ 2.0 Å across Astex Diverse, Self-Docking, Cross-Docking"
echo ""

POSEX_LOG="${LOG_DIR}/posex_cd_corrected_${TIMESTAMP}.log"
bash "$LAUNCHER" \
  --datasets posex_cd \
  --ga-population 900 \
  --run-id "posex_cd_corrected_${TIMESTAMP}" \
  2>&1 | tee "$POSEX_LOG"

POSEX_STATUS=$?
if [ $POSEX_STATUS -eq 0 ]; then
  echo "✓ PoseX campaign PASSED"
else
  echo "✗ PoseX campaign FAILED (exit code $POSEX_STATUS)"
  echo "Check log: $POSEX_LOG"
  exit $POSEX_STATUS
fi

echo ""
echo "=========================================="
echo "CAMPAIGN COMPLETE"
echo "  Phase 1 (Astex): ✓"
echo "  Phase 2 (PoseX): ✓"
echo "=========================================="
echo ""
echo "Next: Analyze results and promote to Tier-2 suite"
