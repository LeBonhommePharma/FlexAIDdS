#!/usr/bin/env bash
#
# Small Real Benchmark Example: 1STP (Biotin / Streptavidin) Self-Docking
#
# This example demonstrates the modern FlexAIDδS DatasetRunner features
# introduced in 2026-05:
#   - Per-entry checkpointing (individual target results)
#   - --resume with automatic cost-hint loading (CostHistory + EMA)
#   - --package for professional reproducibility artifacts
#   - Rich cost/timing reporting in final Markdown
#
# It uses a single well-known complex (PDB 1STP) as a "tiny real benchmark".
#
# Usage (safe dry-run by default):
#   bash .grok/skills/flexaidds/examples/small_real_benchmark_1stp.sh
#
# For a real run (requires working FlexAIDδS binary + data):
#   FLEXAIDDS_BINARY=/path/to/FlexAIDδS \
#   bash .grok/skills/flexaidds/examples/small_real_benchmark_1stp.sh --real
#
# Recommended for publications / audits:
#   Always end with --package so you get the full validation zip.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SKILL_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

DATASET_SLUG="1stp_biotin_demo"
RESULTS_DIR="results/benchmarks/1stp_demo"
TIER=1

echo "=== FlexAIDδS Small Real Benchmark Example ==="
echo "Complex: 1STP (Streptavidin + Biotin)"
echo "Features demonstrated: per-entry saving, resume + cost history, --package"
echo ""

# Allow overriding the binary
BINARY="${FLEXAIDDS_BINARY:-}"

# Parse simple flag
REAL_RUN=false
if [[ "${1:-}" == "--real" ]]; then
    REAL_RUN=true
fi

if $REAL_RUN && [[ -z "$BINARY" ]]; then
    echo "ERROR: --real requested but FLEXAIDDS_BINARY is not set."
    echo "Example: FLEXAIDDS_BINARY=/path/to/your/FlexAIDδS bash $0 --real"
    exit 1
fi

# Step 0: Ensure runtime data (matrices + *.def files)
echo "[0/4] Ensuring critical docking data..."
python3 "$SKILL_ROOT/scripts/ensure_docking_data.py" --info

# Step 1: (Optional) Redock 1STP to produce real inputs
# For this tiny demo we use a synthetic tiny dataset defined inline.
# In a true "real" workflow you would run:
#   python3 "$SKILL_ROOT/scripts/redock_from_pdb.py" 1STP --include-modified-residues
# and then point a custom YAML at the produced files.

echo "[1/4] Preparing tiny 1STP demo dataset (synthetic for demo safety)..."

mkdir -p "$SKILL_ROOT/examples/data/1stp_demo"

# Create a minimal dataset YAML for this single-complex demo
cat > "$SKILL_ROOT/examples/data/1stp_demo.yaml" << 'YAML'
slug: 1stp_biotin_demo
name: "1STP Biotin Self-Docking Demo (Tiny)"
description: "Minimal single-complex benchmark for testing per-entry resume + cost tracking features."
tier: 1
tier1_subset_size: 1
benchmark_order: 0
structural_states:
  - holo
targets:
  - 1stp
metrics:
  - docking_power_top1
  - entropy_rescue_rate
expected_baselines:
  docking_power_top1: 1.0
  entropy_rescue_rate: 0.0
YAML

echo "  Created: $SKILL_ROOT/examples/data/1stp_demo.yaml"

# Step 2: First run (or partial run)
echo "[2/4] Running first pass (may be partial)..."
CMD=(
    python3 "$SKILL_ROOT/scripts/dataset_runner.py"
    --dataset "$DATASET_SLUG"
    --datasets-dir "$SKILL_ROOT/examples/data"
    --tier "$TIER"
    --results-dir "$RESULTS_DIR"
    --workers 2
    --package
)

if $REAL_RUN; then
    CMD+=(--binary "$BINARY")
else
    CMD+=(--dry-run)
    echo "  (Running in --dry-run mode for safety. Use --real for actual docking.)"
fi

"${CMD[@]}"

echo ""
echo "[3/4] Demonstrating --resume + cost-aware scheduling..."
# Second run with --resume (will skip already-completed entries and use CostHistory)
RESUME_CMD=(
    python3 "$SKILL_ROOT/scripts/dataset_runner.py"
    --dataset "$DATASET_SLUG"
    --datasets-dir "$SKILL_ROOT/examples/data"
    --tier "$TIER"
    --results-dir "$RESULTS_DIR"
    --workers 2
    --resume
    --package
)

if $REAL_RUN; then
    RESUME_CMD+=(--binary "$BINARY")
else
    RESUME_CMD+=(--dry-run)
fi

"${RESUME_CMD[@]}"

echo ""
echo "[4/4] Verification"
echo "  Per-entry results: $RESULTS_DIR/1stp_biotin_demo/tier1/"
echo "  Rich manifest with costs: $RESULTS_DIR/1stp_biotin_demo/tier1/_entry_manifest.json"
echo "  Cost history (EMA):       $RESULTS_DIR/1stp_biotin_demo/tier1/.cost_history.json"
echo "  Validation package:       Look for flexaidds_validation_package_*.zip in parent of results"
echo ""
echo "=== Example complete ==="
echo "Open the generated Markdown report and VALIDATION_SUMMARY.md to see the new rich cost tables."
echo "This pattern (tiny focused benchmark + --resume + --package) is excellent for development and publications."