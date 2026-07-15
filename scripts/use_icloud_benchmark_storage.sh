#!/usr/bin/env bash
# Pin all FlexAIDdS benchmark I/O to iCloud Drive for this shell (+ ~/.flexaidds_env).
# No machine-specific usernames: uses $HOME only.
set -euo pipefail

ICLOUD_ROOT="${HOME}/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$ICLOUD_ROOT}"
export FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$FLEXAIDDS_ICLOUD/results}"
export FLEXAIDDS_QUEUE_ROOT="${FLEXAIDDS_QUEUE_ROOT:-$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1}"

mkdir -p \
  "$FLEXAIDDS_ICLOUD/results/working" \
  "$FLEXAIDDS_ICLOUD/results/campaigns" \
  "$FLEXAIDDS_ICLOUD/results/oracle_ceiling" \
  "$FLEXAIDDS_ICLOUD/logs" \
  "$FLEXAIDDS_ICLOUD/queues" \
  "$FLEXAIDDS_ICLOUD/provenance"

ENVF="${HOME}/.flexaidds_env"
touch "$ENVF"
python3 - <<'PY'
from pathlib import Path
import os
env = Path.home() / ".flexaidds_env"
lines = env.read_text().splitlines() if env.exists() else []
drop = {"export FLEXAIDDS_ICLOUD=", "export FLEXAIDDS_RESULTS=", "export FLEXAIDDS_QUEUE_ROOT="}
keep = [l for l in lines if not any(l.startswith(p) for p in drop)]
icloud = Path.home() / "Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"
keep += [
    f"export FLEXAIDDS_ICLOUD='{icloud}'",
    f"export FLEXAIDDS_RESULTS='{icloud}/results'",
    f"export FLEXAIDDS_QUEUE_ROOT='{icloud}/queues/three_engine_entropy_q1'",
]
env.write_text("\n".join(keep) + "\n")
print("Updated", env)
print("FLEXAIDDS_ICLOUD=", icloud)
print("FLEXAIDDS_RESULTS=", icloud / "results")
PY

echo "Source this file or open a new shell after: source ~/.flexaidds_env"
echo "Queue scripts: $FLEXAIDDS_QUEUE_ROOT/scripts/"
