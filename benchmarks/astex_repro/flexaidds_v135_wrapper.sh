#!/bin/bash
# v135 FlexAIDdS wrapper — patches GA params in rN/dock_config.json before exec.
#
# The benchmark_datasets oracle-ceiling runner hardcodes bad GA params:
#   sharing_alpha=2.28571 (correct: 4)
#   num_chromosomes=1750  (correct: 1000)
#   num_generations=2000  (correct: 875)
#
# This wrapper intercepts the --config <path> argument, patches those fields
# in-place, then execs the real FlexAIDdS binary with the original argument list.

REAL_BINARY="/Users/lp.more/Projects/FlexAIDdS/build_v135/FlexAIDdS"

# Extract --config <path> from argv
CONFIG_PATH=""
NEXT_IS_CONFIG=0
for arg in "$@"; do
    if [ "$NEXT_IS_CONFIG" = "1" ]; then
        CONFIG_PATH="$arg"
        break
    fi
    [ "$arg" = "--config" ] && NEXT_IS_CONFIG=1
done

# Patch GA params in-place before FlexAIDdS reads the config
if [ -n "$CONFIG_PATH" ] && [ -f "$CONFIG_PATH" ]; then
    python3 - "$CONFIG_PATH" <<'PYEOF'
import json, sys
path = sys.argv[1]
with open(path) as f:
    d = json.load(f)
if 'ga' in d:
    d['ga']['sharing_alpha']   = 4      # was 2.28571
    d['ga']['num_chromosomes'] = 1000   # was 1750
    d['ga']['num_generations'] = 875
with open(path, 'w') as f:
    json.dump(d, f, indent=2)
PYEOF
fi

exec "$REAL_BINARY" "$@"
