#!/bin/bash
# Wrapper: appends -c mif_off.json to the harness's fixed argv (no hook exists).
exec /Users/lp.more/Projects/FlexAIDdS/ab_mac_20260806T133329/wt_post_cpu/build/FlexAID "$@" -c /Users/lp.more/Projects/FlexAIDdS/interventions/mif_off.json
