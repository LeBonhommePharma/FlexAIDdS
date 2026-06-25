# Multi-Cleft Restoration Comparison Snapshot

Snapshot generated on 2026-06-25 from the completed revived multi-cleft campaign.

Regeneration command:

```bash
python3 summarize_multicleft_astex.py \
  --result-dir /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/multicleft_full_top3_dirty_transport/results \
  --output-prefix /private/tmp/multicleft_astex_audit
```

| Campaign | Evidence | Cleft handling | Native/site handling | N | Successes <2 A | Success % | Mean RMSD | Median RMSD | Notes |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| Old 2017-2019 multi-cleft | evidence pending | independent GA per major cleft | evidence pending | TBD | TBD | TBD | TBD | TBD | Do not fill until old logs/table are recovered. |
| Current single-GA v112 | /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v112_20260624_2130_oracle_full85/astex_diverse_results.csv | single GA | self-dock/oracle reference protocol | 85 | 38 | 44.7 | 4.21 | 2.35 | Local current comparator. |
| Revived multi-cleft | /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/multicleft_full_top3_dirty_transport/results/*/result.csv (242 files) | one DatasetRunner entry per Get_Cleft sphere file; selected by best_score | autonomous; no oracle-site env | 85 | 4 | 4.7 | 19.76 | 13.29 | Oracle-best diagnostic successes across clefts: 5. Clean rerun remains live for provenance. |

The completed dirty-transport run provides the full 85-target revived multi-cleft comparison row. The clean rerun remains live as a provenance check, but it is not required for the comparison table.
