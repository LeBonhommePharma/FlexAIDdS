# Multi-Cleft Restoration Comparison Snapshot

Snapshot generated on 2026-06-25 from the live restored multi-cleft campaign.

Regeneration command:

```bash
python3 summarize_multicleft_astex.py \
  --result-dir /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/multicleft_full_top3_clean_ef450149_reuseclefts/results \
  --output-prefix /private/tmp/multicleft_astex_audit
```

| Campaign | Evidence | Cleft handling | Native/site handling | N | Successes <2 A | Success % | Mean RMSD | Median RMSD | Notes |
|---|---|---|---|---:|---:|---:|---:|---:|---|
| Old 2017-2019 multi-cleft | evidence pending | independent GA per major cleft | evidence pending | TBD | TBD | TBD | TBD | TBD | Do not fill until old logs/table are recovered. |
| Current single-GA v112 | /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/v112_20260624_2130_oracle_full85/astex_diverse_results.csv | single GA | self-dock/oracle reference protocol | 85 | 38 | 44.7 | 4.21 | 2.35 | Local current comparator. |
| Revived multi-cleft | /Users/lp.more/Documents/PhD/Programs/FlexAIDdS/results/multicleft_full_top3_clean_ef450149_reuseclefts/results/*/result.csv (187 files) | one DatasetRunner entry per Get_Cleft sphere file; selected by best_score | autonomous; no oracle-site env | 67 | 3 | 4.5 | 20.55 | 19.23 | Oracle-best diagnostic successes across clefts: 4. |

The revived multi-cleft run is still live. This snapshot reflects the current partial campaign state, not a final 85/85 closeout.
