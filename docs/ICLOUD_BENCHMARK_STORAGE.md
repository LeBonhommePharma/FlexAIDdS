# iCloud benchmark storage (production)

**Non-negotiable:** every production / claim / three-engine / residual campaign writes **final artifacts to iCloud Drive**.

## Capacity

- **iCloud Drive quota:** ~2 TB (operator: ~1.5 TB free).  
- **Local APFS cache:** may be tight even when cloud free space is large — prefer streaming results to iCloud OUT, keep binaries on local disk.

## Canonical paths (environment)

| Variable | Purpose | Typical value |
|----------|---------|----------------|
| `FLEXAIDDS_ICLOUD` | Benchmark root | `…/CloudDocs/FlexAIDdS_benchmarks` |
| `FLEXAIDDS_RESULTS` | Campaign results | `$FLEXAIDDS_ICLOUD/results` |
| `FLEXAIDDS_QUEUE_ROOT` | Three-engine queue | `$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1` |
| `FLEXAIDDS_WORKING` | Scratch under results | `$FLEXAIDDS_RESULTS/working` |

```bash
source ~/.flexaidds_env   # optional build pins
source "$FLEXAIDDS_ROOT/scripts/use_icloud_benchmark_storage.sh"
```

## What goes where

| Content | Location |
|---------|----------|
| `result.csv`, elected poses, RUN_RECEIPT, claim aggregates | **iCloud** `$FLEXAIDDS_RESULTS/campaigns/…` |
| Queue logs, work trees, STATUS | **iCloud** `$FLEXAIDDS_QUEUE_ROOT/` |
| FlexAID / FlexAIDdS / `benchmark_datasets` Mach-O | **Local** staging (queue `bin/` → local symlink) |
| Oracle-ceiling archive (if still under `~/flexaidds_results`) | **Migrate** to `$FLEXAIDDS_RESULTS/archive/` when local free allows |

## Enforcement

- `scripts/run_C0_full85.sh` — refuses non-iCloud `OUT`  
- `scripts/run_flexaid_arm_pilot8.sh` — refuses non-iCloud `OUT`  
- `scripts/require_icloud_out.sh` — shared helper  
- `scripts/monitor_all_benchmarks.py` — reports paths under iCloud first  

## Operator checklist

1. Before launch: `source scripts/use_icloud_benchmark_storage.sh`  
2. Confirm: `echo $FLEXAIDDS_RESULTS` contains `CloudDocs`  
3. After run: `ls "$FLEXAIDDS_RESULTS/campaigns/<name>"`  
4. Never point claim campaigns at `~/flexaidds_results/…` for new work  

## Note on local disk pressure

If `df` shows local Data volume nearly full while iCloud still has free cloud quota, **do not** fill `~/flexaidds_results`. Keep writing campaign OUT under CloudDocs (iCloud evicts local copies as needed). Stage large archives with `brctl` / wait for upload rather than dual local copies.
