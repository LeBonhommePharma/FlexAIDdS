# iCloud benchmark storage (production)

> **Source of truth:** `AGENTS.md` § *Benchmark storage (Non-Negotiable) — local-first / thin-iCloud*. This page expands paths and operator runbooks.  
> **Comparative science hub:** [`implementation/COMPARATIVE_SCIENCE_README.md`](implementation/COMPARATIVE_SCIENCE_README.md).

**Architecture (enforced default):**

| Layer | Rule |
|-------|------|
| **Live compute** | Local APFS only — `$FLEXAIDDS_LOCAL_ROOT` (default `~/flexaidds_results`) |
| **iCloud** | Thin durable mirror (`result.csv`, receipts, thin OPS) — **not** a working FS |
| **Any CloudDocs read** | Timeout-bounded via `scripts/icloud_safe_io.py`, **or** materialize → local pin-cache then hash/read |
| **Ops / cron** | Local-first scan; wall-clock cap; optional walker reaper |
| **Dockers** | **Never** killed by reaper (`FlexAIDdS`, `benchmark_datasets`, claim `caffeinate`) |

**Non-negotiable (durability):** final claim artifacts (`result.csv`, RUN_RECEIPT, thin status) must land under **iCloud Drive** as a **mirror** of local work.

**Non-negotiable (anti-hang, 2026-07-15+):** **live GA I/O must not hit CloudDocs fileprovider.** Writing / re-reading thousands of pose PDBs under `Mobile Documents/…` stalls DatasetRunner (0% CPU for hours, no FlexAIDdS child). Fix = **local live OUT + periodic thin sync to iCloud**.

## Capacity

- **iCloud Drive quota:** ~2 TB (operator: ~1.5 TB free).  
- **Local APFS:** use `~/flexaidds_results` for live claim OUT, binaries, and logs; free space must allow full pose trees.

## Canonical paths (environment)

| Variable | Purpose | Typical value |
|----------|---------|----------------|
| `FLEXAIDDS_LOCAL_ROOT` | **Live work root** (GA OUT, logs, bins, pin-cache) | `~/flexaidds_results` |
| `FLEXAIDDS_ICLOUD` | Durable thin-mirror root | `…/CloudDocs/FlexAIDdS_benchmarks` |
| `FLEXAIDDS_RESULTS` | iCloud campaign mirror | `$FLEXAIDDS_ICLOUD/results` |
| `FLEXAIDDS_QUEUE_ROOT` | Three-engine queue (inputs; optional logs) | `$FLEXAIDDS_ICLOUD/queues/three_engine_entropy_q1` |
| `C0_CLAIM_LOCAL_OUT` | Live claim OUT | `$FLEXAIDDS_LOCAL_ROOT/campaigns/C0_full85_claim_…` |
| `C0_CLAIM_ICLOUD_OUT` | Durable mirror | `$FLEXAIDDS_ICLOUD/results/campaigns/C0_full85_claim_…` |
| `FLEXAIDDS_PIN_CACHE` | Materialize pin-cache (optional override) | `$FLEXAIDDS_LOCAL_ROOT/pins/materialize` |

```bash
source ~/.flexaidds_env   # optional build pins
bash "$FLEXAIDDS_ROOT/scripts/ensure_local_first_layout.sh"
source "$FLEXAIDDS_ROOT/scripts/claim_local_staging_paths.sh"
# default claim launch (local OUT):
bash "$FLEXAIDDS_ROOT/scripts/run_C0_claim_clean.sh"
# optional background mirror:
nohup bash "$FLEXAIDDS_ROOT/scripts/claim_icloud_sync_loop.sh" &
```

## What goes where (anti-hang)

| Content | Location |
|---------|----------|
| **Live** GA poses, `dock_config`, stdout/stderr, elected PDBs | **Local** `$C0_CLAIM_LOCAL_OUT` |
| Live claim logs / pid / ops status | **Local** `$C0_CLAIM_LOCAL_LOGDIR` |
| FlexAIDdS + `benchmark_datasets` Mach-O + `data/` | **Local** `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/C` |
| Pin-cache (materialized CloudDocs files) | **Local** `$FLEXAIDDS_LOCAL_ROOT/pins/materialize` |
| Manifest / site inputs (read-only) | Prefer **local copy** of inputs; fall back to queue |
| Durable `result.csv` + RUN_RECEIPT + thin OPS status | **iCloud mirror** via `sync_claim_local_to_icloud.sh` (every N min) |
| Full pose tree on iCloud | Optional `--with-poses` only; not required for resume |

## Why direct iCloud OUT hangs

1. DatasetRunner **re-pools** completed targets (Fix B / 3Dsig re-elect) by walking pose trees.  
2. CloudDocs coordinates every open/stat; large pose dirs → **fileprovider stall**.  
3. Parent `benchmark_datasets` sits at **0% CPU** with **no FlexAIDdS child** for hours.  
4. Ops tools that `find` / `rglob` / bulk `md5` under CloudDocs also hang — use local logs + glob only.

## Production anti-hang toolkit (required for agents & cron)

| Tool | Role |
|------|------|
| `scripts/ensure_local_first_layout.sh` | Idempotent local dirs: `campaigns/`, `logs/{ops,ops_monitor,C0_claim}`, `pins/materialize`, `three_engine_entropy_q1/{bin/C,data,inputs}` |
| `scripts/icloud_safe_io.py` | Detect CloudDocs; **timeout-bounded** read/md5; **materialize** files to `~/flexaidds_results/pins/materialize/` before hashing |
| `scripts/benchmark_ops_monitor.py` | **Local-first** campaign scan; one-level `*/result.csv` only; no `**/` rglob on CloudDocs |
| `scripts/run_benchmark_ops_monitor.sh` | Wall-clock timeout (default 90s); optional `--reap-walkers` |
| `scripts/reap_hung_icloud_walkers.sh` | Kill stuck `find`/`md5`/`rglob` on CloudDocs **only** — **never** FlexAIDdS / `benchmark_datasets` / claim caffeinate / sync loop |
| `scripts/claim_icloud_sync_loop.sh` | Thin local → iCloud mirror every N min; status written **local first** |

### Rules for AI agents (hard)

1. **Never** `find ~` or `Path.rglob` under `Mobile Documents/`.  
2. Hash CloudDocs files only via:  
   `python3 scripts/icloud_safe_io.py md5 <path>` (materializes then hashes locally).  
3. Prefer paths under `$FLEXAIDDS_LOCAL_ROOT` (`~/flexaidds_results`).  
4. If ops cron stalls, run:  
   `bash scripts/reap_hung_icloud_walkers.sh` then re-run the monitor.  
5. Do **not** set `HOMEBREW_NO_REQUIRE_TAP_TRUST=1` as a hang “fix” — unrelated.

### CLI examples

```bash
# Safe hash (works even if path is under CloudDocs)
python3 scripts/icloud_safe_io.py md5 \
  "$HOME/Documents/PhD/AtomTypes/MC_st0r5.2_6.dat"

# Materialize a CloudDocs file into local pin-cache
python3 scripts/icloud_safe_io.py materialize "/path/under/CloudDocs/file.dat"

# Ops tick (scheduler)
bash scripts/run_benchmark_ops_monitor.sh --reap-walkers
```

## Enforcement

- `scripts/ensure_local_first_layout.sh` — called by `run_C0_claim_clean.sh` when `USE_LOCAL=1`.  
- `scripts/run_C0_claim_clean.sh` — **default `FLEXAIDDS_CLAIM_LOCAL=1`** (local OUT). Use `--icloud-out` only for experiments.  
- `scripts/sync_claim_local_to_icloud.sh` — thin rsync of `*/result.csv` (+ optional extras), timeout-bound.  
- `scripts/claim_icloud_sync_loop.sh` — every N min: sync + hang-safe ops status (no deep iCloud find).  
- `scripts/require_icloud_out.sh` — still used for pure-iCloud arms / legacy `--icloud-out`.  
- Unit tests: `tests/test_icloud_safe_io.py` (`python3 -m pytest tests/test_icloud_safe_io.py -q`).  


## Operator checklist (claim C0)

1. Kill any 0% CPU `benchmark_datasets` with no FlexAIDdS child.  
2. Seed local OUT with completed `result.csv` only (not full pose trees).  
3. `bash scripts/run_C0_claim_clean.sh` → OUT under `~/flexaidds_results/campaigns/…`.  
4. Confirm FlexAIDdS child **CPU ≫ 0** within ~30 s.  
5. Keep `claim_icloud_sync_loop.sh` running for durable mirror.  
6. Never dual-launch two claim workers on the same campaign.  

## Note on local disk pressure

If local Data volume is nearly full, free space before claim-scale GA (pose trees are large). Prefer thin iCloud mirror (result.csv only) rather than full dual trees.
