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
- `scripts/sync_three_engine_local_to_icloud.sh` — thin A/B0/B OUT mirror (timeout-bound).  
- `scripts/require_icloud_out.sh` — still used for pure-iCloud arms / legacy `--icloud-out`.  
- `scripts/sync_agent_homes_to_icloud.sh` + `scripts/agent_icloud_paths.py` — agent-home backup/restore (below).  
- Unit tests: `tests/test_icloud_safe_io.py`, `tests/test_agent_icloud_paths.py`.  

## Save benchmarks to iCloud (local → durable mirror)

Live claim / three-engine OUT stays under `$FLEXAIDDS_LOCAL_ROOT`. Push thin artifacts to CloudDocs:

```bash
# Claim campaign (result.csv + receipts; hang-safe)
bash scripts/sync_claim_local_to_icloud.sh --dry-run
bash scripts/sync_claim_local_to_icloud.sh

# Optional periodic loop (every N min)
nohup bash scripts/claim_icloud_sync_loop.sh &

# Three-engine A/B0/B arms
bash scripts/sync_three_engine_local_to_icloud.sh --dry-run
bash scripts/sync_three_engine_local_to_icloud.sh
```

Resolved paths (from `claim_local_staging_paths.sh` / `use_local_first_benchmark_storage.sh`):

| Role | Path |
|------|------|
| LOCAL claim OUT | `$FLEXAIDDS_LOCAL_ROOT/campaigns/<C0_CAMPAIGN_ID>` |
| REMOTE claim mirror | `$FLEXAIDDS_ICLOUD/results/campaigns/<C0_CAMPAIGN_ID>` |
| LOCAL three-engine | `$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/…` |
| REMOTE three-engine | `$FLEXAIDDS_ICLOUD/results/campaigns/three_engine/…` |

Never run live GA with OUT on CloudDocs. Full pose trees on iCloud only with `--with-poses`.

## Agent homes → iCloud (Claude, Claude Science, Codex, Grok)

**Live agent homes stay on local APFS** (`~/.claude`, Application Support Claude, `~/.claude-science`, `~/.codex`, `~/.grok`).  
**iCloud holds a durable mirror only** under `$FLEXAIDDS_ICLOUD/agent_homes/` — never replace homes with CloudDocs symlinks.

| Agent | Local source | Remote name under `agent_homes/` |
|-------|--------------|----------------------------------|
| Claude Code | `~/.claude` | `dot_claude` |
| Claude Desktop | `~/Library/Application Support/Claude` | `Application_Support_Claude` |
| Claude Science | `~/.claude-science` (selective) | `dot_claude_science` |
| Codex | `~/.codex` | `dot_codex` |
| Grok | `~/.grok` | `dot_grok` |

Default backup is **thin**: excludes caches, `vm_bundles/`, full `conda/` / `runtime/` under Claude Science (~11G reinstallable), large lock/sqlite-wal noise. Use `--full` only when you intentionally want a heavy archive.

```bash
# Map + dry-run (no I/O to CloudDocs)
bash scripts/sync_agent_homes_to_icloud.sh --dry-run
bash scripts/sync_agent_homes_to_icloud.sh --print-map

# Backup local → iCloud (timeout-wrapped rsync)
bash scripts/sync_agent_homes_to_icloud.sh --backup
bash scripts/sync_agent_homes_to_icloud.sh --backup --agents claude,codex,grok

# Pure path helpers (testable)
python3 scripts/agent_icloud_paths.py --print-map
python3 scripts/agent_icloud_paths.py --print-excludes claude_science
```

### Restore (iCloud archive → local) — seed rsync direction

Archive batches under  
`$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_<UTC>/`  
use the same remote names (`dot_claude`, `Application_Support_Claude`, …).

**Seed pattern (Claude only)** — matches the operator restore used after SSD archive:

```bash
A="$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/archived_from_ssd/archive_batch_20260725T095624Z"
rsync -a "$A/dot_claude/" "$HOME/.claude/"
rsync -a "$A/Application_Support_Claude/" "$HOME/Library/Application Support/Claude/"
```

**Shipped restore helpers** (print or run):

```bash
# Print exact seed block
bash scripts/sync_agent_homes_to_icloud.sh --print-seed-restore \
  --archive-batch "$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_20260725T095624Z"

# Print rsync lines for all agents in a batch
bash scripts/sync_agent_homes_to_icloud.sh --print-restore-cmds \
  --archive-batch "$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_20260725T095624Z"

# Execute restore (timeout-wrapped; does not --delete local extras)
bash scripts/sync_agent_homes_to_icloud.sh --restore --dry-run \
  --archive-batch "$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_20260725T095624Z"
bash scripts/sync_agent_homes_to_icloud.sh --restore \
  --archive-batch "$FLEXAIDDS_ICLOUD/archived_from_ssd/archive_batch_20260725T095624Z" \
  --agents claude,claude_app
```

Restore never makes live `$HOME` agent dirs depend on CloudDocs for runtime.

## Operator checklist (claim C0)

1. Kill any 0% CPU `benchmark_datasets` with no FlexAIDdS child.  
2. Seed local OUT with completed `result.csv` only (not full pose trees).  
3. `bash scripts/run_C0_claim_clean.sh` → OUT under `~/flexaidds_results/campaigns/…`.  
4. Confirm FlexAIDdS child **CPU ≫ 0** within ~30 s.  
5. Keep `claim_icloud_sync_loop.sh` running for durable mirror.  
6. Never dual-launch two claim workers on the same campaign.  
7. Periodically: `bash scripts/sync_agent_homes_to_icloud.sh --backup` for agent configs/sessions.  

## Note on local disk pressure

If local Data volume is nearly full, free space before claim-scale GA (pose trees are large). Prefer thin iCloud mirror (result.csv only) rather than full dual trees.
