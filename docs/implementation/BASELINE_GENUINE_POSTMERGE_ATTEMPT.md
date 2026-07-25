# Post-merge genuine baseline re-aggregation attempt (2026-07-25)

> **PROCEDURE ONLY — no measured post-merge rates.**  
> `data_missing=true` · `reaggregate_ran=false` · Do **not** invent N / genuine % / BCR from this file.

**Hub:** [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md)  
**Compared against:** pre-merge genuine baseline in [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md)  
**Campaign baseline id:** `v_autonomous_20260724_160919`  
**Matrix pin (claim comparisons):** **`9dc93717dfed0698006d88dd6a9627bc`** (9dc9)

This document records **how** to re-aggregate once campaign trees are local. It is **not** a results table.

---

## Outcome of this session

| Check | Result |
|-------|--------|
| Local root | `$HOME/flexaidds_results` (`FLEXAIDDS_LOCAL_ROOT` default) |
| Search | `find … -maxdepth 4` / depth-limited local walk only — **no** `Path.rglob` / `find` under `Mobile Documents` / CloudDocs |
| `result.csv` / summary CSV on local APFS | **0 hits** |
| Campaign dirs present | Stubs only: `campaigns/C0_full85_9dc9_v{3,4,5,6}_20260724/` each with `MOVED_TO_ICLOUD.txt` |
| Bulk archive note | `~/flexaidds_results/MOVED_TO_ICLOUD_20260725T095624Z.txt` — local bulk removed after rsync to iCloud archive batch `20260725T095624Z` |
| Post-PR#300/#301 full result tree locally | **Not available** |
| Aggregator invoked with real rows | **No** |
| Post-merge N / genuine % / BCR / seed_echo | **Unknown** (do not invent) |

**Symlinks under `~/flexaidds_results/` that point into CloudDocs** (do not walk):

- `campaigns_icloud_archive` → `…/CloudDocs/FlexAIDdS_benchmarks/archived_from_ssd/flexaidds_results__campaigns`
- `icloud_archive_batch_20260725T095624Z` → `…/archive_batch_20260725T095624Z/flexaidds_results`
- `icloud_archive_batch_root_20260725T095624Z` → archive batch root

---

## Pre-merge baseline (do not invent; fixed reference)

Source: `docs/implementation/BASELINE_GENUINE_2026-07-24.md`

| Metric | Value |
|--------|------:|
| **GENUINE top-1 ≤2 Å** | **20 / 79 = 25.3%** |
| Best-cluster ≤2 Å (BCR) | 22 / 79 = **27.8%** |
| Seed-echo contamination | **0** |
| Targets with scores in denominator | **79** (80/85 finished; 5 never scored) |
| Campaign | `v_autonomous_20260724_160919` |
| Mode | Autonomous blind (not JCIM defined-cleft FLRP) |
| Matrix era | 9dc9 |

### Success PDBs (20) — baseline genuine set

```
1HNN  1HQ2  1OPK  1P62  1Q1G  1Q41  1R1H  1T46  1TZ8  1U4D
1UML  1V4S  1W1P  1Y6B  1Y6R  1YQY  1YWR  2BM2  2BSM  2HB1
```

### Overlap methodology (for future A/B vs this baseline)

When a post-merge (or new) campaign is aggregated:

1. Build set **G_base** = the 20 PDBs above (uppercase).  
2. Build set **G_new** = PDBs with **genuine success** under the same contract (below).  
3. Report:
   - `|G_new| / N_scored` and rate  
   - `|G_base ∩ G_new|` (retained successes)  
   - `|G_new − G_base|` (new wins)  
   - `|G_base − G_new|` (regressions)  
4. Do **not** mix contracts: genuine S1-style ≠ STRICT `claim_ready` ≠ 3Dsig **S_top10** bootstrap.  
5. Label mode (autonomous vs defined-cleft) and matrix MD5 on every table.  
6. Prefer **same denominator policy** when comparing to 25.3%: that baseline used **scored N=79**, not the frozen Astex-85 manifest (85). If using `aggregate_claim_metrics.py` default fixed-85 denominator, report both **rate@85** and **rate@scored** (or `N_claim` / scored rows) so numbers are not inflated/deflated by admission drops.

---

## Metric contract (fail-closed)

**Genuine success (matches baseline headline language):**

```
seed_echo == 0
AND elected rank-0 ordered RMSD ≤ 2.0 Å
```

- **RMSD column:** `rmsd_to_crystal` only (ordered direct).  
- **Never** `rmsd_hungarian` for S1 / genuine claim rates.  
- Fallback only if legacy rows lack `rmsd_to_crystal`: ordered `rmsd_top1` (see `elected_rmsd()` in aggregator).  
- Missing/blank `seed_echo` **fails closed** (not treated as success).  
- **BCR / sampling ceiling:** `best_cluster_rmsd` / `rmsd_bcr` / `conditional_scanned_pool_ceiling` ≤ 2.0 Å (aggregator **S3**, diagnostic only).  
- Matrix for claim comparisons: **9dc9** `9dc93717dfed0698006d88dd6a9627bc`.  
- Do **not** count seed-echo poses as successes even if RMSD ≤ 2.0.

Normative code: `scripts/aggregate_claim_metrics.py`  
Normative protocol: `benchmarks/protocols/admission_metrics_contract.md` (when present).

| Baseline term | Aggregator metric | Notes |
|---------------|-------------------|--------|
| Genuine top-1 ≤2 Å | **S1** (`--headline s1 --diagnostic-only`) | `is_s1`: seed_echo false + `rmsd_to_crystal` ≤ 2.0 |
| BCR | **S3** (`--headline s3 --diagnostic-only`) | Pool ceiling; never headline claim success |
| STRICT claim_ready | **STRICT** (default headline) | Stricter than baseline “genuine”; needs PB/tENCoM/etc. when columns present |
| seed_echo count | Dropped rows with `seed_echo!=0` + any non-zero flags in raw rows | Fail-closed admission |

---

## Exact aggregation commands (after local materialize)

### 1) Layout + env (local-first)

```bash
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_ROOT="$(git -C /path/to/FlexAIDdS rev-parse --show-toplevel)"
cd "$FLEXAIDDS_ROOT"
bash scripts/ensure_local_first_layout.sh
# optional: source claim staging helpers
# bash scripts/claim_local_staging_paths.sh
```

### 2) Materialize **one** campaign tree from iCloud archive → local (no CloudDocs walk)

Archive batch root (from local stub note; path is informational — use **rsync of known leaf names**, never `find` under CloudDocs):

```text
$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/archived_from_ssd/archive_batch_20260725T095624Z/flexaidds_results/
```

Restore candidates (pick the post-merge / claim campaign of interest; names from local stubs + baseline id):

```bash
ARCHIVE_ROOT="$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks/archived_from_ssd/archive_batch_20260725T095624Z/flexaidds_results"
LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"

# Example: known C0 stub names (post-9dc9 claim series)
for name in \
  C0_full85_9dc9_v6_20260724 \
  C0_full85_9dc9_v5_20260724 \
  C0_full85_9dc9_v4_20260724 \
  C0_full85_9dc9_v3_20260724 \
  v_autonomous_20260724_160919
do
  src="$ARCHIVE_ROOT/campaigns/$name"
  # also try run-dir at archive flexaidds_results root
  if [[ ! -d "$src" ]]; then src="$ARCHIVE_ROOT/$name"; fi
  # Prefer timeout-bounded single-file probes via icloud_safe_io before bulk rsync:
  #   python3 scripts/icloud_safe_io.py is-cloud "$src"
  # Full tree restore only for the chosen campaign (operator-driven):
  #   rsync -a "$src/" "$LOCAL_ROOT/campaigns/$name/"
done
```

Thin CSV-only seed (preferred when full pose trees are huge):

```bash
# After a single known path is available, materialize individual CSVs:
python3 scripts/icloud_safe_io.py materialize "/path/under/CloudDocs/.../result.csv"
# Hashes only via safe IO:
python3 scripts/icloud_safe_io.py md5 "/path/under/CloudDocs/.../result.csv"
```

**Forbidden:** `find` / `Path.rglob` under `Mobile Documents/` or any CloudDocs path.

### 3) Verify local CSV inventory (local only)

```bash
# maxdepth-limited; stays under local APFS
find "$LOCAL_ROOT/campaigns" -maxdepth 4 -name 'result.csv' 2>/dev/null | head -100
# or summary flat:
find "$LOCAL_ROOT/campaigns" -maxdepth 3 \( -name 'summary.csv' -o -name 'results.csv' -o -name 'claim_summary.csv' \) 2>/dev/null
```

Aggregator discovery rule (`load_campaign_rows`):

1. Prefer `$CAMP/*/result.csv` (one authoritative row per target dir).  
2. Else flat summary names: `astex_diverse_results.csv`, `astex_crossdock_85_results.csv`, `results.csv`, `summary.csv`, `claim_summary.csv`.

### 4) Run genuine (S1) + BCR (S3) aggregation — matrix 9dc9

```bash
CAMP="$LOCAL_ROOT/campaigns/<POSTMERGE_OR_AUTONOMOUS_CAMPAIGN>"
PIN=9dc93717dfed0698006d88dd6a9627bc
OUT_DIR="$LOCAL_ROOT/campaigns/analysis_postmerge_genuine"
mkdir -p "$OUT_DIR"

# Genuine-style headline (RMSD-only diagnostic; seed_echo fail-closed inside aggregator)
python3 scripts/aggregate_claim_metrics.py "$CAMP" \
  --matrix-md5 "$PIN" \
  --headline s1 --diagnostic-only \
  --json "$OUT_DIR/s1_genuine.json"

# BCR / pool ceiling (diagnostic only)
python3 scripts/aggregate_claim_metrics.py "$CAMP" \
  --matrix-md5 "$PIN" \
  --headline s3 --diagnostic-only \
  --json "$OUT_DIR/s3_bcr.json"

# Optional: STRICT claim_ready (stricter than baseline genuine)
python3 scripts/aggregate_claim_metrics.py "$CAMP" \
  --matrix-md5 "$PIN" \
  --headline strict \
  --json "$OUT_DIR/strict.json"
```

Single flat CSV:

```bash
python3 scripts/aggregate_claim_metrics.py --csv "$CAMP/summary.csv" \
  --matrix-md5 "$PIN" --headline s1 --diagnostic-only --json "$OUT_DIR/s1.json"
```

### 5) Seed-echo count (fail-closed audit)

From the JSON report:

- `N_dropped` + `dropped_rows` reasons containing `seed_echo!=0`  
- Or scan local CSVs only:

```bash
# Local tree only
python3 - <<'PY'
import csv, sys
from pathlib import Path
camp = Path(sys.argv[1])
n_echo = n_rows = 0
for p in sorted(camp.glob("*/result.csv")):
    with p.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        continue
    r = rows[0]
    n_rows += 1
    se = str(r.get("seed_echo", "")).strip()
    if se not in ("0", "0.0", "False", "false", "NO", "no"):
        n_echo += 1
        print(r.get("pdb_id") or r.get("pdb") or p.parent.name, "seed_echo=", se)
print(f"rows={n_rows} seed_echo_nonzero_or_missing={n_echo}")
PY
"$CAMP"
```

### 6) Fields to record when re-aggregation succeeds

Fill and promote this section (replace procedure-only status):

| Field | Value |
|-------|--------|
| campaign_id / path | |
| N_raw / N_claim / N_scored | |
| N_denominator (manifest 85 vs scored) | |
| genuine S1 n / rate | |
| BCR S3 n / rate | |
| seed_echo non-zero count | |
| matrix_md5 | must be 9dc9 for claim compare |
| binary path + sha256 | |
| git_commit of binary | |
| mode | autonomous vs defined-cleft |
| post-PR / merge note | e.g. after PR #300/#301 |
| \|G_base ∩ G_new\| / regressions / new wins | |

---

## Relation to comparative pipeline

- Hub: [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md).  
- This genuine table is **not** the three-engine **S_top10** comparative table (P5).  
- Comparative gates remain in [`CAMPAIGN_STATUS_2026-07-25.md`](CAMPAIGN_STATUS_2026-07-25.md) (P2 hold on oracle).  
- Matrix for claim/genuine comparisons stays **9dc9**; do not compare to **72d7** packing-fork rates without labeling confounds.

---

## Summary for agents

```
data_missing=true
reaggregate_ran=false
mode=procedure_only
baseline_genuine=20/79=25.3%
baseline_bcr=22/79=27.8%
baseline_seed_echo=0
baseline_campaign=v_autonomous_20260724_160919
matrix_pin=9dc93717dfed0698006d88dd6a9627bc
aggregator=scripts/aggregate_claim_metrics.py
genuine_cmd=python3 scripts/aggregate_claim_metrics.py <camp> --matrix-md5 9dc93717dfed0698006d88dd6a9627bc --headline s1 --diagnostic-only
bcr_cmd=... --headline s3 --diagnostic-only
local_result_csv_hits=0
archive_batch=20260725T095624Z
```
