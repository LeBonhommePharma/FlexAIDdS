# 3Dsig S_top10 — 10k bootstrap runbook

**Metric:** S_top10 = success if min(RMSD among top-10 ranked modes) **&lt; 2.0 Å**  
**Headline:** median success rate over **10 000** bootstrap resamples of N cases  
**Deck targets (Astex Diverse):** FlexAID **0.66** · FlexAIDdS **0.69**  
**Protocol:** `docs/implementation/3dsig_red_pair_protocol.md`

Live OUT (local-first)::

```text
~/flexaidds_results/campaigns/three_engine/{A,B0,B}/3dsig_r10/
```

Do **not** dual-launch arms or thrash iCloud while docks run.

---

## 1. Schema gaps (audit)

| Source | Top-10 RMSDs? | Notes |
|--------|---------------|--------|
| `parse_flexaid_arm_results.py` | **Yes (now)** | Emits `mode_rmsd_0..9`, `success_s_top10`, `n_top10_rmsds`, `rmsd_extract_gap` |
| Older `result.csv` | **No** | Only `rmsd_top1` / `rmsd_bcr` — incomplete for S_top10 |
| Ranked pose PDBs | **Required** | `{PDB}_r{R}_{rank}.pdb` or FO `{PDB}_r{R}_{minPts}_{rank}.pdb` with REMARK RMSD |
| `.cad` only | **No crystal RMSD** | Inter-cluster RMSDs only; cannot recover S_top10 |
| `.rrd` | **Yes (fallback)** | Written after pose PDBs; heads usable if present |

**Blocking gap (2026-07-15 pilot):** FlexAID A binary **segfaults at clustering** after writing `.cad` / `_INI.pdb` / `_par.res`, so **no ranked pose PDBs** and empty `mode_rmsd_*`. Logs show:

```text
clustering all individuals in GA... Segmentation fault: 11
```

Until pose PDBs (or `.rrd`) appear, bootstrap correctly reports **NA**.

---

## 2. End-to-end recipe

### Extract A → bootstrap

```bash
ROOT="${FLEXAIDDS_ROOT:-$HOME/Projects/FlexAIDdS}"
OUTA="$HOME/flexaidds_results/campaigns/three_engine/A/3dsig_r10"
mkdir -p "$HOME/flexaidds_results/campaigns/three_engine/bootstrap"

python3 "$ROOT/scripts/extract_3dsig_s_top10_from_arm.py" \
  --arm-dir "$OUTA" --score cf --strategy auto \
  --json-out "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/A_cases.json"

python3 "$ROOT/scripts/bootstrap_3dsig_s_top10.py" \
  --cases "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/A_cases.json" \
  --bootstraps 10000 --label A \
  --json-out "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/A_s_top10.json"
```

### Extract B (entropy) → bootstrap

```bash
OUTB="$HOME/flexaidds_results/campaigns/three_engine/B/3dsig_r10"

python3 "$ROOT/scripts/extract_3dsig_s_top10_from_arm.py" \
  --arm-dir "$OUTB" --score acf --strategy auto \
  --json-out "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/B_cases.json"

python3 "$ROOT/scripts/bootstrap_3dsig_s_top10.py" \
  --cases "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/B_cases.json" \
  --bootstraps 10000 --label B \
  --json-out "$HOME/flexaidds_results/campaigns/three_engine/bootstrap/B_s_top10.json"
```

### Optional B0 (master CF control)

Same as A with `--label B0` and `B0/3dsig_r10`.

### Compare to deck

Bootstrap stdout includes:

```text
deck compare: live_median=… vs FlexAID=0.66 / FlexAIDdS=0.69
```

Or re-parse JSON:

```bash
python3 - <<'PY'
import json
from pathlib import Path
base = Path.home()/"flexaidds_results/campaigns/three_engine/bootstrap"
for lab in ("A","B0","B"):
    p = base/f"{lab}_s_top10.json"
    if not p.is_file():
        print(lab, "missing"); continue
    d = json.loads(p.read_text())
    print(lab, "median=", d.get("median"), "n=", d.get("n_cases"), "status=", d.get("status"))
print("deck: FlexAID=0.66 FlexAIDdS=0.69")
PY
```

### Alternative: arm-dir via result.csv

After re-parse with the extended parser:

```bash
# per finished case (safe while other docks run)
python3 scripts/parse_flexaid_arm_results.py \
  --arm A --pdb 1G9V \
  --out-dir ~/flexaidds_results/campaigns/three_engine/A/3dsig_r10/1G9V \
  --matrix-md5 72d7c7396702331d96ff12d18f831796

python3 scripts/bootstrap_3dsig_s_top10.py \
  --arm-dir ~/flexaidds_results/campaigns/three_engine/A/3dsig_r10 \
  --label A --bootstraps 10000
```

Prefer **extract_*** when pose PDBs exist but `result.csv` is stale.

---

## 3. Ranking strategies (`extract_3dsig_s_top10_from_arm.py`)

| `--strategy` | Behavior |
|--------------|----------|
| `auto` (default) | `global` if any multi-rank pose exists, else `restart_heads` |
| `global` | Pool emission ranks &lt; 10 across restarts, sort by score, keep 10 |
| `restart_heads` | Rank-0 from each restart, sort by score (≤10) — natural for R=10 |

| `--score` | Arm |
|-----------|-----|
| `cf` | A / B0 — CF.app ascending |
| `acf` | B — soft-β \(\tilde G\) when REMARK present, else CF |

---

## 4. Partial data policy

- **0 cases with RMSDs** → bootstrap prints **NA**, exit 0, JSON `status: "NA"`.
- **Partial N** (e.g. 3/8 pilot) → bootstrap runs on available cases only; label as pilot, not Astex-85 claim.
- **cad_only** gap → do not invent RMSD from `.cad` inter-cluster matrix.
- Never claim deck reproduction until N matches the claimed set and pose RMSDs exist.

---

## 5. Re-run when docks finish

Safe, non-destructive polling (no kill, no dual-launch):

```bash
# once more cases show *_0.pdb or *.rrd:
python3 scripts/extract_3dsig_s_top10_from_arm.py --arm-dir "$OUTA" --json-out .../A_cases.json
python3 scripts/bootstrap_3dsig_s_top10.py --cases .../A_cases.json --label A --bootstraps 10000
```

iCloud sync only after arms complete: `bash scripts/sync_three_engine_local_to_icloud.sh`.
