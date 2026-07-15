# Matrix Pin — JCIM 2015 (`MC_st0r5.2_6.dat`)

**Canonical claim MD5**: `72d7c7396702331d96ff12d18f831796`  
**File name**: `MC_st0r5.2_6.dat`  
**Size (claim pin bytes)**: 18043  
**Provenance lineage**: Sippl reference-ratio VCT contact potential used with FlexAID (Gaudreault & Najmanovich, *J. Chem. Inf. Model.* 2015).  
**Validation timestamp (this audit)**: **2026-07-15T12:05:12Z** (local 2026-07-15 08:05:12 EDT)

> Related: default pin also recorded in `benchmarks/protocols/admission_metrics_contract.md`  
> and `benchmarks/protocols/three_engine_entropy_comparison.md`.  
> Scientific entry audit (values, anomalies): `docs/VCT_MATRIX_AUDIT.md` (audited against a `build_lto` tree copy — re-check MD5 before treating that path as claim).

---

## 1. Decision (authoritative)

| Role | MD5 (full) | Status |
|------|------------|--------|
| **JCIM / three-engine / C0 claim pin** | `72d7c7396702331d96ff12d18f831796` | **ACTIVE PIN** — use for all claim docks, admission, aggregation |
| **Repo / WRK / CMake build copies** | `9dc93717dfed0698006d88dd6a9627bc` | **STALE for claims** — experimental softcore / 6-entry edit (`git` commit `8c0c840ff`, 2026-06-12) |
| **PhD AtomTypes archive** | `204b75ef31b69e4a14deecf8a48c3f71` | **Semantic twin of 72d7** (820 pairwise values identical); **not** the claim pin (whitespace / fixed-width formatting differs on 4 lines). Dated mtime 2015-05-11 |

**Do not mix MD5s within a campaign.** Admission requires exact byte pin `72d7…` (see `matrix_md5` contract).

---

## 2. Triple-validation (on-disk, this session)

Claim pin file (local durable pin, not iCloud walk):

`~/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin`

| Check | Path | MD5 | Match pin? |
|-------|------|-----|------------|
| Pin artifact | `/Users/lp.more/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin` | `72d7c7396702331d96ff12d18f831796` | **yes** (identity) |
| Campaign data | `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` | `72d7c7396702331d96ff12d18f831796` | **yes** (`cmp` OK) |
| Engine B data | `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/bin/B/data/MC_st0r5.2_6.dat` | `72d7c7396702331d96ff12d18f831796` | **yes** (`cmp` OK) |
| Engine C copy | `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/bin/C/MC_st0r5.2_6.dat` | `72d7c7396702331d96ff12d18f831796` | **yes** (`cmp` OK) |

Triple (pin × campaign data × bin mirrors) is **PASS** for claim arm work.

---

## 3. Full path → MD5 inventory (known local paths only)

No deep iCloud / CloudDocs tree walks. Paths verified with `md5 -q` on 2026-07-15.

### 3.1 Claim pin family — `72d7c7396702331d96ff12d18f831796` (size 18043)

| Path | MD5 | Notes |
|------|-----|--------|
| `/Users/lp.more/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin` | `72d7c739…` | Named claim pin |
| `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` | `72d7c739…` | **IMATRX target for A/B/B0** |
| `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/bin/B/data/MC_st0r5.2_6.dat` | `72d7c739…` | B DEPSPA mirror |
| `/Users/lp.more/flexaidds_results/three_engine_entropy_q1/bin/C/MC_st0r5.2_6.dat` | `72d7c739…` | C binary-dir copy |

### 3.2 Stale experimental family — `9dc93717dfed0698006d88dd6a9627bc` (size 18040)

| Path | MD5 | Notes |
|------|-----|--------|
| `/Users/lp.more/Projects/FlexAIDdS/MC_st0r5.2_6.dat` | `9dc93717…` | Git-tracked repo root; CMake copies this into builds |
| `/Users/lp.more/Projects/FlexAIDdS/WRK/MC_st0r5.2_6.dat` | `9dc93717…` | WRK; **byte-identical to repo root** |
| `/Users/lp.more/Projects/FlexAIDdS/build/MC_st0r5.2_6.dat` | `9dc93717…` | CMake copy from source |
| `/Users/lp.more/Projects/FlexAIDdS/build_lto/MC_st0r5.2_6.dat` | `9dc93717…` | CMake copy |
| `/Users/lp.more/Projects/FlexAIDdS/build_claim/MC_st0r5.2_6.dat` | `9dc93717…` | CMake copy |
| `/Users/lp.more/Projects/FlexAIDdS/build_metal_fix/MC_st0r5.2_6.dat` | `9dc93717…` | CMake copy |
| `/Users/lp.more/Projects/FlexAIDdS/build_v137/MC_st0r5.2_6.dat` | `9dc93717…` | CMake copy |
| `/Users/lp.more/Projects/FlexAIDdS/benchmarks/astex_repro/engine/MC_st0r5.2_6.dat` | `9dc93717…` | Astex repro engine pin (legacy handoff: 9dc9) |

### 3.3 Lab archive (semantic JCIM values, different bytes) — `204b75ef31b69e4a14deecf8a48c3f71` (size 18040)

| Path | MD5 | Notes |
|------|-----|--------|
| `/Users/lp.more/Documents/PhD/AtomTypes/MC_st0r5.2_6.dat` | `204b75ef…` | 2015-era archive; **numeric values match 72d7** |
| `/Users/lp.more/Library/Mobile Documents/com~apple~CloudDocs/Documents/PhD/AtomTypes/MC_st0r5.2_6.dat` | `204b75ef…` | Same file via CloudDocs mount (no deep find) |

### 3.4 Not found (checked)

| Location | Result |
|----------|--------|
| `.grok/skills/flexaidds/data/MC_st0r5.2_6.dat` | **Absent** (skill data has defs / other matrices only) |
| `$FLEXAIDDS_QUEUE_ROOT` (iCloud queues path in `~/.flexaidds_env`) | **Not scanned** (CloudDocs policy); local campaign OUT uses `three_engine_entropy_q1/data` |
| `~/flexaidds_results/campaigns/three_engine/{A,B,B0}` | No matrix file; thin campaign OUT only |
| `bin/A/` under three_engine | Binary only; IMATRX via shared `data/` |

---

## 4. A / B / B0 work CONFIG → IMATRX (claim run)

Campaign work root:

`/Users/lp.more/flexaidds_results/three_engine_entropy_q1/work/{A,B,B0}/`

| Arm | CONFIG.inp with IMATRX | Unique IMATRX path | Resolved MD5 | Non-claim count |
|-----|------------------------|--------------------|--------------|-----------------|
| **A** | 88 | `…/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` | `72d7c739…` | **0** |
| **B** | 88 | same | `72d7c739…` | **0** |
| **B0** | 88 | same | `72d7c739…` | **0** |

Sample (all arms equivalent):

```text
IMATRX /Users/lp.more/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat
DEPSPA /Users/lp.more/flexaidds_results/three_engine_entropy_q1/data
```

**Live dock check (2026-07-15):** arm A pilot running FlexAID on `work/A/1GPK/restart_0/CONFIG.inp` — IMATRX already `72d7`.  
**Action taken:** **no overwrite** of live or work CONFIG files (already correct).

---

## 5. Scientific delta: 72d7 (claim) vs 9dc9 (repo/WRK)

Same 40-type triangular layout (820 keys). **7 entries differ** (repo = softcore-style experiment):

| Entry | Claim `72d7` | Repo/WRK `9dc9` |
|-------|--------------|-----------------|
| `[2-4]` | −149.40 | −86.6 |
| `[10-13]` | −86.64 | −15.0 |
| `[10-35]` | 0 | −175.0 |
| `[12-35]` | 0 | −110.0 |
| `[13-40]` | 33.99 | 90.0 |
| `[14-40]` | 43.24 | 90.0 |
| `[15-40]` | 29.56 | 90.0 |

`docs/VCT_MATRIX_AUDIT.md` quotes the **claim-family** values for these keys (e.g. `[2-4] = -149.40`). The 9dc9 set matches the June 2026 softcore / matrix-edit experiment, **not** the frozen JCIM claim pin.

### 72d7 vs AtomTypes `204b`

- Pairwise **values: 0 differences** (semantic JCIM archive).
- **Bytes differ** on 4 fixed-width lines only (e.g. `-149.4` vs `-149.40`, spacing on `13-40` / `14-40` / `15-40`).
- **Never substitute 204b for claim pin checks** — admission is MD5-exact.

---

## 6. Restore recommendation

| Target | Current | Recommendation | Safe now? |
|--------|---------|----------------|-----------|
| Live / work `CONFIG.inp` A·B·B0 | Point to claim `data/` @ 72d7 | **Leave alone** | N/A (already correct) |
| `three_engine_entropy_q1/data/MC_st0r5.2_6.dat` | 72d7 | **Keep** (claim source of truth for campaign) | yes |
| `pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin` | 72d7 | **Keep** | yes |
| Repo root + `WRK/` + `build*/` | 9dc9 | **Restore from pin for future claim-aligned local builds** when no process is reading those paths as IMATRX; prefer pointing docks at `three_engine…/data` rather than silent mid-run swaps | Restore of **WRK** is safe if unused as live IMATRX; **repo root** is git-tracked experimental matrix — restore only via intentional commit after campaign freeze, or leave 9dc9 as labeled experiment and never use it for claim IMATRX |
| `benchmarks/astex_repro/engine/` | 9dc9 | Treat as **legacy repro pin**; do not auto-merge into JCIM claim without re-baselining that campaign | defer |
| iCloud queue root | not scanned | Prefer local `three_engine_entropy_q1/data` + pin artifact; if queue stages matrix copies, materialize with `scripts/icloud_safe_io.py` and assert MD5 == 72d7 | no deep find |

### Restore actions taken this session

| Action | Result |
|--------|--------|
| Overwrite live dock CONFIG | **Not done** (already 72d7; arm A live) |
| Overwrite campaign `data/` matrix | **Not done** (already 72d7) |
| Overwrite repo / WRK / build matrices | **Not done** — deferred; recommendation only (9dc9 is intentional experimental tree; claim path is separate) |

**How to restore WRK (when safe, offline):**

```bash
PIN="$HOME/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin"
# verify
test "$(md5 -q "$PIN")" = "72d7c7396702331d96ff12d18f831796"
cp -p "$PIN" "$HOME/Projects/FlexAIDdS/WRK/MC_st0r5.2_6.dat"
# optional: repo root (creates dirty git tree — commit deliberately)
# cp -p "$PIN" "$HOME/Projects/FlexAIDdS/MC_st0r5.2_6.dat"
```

**How to verify any matrix:**

```bash
md5 -q path/to/MC_st0r5.2_6.dat
# must print: 72d7c7396702331d96ff12d18f831796  for claim eligibility
```

---

## 7. Operational rules (claim campaigns)

1. **IMATRX** must resolve to bytes MD5 `72d7c7396702331d96ff12d18f831796`.
2. Prefer absolute path under local APFS:  
   `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat`  
   (default local root `~/flexaidds_results`).
3. Never point claim docks at repo root / `build_lto` matrix while those remain 9dc9.
4. Do not “fix mid-run” by swapping the file under a live `IMATRX` path; stage new trees instead.
5. Aggregators:  
   `python3 scripts/aggregate_claim_metrics.py <dir> --matrix-md5 72d7c7396702331d96ff12d18f831796`

---

## 8. Audit command log (reproducible)

```bash
# pin + campaign data
md5 -q ~/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin
md5 -q ~/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat
cmp ~/flexaidds_results/pins/MC_st0r5.2_6.dat.JCIM2015_claim_pin \
    ~/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat

# repo / WRK
md5 -q /Users/lp.more/Projects/FlexAIDdS/MC_st0r5.2_6.dat
md5 -q /Users/lp.more/Projects/FlexAIDdS/WRK/MC_st0r5.2_6.dat

# AtomTypes archive (formatting twin)
md5 -q ~/Documents/PhD/AtomTypes/MC_st0r5.2_6.dat
```

---

*Last updated: 2026-07-15T12:05:12Z — matrix triple-validate / restore pin decision for MC_st0r5.2_6.dat.*
