# Comparative science — hub

**Entry point** for FlexAID / FlexAIDdS comparative methodology, the pre-merge genuine baseline, the P0–P5 pipeline, and Wave 3 sampling work.

| Field | Value |
|-------|--------|
| **Status** | Operator hub (docs only). Not a published success rate — unverified / no METHODOLOGY.md §0 receipt |
| **Last science snapshot** | 2026-07-25 |
| **Branch for pipeline work** | `feat/comparative-p0-p5-pipeline` |
| **Machine pins** | [`arm_pins.json`](arm_pins.json) |
| **Methodology parents** | `METHODOLOGY.md`, [`docs/ICLOUD_BENCHMARK_STORAGE.md`](../ICLOUD_BENCHMARK_STORAGE.md) |

If you only open one file, open this one. Deep specs live in the linked docs below — do not duplicate numbers elsewhere without linking back here or to the baseline record.

---

## 1. Goal (three arms)

Measure docking success on the **same** Astex Diverse native redock task under **frozen fairness axes**, varying only the intended science identity of each arm:

| Arm | Science identity | TEMPER | Clustering | Source pin (short) |
|-----|------------------|--------|------------|--------------------|
| **A** | JCIM 2015-era **CF-only** FlexAID | 0 | CF | FlexAID `f766a14` |
| **B** | **First entropy** FlexAID (soft-β / FO lineage) | 21 | FO | FlexAID `1a6ae0b` |
| **C** | Current **FlexAIDdS** (same fairness axes) | 21 | FO | FlexAIDdS build commit at binary time |

**Fairness axes (identical across arms):** N=85 Astex native · 10 sims · budget 1000×2000 · matrix MD5 **9dc9** · seed-off · headline **S_top10** (10k bootstrap median) · serial host · local-first I/O.

Full design: [`COMPARATIVE_GOAL_METHODOLOGY.md`](COMPARATIVE_GOAL_METHODOLOGY.md) (G1–G9).  
Arm tables + protocol: [`COMPARATIVE_BENCHMARK_METHODOLOGY.md`](COMPARATIVE_BENCHMARK_METHODOLOGY.md).  
Pins: [`arm_pins.json`](arm_pins.json).

---

## 2. Where we are (honest snapshot)

### 2.1 Pre-merge genuine baseline (OPS session record — not publishable)

Campaign **`v_autonomous_20260724_160919`** — recorded in [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md). **Unverified / no receipt. Not a current docking-power rate.**

| Metric | Value | Role |
|--------|------:|------|
| **Genuine top-1 ≤2 Å** | **20 / 79 = 25.3%** | OPS session record — not a published rate |
| **BCR** (best-cluster ≤2 Å) | **22 / 79 = 27.8%** | Sampling ceiling |
| Election gap (BCR − genuine) | **~2 targets** | Election wall essentially closed |
| Seed-echo | **0** | Clean multi-target number |
| Denominator | **79 scored** (80/85 finished; 5 never scored) | Not fixed-85 unless labeled |

**Science conclusion (2026-07-25):** election / `free_energy_strict` gap is closed on this baseline. **Sampling is the bottleneck** (BCR ~28%). Route: **Wave 3** BCR raisers — see [`WAVE3_SAMPLING_BCR_PLAN.md`](WAVE3_SAMPLING_BCR_PLAN.md). Softβ S1 reorders heads; it cannot invent near-natives when BCR = 0.

This 25.3% figure is an **OPS session record**, not a published rate — **not** JCIM top-10, **not** 3Dsig S_top10, **not** a substitute for the three-arm comparative table.

### 2.2 Post-merge re-aggregation

[`BASELINE_GENUINE_POSTMERGE_ATTEMPT.md`](BASELINE_GENUINE_POSTMERGE_ATTEMPT.md) is **PROCEDURE ONLY** (`data_missing=true`). Local `$FLEXAIDDS_LOCAL_ROOT` had **0** `result.csv` hits; bulk trees moved to iCloud archive. **Do not invent** post-merge N / genuine % / BCR.

### 2.3 Comparative P0–P5 pipeline

Code: `python/flexaidds/comparative_phases/` · CLIs below · **17** unit tests in `python/tests/test_comparative_phases.py`.

| Phase | Status | Note |
|-------|--------|------|
| P0–P1 | Pass (layout + reconstruction receipts if bins empty) | Matrix **9dc9** |
| **P2** | **HOLD — live blocker** | Needs real `native_cf_oracle_gate` JSON (`ok` / `exit_code` / `ranking_forbidden`), not empty/deferred |
| P3–P4 | Pending behind P2 | Fail-closed serial |
| P5 | Scaffolding OK | Needs arm `result.csv` for real S_top10 |
| Unit tests | **17/17** | `test_comparative_phases.py` |

Live campaign ops: [`CAMPAIGN_STATUS_2026-07-25.md`](CAMPAIGN_STATUS_2026-07-25.md). **No live science dock** until pre-gates pass. Binaries A/B still source-pinned, binary missing.

### 2.4 Matrix pin

Claim / comparative / baseline-era matrix:

```text
MC_st0r5.2_6.dat
MD5 9dc93717dfed0698006d88dd6a9627bc   # "9dc9"
```

**Not** the **72d7** packing-sweetened fork. Live default path:  
`$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat`

---

## 3. Doc map

| Doc | What it is |
|-----|------------|
| **This hub** | Entry point, snapshot, quickstart, metric rules |
| [`COMPARATIVE_GOAL_METHODOLOGY.md`](COMPARATIVE_GOAL_METHODOLOGY.md) | G1–G9 acceptance, layer model, fairness axes |
| [`COMPARATIVE_BENCHMARK_METHODOLOGY.md`](COMPARATIVE_BENCHMARK_METHODOLOGY.md) | Arm A/B/C protocol detail, frozen axes |
| [`arm_pins.json`](arm_pins.json) | Machine-readable commits, matrix, paths |
| [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md) | OPS session 25.3% / 27.8% / seed_echo=0 — **not publishable** |
| [`BASELINE_GENUINE_POSTMERGE_ATTEMPT.md`](BASELINE_GENUINE_POSTMERGE_ATTEMPT.md) | Procedure-only re-agg (no numbers yet) |
| [`CAMPAIGN_STATUS_2026-07-25.md`](CAMPAIGN_STATUS_2026-07-25.md) | Live ops: bins, iCloud, next steps, blockers |
| [`FORWARD_SUCCESS_RATE_PLAN.md`](FORWARD_SUCCESS_RATE_PLAN.md) | Normative KEEP/DEFER/REJECT sequencing |
| [`WAVE3_SAMPLING_BCR_PLAN.md`](WAVE3_SAMPLING_BCR_PLAN.md) | Sampling / BCR implementation plan |
| [`softbeta_election_policy.md`](softbeta_election_policy.md) | Softβ S1 default OFF; election vs sampling |
| [`3dsig_red_pair_protocol.md`](3dsig_red_pair_protocol.md) | Red-pair / 3Dsig family protocol |
| [`3dsig_shannon_ranking.md`](3dsig_shannon_ranking.md) | Shannon / soft-β ranking notes |
| `METHODOLOGY.md` (repo root) | Global build / parity / claim gates |
| `benchmarks/protocols/admission_metrics_contract.md` | Aggregator metric contract |
| `docs/ICLOUD_BENCHMARK_STORAGE.md` | Local-first vs thin iCloud mirror |

**Code / CLI (not docs):**

| Path | Role |
|------|------|
| `python/flexaidds/comparative_phases/` | P0–P5 gate modules |
| `scripts/run_comparative_phases.py` | Full pipeline driver (`--pipeline-dry`) |
| `scripts/comparative_phase_gate.py` | Per-phase gate (`--dry-run`) |
| `scripts/aggregate_claim_metrics.py` | S1 / S3 / STRICT aggregation |
| `scripts/bootstrap_3dsig_s_top10.py` | Comparative S_top10 bootstrap |
| `scripts/wave3_preflight.sh` | Matrix + seed-off preflight (no dock) |
| `python/tests/test_comparative_phases.py` | 17 unit tests |

---

## 4. Operator quickstart

### Env

```bash
export FLEXAIDDS_ROOT="$(git rev-parse --show-toplevel)"
cd "$FLEXAIDDS_ROOT"
export PYTHONPATH="$FLEXAIDDS_ROOT/python"
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"

# Prefer conda envs when present:
#   python  → numpy for phase CLIs
#   cpp-python-core → pytest
```

Generic form (if conda science envs are installed):

```bash
# Phase dry-run (numpy):
~/.claude-science/conda/envs/python/bin/python scripts/run_comparative_phases.py --pipeline-dry
~/.claude-science/conda/envs/python/bin/python scripts/comparative_phase_gate.py --dry-run

# Unit tests (pytest):
PYTHONPATH=$PWD/python ~/.claude-science/conda/envs/cpp-python-core/bin/python -m pytest \
  python/tests/test_comparative_phases.py -q
```

Flag on the pipeline driver is **`--pipeline-dry`**, not `--dry-run`. Gate helper uses `--dry-run`.

### Wave 3 preflight (no dock)

```bash
bash scripts/wave3_preflight.sh
# Expect matrix MD5 9dc93717… and seed-off env echo
```

### Aggregate genuine / BCR (only after local CSVs exist)

```bash
# Fail closed if trees missing — do not invent rates.
CAMP="$FLEXAIDDS_LOCAL_ROOT/campaigns/<campaign_id>"
PIN=9dc93717dfed0698006d88dd6a9627bc

python3 scripts/aggregate_claim_metrics.py "$CAMP" \
  --matrix-md5 "$PIN" --headline s1 --diagnostic-only --json /tmp/s1_genuine.json

python3 scripts/aggregate_claim_metrics.py "$CAMP" \
  --matrix-md5 "$PIN" --headline s3 --diagnostic-only --json /tmp/s3_bcr.json
```

Full materialize procedure: [`BASELINE_GENUINE_POSTMERGE_ATTEMPT.md`](BASELINE_GENUINE_POSTMERGE_ATTEMPT.md).  
**Never** `find` / `Path.rglob` under `Mobile Documents/` / CloudDocs; use `scripts/icloud_safe_io.py` for any CloudDocs path.

### Published anchors (not our measured rates)

| Contract | Value | Label carefully as |
|----------|------:|--------------------|
| JCIM 2015 Astex native FLRP **top-1** | 45.2% | Published top-1 |
| JCIM 2015 Astex native FLRP **top-10** | 66.7% | Published top-10 |
| 3Dsig 2017 S_top10-style medians | ~0.66 / ~0.69 | Bootstrap median family |

---

## 5. Metric label rules (fail-closed language)

| Label | Definition | Use as |
|-------|------------|--------|
| **Genuine / S1** | Rank-0 ordered RMSD ≤ 2.0 Å **and** `seed_echo=0` | Modern product KPI; OPS 25.3% is not a published rate |
| **Top-1** | Elected first pose ≤ 2.0 Å (publish/JCIM style) | Published anchors; label mode |
| **Top-10** | Any of ranks 0..9 ≤ 2.0 Å (case success) | JCIM Table 2 style |
| **S_top10** | Case success (top-10 style) → **median** of 10k bootstrap resamples | **Comparative headline** (3Dsig family) |
| **BCR** | Best cluster-head (or pool ceiling) RMSD ≤ 2.0 Å | **Diagnostic** sampling ceiling only |
| **STRICT / claim_ready** | Aggregator default when PB/tENCoM columns present | Stricter than baseline “genuine” |

**Rules:**

1. Never mix bare percentages across contracts without the label.  
2. Never count seed-echo poses as success.  
3. Prefer `rmsd_to_crystal` (ordered); never Hungarian for genuine claim rates.  
4. Report **denominator** (scored N vs fixed-85) on every table.  
5. Label **mode**: autonomous blind vs defined-cleft FLRP.  
6. Matrix MD5 on every claim table (**9dc9** for this program).  
7. CF / soft-β on CF are **scoring proxies**, not thermodynamic ΔG.

---

## 6. Explicit non-goals

- Matching published 45.2% / 66.7% / 0.66 as a “done” gate for methodology.  
- Claiming post-merge rates while `data_missing=true`.  
- Softβ DatasetRunner S1 as a substitute for arm B entropy identity.  
- Re-pinning claim production matrix to **72d7**.  
- Dual full-85 launch on the ~18 GiB science host.  
- Walking CloudDocs with `find` / `rglob`.  
- Physical ΔG / vib-Shannon as default docking ranker without separate labels.  
- Inventing binary SHA256 for missing A/B builds.

---

## 7. Recommended reading order

1. **This hub** (you are here)  
2. Baseline numbers → [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md)  
3. Goal design → [`COMPARATIVE_GOAL_METHODOLOGY.md`](COMPARATIVE_GOAL_METHODOLOGY.md)  
4. Live blockers → [`CAMPAIGN_STATUS_2026-07-25.md`](CAMPAIGN_STATUS_2026-07-25.md)  
5. Rate-raise sequence → [`FORWARD_SUCCESS_RATE_PLAN.md`](FORWARD_SUCCESS_RATE_PLAN.md) §0–§2  
6. Next engineering wave → [`WAVE3_SAMPLING_BCR_PLAN.md`](WAVE3_SAMPLING_BCR_PLAN.md)  
7. When trees return → [`BASELINE_GENUINE_POSTMERGE_ATTEMPT.md`](BASELINE_GENUINE_POSTMERGE_ATTEMPT.md)

---

## 8. Single next scientific steps (pointers only)

| Track | Next step | Detail doc |
|-------|-----------|------------|
| **Comparative A/B/C** | Build pin A+B binaries; real P2 oracle JSON; pilot8 serial | Campaign status §6 |
| **Genuine rates** | Materialize one local campaign tree → re-aggregate S1/S3 | Post-merge attempt |
| **Product KPI** | Wave 3 sampling / BCR raisers (flag-gated) | Wave 3 plan |

**Forbidden as next step:** dual full-85, Softβ-as-sampling-fix, inventing post-merge %, quoting incomplete trees as full-85.
