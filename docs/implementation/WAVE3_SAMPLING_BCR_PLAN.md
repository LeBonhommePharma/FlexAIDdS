# Wave 3 — Sampling / BCR-raiser implementation plan

**TL;DR:** On the pre-merge genuine baseline (**20/79 = 25.3%**, BCR **27.8%**, election gap ~2), scoring/election is no longer the wall — **sampling is**. Wave 3 raises **BCR** first (near-natives among heads), then genuine S1, via flag-gated search levers (BOOM/share, coarse-init, niche metric). Softβ S1 is not a sampling fix. Hub: [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md).

**Status:** concrete engineering plan (no multi-hour GA claim run required to accept this doc).  
**Date:** 2026-07-25  
**Parents:**  
[`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md) (hub) ·  
[`FORWARD_SUCCESS_RATE_PLAN.md`](FORWARD_SUCCESS_RATE_PLAN.md) §Wave 3 ·  
[`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md) ·  
[`COMPARATIVE_GOAL_METHODOLOGY.md`](COMPARATIVE_GOAL_METHODOLOGY.md) ·  
[`softbeta_election_policy.md`](softbeta_election_policy.md) · `METHODOLOGY.md`

**Does not invent rates.** All percentages below are documented baselines or published anchors.

---

## 0. Goal

Raise **BCR** (sampling ceiling: min cluster-head RMSD ≤ 2.0 Å), then **genuine** (S1 ∧ seed_echo=0), toward the **JCIM 2015 top-1 floor of 45.2%** as a *goal metric* — **without** reopening Softβ DatasetRunner S1 as the primary lever.

| Metric | Load-bearing baseline | Role in Wave 3 |
|--------|----------------------:|----------------|
| Genuine top-1 | **20/79 = 25.3%** | Product KPI after sampling lifts |
| BCR | **22/79 = 27.8%** | **Primary Wave 3 objective** |
| Election gap | **~2 targets** | Already closed enough — do not re-litigate Softβ |
| Seed-echo | **0** | Must stay 0 on claim paths |
| Campaign | `v_autonomous_20260724_160919` (pre-merge) | Reference only; re-baseline after binary pin |

**Science slogan (2026-07-25):** election / `free_energy_strict` wall is closed on this baseline; **sampling is the bottleneck** (BCR ceiling ~28%). Softβ reorders heads; it cannot create near-natives when BCR = 0.

**Published anchors (do not mix bare %):** JCIM top-1 **45.2%** · JCIM top-10 **66.7%** · 3Dsig S_top10-style ~0.66 / ~0.69.

**Claim matrix pin:** `MC_st0r5.2_6.dat` MD5 **`9dc93717dfed0698006d88dd6a9627bc`** (NOT 72d7 packing fork).

---

## 1. Why Wave 3 (not Softβ, not CAP, not dual full85)

```text
BCR floor (sampling)  ──W3──►  near-natives appear among heads
         │
         ▼
S1 / genuine (election already ~closed on baseline)
         │
         ▼
Wall / scoring (W2) gates memetic polish only
         │
         ▼
Comparative P0–P5 + claim full85 (after P2 oracle PASS)
```

| Layer | Can raise | Cannot |
|-------|-----------|--------|
| **Sampling (this wave)** | BCR | Invent ΔG |
| Softβ S1 / ACF election | Convert BCR>0 → S1 | Fix BCR=0 |
| COM_BURIAL_CAP | Mask burial free-lunch (UNCITABLE @ −130) | Substitute for wall redesign |
| Bigger pop alone | Marginal budget | Visit new basins if niche/BOOM collapses |

---

## 2. Prioritized KEEP levers (file paths + knobs)

Implement in order. Prefer **flag-gated / default-OFF or pilot-only JSON** for any behavioral change until serial pilot PASS.

### P0 — Preflight & measurement hygiene (no dock)

| Item | Path / command | Notes |
|------|----------------|-------|
| Matrix pin | `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` | MD5 must be **9dc9…** |
| Preflight | `scripts/wave3_preflight.sh` | Matrix + seed-off env echo only |
| Baseline doc | `docs/implementation/BASELINE_GENUINE_2026-07-24.md` | Denominator 79; pre-merge caveat |
| Metrics | `scripts/aggregate_claim_metrics.py`, `scripts/bootstrap_3dsig_s_top10.py` | BCR / genuine / S_top10 labeled |

### K1 — S1 search coverage: BOOM interval × sigma/share coupling  (**W3.1**)

**Role:** anti-collapse diversity so near-native basins are visited **more than once** (freq>1), not wiped or starved.

| Component | File | Current claim-path behavior |
|-----------|------|-----------------------------|
| Periodic BOOM injection | `LIB/gaboom.cpp` (~974–1016) | Replaces worst `(boom_inject_fraction × half)` chroms every `boom_inject_interval` gens with **fresh random** |
| BOOM knobs parse | `LIB/config_parser.cpp` (`ga.boom_inject_interval` default 100; `ga.boom_inject_fraction`) · `LIB/config_defaults.h` | Defaults interval=100, fraction=1.0 |
| Env → boom fraction | `FLEXAIDDS_BOOM_FRAC` → `ProtocolConfig::boom_frac` (`LIB/ProtocolConfig.cpp`) | Applied when set |
| **Claim DatasetRunner** | `LIB/DatasetRunner.cpp` (~6044–6058) | **Hardcodes `boom_inject_fraction: 0.0`** for no-seed modes — full 1.0 every 100 gens was catastrophic (wipes progress; CF≈0 stagnation) |
| Niche sigma (SIGMA_SHARE) | `LIB/gaboom.cpp` (~473–483): `sig_share` from gene bounds / peaks / scale | Logged as `SIGMA_SHARE=` |
| Share knobs | `ga.sharing_alpha` (env `FLEXAIDDS_SHARING_ALPHA`); `ga.sharing_peaks` / `ga.sharing_scale` (legacy SHAREPEK / SHARESCL in `gaboom.cpp` INI parse) | Fairness freeze: SHARESCL **10**, SHAREPEK **5**, alpha **4** baseline |
| Niche distance | `LIB/calc_rmsp.cpp` + `gaboom.cpp` pshare (~2741+) | Gene-space RMSP vs `sig_share` |

**Implementation sequence (code, still pilot-gated):**

1. **Do not** restore `boom_inject_fraction=1.0` on claim path.  
2. Add a **pilot-only** path (JSON and/or env) for a *mild* factorial, e.g.:
   - `boom_inject_interval` ∈ {100, 200, 400}
   - `boom_inject_fraction` ∈ {0.0 (control), 0.05, 0.10, 0.25} of worst **half** only  
3. Couple with sigma/share scale A/B: hold SHARESCL/SHAREPEK at fairness freeze first; optional `FLEXAIDDS_SHARING_ALPHA` 2.0 vs 4.0 only after boom control is stable.  
4. Acceptance (single-run diversity bar, restarts=1): near-native **freq>1** on tight failures (FORWARD cites 1J3J/1K3U); BCR improves vs control on discriminating set without SEC thrash.  
5. Receipt must log `boom_inject_interval`, `boom_inject_fraction`, `sharing_alpha`, `SIGMA_SHARE` (or peaks/scale).

### K2 — Coarse-init CF-rank seeds + orientation budget  (**W3.3**, S2)

**Role:** gen-0 starts with VCT-scored pocket placements, not all-clash randoms.

| Component | File | Knobs |
|-----------|------|-------|
| Coarse pocket scan | `LIB/coarse_init.cpp`, `LIB/coarse_init.h` | `n_orientations`, `n_seeds`, `grid_step` |
| Parse | `LIB/config_parser.cpp` (`coarse_init.enabled`, `grid_step`, `n_seeds`, `n_orientations` default **64**) | |
| Claim emission | `LIB/DatasetRunner.cpp` (~6000–6006) | **Always `enabled: true`**, `n_seeds: 25`, `n_orientations: 64`, grid_step from env/default |
| Grid step env | `FLEXAIDDS_COARSE_GRID_STEP` (`DatasetRunner.cpp` ~5867) | Pilot A/B step |
| Seed filter | `coarse_init.cpp`: keep CF **&lt; CLASH_THRESHOLD** (CF-rank, not CF&lt;0) | Regression: `tests/test_coarse_init_claim_default.py` |
| GA inject | `LIB/gaboom.cpp` populate / coarse seed block (~3072+) | Do not ring-randomise coarse seeds after score |

**Pilot A/B:** `n_orientations` **64 vs 256** on tight sites (FORWARD: 1OF1/1J3J/1K3U); cost ~linear — **pilot only**. Confirm `coarse_init.enabled=true` and seed_fraction=0 in dock JSON.

### K3 — Gene niche (E4) after boom pilot  (**W3.2**)

| Component | File | Plan |
|-----------|------|------|
| Gene-space niche | `LIB/calc_rmsp.cpp`, `gaboom.cpp` pshare | Basin starvation / BOOM wasted if gene0+angles dominate Å proximity |
| Cartesian niche | **not shipped as product default** | KEEP as search wave: niche distance in Å (or gene0-decoupled share) so multi-modal burial vs native can co-exist |

Gate: only after K1 pilot shows BOOM/share do not thrash SEC; default OFF flag until PASS.

### K4 — Wall before memetic interlock  (**W2.1 → W3.4**, not a Softβ substitute)

| Component | File | Env / flag |
|-----------|------|------------|
| Soft wall | `LIB/soft_wall.h` | Product soft-core / cap design surface |
| WAL knobs | `LIB/vcfunction.cpp` | `FLEXAIDDS_SOFTCORE_WAL`, `FLEXAIDDS_WAL_COERCIVE`, `FLEXAIDDS_WAL_STIFF` |
| COM hacks | — | `FLEXAIDDS_COM_BURIAL_CAP` **REJECT** as product default (see §3) |
| Memetic (E5) | FORWARD W3.4 | **Only if wall pilot PASS** |

**Ops interlock name:** treat **`WALL_PILOT_PASS`** as a *campaign gate label* (workorder / receipt field / science checklist), not a required engine env today. No memetic local-refine enablement without:

1. Frozen-pose / probe_cf wall A/B (W2.1): reduce com-over-burial invert on 1G9V-class **without** flipping clean native-as-CF-min probes.  
2. Explicit **PASS** recorded in campaign notes + receipt `protocol_config` extras if/when a flag is wired.  
3. Default OFF until that PASS.

If wall redesign is still OPEN, Wave 3 still runs **K1–K3 only** (BCR sampling); defer memetic.

### K5 — Already-landed hygiene (do not regress)

| ID | Status | Path |
|----|--------|------|
| Hard-clash severity (S3) | KEEP landed | clash path / CF gradient (see FORWARD) |
| Water O.3 typing (S4) | KEEP landed | typing lineage |
| free_energy_strict election | Baseline election gap ~closed | `LIB/SoftBetaFreeEnergy.h`, `LIB/cluster.cpp` |
| Softβ S1 product | **OFF** | `FLEXAIDDS_SOFTBETA_ELECTION=0` (see softbeta policy) |

---

## 3. Explicit REJECT list (do not reopen without new evidence)

| Item | Why |
|------|-----|
| **Bigger pop / more gens as sole lever** (FORWARD S5) | Basin visited once / not developed; budget not main bottleneck |
| **Softβ DatasetRunner S1 default ON** | Reorders heads only; BCR=0 ⇒ no S1; not a sampling fix |
| **COM_BURIAL_CAP as product default** (esp. −130 full85) | OOM, incomplete, single-target-tuned, UNCITABLE |
| **Dual full85** (or workers ≥6 on ~18 GiB) | OOM proven; serial science only |
| **72d7 matrix re-pin** | Packing fork; confounds claim |
| **Explicit population election default ON** | Elects over-burial basins (1G9V-class) |
| **Memetic before wall PASS** | Expected to deepen burial free-lunch |
| **CMA-ES primary as-is** (S6) | Invalid CF=10000 probes until gene-0 discrete decode |
| Mixing S1 single-run % with JCIM top-10 66.7% or 3Dsig 0.66 | Contract mismatch |

---

## 4. Pilot protocol (serial, seed-off, 9dc9)

### 4.1 Box rules

- Host: **one** science owner; workers **2–4** max; `OMP_NUM_THREADS=1` per worker.  
- **No dual full85.** No mmap rebuild of a live claim binary mid-run.  
- Budget when citing rates: **pop×gen = 1000×2000**, **R=10** (or R=1 for diversity-bar micro-probes, labeled).  
- Matrix MD5 **`9dc93717dfed0698006d88dd6a9627bc`**.  
- Seed: **off** — `FLEXAIDDS_SEED_ELITISM=0`, `FLEXAIDDS_NATIVE_SEED_FRAC=0` (or unset / empty reflig seed_fraction).  
- Softβ: **`FLEXAIDDS_SOFTBETA_ELECTION=0`**.  
- Local-first OUT under `$FLEXAIDDS_LOCAL_ROOT` (not CloudDocs live I/O).

### 4.2 Discriminating micro-set (2–3 targets)

Use for **mechanism + sampling factorial** before pilot8:

| PDB | Why |
|-----|-----|
| **1P62** | pilot8 / canary gate class; integrity + native-CF gate path (`scripts/run_pilot8_canary_gates.sh`) |
| **1T40** | same canary pair as comparative Phase 2 |
| **1G9V** (class) | burial / election pathology class; BCR vs elected gap historically informative — **sampling still required** if BCR high-Å |

Optional tight-site add-ons for orientation BOOM (not required for first gate): 1J3J, 1K3U (FORWARD W3.1/W3.3).

### 4.3 Pilot8 gate set

```text
{1G9V, 1GPK, 1MEH, 1P62, 1Q4G, 1R9O, 1T40, 2BYS}
```

Per `COMPARATIVE_GOAL_METHODOLOGY.md` fairness axes. Run **serial** after micro-set shows a BCR signal (or documents null with receipts).

### 4.4 Suggested factorial (micro-set, R=1 labeled OR R=10 claim-style)

| Arm | boom_inject_fraction | boom_inject_interval | n_orientations | Notes |
|-----|---------------------:|---------------------:|---------------:|-------|
| C0 | 0.0 | 100 | 64 | Claim-path control (DatasetRunner today) |
| B05 | 0.05 | 200 | 64 | Mild diversity |
| B10 | 0.10 | 200 | 64 | Mild+ |
| O256 | 0.0 | 100 | 256 | Orientation-only cost A/B |

One variable at a time when possible. Log binary SHA + git commit on every arm.

### 4.5 Pre-dock commands (no multi-hour GA)

```bash
export FLEXAIDDS_ROOT="$(git rev-parse --show-toplevel)"
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
bash scripts/wave3_preflight.sh
# Optional canary (no dock):
bash scripts/run_pilot8_canary_gates.sh --pdb 1P62,1T40 --report-only
# Comparative pipeline dry (does not dock):
PYTHONPATH=$PWD/python python3 scripts/run_comparative_phases.py --pipeline-dry
```

### 4.6 Launch shape (when operator authorizes live pilot)

Prefer existing claim/local-first launchers with **subset PDB list**, serial, seed-off — e.g. patterns in `scripts/run_C0_claim_clean.sh`, `scripts/run_pilot8_canary_gates.sh`, or DatasetRunner with explicit target list. **Do not** invent a second full85 campaign for Wave 3.

---

## 5. Success metrics & receipt fields

### 5.1 Metrics (report all; invent none)

| ID | Definition | Wave 3 use |
|----|------------|------------|
| **BCR** | Min cluster-head RMSD ≤ 2.0 Å | **Primary** — sampling success |
| **genuine** | rank-0 ≤ 2.0 Å ∧ seed_echo=0 | Secondary product KPI |
| **S1** | rank-0 ≤ 2.0 Å | Modern top-1 (may equal genuine if seed clean) |
| **S_top10** | any rank 0..9 ≤ 2.0 Å | Comparative headline when ranks present |
| **election gap** | BCR successes − genuine successes | Must not reopen large gap |
| **seed_echo** | contamination flag | Must stay **0** |

**Pass signals (pilot, not full85 promise):**

- Micro-set: BCR improves vs C0 control on ≥1 discriminating fail **or** near-native head frequency ↑ with stable SEC; no seed_echo.  
- Pilot8: BCR (and if heads exist, S1/genuine) non-regression on known goods; document nulls with receipts.  
- **Do not** claim ≥45.2% until full-85 claim protocol after comparative gates.

### 5.2 Receipt fields (every arm)

Write / verify in `RUN_RECEIPT.json` (`LIB/RunReceipt.{h,cpp}`) + campaign notes:

| Field | Required value / note |
|-------|------------------------|
| `matrix_md5` | `9dc93717dfed0698006d88dd6a9627bc` |
| `matrix_path` | live local path (not iCloud hash walk) |
| `binary_sha256` / `binary_path` | claim binary |
| `git_commit` | build commit |
| `pop` / `gen` / `restarts` | 1000 / 2000 / 10 (or labeled micro) |
| `seed_elitism` | **0** / false |
| `protocol_config` | full ProtocolConfig snapshot |
| Softβ env | `FLEXAIDDS_SOFTBETA_ELECTION=0` |
| Native seed | `FLEXAIDDS_NATIVE_SEED_FRAC=0` |
| Boom knobs | `boom_inject_interval`, `boom_inject_fraction` (JSON + receipt extras) |
| Coarse | `coarse_init.enabled`, `n_orientations`, `grid_step` |
| Share | `sharing_alpha` (+ peaks/scale if non-default) |
| Wall gate | `WALL_PILOT_PASS` yes/no/unknown (ops label) |
| Campaign | e.g. wave3 pilot id; baseline ref `v_autonomous_20260724_160919` |

---

## 6. Interlocks

| Gate | Meaning | Blocked work if FAIL |
|------|---------|----------------------|
| **Matrix 9dc9** | Preflight MD5 | Any claim-style pilot |
| **Seed-off** | SEED_ELITISM=0, NATIVE_SEED_FRAC=0 | Genuine claims |
| **Softβ OFF** | Primary lever remains sampling | Softβ default-ON merges as “fix” |
| **Election closed** | free_energy_strict path; gap ~2 on baseline | Reopening Softβ S1 as Wave 3 primary |
| **WALL_PILOT_PASS / wall-before-memetic** | W2.1 probe_cf / wall redesign PASS | **E5 memetic** only |
| **No dual full85** | Serial one heavy campaign | Parallel full85 |
| **Comparative P2 oracle** | Real `native_cf_oracle_gate` JSON PASS | Claim full85 / P3–P4 advance |

---

## 7. Interaction with comparative P0–P5

Wave 3 sampling work **feeds** arm C quality but **does not replace** the comparative pipeline.

| Phase | Comparative role | Wave 3 relation |
|-------|------------------|-----------------|
| **P0** layout + matrix | `scripts/comparative_p0_layout.sh` | Same matrix pin as Wave 3 |
| **P1** pin binaries A/B/C | reconstruction labels allowed | Rebuild C after Wave 3 code lands |
| **P2** native CF oracle | **Live blocker** (CAMPAIGN_STATUS): needs real gate JSON (`ok` / `exit_code` / `ranking_forbidden`) | **Must PASS before claim full85**; Wave 3 pilots on 1P62/1T40 can *feed* oracle evidence but P2 is fail-closed |
| **P3** pilot8 serial A→B→C | Fairness axes frozen | Run Wave 3 BCR pilots **on FlexAIDdS (C-class)** first; do not dual-launch full comparative full85 |
| **P4** full85 serial | Only after P2 + pilot interpretability | Wave 3 full85 genuine claim only after this gate stack |
| **P5** S_top10 bootstrap | Headline comparative | Report separately from genuine/BCR |

**Hard rule:** **P2 oracle must pass before claim full85.** Sampling improvements that raise BCR on micro-sets do not authorize a dual full85 or a Softβ-confounded “entropy” arm.

CLI dry checks:

```bash
PYTHONPATH=$PWD/python python3 scripts/run_comparative_phases.py --pipeline-dry
PYTHONPATH=$PWD/python python3 scripts/comparative_phase_gate.py --dry-run
```

---

## 8. Chunked delivery (implementation order)

| Chunk | Deliverable | Verification gate |
|-------|-------------|-------------------|
| **C0** | This plan + `scripts/wave3_preflight.sh` | Preflight exit 0 on matrix pin; no dock |
| **C1** | Mild BOOM fraction path (not 1.0) + receipt fields | Unit/config test; default claim path still fraction=0 unless pilot flag |
| **C2** | Coarse `n_orientations` pilot override (64/256) | JSON emission test; cost note |
| **C3** | Serial micro-set A/B (1P62, 1T40, ±1G9V) | BCR + seed_echo=0 receipts |
| **C4** | Gene niche / Cartesian share design (flag OFF) | Unit tests only until C3 PASS |
| **C5** | Memetic only if `WALL_PILOT_PASS` | Blocked until W2.1 |
| **C6** | Re-baseline genuine/BCR vs `v_autonomous_20260724_160919` after pin | Document rates; no invented uplift |

---

## 9. What “done” means for Wave 3 (this doc’s scope)

| Done now (doc wave) | Not done (requires authorized live GA) |
|---------------------|----------------------------------------|
| KEEP/REJECT inventory with paths | Multi-hour full85 dock |
| Pilot protocol + preflight script | Measured BCR ≥ JCIM 45.2% |
| Interlocks + comparative P2 rule | WALL_PILOT_PASS science experiment |
| Receipt field checklist | Production default change for BOOM fraction |

---

## 10. References (on disk)

- `docs/implementation/FORWARD_SUCCESS_RATE_PLAN.md`  
- `docs/implementation/BASELINE_GENUINE_2026-07-24.md`  
- `docs/implementation/COMPARATIVE_GOAL_METHODOLOGY.md`  
- `docs/implementation/CAMPAIGN_STATUS_2026-07-25.md`  
- `docs/implementation/softbeta_election_policy.md`  
- `docs/implementation/protocol-config.md`  
- `LIB/gaboom.cpp`, `LIB/coarse_init.cpp`, `LIB/DatasetRunner.cpp`, `LIB/soft_wall.h`, `LIB/SoftBetaFreeEnergy.h`, `LIB/calc_rmsp.cpp`  
- `scripts/wave3_preflight.sh`, `scripts/run_pilot8_canary_gates.sh`, `scripts/comparative_p0_layout.sh`  
- `tests/test_coarse_init_claim_default.py`, `tests/test_forward_success_rate_plan.py`
