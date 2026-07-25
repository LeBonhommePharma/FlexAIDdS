# Campaign status — comparative FlexAID / FlexAIDdS (2026-07-25)

**Status:** **NO LIVE SCIENCE DOCK.** Comparative N=85 three-arm is **blocked** until pre-gates pass.  
**FlexAIDdS HEAD:** `99d17e4f` (`main`)  
**Goal design (phases G1–G9):** [`COMPARATIVE_GOAL_METHODOLOGY.md`](COMPARATIVE_GOAL_METHODOLOGY.md)  
**Arm specs:** [`COMPARATIVE_BENCHMARK_METHODOLOGY.md`](COMPARATIVE_BENCHMARK_METHODOLOGY.md) · pins [`arm_pins.json`](arm_pins.json)  
**Parents:** `3dsig_red_pair_protocol.md`, `METHODOLOGY.md`, Claude Science audit stack, `docs/ICLOUD_BENCHMARK_STORAGE.md`

---

## 1. Where the campaign is

| Track | State | Evidence |
|-------|--------|----------|
| **Live GA / FlexAID / FlexAIDdS dock** | **None** | No active claim/red-pair process at check time |
| **three_engine bin A / B** | **Empty** | `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/{A,B}/` have no `FlexAID` Mach-O |
| **Matrix pin (local)** | **OK** | `~/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` MD5 **`9dc93717dfed0698006d88dd6a9627bc`** |
| **C0 full-85 claim arms** | **SUSPENDED / archived** | Local stubs only (`MOVED_TO_ICLOUD`); `C0_SUSPENDED.md` — do not relaunch C0 alongside red-pair |
| **Pilot8 red-pair (historical)** | **Completed, SCIENCE FAIL** | `PILOT8_ANALYSIS_LATEST.md`: B0 and B both **S1 = S_top10 = BCR = 0/8** at ≤2 Å; matrix was **72d7** (wrong pin for claim); arm A cad-only / unusable |
| **Full-85 red-pair prep** | **Prepared, not citable** | Manifests under iCloud archive; SCIENCE_HOLD (B0 not control if A≡B; oracle gate for FO@298); prior scratch frozen S10=0 pathology |
| **resolve_build (FlexAIDdS)** | **Fail** | Pinned engine SHA missing / mismatch — rebuild before FlexAIDdS arm C |
| **Diagnostic (Science-grade)** | **Mechanism only** | `diagnostic/truth.ndjson` genuine_6=0; 1G9V ΔCF(false−native)≈−70 — CF prefers burial false min |

**Published anchors (not our measured rates):**

| Contract | Value |
|----------|--------|
| JCIM 2015 Astex native FLRP **top-1** | 45.2% |
| JCIM 2015 Astex native FLRP **top-10** | 66.7% |
| 3Dsig 2017 S_top10-style medians | FlexAID ~0.66 / FlexAID+entropy ~0.69 |

Do **not** mix top-1 with top-10 / S_top10 without labels.

---

## 2. Comparative arm pin (source of truth)

Repo: **`/Users/lp.more/Projects/FlexAID`** (`https://github.com/LeBonhommePharma/FlexAID.git`)

| Arm | Science identity | Source | Commit (full) | TEMPER | CLUSTA | Binary path (local live) |
|-----|------------------|--------|---------------|--------|--------|---------------------------|
| **A** | JCIM 2015-era **CF-only** | FlexAID `master` | **`f766a14e256c4b0ca45df77f28db2bfcad82a3b2`** (2015-12-16 TEMPER0 CF fix) | **0** | **CF** | `…/bin/A/FlexAID` |
| **B** | **First entropy** FlexAID (soft-β / FO lineage) | FlexAID `entropy` (= `origin/Entropy`) | **`1a6ae0b074084eadbaeee5c2c7973777a5cacf5e`** (2020-04-29 tip) | **21** | **FO** (single literature MinPts) | `…/bin/B/FlexAID` |
| **B0** | Binary control only | Same binary as **B** | (SHA of B) | 0 | CF | same as B |
| **C** | Current **FlexAIDdS** | FlexAIDdS `main` | build commit at binary time (doc snapshot `99d17e4f`) | **21** (entropy parity) | **FO** | `…/bin/C/FlexAIDdS` + `benchmark_datasets` |

**Notes**

- Master tip `9aa7995` is README-only (2026) — **not** a science binary pin.  
- Entropy birth: `f7e18d4` (2014-09-04 TEMPER); claim arm B uses **full entropy branch tip** `1a6ae0b`.  
- Alternate A code tip if rebuild needs later CF fix: `53771bd` (2016-01-05).  
- If historical A cannot build: label **“CF reconstruction”** on current `--legacy` TEMPER0 — never claim historical A SHA.  
- Softβ DatasetRunner S1 remains **OFF** for arm B identity (engine TEMPER+FO path).

Full JSON: `docs/implementation/arm_pins.json`.

---

## 3. Frozen fairness axes (all arms)

Identical across A / B / C:

- Dataset: Astex Diverse **N=85** (pilot gate: **pilot8** first)  
- Sims: **10** · Budget: **1000×2000 = 2e6** evals/sim  
- Matrix MD5: **`9dc93717dfed0698006d88dd6a9627bc`**  
- Seed: **off** (claim)  
- PSHARE: SHARESCL **10**, SHAREPEK **5**  
- Success headline: **S_top10** (any of ranks 0..9 ≤ 2.0 Å), median of **10k** bootstrap  
- Secondary: S1, BCR, (C only) PoseBusters / genuine  
- Host: **serial** one heavy arm; ~18 GiB — never dual full85  

---

## 4. iCloud Drive layout (load / save)

**Architecture:** local-first compute · iCloud = thin durable mirror (see `docs/ICLOUD_BENCHMARK_STORAGE.md`).

| Role | Variable | Default path |
|------|----------|--------------|
| Live root | `FLEXAIDDS_LOCAL_ROOT` | `~/flexaidds_results` |
| iCloud root | `FLEXAIDDS_ICLOUD` | `~/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks` |
| Thin results mirror | `FLEXAIDDS_RESULTS` | `$FLEXAIDDS_ICLOUD/results` (or legacy `…/CloudDocs/flexaidds/results`) |

```text
$FLEXAIDDS_LOCAL_ROOT/                          # LIVE (APFS only)
  three_engine_entropy_q1/
    bin/{A,B,C}/FlexAID|FlexAIDdS               # Mach-O — never sync as working FS
    data/MC_st0r5.2_6.dat                        # matrix pin
    inputs/                                      # prepared INP/GA
  campaigns/three_engine/{A,B0,B,C}/<campaign>/  # live GA OUT, poses, logs
  logs/
  pins/materialize/                              # CloudDocs materialize cache

$FLEXAIDDS_ICLOUD/                              # THIN MIRROR only
  results/campaigns/three_engine/{A,B0,B,C}/<campaign>/
    */result.csv
    RUN_RECEIPT*
    *oracle_status.json
    bootstrap/*.json
  pins/                                          # optional receipt copies of SHA/matrix
  archived_from_ssd/                             # cold archive of old campaigns
```

**Operator rules**

1. Export env: `export FLEXAIDDS_LOCAL_ROOT=$HOME/flexaidds_results`  
   `export FLEXAIDDS_ICLOUD="$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks"`  
2. `bash scripts/ensure_local_first_layout.sh`  
3. Live docks → local OUT only.  
4. Sync: `bash scripts/sync_three_engine_local_to_icloud.sh --campaign <id>` (result.csv + receipts).  
5. **Never** `find` / `Path.rglob` under `Mobile Documents/`. Hash via `python3 scripts/icloud_safe_io.py md5 <path>`.  
6. Load benchmarks for analysis: prefer local pin-cache; else materialize thin CSV only.

---

## 5. Why not N=85 three-arm yet (Science order)

Claude Science methodology: **mechanism / measurement / binary identity before claim-scale**.

| Blocker | Why it blocks comparative science |
|---------|-----------------------------------|
| **bin A / B missing** | Cannot attribute rates to JCIM vs first-entropy SHAs |
| **Pilot8 BCR=0 on both CF and FO** (prior run) | SCIENCE GATE FAIL for docking recovery — ranking cannot invent ≤2 Å poses; prior run also used **72d7** matrix |
| **1G9V CF false-min preference** (Δ≈−70) | CF landscape issue / prep / search — not Softβ-first |
| **RMSD truncation still in Python DatasetRunner** | Inflates/confounds modern C rates until fail-closed |
| **C0 dual-launch risk** | Forbidden with red-pair on this host |

---

## 6. Single next scientific step (explicit)

### **NEXT:** Build and pin FlexAID **A** + **B** binaries from the commits above; run a **2-target native CF oracle + pilot8 preflight** (not full 85).

**Do this now (ordered):**

1. **Build A** from FlexAID `@f766a14` → install  
   `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/A/FlexAID`  
   Record `shasum -a 256` into campaign receipt + update `arm_pins.json` `binary_sha256`.

2. **Build B** from FlexAID `@1a6ae0b` (branch `entropy`) →  
   `…/bin/B/FlexAID` · record SHA256 (must **≠** A).

3. **Confirm matrix** on live path MD5 `9dc93717…` (already OK locally).

4. **Mechanism gate (2 targets, e.g. 1P62 + 1T40):**  
   `bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40 …`  
   **Native CF oracle must not fail closed** (native CF ≫ decoy CF on healthy targets).  
   If oracle fails → **SCIENCE HOLD** — fix prep/ligand/cleft, **do not** launch N=85.

5. Only after oracle PASS: **serial pilot8** A → (optional B0) → B with seed-off, 1000×2000, R=10, matrix 9dc9.  
   Analyze S_top10 / S1 / BCR. Full 85 **only if** pilot8 shows non-zero BCR and schema `mode_rmsd_0..9`.

6. **Arm C** (FlexAIDdS) only after A/B pilot is interpretable; same fairness axes; rebuild engine so `resolve_build.py --check` passes.

**Forbidden as next step:** dual full-85, Softβ-as-sampling-fix, re-enabling C0 claim thrash, quoting pilot rates as deck 0.66/0.69.

---

## 7. Verification plan (this document)

| Check | Observation expected |
|-------|----------------------|
| `docs/implementation/arm_pins.json` exists | A=`f766a14…`, B=`1a6ae0b…`, matrix `9dc93717…` |
| Local matrix MD5 | equals pin |
| `bin/A` and `bin/B` | empty **until** build step done (`SOURCE_PINNED_BINARY_MISSING`) |
| No live dual campaign | no concurrent full85 owners |
| iCloud section | local live vs thin mirror paths defined |

---

## 8. References

- Gaudreault & Najmanovich, *JCIM* 2015, 55:1323–1336  
- Morency, 3Dsig 2017  
- `docs/implementation/COMPARATIVE_BENCHMARK_METHODOLOGY.md`  
- `docs/implementation/3dsig_red_pair_protocol.md`  
- `docs/ICLOUD_BENCHMARK_STORAGE.md`  
- Archive notes: pilot8 analysis (BCR=0), SCIENCE_HOLD, C0_SUSPENDED  
