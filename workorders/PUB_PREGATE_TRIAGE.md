# Publication pre-gate triage (P2 / A·B·C binaries)

**Date:** 2026-07-28  
**After:** S4 PHENOTYPE_UNIQUE closed **PASS_LIVENESS** (null magnitude)  
**Local root:** `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1` → `~/flexaidds_results/three_engine_entropy_q1`  
**Explicit:** **full-85 is NOT authorized** by this memo.

---

## 0. Prerequisite completed: S4 close-out

| Item | Result |
|------|--------|
| OUT | `s4_pheno_unique_near_miss_20260727_211213` |
| Status | **PASS_LIVENESS** · ACCEPT_S4_PHENO=**False** |
| mean_dBCR | **−0.057** (floor −0.5 not met) |
| L4 | control 0 / tx 4 `[NEW-SEARCH-ARCH]` |
| Pins | **PINS_OK** (shared SHA `afd5cf42…`) |
| Sol #9 | **Released** · `may_dock=true` |
| Workorder | `S4_PHENOTYPE_UNIQUE_POSTERIORI.md` |

S4 does **not** unlock full-85 or comparative N=85 science dock.

---

## 1. Matrix pin

| Check | Value |
|-------|--------|
| File | `~/flexaidds_results/three_engine_entropy_q1/data/MC_st0r5.2_6.dat` |
| MD5 | **`9dc93717dfed0698006d88dd6a9627bc`** |
| Claim pin | **OK** |

---

## 2. Binary staging (A / B / C)

Paths under `~/flexaidds_results/three_engine_entropy_q1/bin/{A,B,C}/`.

| Arm | Science identity (pin) | Staged path | SHA256 (staged) | Status |
|-----|------------------------|-------------|-----------------|--------|
| **A** | JCIM-era CF-only FlexAID @ `f766a14e…` | `bin/A/FlexAID` | `62bdab8da7ef3250…` | **STAGED_RECONSTRUCTION** — **not** historical SHA |
| **B** | First entropy FlexAID @ `1a6ae0b0…` | `bin/B/FlexAID` | `c37d169f81146623…` | **STAGED_RECONSTRUCTION** — **not** historical SHA |
| **C** | FlexAIDdS current (worktree build) | `bin/C/FlexAIDdS` | `afd5cf42d8cb726d…` | **STAGED** (matches S4 pilot binary) |

### Identity labels (on disk)

- `bin/A/IDENTITY.txt` — reconstruction label; historical pin commit recorded  
- `bin/B/IDENTITY.txt` — same  
- `bin/C/SHA256.txt` — current FlexAIDdS  

### What this means for comparative claims

| Claim language | Allowed? |
|----------------|----------|
| “Arm C = this FlexAIDdS binary SHA” | Yes, for C |
| “Arm A = JCIM 2015 binary f766a14” | **No** until rebuild from that commit and re-pin SHA |
| “Arm B = entropy tip 1a6ae0b” | **No** until rebuild from that commit and re-pin SHA |
| “CF reconstruction / modern FlexAID labeled TEMPER0” | Yes if labeled **reconstruction** (arm_pins.json already allows this) |

**Remaining work for true historical A/B:** build FlexAID at pinned commits on this host (macOS; only Linux makefiles in `Projects/FlexAID/BIN` today — may need cmake path from FlexAID README or Linux builder), then replace staged files + update `docs/implementation/arm_pins.json` `binary_sha256` fields.

---

## 3. P2 native-CF oracle

### Contract

`scripts/native_cf_oracle_gate.py`: PASS when `CF_native ≤ best_ga_cf + tol` (lower CF better).  
Missing/sentinel `CF_native` → **exit 3**, `ranking_forbidden=true`.

### What we ran

Targets from S4 control + G4.3 control leaves (`1N1M`, `1L7F`):

| Leaf | CF_native | best_ga_cf | Gate exit | ranking_forbidden |
|------|-----------|------------|-----------|-------------------|
| S4 1N1M | **missing** (INI sentinel) | −99.314 | 3 | true |
| S4 1L7F | **missing** | −157.729 | 3 | true |
| G4.3 1N1M | **missing** | −99.314 | 3 | true |
| G4.3 1L7F | **missing** | −157.729 | 3 | true |

Aggregate:  
`~/flexaidds_results/three_engine_entropy_q1/oracle/oracle_status.json`  
→ **`status: hold`**, `ok: false`, `full85_authorized: false`.

### Pipeline dry

`run_comparative_phases.py --pipeline-dry --oracle-json …` reports phase state:

| Phase | State |
|-------|--------|
| P0 | pass (matrix 9dc9) |
| P1 | pass (reconstruction receipts labeled) |
| **P2** | **hold** |
| P3–P4 | pending behind P2 |
| P5 | scaffolding/dry only |

**next_allowed: P2** — comparative science dock still blocked.

### Why CF_native is missing

INI PDBs in these OUT trees do not carry a usable `REMARK CF=` for the crystal reference (sentinel / null).  
`result.csv` `cf_native` is 0.0000 on these leaves — not a valid competitive score.

**P2 clear requires (choose one path):**

1. **Score-native receipt:** run `probe_cf` (or engine path that writes CF on crystal pose) with production `--config` + matrix 9dc9; write `cf_native` into oracle inputs; re-run gate.  
2. **Dedicated canary work dirs** with valid INI REMARK CF (pilot8-style prep) and unseeded GA pose REMARK CF under the **same** prep/matrix.  
3. Do **not** `--force-p2-pass` for any claim comparative table.

---

## 4. Gate board (honest)

| Gate | Status | Blocks full-85 / N=85 comparative? |
|------|--------|-------------------------------------|
| Sol #9 free | **yes** | No (cleared) |
| Matrix 9dc9 | **pass** | No |
| Arm C binary | **staged** | No for C-only work |
| Arm A historical SHA | **fail** (reconstruction only) | **Yes** for JCIM-identity claim |
| Arm B historical SHA | **fail** (reconstruction only) | **Yes** for 3Dsig-identity claim |
| P2 native CF oracle | **hold** | **Yes** |
| Phase-4 sampling ACCEPT | **not met** (null stack + S4 liveness) | **Yes** for claim full-85 |
| S4 PHENOTYPE_UNIQUE | **PASS_LIVENESS only** | Does not clear above |

**full-85: NOT authorized.**

---

## 5. Ordered next actions (after this memo)

### 5a. Unblock P2 (priority)

1. Produce crystal-pose CF under production LOCCLF + 9dc9 for canary set (at least pilot8 / SEARCH-MISS panel).  
2. Ensure unseeded GA poses on same prep have REMARK CF.  
3. Re-run `native_cf_oracle_gate.py` → real PASS/FAIL (not missing).  
4. Feed `oracle_status.json` with `status=pass` only if contract satisfied.  
5. Re-run `run_comparative_phases.py --pipeline-dry` until P2 ≠ hold.

### 5b. Historical A/B binaries (if comparative identity required)

1. Build FlexAID @ `f766a14e` → install to `bin/A`, record SHA in arm_pins.json.  
2. Build FlexAID @ `1a6ae0b0` (entropy) → `bin/B`, record SHA.  
3. If unbuildable on macOS: keep **reconstruction** labels and never claim historical identity in the paper.

### 5c. Optional residual science (orthogonal)

- One-var **B** pilot: `FLEXAIDDS_BASIN_REINJECT=1` only (near-miss), after a priori.  
- Do **not** bundle with PHENOTYPE_UNIQUE for claim.

### 5d. Only after 5a (+ 5b if claiming historical A/B)

- Pilot8 fairness smoke under 9dc9 + seed-off  
- Then multi-seed N=85 S_top10 under admission contract  
- Report S1 / S2 / STRICT / `claim_ready` separately  

---

## 6. Artifacts written this triage

| Path | Role |
|------|------|
| `three_engine_entropy_q1/bin/{A,B,C}/` | Staged binaries + SHA/IDENTITY |
| `three_engine_entropy_q1/oracle/s4_panel_native_cf_oracle/*` | Per-target gate JSON |
| `three_engine_entropy_q1/oracle/oracle_status.json` | Aggregate P2 hold |
| `three_engine_entropy_q1/oracle/pipeline_dry.json` | Dry pipeline state |
| `workorders/S4_PHENOTYPE_UNIQUE_POSTERIORI.md` | S4 close |
| `workorders/PUB_PREGATE_TRIAGE.md` | This memo |

---

## 7. One-line status

**S4 A finished PASS_LIVENESS (null mag); lock free; C staged; A/B reconstruction-only; P2 HOLD (no CF_native) — full-85 still forbidden.**
