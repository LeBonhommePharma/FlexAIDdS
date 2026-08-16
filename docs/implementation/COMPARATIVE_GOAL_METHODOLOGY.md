# Comparative goal methodology — how we fairly compare FlexAID (JCIM 2015), first entropy FlexAID, and FlexAIDdS

**Hub (start here):** [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md)  
**Status:** normative **goal-fulfillment design** (methodology only; does not claim new rates).  
**Audience:** every agent (Claude Science, Claude Code, Grok, Codex) and human operators.  
**Parents:** `METHODOLOGY.md` §0/§3; [`COMPARATIVE_BENCHMARK_METHODOLOGY.md`](COMPARATIVE_BENCHMARK_METHODOLOGY.md); `3dsig_red_pair_protocol.md`; `softbeta_election_policy.md`; `docs/ICLOUD_BENCHMARK_STORAGE.md`.  
**Machine pins:** [`arm_pins.json`](arm_pins.json).  
**Live status:** [`CAMPAIGN_STATUS_2026-07-25.md`](CAMPAIGN_STATUS_2026-07-25.md).  
**Genuine baseline (OPS session, not publishable):** [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md) (20/79 = 25.3%; unverified / no receipt).

This document answers: *What does “done” mean for the comparative goal, and what procedure guarantees the comparison is scientifically valid?*

---

## 0. Goal (one sentence)

**Measure docking success on the same Astex Diverse native redock task under frozen fairness axes, varying only the intended science identity of each arm** — so that differences between **JCIM-era FlexAID (CF-only)**, **first FlexAID entropy (soft-β / FO)**, and **current FlexAIDdS** are attributable to ranking/search architecture, not confounds.

### 0.1 What “fulfill the goal” means (acceptance of the *methodology*, not of a success rate)

The goal is **fulfilled** when all of the following are true:

| # | Criterion | Evidence artifact |
|---|-----------|-------------------|
| G1 | Three science arms run under **identical fairness axes** (§2) | Per-arm `RUN_RECEIPT` + `arm_pins.json` SHAs |
| G2 | Arm **A** binary is FlexAID **JCIM-era CF** pin (or labeled reconstruction) | `binary_sha256` + `git_commit` = `f766a14…` or `reconstruction=true` |
| G3 | Arm **B** binary is FlexAID **first-entropy** pin (or labeled reconstruction) | SHA ≠ A; commit = `1a6ae0b…` or `reconstruction=true` |
| G4 | Arm **C** is current FlexAIDdS with **same axes**; entropy path labeled | build commit + soft-β/FO receipt fields |
| G5 | Headline metric is **S_top10** (≤2.0 Å), **10k bootstrap median**, N=85 (or pilot with N labeled) | `bootstrap_3dsig_s_top10.py` JSON per arm |
| G6 | **No seed-echo / native INI** counted as claim success | `seed_elitism=0`, `seed_echo` flags |
| G7 | Matrix MD5 **`9dc93717dfed0698006d88dd6a9627bc`** on every arm | receipt `matrix_md5` |
| G8 | Results durable: local live OUT + **thin iCloud mirror** of CSV/receipts | sync script log |
| G9 | Report table states **anchors vs measured** without metric swap (top-1 ≠ top-10 ≠ S_top10) | final `COMPARATIVE_TABLE.md` |

**Non-goals of “fulfilled”:** matching published 45% / 66% / 0.66 numbers; physical ΔG claims; Softβ DatasetRunner S1 as a substitute for arm B.

---

## 1. Scientific layer model (must not conflate)

| Layer | Object | Units | Elects modes? | Role in this goal |
|-------|--------|-------|---------------|-------------------|
| **L1 Search** | GA fitness = Voronoi **CF** scoring proxy | a.u. | No | Exploration (all arms) |
| **L2 Cluster** | CF clusters **or** FO density modes | geometry + CF | Groups | A: CF · B/C: FO |
| **L3 Rank** | Soft-β \(\tilde G=\tilde H-T\tilde S\) on **CF** (engine ACF when T>0) | CF a.u., T score-temp | **Yes (B, C)** | Entropy arm identity |
| **L4 Thermo ledger** | StatMech / vib / ShannonThermoStack | kcal when calibrated | **No** (unless separately labeled) | **Out of band** for this goal |

**Language rule (Claude Science / AGENTS):** never call L1 or L3 “binding free energy ΔG.” Use “CF scoring proxy” and “soft-β ranking on CF.”

Identity for L3:

\[
\tilde G = E_{\min} - T \ln Z_{\mathrm{local}},\quad
p_i \propto e^{-\mathrm{CF}_i / T}
\]

with **local** \(Z\) inside each mode. \(T\) is a **score temperature** (arm B default **21**), not \(k_B T\) in kcal.

---

## 2. Fairness axes (frozen — identical for A, B, C)

If any axis differs, the arm pair is **not** a valid comparative cell (label as exploratory).

| Axis | Frozen value |
|------|----------------|
| Dataset | Astex Diverse **native** cognate-pocket redock, **N = 85** (primary) |
| Pilot gate | **pilot8** = `{1G9V,1GPK,1MEH,1P62,1Q4G,1R9O,1T40,2BYS}` before full 85 |
| Sims / case | **10** independent restarts |
| Budget / sim | **2 000 000** evals = pop **1000** × gen **2000** (fixed gen; DoF pop-scale only if **documented on all arms**) |
| Matrix | `MC_st0r5.2_6.dat` MD5 **`9dc93717dfed0698006d88dd6a9627bc`** |
| Niche (PSHARE) | SHARESCL **10**, SHAREPEK **5**, SHAREALF **4** |
| Seed (claim) | **Off** — no native pose seed / no seed_elitism |
| Site prep | **Defined-cleft / cognate pocket** class identical across arms (same apo, same sphere/cleft inputs) |
| Ligand | Same cognate ligand file + integrity gate (latm / atom count) |
| RMSD cutoff | Heavy-atom, symmetry-aware preferred; success **≤ 2.0 Å** inclusive |
| RMSD engine | spyrmsd 0.9.0 preferred; Hungarian only if logged |
| Headline | Case-level **S_top10** → **median** of **10 000** bootstrap resamples of cases |
| Host | **Serial** one heavy arm; no dual full-85; local-first I/O |

**Prep identity rule:** generate inputs **once** per target (apo, ligand, cleft/sphere, matrix path); arms differ only in **binary + TEMPER + CLUSTA + election path**. Same `SEED_BASE` for restart seeds across arms if restarts are seed-indexed.

---

## 3. Arm identities (what may differ)

| Arm | Science identity | Source pin | Binary install | TEMPER | CLUSTA | Election |
|-----|------------------|------------|----------------|--------|--------|----------|
| **A** | JCIM 2015-era **CF-only** | FlexAID `master` @ **`f766a14e256c4b0ca45df77f28db2bfcad82a3b2`** | `…/bin/A/FlexAID` | **0** | **CF** | min CF |
| **B** | **First entropy** FlexAID (3Dsig lineage) | FlexAID `entropy` @ **`1a6ae0b074084eadbaeee5c2c7973777a5cacf5e`** | `…/bin/B/FlexAID` | **21** | **FO** (single literature MinPts) | engine soft-β / ACF on CF |
| **B0** | Binary control (optional) | same Mach-O as B | same as B | 0 | CF | min CF |
| **C** | **Current FlexAIDdS** | FlexAIDdS build commit + SHA256 | `…/bin/C/FlexAIDdS` | **21** (parity with B) | **FO** | engine soft-β path; Softβ DatasetRunner S1 **OFF** unless labeled |

### 3.1 Reconstruction labels (when historical build fails)

| Situation | Allowed label | Forbidden claim |
|-----------|---------------|-----------------|
| A binary cannot build from `f766a14` | **“CF reconstruction”** — current `--legacy` TEMPER0 CF, same axes | “Historical JCIM binary” / “byte-identical 2015” |
| B binary cannot build from `1a6ae0b` | **“Entropy reconstruction”** — modern engine TEMPER21+FO, same axes | “First-entropy SHA” / “3Dsig binary replica” |
| A SHA == B SHA | Report **B vs B0 only** as entropy ON/OFF; A is not independent | “Three independent engines” |

### 3.2 Forbidden confounds (invalidate comparative cell)

- Matrix **72d7** packing fork vs **9dc9**
- Softβ + `COM_BURIAL_CAP` + autonomous multi-knob as “entropy”
- DatasetRunner Softβ S1 alone as arm B (reorders heads ≠ engine FO path)
- Seed-echo / INI RMSD as success
- Physical StatMech / tENCoM ΔG as the red-bar objective without a **separate** column
- Different pop×gen, restarts, or cleft class across arms
- Dual full-85 launch on ~18 GiB host
- Prefix-truncated RMSD (fail-closed only: sentinel −1)

---

## 4. Metrics hierarchy (report all; headline one)

### 4.1 Primary (goal contract)

| ID | Definition | Use |
|----|------------|-----|
| **S_top10** | Any of emitted ranks **0..9** has RMSD ≤ 2.0 Å | **Headline comparative** (3Dsig family) |
| **S_top10 median** | Median of 10 000 bootstrap resamples of case-level S_top10 | Published-style bar |

Requires CSV columns `mode_rmsd_0`…`mode_rmsd_9` (or documented rank table). **Fail-closed:** missing columns ⇒ S_top10 = N/A, not imputed from BCR.

### 4.2 Secondary (always log; never replace headline without label)

| ID | Definition | Use |
|----|------------|-----|
| **S1** | Rank-0 RMSD ≤ 2.0 Å | Modern top-1 KPI |
| **BCR** | Min cluster-head RMSD ≤ 2.0 Å | Sampling ceiling |
| **genuine** | S1 ∧ seed_echo=0 ∧ non-seed pose_source | Fail-closed claim hygiene |
| **S2** (C optional) | S1 ∧ PoseBusters pass | Modern package only |

### 4.3 Published anchors (protocol labels only — do not swap)

| Anchor | Value | Label when citing |
|--------|-------|-------------------|
| JCIM 2015 Astex native FLRP **top-1** | **45.2%** | Top-1 only — **not** S_top10 |
| JCIM 2015 Astex native FLRP **top-10** | **66.7%** | Top-10 only — **not** top-1 |
| 3Dsig 2017 red medians | FlexAID **~0.66** / entropy **~0.69** | S_top10-family bootstrap medians |

**Comparative interpretation rules**

1. Arm A S_top10 ↔ 3Dsig FlexAID ~0.66 and JCIM top-10 family (not JCIM top-1).  
2. Arm A S1 ↔ JCIM top-1 **only if** multi-run ranking protocol is documented as equivalent.  
3. Arm B S_top10 ↔ 3Dsig FlexAID+entropy ~0.69.  
4. Arm C vs B: packaging/modern stack under same L3 objective — not “new physics” unless L4 is separately validated.

---

## 5. Phased execution (ensures we *can* compare)

Science order: **mechanism → pilot → full N**. Skipping a gate invalidates later claim tables.

```text
Phase 0  Pins & storage
Phase 1  Build A + B (+ C later)
Phase 2  Mechanism / oracle (2–3 targets)
Phase 3  Pilot8 serial A → B0? → B  [then C]
Phase 4  Full85 serial (same order)
Phase 5  Bootstrap + comparative table + iCloud thin mirror
```

### Phase 0 — Pins and storage

```bash
export FLEXAIDDS_ROOT="$(git rev-parse --show-toplevel)"
export FLEXAIDDS_LOCAL_ROOT="${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}"
export FLEXAIDDS_ICLOUD="${FLEXAIDDS_ICLOUD:-$HOME/Library/Mobile Documents/com~apple~CloudDocs/FlexAIDdS_benchmarks}"
export FLEXAIDDS_RESULTS="${FLEXAIDDS_RESULTS:-$FLEXAIDDS_ICLOUD/results}"
bash scripts/ensure_local_first_layout.sh
# Matrix on live path
md5 -q "$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/data/MC_st0r5.2_6.dat"
# must equal 9dc93717dfed0698006d88dd6a9627bc
```

**Exit:** matrix pin OK; local dirs exist; env vars set.

### Phase 1 — Build source-pinned binaries

| Step | Action | Exit criterion |
|------|--------|----------------|
| 1.1 | Build FlexAID @ `f766a14` → `bin/A/FlexAID` | Mach-O exists; `shasum -a 256` recorded in receipt |
| 1.2 | Build FlexAID @ `1a6ae0b` (branch `entropy`) → `bin/B/FlexAID` | SHA **≠** A; receipt updated |
| 1.3 | Build FlexAIDdS → `bin/C/` when entering arm C | `resolve_build.py --check` or documented SHA |

If build fails → reconstruction path (§3.1) with explicit label; still record SHA.

### Phase 2 — Mechanism / native CF oracle (not ranking)

**Purpose:** prove prep + CF landscape can place near-native among searchable poses before spending pilot8×3 arms.

```bash
# Example (adjust paths): 2 discriminating targets
bash scripts/run_pilot8_canary_gates.sh --arm B0 --pdb 1P62,1T40 \
  --work-root "$FLEXAID_WORK_ROOT" --results-root "<LOCAL_OUT>"
```

| Outcome | Decision |
|---------|----------|
| Native CF competitive vs decoy / oracle panel **PASS** | Proceed Phase 3 |
| BCR = 0 and native CF ≫ elected false min (e.g. 1G9V-class ΔCF) | **SCIENCE HOLD** — fix prep/ligand/cleft/search; **do not** full85 |
| Oracle FAIL | Fail-closed; no Softβ experiments as “fix” |

### Phase 3 — Pilot8 (N=8) serial

```bash
export FLEXAID_POP=1000 FLEXAID_GEN=2000 FLEXAID_RESTARTS=10
export FLEXAID_CAMPAIGN="comparative_pilot8_$(date +%Y%m%d)"
# Serial only:
bash scripts/run_3dsig_red_pair_serial.sh --only A
# after A complete:
bash scripts/run_3dsig_red_pair_serial.sh --only B   # or --from B0 if B0 wanted
# C separately with same campaign fairness axes
```

| Exit | Proceed to full85? |
|------|---------------------|
| Schema `mode_rmsd_0..9` present; BCR or S_top10 **interpretable** (even if low) | Yes if not SCIENCE HOLD |
| BCR = 0/8 **and** S_top10 = 0/8 on **both** A and B under **9dc9** matrix | **HOLD** — treat as sampling/prep failure (see prior 72d7 pilot history); do not claim entropy ranking |
| Arm A incomplete / cad-only | Fix A binary; do not report A |

### Phase 4 — Full N=85 serial

Same knobs as Phase 3; campaign id distinct (`comparative_full85_*`).  
Order: **A → (optional B0) → B → C**. One arm at a time.

### Phase 5 — Analysis and durability

```bash
for ARM in A B C; do
  python3 scripts/bootstrap_3dsig_s_top10.py \
    --arm-dir "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/$ARM/${FLEXAID_CAMPAIGN}" \
    --bootstraps 10000 \
    --json-out "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/analysis/${ARM}_s_top10.json"
done
bash scripts/sync_three_engine_local_to_icloud.sh --campaign "$FLEXAID_CAMPAIGN"
```

Emit `COMPARATIVE_TABLE.md` with columns:

```text
arm | binary_sha | commit | temper | clusta | N | S_top10_median [CI] | S1 | BCR | matrix_md5 | reconstruction?
```

---

## 6. Storage design (compare later without data loss)

### 6.1 Architecture

| Layer | Where | Contents |
|-------|--------|----------|
| **Live** | `$FLEXAIDDS_LOCAL_ROOT` (default `~/flexaidds_results`) | GA poses, logs, Mach-O bins, matrix, full OUT trees |
| **Thin durable** | `$FLEXAIDDS_ICLOUD/results/…` | `result.csv`, RUN_RECEIPT, bootstrap JSON, oracle status |
| **Cold** | `$FLEXAIDDS_ICLOUD/archived_from_ssd/` | Old campaigns |

### 6.2 Layout

```text
$FLEXAIDDS_LOCAL_ROOT/
  three_engine_entropy_q1/
    bin/{A,B,C}/FlexAID|FlexAIDdS
    data/MC_st0r5.2_6.dat
    inputs/
  campaigns/three_engine/{A,B0,B,C}/<campaign_id>/
  campaigns/three_engine/analysis/
  logs/
  pins/materialize/

$FLEXAIDDS_ICLOUD/   # …/CloudDocs/FlexAIDdS_benchmarks
  results/campaigns/three_engine/{A,B0,B,C}/<campaign_id>/
  pins/
  archived_from_ssd/
```

### 6.3 Hard rules

1. **Never** live-write GA under `Mobile Documents/`.  
2. **Never** `find` / `Path.rglob` CloudDocs.  
3. Hash CloudDocs only via `scripts/icloud_safe_io.py`.  
4. Thin sync after each arm completes.  
5. Comparative analysis reads **local** CSV first; iCloud is backup.

---

## 7. Receipt schema (every arm run)

Minimum fields (JSON or RUN_RECEIPT):

```json
{
  "campaign_id": "comparative_pilot8_YYYYMMDD",
  "arm": "A|B0|B|C",
  "science_identity": "jcim_cf|first_entropy|flexaidds",
  "source_repo": "FlexAID|FlexAIDdS",
  "git_commit": "<full sha of binary source>",
  "binary_sha256": "<sha256 of Mach-O>",
  "reconstruction": false,
  "matrix_md5": "9dc93717dfed0698006d88dd6a9627bc",
  "pop": 1000,
  "gen": 2000,
  "restarts": 10,
  "seed_elitism": 0,
  "temper": 0,
  "clusta": "CF",
  "fo_minpts_policy": null,
  "site_class": "defined_cleft",
  "rmsd_engine": "spyrmsd|hungarian",
  "n_targets_attempted": 8,
  "n_targets_with_result_csv": 8
}
```

Without receipt fields G1–G4, the arm **cannot** enter the comparative table.

---

## 8. Decision tree (when numbers appear)

```text
                    ┌─ BCR≈0 all arms ──► sampling/prep problem; ranking A/B meaningless
S_top10 table ──────┤
                    ├─ BCR>0, S1≪BCR ──► election/ranking lever (entropy may help)
                    └─ S1≈BCR high ────► search+rank both working; compare to anchors
```

| Pattern | Scientific conclusion allowed |
|---------|-------------------------------|
| A ≈ published family, B > A on S_top10 | Entropy ranking helps under this protocol |
| A and B both ~0 under 9dc9 | **Cannot** claim entropy hurts/helps; fix mechanism first |
| C ≈ B, both high | FlexAIDdS packages first-entropy objective faithfully |
| C ≪ B under same axes | Packaging / election identity bug (FO sidecars, T mismatch, etc.) |
| C ≫ B multi-knob | Confounded — not attributable to “entropy” alone |

---

## 9. Relation to other documents

| Document | Role |
|----------|------|
| **This file** | Goal definition, phases, decision rules, fulfillment criteria |
| `COMPARATIVE_BENCHMARK_METHODOLOGY.md` | Arm specs, fairness table, confound checklist |
| `arm_pins.json` | Machine-readable commits/paths |
| `3dsig_red_pair_protocol.md` | Deck metric freeze, FO MinPts, launch order |
| `METHODOLOGY.md` | Global env invariants, determinism, CI accuracy gate |
| `ICLOUD_BENCHMARK_STORAGE.md` | Anti-hang storage ops |
| `CAMPAIGN_STATUS_*.md` | Point-in-time campaign state |

When documents conflict on **claim rates**, prefer **this file + METHODOLOGY.md**. When conflict on **red-pair knobs**, prefer `3dsig_red_pair_protocol.md` if this file is silent.

---

## 10. Operator checklist (print and tick)

- [ ] Env: `FLEXAIDDS_LOCAL_ROOT`, `FLEXAIDDS_ICLOUD`, matrix MD5 verified  
- [ ] Bin A from `f766a14` (or reconstruction labeled)  
- [ ] Bin B from `1a6ae0b` (or reconstruction labeled); SHA ≠ A  
- [ ] Phase 2 oracle not SCIENCE HOLD  
- [ ] Pilot8 serial complete; `mode_rmsd_*` present  
- [ ] Full85 only after pilot interpretable  
- [ ] Bootstrap 10k per arm  
- [ ] Thin iCloud sync  
- [ ] Comparative table with anchors correctly labeled  
- [ ] No dual full85 / no C0 thrash during red-pair  

---

## 11. Single next step (current machine state)

As of campaign status 2026-07-25: **bin A/B empty** → execute **Phase 1 then Phase 2** (build pins → 2-target oracle). Do **not** open Phase 4 (N=85) until Phase 3 exit criteria pass under matrix **9dc9**.

---

## 12. References

1. Gaudreault F, Najmanovich RJ. *J. Chem. Inf. Model.* 2015;55:1323–1336. doi:10.1021/acs.jcim.5b00078  
2. Morency L-P. 3Dsig / ISMB-ECCB 2017 (entropy ranking deck).  
3. In-repo: `METHODOLOGY.md`; `arm_pins.json`; `scripts/run_3dsig_red_pair_serial.sh`; `scripts/bootstrap_3dsig_s_top10.py`; `scripts/generate_flexaid_inp.py` (`ARM_SPEC`).  
