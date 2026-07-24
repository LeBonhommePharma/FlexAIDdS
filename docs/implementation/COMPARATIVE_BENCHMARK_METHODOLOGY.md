# Comparative benchmarking methodology — FlexAID 2015 vs first Shannon/soft-β entropy vs FlexAIDdS

**Status:** normative protocol for fair, objective three-arm comparison (methodology only; does not claim new dock rates).  
**Parents:** `METHODOLOGY.md` §0/§3; `docs/implementation/3dsig_red_pair_protocol.md` §§1–2; `docs/implementation/3dsig_shannon_ranking.md` §§1–2; `docs/implementation/softbeta_election_policy.md`; `docs/classic_entropy_ranking.md`.  
**Published anchors:** Gaudreault & Najmanovich, *J. Chem. Inf. Model.* **55**, 1323–1336 (2015), doi:10.1021/acs.jcim.5b00078; Morency, 3Dsig/ISMB-ECCB 2017 deck (`Morency_LP_3Dsig_2017.pdf`).  
**Environment:** single serial science owner on ~18 GiB hosts; never dual-launch full Astex campaigns.

---

## 1. Purpose

Reproduce a **comparative methodology** that places three scientific arms on **identical fairness axes**, so that measured differences are attributable only to the intended ranking/search variable:

| Arm | Science identity | What may change vs other arms |
|-----|------------------|-------------------------------|
| **A** | **FlexAID 2015 JCIM-era CF-only** | TEMPER 0, CF clustering/election (no entropy ranking) |
| **B** | **First Shannon / soft-β entropy concept** in the FlexAID lineage | TEMPER >0, density-mode (FO) clustering, soft-β \(\tilde G=\tilde H-T\tilde S\) / ACF on **CF a.u.** |
| **C** | **Current FlexAIDdS** | Modern engine packaging (DatasetRunner / FlexAIDdS stack) under the **same** fairness axes; entropy path must still be the soft-β CF objective unless explicitly labeled otherwise |

Arms must **not** be confounded packages (e.g. Softβ + `COM_BURIAL_CAP` + matrix **72d7** + autonomous cleft as a silent “entropy” arm).

---

## 2. Fairness axes (frozen — must match across arms)

These axes are **identical** for A, B, and C. Only the science variables in §3 differ.

| Axis | Frozen value | Source |
|------|--------------|--------|
| Primary dataset | Astex Diverse **N = 85** (native cognate-pocket redock story) | JCIM 2015 / 3Dsig / `3dsig_shannon_ranking.md` §2 |
| Sims / case | **10** independent simulations | 3Dsig deck; `FLEXAID_RESTARTS=10` / `FLEXAIDDS_RESTARTS=10` |
| Budget / sim | **2 000 000** energy evaluations | pop × gen = **1000 × 2000** (fixed gen; pop×DoF only if **documented**) |
| Matrix | `MC_st0r5.2_6.dat` MD5 **`9dc93717dfed0698006d88dd6a9627bc`** | JCIM production pin; **not** the 72d7 packing fork |
| Seed policy (claim-style) | **No native-pose seed** (`seed_elitism=0` / no INI crystal pose as claim success) | METHODOLOGY / red-pair |
| RMSD | Heavy-atom, symmetry-aware; success cutoff **≤ 2.0 Å** (inclusive) | JCIM / 3Dsig |
| RMSD engine | Prefer spyrmsd 0.9.0; Hungarian fallback only if logged | `METHODOLOGY.md` §0 |
| Headline comparative statistic | **Median** case success over **10 000** bootstrap resamples (with replacement) | 3Dsig |
| SEC / security channel | Off for accuracy benchmarks (`FLEXAIDDS_NO_SEC=1` when applicable) | `METHODOLOGY.md` §0 |
| Host concurrency | Serial science path on 18 GiB; workers × OMP ≤ P-cores | red-pair / AGENTS |

**Secondary datasets** (Astex non-native, HAP2) may be reported **with the same axes** but are not required for the primary three-arm table.

---

## 3. Three arms (protocol pin)

### 3.1 Arm A — FlexAID 2015 JCIM-era CF-only

| Field | Specification |
|-------|----------------|
| **Intent** | CF-only docking / ranking as in FlexAID 2015 comparative design |
| **TEMPER** | **0** |
| **Clustering** | `CLUSTA CF` |
| **Election** | Lowest CF / CF-path emission (no soft-β \(\tilde G\)) |
| **Binary preference** | Historical FlexAID pin under `three_engine_entropy_q1/bin/A/FlexAID` (SHA **must differ** from master B when claiming independent A) |
| **CONFIG generators** | `scripts/generate_flexaid_inp.py` arm `A` (`ARM_SPEC["A"]`) |
| **Launcher** | `scripts/run_flexaid_arm_pilot8.sh A` → `scripts/run_3dsig_red_pair_serial.sh --only A` |
| **Protocol-equivalent reconstruction** (if binary A missing) | Same matrix/budget/seed-off; CONFIG `TEMPER 0` + `CLUSTA CF`; current master binary as **B0-style CF control only** — label as **“CF reconstruction, not historical A SHA”** and do not claim byte-identical 2015 executable |

### 3.2 Arm B — first Shannon / soft-β entropy concept

| Field | Specification |
|-------|----------------|
| **Intent** | **First entropy-in-election** concept as in 3Dsig 2017 / FlexAID lineage: soft-β free energy on the **CF scoring proxy**, not ShannonThermoStack physical ΔG |
| **Ranking objective** | \(\displaystyle \tilde G=\tilde H-T\tilde S,\quad p_i\propto e^{-\mathrm{CF}_i/T},\quad \tilde G\equiv E_{\min}-T\ln Z_{\mathrm{local}}\) (cluster **ACF**) |
| **Implementation pin** | Single math: `LIB/SoftBetaFreeEnergy.h` (`flexaids::soft_beta::free_energy` / `acf`). Shared by `cluster.cpp` ACF, BindingMode classic F, and (only if flagged) DatasetRunner Softβ S1 |
| **First shared SoftBeta identity commit (repo)** | `ee7c3203` — *Add: SoftBeta ranking identity (G̃=H̃−T·S̃ ≡ ACF)*; DatasetRunner Shannon elect earlier at `c82e6fc2` |
| **TEMPER** | **21** (LP-optimized engine soft-T for FO/ACF; **not** physical \(k_B T\) kcal). Override `--temper 298` only with receipt note |
| **Clustering** | **`CLUSTA FO`**, **exactly one** literature MinPts pass (`fo_choose_minpts`) — **not** DP, not triple MinPts ladder |
| **Election path** | **Engine** ACF / soft-β emission when \(T>0\) (red-pair science path). DatasetRunner `FLEXAIDDS_SOFTBETA_ELECTION` is a **separate** S1 rescoring flag (default **OFF**) and is **not** required for arm B engine path |
| **Binary preference** | Master FlexAID / FlexAIDdS entropy-capable binary under `bin/B/FlexAID` (or staged equivalent) |
| **CONFIG** | `scripts/generate_flexaid_inp.py` arm `B` |
| **Forbidden confounds** | Do not package arm B as Softβ + COM_BURIAL_CAP + 72d7 + autonomous multi-knob “entropy” |

### 3.3 Arm C — current FlexAIDdS

| Field | Specification |
|-------|----------------|
| **Intent** | Current FlexAIDdS stack under the **same fairness axes** as A/B |
| **Binary preference** | Current `FlexAIDdS` / `benchmark_datasets` (e.g. local `~/flexaidds_results/baseline_engine/` when used for crash-proof C0) — record **binary_sha256** + build commit in receipt |
| **Entropy path for parity with B** | Soft-β / ACF on CF (engine TEMPER>0 + FO **or** explicit Softβ S1 with same \(T\)); **do not** report StatMech / ShannonThermoStack / tENCoM vib as the 3Dsig red-bar objective without a separate labeled column |
| **Modern extras (optional columns only)** | PoseBusters (S2), BCR, STRICT claim_ready, vib REMARKs — never replace **S_top10** headline |
| **Default production flags** | Softβ S1 DatasetRunner **OFF** unless opted in (`softbeta_election_policy.md`). Arm C comparative runs that intend 3Dsig entropy must **explicitly** enable the engine entropy path (TEMPER/FO) or Softβ S1 and log it |
| **Cleft / mode** | For JCIM/3Dsig-comparable fair redock prefer **defined-cleft / cognate pocket** consistent with A/B prep; if `autonomous` is used, label rates as **autonomous**, not JCIM-native FLRP |

### 3.4 Optional control B0 (not a fourth science arm)

`B0` = master binary, TEMPER 0, CLUSTA CF (`ARM_SPEC["B0"]`). Use only to separate **binary drift** from **entropy ranking**. If `bin/A` SHA equals master, B0 is a twin — not an independent historical control.

---

## 4. Metrics (never mix without labels)

| ID | Definition | Role |
|----|------------|------|
| **S_top10** | Any of the **top-10 ranked modes** (emitted rank order) has RMSD **≤ 2.0 Å** | **Primary comparative / 3Dsig claim contract** |
| **S1** | Rank-0 elected pose RMSD ≤ 2.0 Å | Modern stricter KPI (secondary) |
| **S2** | S1 ∧ PoseBusters | Modern secondary (not in 2017 deck) |
| **BCR** | Best cluster-head RMSD ≤ 2.0 Å | Sampling ceiling diagnostic only |
| **genuine** (audit) | S1 ∧ seed_echo=0 ∧ non-seed pose_source | Fail-closed claim hygiene |

**Headline:** for each arm, case-level S_top10 → **median** of **10 000** bootstrap resamples of the case set (`scripts/bootstrap_3dsig_s_top10.py`).

**Fail-closed:** bootstrap requires `mode_rmsd_0..9` (or documented rank table). Never treat `rmsd_bcr` / raw `success` as S_top10.

### 4.1 Published anchors — protocol labels only

JCIM 2015 (Gaudreault & Najmanovich, Table 2 — Astex native **FLRP**): report **each row with its own metric**. Do not swap labels.

| Anchor | Value | Protocol label (required when citing) |
|--------|-------|----------------------------------------|
| JCIM 2015 Astex native FLRP **top-1** (rank-1) | **45.2%** | **Top-1** success among published multi-run / ranking protocol in Table 2 — **not** top-10 |
| JCIM 2015 Astex native FLRP **top-10** | **66.7%** | **Top-10** success (Table 2) — **not** top-1 |
| 3Dsig 2017 Astex Diverse red bars | FlexAID **~0.66** / FlexAIDdS **~0.69** | **S_top10-style** median success over **10 000** bootstrap resamples of cases (deck); same family as top-10-ish rates, **not** JCIM top-1 45.2% |
| Morency poster (context) | FlexAID 66% / entropy 69% | Deck/poster S_top10-family; cite as such, not as JCIM top-1 |

**Forbidden:**
- Labeling **45.2%** as top-10 (it is **top-1** in Table 2).
- Labeling **66.7%** as top-1 (it is **top-10** in Table 2).
- Comparing modern **S1 single-run** claims to **JCIM top-10 (66.7%)** or **3Dsig S_top10 medians (~0.66)** without stating they are different contracts from **JCIM top-1 (45.2%)**.

---

## 5. Analysis pipeline (after arms complete)

```bash
# 1) Prepare + run arms serially (pilot8 first, then full85)
bash scripts/run_3dsig_red_pair_serial.sh --dry-run   # matrix + path smoke
bash scripts/run_3dsig_red_pair_serial.sh             # A → B0 → B when binaries present

# 2) Per arm, after mode_rmsd_* are filled:
python3 scripts/bootstrap_3dsig_s_top10.py \
  --arm-dir "$FLEXAIDDS_LOCAL_ROOT/campaigns/three_engine/A/3dsig_r10" \
  --bootstraps 10000 --json-out A_s_top10.json

# 3) Report table: arm × (S_top10 median [CI], S1 rate, BCR rate) — labels mandatory
```

Arm **C** (current FlexAIDdS) uses the same bootstrap entry point on a result tree that emits `mode_rmsd_0..9` in emitted rank order. If only rank-0 is available, report **S1 only** and mark S_top10 as **N/A (schema incomplete)** — do not impute.

---

## 6. Confound checklist

### 6.1 Must **match** across A / B / C

- [ ] Dataset N and target list (Astex Diverse 85 primary)
- [ ] pop × gen = 1000 × 2000 (2e6) per sim
- [ ] Restarts = 10 for claim rates (or labeled lower for pilot)
- [ ] Matrix MD5 `9dc93717…`
- [ ] No native-pose seed for claim-style rates
- [ ] RMSD engine + ≤2.0 Å cutoff
- [ ] Bootstrap 10k median for headline S_top10
- [ ] Same cleft/site definition class (cognate defined-cleft vs autonomous labeled)

### 6.2 Must **differ** only by science variable

| Arm | TEMPER | CLUSTA | Ranking objective |
|-----|--------|--------|-------------------|
| A | 0 | CF | min CF |
| B | 21 (default) | FO (single MinPts) | soft-β \(\tilde G\) / ACF on CF |
| C | documented (prefer match B for entropy parity) | documented | soft-β CF **or** CF-only control; never unlabeled multi-knob |

### 6.3 Explicit non-comparability (do not mix)

| Situation | Why invalid as fair arm comparison |
|-----------|-------------------------------------|
| 72d7 matrix vs 9dc9 | Different VCT landscape (packing fork) |
| Softβ + COM_BURIAL_CAP + autonomous as “entropy” | Multiple science variables |
| S1 / JCIM top-1 (45.2%) vs JCIM top-10 (66.7%) or 3Dsig S_top10 (~0.66) | Different success contracts (Table 2 + deck) |
| Softβ S1 DatasetRunner alone as arm B | Reorders heads only; **≠** engine TEMPER+FO sampling/election path unless documented |
| Physical StatMech / tENCoM ΔG vs soft-β \(\tilde G\) | Different physics claims |
| Seed-echo / INI RMSD as success | Violates claim seed-off |

---

## 7. Binary / commit pins (fill at run time)

See companion machine-readable map: campaign receipts + `{scratch}/arm_pins.json` template fields below.

| Arm | Preferred binary path | Commit / identity | Status if missing |
|-----|----------------------|-------------------|-------------------|
| A | `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/A/FlexAID` | Historical pin SHA (must ≠ B when claiming independent A) | **Protocol-equivalent reconstruction** (§3.1) |
| B | `…/bin/B/FlexAID` or entropy-capable master | Soft-β identity ≥ `ee7c3203` / FO single MinPts | Reconstruct TEMPER21+FO on current tree; label reconstruction |
| C | `FlexAIDdS` + `benchmark_datasets` (baseline_engine or `build/`) | `git rev-parse` of **binary build** + `binary_sha256` | Always available from current tree |

**Receipt fields (every arm):** `matrix_md5`, `binary_sha256`, `git_commit` (of binary), `pop`, `gen`, `restarts`, `seed_elitism`, `temper`, `clusta`, `mode` (defined-cleft vs autonomous).

---

## 8. Relation to existing three-engine red-pair

This document **extends** the two red bars (FlexAID vs FlexAID+entropy) into an explicit **three-arm** comparative table:

| This document | Red-pair map |
|---------------|--------------|
| Arm A | Red FlexAID / three-engine **A** |
| Arm B | Red FlexAIDdS entropy / three-engine **B** (TEMPER21+FO) |
| Arm C | Current FlexAIDdS packaging under same fairness axes (modern KPI columns allowed) |
| B0 | Optional binary-control (not a science arm) |

Primary launch path remains:

```bash
bash scripts/run_3dsig_red_pair_serial.sh   # A → B0 → B
```

Arm C is launched separately with the same matrix/budget/seed-off and must emit comparable `mode_rmsd_*` for S_top10.

---

## 9. Non-goals

- Claiming restored docking rates without running this protocol end-to-end.
- Full dual-campaign thrash on 18 GiB hosts.
- Equating modern physical free-energy ledgers with 2017 soft-β CF ranking without labels.
- Changing production Softβ S1 defaults on `main` as part of this methodology.

---

## 10. References

1. Gaudreault F, Najmanovich RJ. FlexAID: Revisiting Docking on Non-Native-Complex Structures. *J. Chem. Inf. Model.* 2015;55(7):1323–1336. doi:10.1021/acs.jcim.5b00078  
2. Morency L-P. The Impact of Conformational Entropy on the Accuracy of FlexAID in Binding Mode Prediction. 3Dsig / ISMB-ECCB 2017.  
3. In-repo: `METHODOLOGY.md`; `docs/implementation/3dsig_red_pair_protocol.md`; `docs/implementation/3dsig_shannon_ranking.md`; `scripts/bootstrap_3dsig_s_top10.py`; `scripts/generate_flexaid_inp.py` (`ARM_SPEC`).
