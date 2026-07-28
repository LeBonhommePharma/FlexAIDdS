# S4 PHENOTYPE_UNIQUE near-miss — status report & next execution definition

**Report time:** 2026-07-27 ~21:55 EDT (live process check)  
**OUT:** `/Users/lp.more/flexaidds_results/s4_pheno_unique_near_miss_20260727_211213`  
**Branch tip (binary):** `9569df16` · **a priori commit:** `8f0d8f54`  
**Experiment:** one-var matched A/B — `FLEXAIDDS_PHENOTYPE_UNIQUE=1` vs control unset  

---

## 1. What “last step” was

After Phase-4 near-miss sampling nulls (G4.1 BOOM, ELECTION_V135, G4.3 MUTATION) and packaging residual S1–S5, residual flip order prioritized **new_search_arch**. Options **A+B** were implemented env-gated (`LIB/new_search_arch.h` + `gaboom.cpp`).  

**Last execution step:** launch **pilot 1 = option A only** (not A+B bundled):

| Axis | Value |
|------|--------|
| One variable | `FLEXAIDDS_PHENOTYPE_UNIQUE=1` on treatment only |
| Panel | NEAR_MISS `1N1M`, `1L7F` |
| Restarts | R=2 sequential |
| Budget | pop=1000; gens scaled by DoF (1L7F 8000, 1N1M 4000); CLI base 2000 |
| Matrix | **9dc9** (`md5 9dc93717dfed0698006d88dd6a9627bc`) |
| NO_SEC | 1 |
| Workers | 2 (OMP=1 each) |
| Binary SHA256 | `afd5cf42d8cb726de5b92fb66431095360f6161a088945eba110299cd09e4f57` |
| Sol #9 | Lock held `d20def7b-…` owner `grok-s4-pheno-unique` |
| L4 expected | `[NEW-SEARCH-ARCH] phenotype_unique=1` on **treatment** stderr only; control **0** |
| Magnitude floor | mean_dBCR ≤ −0.5 Å **or** ≥1 BCR&lt;2 **or** elect≤2.5; no wipeout |
| Forbidden | full-85, BOOM/election/MUT_GRAN/basin bundled |

**Not yet run:** option B (`BASIN_REINJECT`) as a separate matched pilot.

---

## 2. Live progress (IN_FLIGHT — not closed)

| Item | State |
|------|--------|
| Launcher | **Alive** (pid 25740, ~42 min wall at check) |
| Evaluator | **Alive**, polling; no `s4_pheno_posteriori.txt` yet |
| Phase | **Control arm only** — treatment `arm_pheno_unique` **not started** |
| ALL_ARMS_DONE | **No** |
| ACCEPT / posteriori | **Not written** (incomplete by design until both arms finish) |

### Control arm generation progress

| Target | Restart | Generations | Notes |
|--------|---------|-------------|--------|
| 1N1M | r0 | 3999 / 4000 done | **result.csv present** |
| 1N1M | r1 | 3999 / 4000 done | Folded into elect result |
| 1L7F | r0 | 7999 / 8000 done | TIMING SUMMARY present |
| 1L7F | r1 | **~6482 / 8000** | **Still running** (~98% CPU) |

### Partial control metrics (already on disk)

| Code | elect RMSD | BCR | wall_s | success_rmsd | pb_pass | claim_ready |
|------|-----------:|----:|-------:|-------------:|--------:|------------:|
| **1N1M** | **6.3999** | **4.1954** | 2477 | 0 | 1* | 0 |
| **1L7F** | — | — | — | — | — | **no result.csv yet** |

\* `pb_pass=1` in this row may reflect row flags / partial path; do **not** promote to STRICT claim without full receipt + bust_cli audit.

### Control L4 (required zero for treatment contrast)

| Target | `[NEW-SEARCH-ARCH]` | `[MUT-GRAN]` | `[BOOM]` | `[BASIN-REINJECT]` |
|--------|--------------------:|-------------:|---------:|-------------------:|
| 1L7F | **0** | 0 | 0 | 0 |
| 1N1M | **0** | 0 | 0 | 0 |

Control has **no** phenotype_unique marker — good so far for matched design.

### Consistency with prior near-miss controls (1N1M)

| Campaign | 1N1M elect | 1N1M BCR |
|----------|-----------:|---------:|
| G4.1 control | 6.3999 | 4.5515 |
| ELECTION control | 6.3999 | 4.0427 |
| G4.3 control | 6.3999 | **4.1954** |
| **S4 control (this run)** | **6.3999** | **4.1954** |

S4 control 1N1M elect matches the multi-campaign **6.3999 attractor**; BCR matches **G4.3 control** exactly on this partial snapshot — useful reproducibility signal, not yet a treatment effect.

### Ops health

| Check | Value |
|-------|--------|
| Free disk | ~40 GiB (≥ 20 floor) |
| may_dock (other experiments) | **false** (lock held) |
| Matrix pin | **9dc9** verified at stamp |
| Engines | 1L7F r1 only active at check |

**ETA (order-of-magnitude):**  
- Finish control 1L7F r1: ~3–8 min  
- Treatment arm (same budget R=2 × 2 targets): ~40–55 min  
- Total remaining: **~45–65 min** from report time if no stall  

---

## 3. Analysis (what we can and cannot conclude yet)

### Can conclude
1. **Launch is valid Sol #9:** lock, matrix 9dc9, stamped binary, workers≤2, local OUT.  
2. **Control path is live and progressing** at expected gen rates (~150 ms/gen 1L7F, ~220 ms/gen 1N1M).  
3. **Control L4 is clean** for NEW-SEARCH-ARCH (0 markers) — necessary for a fair treatment contrast.  
4. **1N1M control elect remains the known false-min attractor (6.40 Å)** — consistent with prior nulls; sampling still hard on this target under production LOCCLF.  
5. **Treatment has not started** — zero scientific conclusion on PHENOTYPE_UNIQUE efficacy.

### Cannot conclude
- Magnitude PASS/FAIL for S4 A  
- L4 pass on treatment (marker not yet observed)  
- Elect regression / improvement vs control  
- Any full-85 or claim_ready implication  

### Risks to watch through finish
| Risk | Mitigation |
|------|------------|
| Evaluator race / premature DONE | Script requires 2+2 valid result.csv and elect≥0 |
| Treatment binary ≠ control | Same stamped SHA for both arms (shared stamp) |
| L4 missing on treatment | Fail L4 if zero `[NEW-SEARCH-ARCH]` after arm complete |
| Confounding with MUT_GRAN / BOOM | Env forbids those on both arms |
| Interpreting 1N1M 6.40 as “new” | It is the established attractor |

---

## 4. How this fits the broader publication picture

Deep-research summary (workflow complete): campaign is **not** peer-review ready.  
This S4 pilot is a **serial near-miss residual**, not the three-arm N=85 claim path.

| Track | Role of this pilot |
|-------|--------------------|
| Phase-4 sampling ACCEPT | Unlikely to unlock alone; tests one architecture residual after null stack |
| Full-85 claim | **Still blocked**; a priori forbids full_85 |
| Comparative A/B/C N=85 | Unblocked only after **P2 + binaries**, not by S4 PASS |
| Methods paper | S4 is one more controlled negative-or-positive **architecture** datum for SI/Methods residual narrative |

---

## 5. Next execution definition (ordered, non-negotiable sequence)

Do **not** launch full-85. Do **not** open a second science dock while Sol #9 holds this OUT.

### Step 1 — Finish S4 (operational)

| ID | Action | Done when |
|----|--------|-----------|
| 1.1 | Let launcher complete control arm (1L7F r1 → aggregate result.csv) | `arm_control/{1L7F,1N1M}/result.csv` both exist, elect ≥ 0 |
| 1.2 | Let launcher run **treatment** `arm_pheno_unique` with env `FLEXAIDDS_PHENOTYPE_UNIQUE=1` only | driver shows `ARM pheno_unique DONE` |
| 1.3 | Wait for `ALL_ARMS_DONE` | line in `driver.log` |
| 1.4 | Do not kill engines unless hang > job-timeout or disk &lt; 20 GiB | — |

**Owner:** current session / caffeinate launcher already running.  
**No new OUT.**

### Step 2 — Score / a posteriori (scientific close of S4)

| ID | Action | Done when |
|----|--------|-----------|
| 2.1 | Confirm evaluator wrote `evidence/s4_pheno_posteriori.txt` + `evidence/accept.txt` | files non-empty |
| 2.2 | Manual cross-check table: elect/BCR control vs pheno for 1L7F/1N1M; mean_dBCR | matches result.csv |
| 2.3 | **L4:** count `[NEW-SEARCH-ARCH]` on treatment ≥1; control = 0 | fail closed if not |
| 2.4 | Apply floors: PASS / PASS_LIVENESS / FAIL / INVALID | status enum only |
| 2.5 | `python3 scripts/benchmark_self_eval.py validate-pins --out $OUT` | **PINS_OK** |
| 2.6 | Release Sol #9: `benchmark_coord.py release --token $TOKEN` | `may_dock=true` |
| 2.7 | Write `workorders/S4_PHENOTYPE_UNIQUE_POSTERIORI.md`; set a priori status CLOSED_* | committed |
| 2.8 | Update `a_posteriori_gate_ledger.md` + freeze pointer if needed | ledger row |

**Decision tree after 2.x:**

| Outcome | Next science residual |
|---------|----------------------|
| PASS magnitude + L4 | Optionally promote A into claim-recipe discussion; still no full-85 without sampling ACCEPT narrative update |
| PASS_LIVENESS / FAIL null | Document; **optional** separate pilot for **B only** (`BASIN_REINJECT=1`) — still one-var |
| INVALID | Fix instrument; re-a priori if re-run |

### Step 3 — Publication pre-gate triage (P2 / binaries) — **after** lock release

Only when Step 2.6 done (lock free). Still **not** full-85.

| ID | Gate | Concrete definition of done |
|----|------|------------------------------|
| **3.1 P2 native-CF oracle** | Clear HOLD on comparative pipeline | Produce real `native_cf_oracle_gate` JSON (or documented receipt) accepted by `scripts/run_comparative_phases.py` / pipeline dry path; P2 status **pass** not HOLD in campaign status |
| **3.2 Stage arm A binary** | JCIM-era CF-only FlexAID | Mach-O at `$FLEXAIDDS_LOCAL_ROOT/three_engine_entropy_q1/bin/A/FlexAID` (or documented path), SHA recorded in `docs/implementation/arm_pins.json` / REPRODUCIBILITY receipt; build from FlexAID pin `f766a14e…` if rebuild required |
| **3.3 Stage arm B binary** | First-entropy FlexAID | Mach-O at `…/bin/B/FlexAID`, SHA recorded; pin `1a6ae0b0…` entropy tip |
| **3.4 Stage arm C binary** | Current FlexAIDdS | Local pin via `resolve_build.py --check` / `--write-pin`; SHA matches receipt; **9dc9** matrix co-located |
| **3.5 Matrix** | Claim pin | Confirm `MC_st0r5.2_6.dat` md5 **9dc93717dfed0698006d88dd6a9627bc** under three_engine data |
| **3.6 Dry pipeline** | Fail-closed serial | `run_comparative_phases.py --pipeline-dry` (or project equivalent) exits clean with P2 no longer HOLD |
| **3.7 Decision memo** | One page | Write `workorders/PUB_PREGATE_TRIAGE.md`: P2 pass/fail, A/B/C SHA table, blockers remaining, **explicit “full-85 NOT authorized”** |

**Out of scope for Step 3:** launching N=85 three-arm science dock; C0 relaunch; dual full-85; memetic unlock.

### Step 4 — (Later goal, after 3.x green) Comparative pilot / claim design

Only if P2 pass **and** A/B/C binaries staged:

- Pilot8 or small fairness smoke under 9dc9 + seed-off  
- Then N=85 multi-seed / S_top10 under admission contract  
- Always report S1 / S2 / STRICT / claim_ready separately  

---

## 6. Explicit non-goals until pre-gates pass

- Full-85 / dual full-85 / C0 claim relaunch  
- Bundling PHENOTYPE_UNIQUE with BOOM, election, MUT_GRAN, or BASIN_REINJECT  
- Claiming true ΔG or claim_ready from S4 alone  
- Reopening VOID levers (WAL, burial weight ladder, memetic flags)  

---

## 7. One-line summary

**S4 A is mid-control (~70% of control budget done; 1N1M control already matches historical 6.40 elect).**  
**Next: finish treatment → full a posteriori + pins + unlock → then P2 oracle + A/B/C binary staging triage — not full-85.**
