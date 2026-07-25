# Docking / Scoring Bug Audit — 2026-07-25

**Repo:** `/Users/lp.more/Projects/FlexAIDdS`  
**Branch audited:** `main` @ `77f911a1` (origin/main)  
**Pilot OUT:** `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740`  
**Binary used by pilot:** `build/FlexAIDdS` (mtime 2026-07-25 02:16; includes `[SEARCH-COVERAGE]` strings from branch work not present on current `main` HEAD)  
**Method:** code path traces + pilot log/result.csv sampling. No full-85, no memetic, no rebuild of long campaigns.  
**Classification key:** **confirmed bug** / **likely bug** / **expected science** / **insufficient evidence**.

---

## 1. Executive summary

- **WAL_COERCIVE is wired and read**, but is a **structural no-op for Voronoi CF.wal** on production soft-wall configs: the cap binds only for per-pair overlap \(o > 1.0\,\text{Å}\) (`soft_wall.h`), and deep interpenetration pairs are **not enumerated** by the Vcontacts surface loop (`vcfunction.cpp` self-comment ~23× undercount). Campaign OFF≡ON (including saturating `cf_wal≥45` sums) is **expected**, not a missing-env bug.
- **BOOM_INTERVAL=50 pilot never exercised BOOM injection.** DatasetRunner claim/autonomous configs hardcode `boom_inject_fraction: 0.0` (`DatasetRunner.cpp:6058`). Injection requires `interval>0 AND fraction>0` (`gaboom.cpp:986-987`). Pilot logs: `boom_fraction 0.000→0.000` and **zero** `[BOOM]` lines. STEP 3 “BOOM” gate measured seed/early-stop noise, not denser BOOM.
- **Early stop is real and common:** many restarts terminate well under 2000 gens via fitness stagnation / entropy convergence (e.g. 1N1M r2 @ gen 310; 1M2Z many ~400–700 gens). Catastrophic-mutation spam does not restore allele entropy on 1N1M/1J3J.
- **Election path looks consistent on this pilot:** Softβ S1 OFF → min finite head CF (rank-0); `seed_echo=0`; BCR is diagnostic pool ceiling and does not replace top-1. 1N1M elect RMSD 5.66 / BCR 4.04 both fail clean-probe bar — not an election-only glitch.
- **Memetic remains correctly blocked** by wall-efficacy FAIL; re-keying the wall pilot gate to **pb_clash** (all-pairs, uncapped by design) is the coherent next scoring experiment — not another WAL_COERCIVE panel.

---

## 2. Confirmed bugs

### B1. `FLEXAIDDS_BOOM_INTERVAL` alone is a no-op on claim/autonomous path  
**Severity:** High (for experiment validity / sampling levers)  
**Class:** Confirmed wiring / product-config interaction bug (not a silent wrong-score bug)

**Trace:**
1. DatasetRunner always emits claim dock_config with  
   `"boom_inject_fraction": 0.0` — intentional after commit `1c142975` (prevent population wipe when seeds off).  
   File: `LIB/DatasetRunner.cpp:6044-6058`.
2. Config parse reads JSON then optional env:  
   `FLEXAIDDS_BOOM_INTERVAL` → `GB->boom_inject_interval`  
   `FLEXAIDDS_BOOM_FRAC` → `GB->boom_inject_fraction`  
   File: `LIB/config_parser.cpp:269-282`.
3. GA injection guard:  
   `if (interval > 0 && fraction > 0.0 && (gen % interval == 0) ...)`  
   File: `LIB/gaboom.cpp:986-987`.

**Repro (from pilot evidence):**
```text
# pilot 1N1M/dock_config.json
"boom_inject_interval": 100,
"boom_inject_fraction": 0,

# 1N1M/stderr.log
[SEARCH-COVERAGE] boom_interval 100→50  boom_fraction 0.000→0.000 ...
# rg '\[BOOM\]' over all restarts → empty
```

**Impact:** STEP 3 W1 gate treated BOOM_INTERVAL=50 as the one variable; actual BOOM inject count stayed 0. 1N1M elect 2.28→5.66 and 1YGC 3.34→1.75 cannot be attributed to denser BOOM. Mean ΔRMSD/BCR is multi-seed noise + early-stop variance under identical BOOM-off policy.

**Note on binary vs tree:** Pilot binary still emits `[SEARCH-COVERAGE]` (from `LIB/SearchCoverage.h` on branch `fix/remediate-campaign-failures`, **not present on current `main` HEAD**). That helper’s own docs state: *interval-only change is a no-op when fraction==0*; fraction restore happens only under `FLEXAIDDS_SEARCH_COVERAGE=1`, which the pilot did not set.

---

### B2. `soft_wall_fitness_energy` ignores `coercive` on legacy `soft_wall_cutoff==0` path  
**Severity:** Medium (production claim uses `soft_wall_cutoff=0.40`, so production soft-core path is unaffected; still a real API lie)  
**Class:** Confirmed code bug

```99:123:LIB/soft_wall.h
inline double soft_wall_fitness_energy(..., bool coercive = false, ...)
{
	if (soft_wall_cutoff > 0.0f) {
		// ... cap skipped only when coercive ...
		if (!coercive && Ewall_sc > WAL_CONTACT_CAP) return WAL_CONTACT_CAP;
		return Ewall_sc;
	}
	const double Ewall_raw = wall_energy_raw_r12(d, cr);
	return (Ewall_raw > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall_raw;  // always caps
}
```

And the `vcfunction` hard-wall branch also ignores `wal_coercive`:

```580:586:LIB/vcfunction.cpp
if (FA->soft_wall_cutoff > 0.0f) {
    Ewall_fitness = soft_wall_fitness_energy(d, cr, FA->soft_wall_cutoff, wal_coercive, wal_stiff);
} else {
    Ewall_fitness = (Ewall > WAL_CONTACT_CAP) ? WAL_CONTACT_CAP : Ewall;  // no wal_coercive
}
```

**Repro:** unit call with `cutoff=0`, deep clash, `coercive=true` vs `false` → identical capped energy. Soft path `cutoff=0.4`, `o=1.5` → OFF=50, ON=112.5 (works).

---

### B3. WAL_COERCIVE cannot affect deep-burial CF even when flag ON (structural scoring blind spot)  
**Severity:** High (for wall-as-burial-opponent thesis / memetic interlock)  
**Class:** Confirmed design defect of using Voronoi wall as the coercive lever (not env-miss)

**Env is read correctly** (not a missing-read bug):

```161:163:LIB/vcfunction.cpp
static const bool wal_coercive =
    (std::getenv("FLEXAIDDS_WAL_COERCIVE") != nullptr &&
     std::getenv("FLEXAIDDS_WAL_COERCIVE")[0] != '0');
```

Passed into soft-core energy at `vcfunction.cpp:582-583`.

**Why OFF≡ON:**
| Fact | Evidence |
|------|----------|
| Cap is **per pair**, not on summed `cf_wal` | `soft_wall.h:118`; `cfs->wal +=` at `vcfunction.cpp:588` |
| Cap binds only for \(o > 1.0\) Å at default k=50, cutoff=0.40 | Numeric: E(0.5)=12.5, E(1.0)=50, E(1.5)=112.5 |
| Summed `cf_wal≥45` can be many sub-cap pairs | saturating panel design |
| Deep clashes **never visit** Voronoi loop | `vcfunction.cpp:589-595` (engine comment; ~23× undercount) |
| Campaign measured byte-identical dCF on 5/5 sat poses | `workorders/WALL_ORACLE.md` |

**Repro:** any production `--config` score-only OFF vs ON on native/falsemin/sat-buried poses with only surface soft contacts → identical `cf_total` / `cf_wal`.

This is **not** “env not applied”; it is “applied flag never binds on reachable contacts.”

---

## 3. Likely bugs / risks

### L1. Catastrophic mutation thrash without diversity recovery  
**Class:** Likely harmful search dynamics (expected under current thresholds)

- 1N1M: ~90 catastrophic-mutation log lines across restarts; allele H stays ~0.25–0.30.
- 1J3J: similar spam from gen ~90; grid ~434k points; elected RMSD ~62 Å (absurd / wrong-basin).
- Guard only fires first half of planned gens (`gaboom.cpp:934`); if SEC/stagnation ends the run early, half-window still allows long thrash sequences.
- BOOM-off + thrash mutation is the only diversity path; it does not re-seed orientation basins.

### L2. Early GA termination vs advertised 2000 gens  
**Class:** Likely risk for eval-budget claims

Pilot TIMING SUMMARY samples (gens timed ≪ 2000):

| Target | Example gens timed |
|--------|--------------------|
| 1N1M r2 | 309 (entropy @ 310) |
| 1N1M r0 | 599 (CF stagnant) |
| 1M2Z r1 | 399 |
| 1J3J r1 | 129 |
| 1YGC r0/r2 | 2000 full |

Stdout patterns: `GA terminated early by fitness stagnation` / `Entropy convergence at generation N`.  
Not a reporting bug — real early stop. Risk: campaigns that claim “2000×pop” evals overstate work when SEC/stagnation fires.

### L3. 1M2Z CF ≈ 0 elected poses  
**Class:** Likely search/site failure, not election math

Elected CF=`-0.4233`, RMSD~13–15 Å, early stagnation with near-zero CF. Indicates basin never found productive contacts. Wall coercive would not help surface CF.wal of productive natives (~26 on crystal score path).

### L4. `SearchCoverage` / `[SEARCH-COVERAGE]` absent from current `main` tree  
**Class:** Likely process/binary drift risk

- Pilot binary strings contain `[SEARCH-COVERAGE]`.
- Current `main` HEAD has **no** `LIB/SearchCoverage.h`, no gaboom call site.
- Feature lives on `fix/remediate-campaign-failures` only.
- Risk: agents re-read main sources and mis-attribute pilot logs; or rebuild main and lose log/behavior silently.

### L5. Election gap metrics (ordered vs Hungarian) can confuse gate tables  
**Class:** Likely documentation footgun (code is intentional)

- `best_cluster_rmsd` / BCR ceiling uses **ordered direct** RMSD (`DatasetRunner.cpp:7123-7129`).
- `rmsd_hungarian` is diagnostic; success_rmsd uses ordered + `!seed_echo` (`DatasetRunner.h:164-165`).
- Gate table mixes elect hungarian-ish columns with BCR ordered — 1N1M elect 5.66 vs BCR 4.04 is a real pool-vs-elect gap, but baseline “elect 2.28 < BCR 4.08” needs same column definitions when comparing campaigns.

---

## 4. Ruled out / expected behavior

| Lead | Verdict | Why |
|------|---------|-----|
| WAL_COERCIVE env not read | **Ruled out** | `vcfunction.cpp:161-163` + soft path passes flag |
| Cap applied before env | **Ruled out** for soft-core path | Cap inside `soft_wall_fitness_energy` after computing E; coercive skips it |
| Production never hits soft-core | **Ruled out** | `soft_wall_cutoff: 0.40` in ops/gates configs + DatasetRunner emit |
| probe_cf omits config in wall scripts | **Ruled out** for oracle scripts | `wall_coercive_oracle.py` requires production config; fail-closed |
| BOOM double-counting / wrong injection path | **Ruled out** for pilot | Injection count = 0; no `[BOOM]` |
| Election bug causing 1N1M 5.66 | **Insufficient as sole cause** | Softβ OFF, rank-0 CF elect; seed_echo=0; BCR also >2; CF identical across 3/5 restarts (−99.3141) → false-min consensus |
| free_energy_strict default flipping rank-0 | **Expected post-E1b on main** | Pilot logs: `Softβ S1 OFF … elect min finite head CF` |
| probe_cf without --config ~200× CF | **Known footgun; not production path** for wall/pilot oracles that pass config | Still dangerous for ad-hoc probes |

**WAL_COERCIVE OFF≡ON on saturating panel:** **expected science/structure**, not measurement error (see B3).

**1YGC genuine elect 1.75:** **expected stochastic win** under BOOM-off; not proof BOOM helped.

---

## 5. Recommended next experiments (one variable each; no full-85)

1. **BOOM efficacy smoke (not a success gate):**  
   One target (e.g. 1N1M), R=1–2, same seed base if possible:  
   - Arm A: default (frac=0)  
   - Arm B: `FLEXAIDDS_BOOM_FRAC=0.5` only (interval default 100)  
   - Arm C: `FLEXAIDDS_BOOM_FRAC=0.5` + `FLEXAIDDS_BOOM_INTERVAL=50`  
   **Pass criterion for the lever (not docking success):** logs show `[BOOM] injection #N` and `boom_fraction` > 0. Compare inject count and SEC gens only.

2. **pb_clash burial oracle (replacement for wall STEP 2):**  
   Score-only native vs sat-buried decoys with production configs:  
   - Arm A: `FLEXAIDDS_PB_CLASH_WEIGHT=0`  
   - Arm B: single weight e.g. `1.0` (or a 3-point ladder as separate arms)  
   Accept if dCF moves toward native-min on ≥4/5 and no clean native CF-min regresses.  
   **Do not** re-run WAL_COERCIVE panels expecting OFF≠ON without constructing per-pair \(o>1\) **and** Voronoi-visible contacts (likely impossible by B3).

3. **Coercive unit fix + micro test (code only):**  
   Make legacy `cutoff==0` and `vcfunction` hard branch honor `coercive`; add gtest: deep clash OFF capped / ON uncapped for both soft and legacy.  
   No docking required.

4. **Early-stop accounting spot check:**  
   Parse TIMING SUMMARY / “terminated early” rates on 8-panel under default knobs; report median gens completed. One variable: e.g. `FLEXAIDDS_NO_SEC` or documented SEC off — only if that env already exists and is intended for A/B (do not invent).

5. **1J3J grid pathology:**  
   Score-only / short R=1 with pocket constraint or finer site definition vs default 434k grid; measure whether elected CF basin leaves the absurd 60+ Å modes. Single site-definition variable.

---

## 6. Files inspected

### Core scoring / wall
- `LIB/soft_wall.h` (WAL_CONTACT_CAP, soft_wall_fitness_energy, coercive)
- `LIB/vcfunction.cpp` (wal_coercive/stiff statics; wall accumulation; Voronoi surface note)
- `LIB/Vcontacts.cpp` (clash_value soft_wall path — no coercive arg)
- `tests/test_soft_wall.cpp`
- `scripts/wall_coercive_oracle.py`, `scripts/wall_saturating_panel.py`
- `workorders/WALL_ORACLE.md`
- `~/flexaidds_results/workorders/WALL_ORACLE_FAIL_EXPLAINED.md` (external OPS note; agrees with B3)
- `docs/SCORING.md` (claims)
- `ops/gates/configs/*_dock_config.json` (soft_wall_cutoff=0.4)

### BOOM / GA / DatasetRunner
- `LIB/config_parser.cpp` (BOOM_INTERVAL / BOOM_FRAC / SIGMA_SCALE)
- `LIB/gaboom.cpp` (boom inject, elitism, diversity catastrophic mutation, SEC/stagnation)
- `LIB/DatasetRunner.cpp` (boom_inject_fraction hard 0.0; election / BCR / seed_echo / result.csv)
- `LIB/DatasetRunner.h` (metric definitions)
- `LIB/gaboom.h` (boom_inject_* fields)
- `LIB/config_defaults.h`
- `docs/implementation/WAVE3_SAMPLING_KNOBS.md`
- Git history: `1c142975` (boom_frac=0 claim path); `LIB/SearchCoverage.h` on `fix/remediate-campaign-failures` (not main)

### Pilot evidence
- `~/flexaidds_results/pilot_w1_boom_interval_20260725_134740/`  
  `pilot.log`, `astex_diverse_results.csv`, `step3_pilot_gate.md`,  
  `1N1M/**/stderr.log|stdout.log|dock_config.json`,  
  `1YGC/**`, `1M2Z/**`, `1J3J/**` samples
- `workorders/STEP3_PILOT_GATE.md`, `workorders/CAMPAIGN_GATE_SUMMARY.md`

### Election / free energy
- `LIB/SoftBetaFreeEnergy.h` (`free_energy_strict`)
- `LIB/cluster.cpp` (rerank hooks)
- Pilot `[SOFTBETA-ELECT]` / `[Z+H]` / `[ELECTED-POSE]` / `[BCR-CEILING]` lines

### Not run
- Full ctest suite (optional; wall unit tests exist, no coercive gtest)
- Full-85 / memetic / rebuild under live benchmark lock

---

## Appendix A — Soft-wall energy vs overlap (production defaults)

`soft_wall_cutoff=0.40`, `k_wal=50`, `WAL_CONTACT_CAP=50`:

| o (Å) | E_sc | Cap binds? | Coercive differs? |
|------:|-----:|:----------:|:-----------------:|
| 0.2 | 4.0 | no | no |
| 0.4 | 8.0 | no | no |
| 0.8 | 32.0 | no | no |
| 1.0 | 50.0 | boundary | no |
| 1.2 | 72.0 | yes | **yes** |
| 2.0 | 200.0 | yes | **yes** |

---

## Appendix B — 1N1M pilot termination (all restarts BOOM-off)

| Restart | End reason (stdout) | ~gens |
|---------|---------------------|------:|
| r0 | CF stagnant 400 gens (best≈−98.56) | 599 |
| r1 | CF stagnant 300 gens (best≈−99.31) | 999 |
| r2 | Entropy convergence gen 310 | 309 |
| r3 | Entropy convergence gen 550 | 549 |
| r4 | CF stagnant 300 gens (best≈−99.31) | 499 |

Elected restart r3 CF=−99.3141, rmsd_hungarian=5.66, best_cluster_rmsd=4.04, seed_echo=0.

---

*Audit only. No gate force-pass. No code fix applied in this session.*
