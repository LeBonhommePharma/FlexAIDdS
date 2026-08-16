# Forward plan — raise FlexAIDdS Astex genuine success rates

**Status:** normative sequencing for OPS / Grok / Claude Code (analysis + next implementation order).  
**Does not claim restored rates.** Incomplete runs are mechanism evidence only when labeled incomplete.  
**Not a published FlexAIDdS success rate** — unverified / no METHODOLOGY.md §0 receipt. Do not cite 25.3% as current docking power.  
**Hub (start here):** [`COMPARATIVE_SCIENCE_README.md`](COMPARATIVE_SCIENCE_README.md).  
**OPS session baseline (not publishable):** [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md) — 20/79 = 25.3% (no receipt).  
**Parents:** `METHODOLOGY.md`; [`COMPARATIVE_BENCHMARK_METHODOLOGY.md`](COMPARATIVE_BENCHMARK_METHODOLOGY.md);  
`$FLEXAIDDS_LOCAL_ROOT/workorders/SYNTHESIS_opus5_audit_and_comcap_verdict_2026-07-24.md` (when present);  
`CAUSAL_ANALYSIS_rate_regression_2026-07-24.md`; `ROOTCAUSE_CORRECTED_2026-07-24.md`; Opus/Science handoffs.

---

## 0. Where we are (honest snapshot — incomplete trees)

| Tree | Genuine (seed_echo=0, rank-0 ≤2 Å) | Label |
|------|----------------------------------|--------|
| **`v_autonomous_20260724_160919` (80/85 scored)** | **20/79 = 25.3%** | OPS session record — **not publishable**; seed-echo **0**; BCR 22/79=27.8%. See [`BASELINE_GENUINE_2026-07-24.md`](BASELINE_GENUINE_2026-07-24.md) |
| C0@9dc9 v4 (partial ~17) | **6/17 (~35%)** | Best partial **defined-cleft** signal; **not** full-85 |
| C0@9dc9 v6 (partial ~11) | 3/11 (~27%) | Same binary as v4 on overlap — **sampling artifact**, not regression |
| control noentropy (partial 17) | 4/17 (~24%) | 9dc9 |
| softβ+CAP@72d7 autonomous (83) | 12/83 (~15%) | **Confounded** (matrix+mode+multi-knob) |
| com-cap fixed CAP=-130 (10 finished) | 1/10 (10%) | **UNCITABLE** — OOM dead run |

**2026-07-25 load-bearing notes:** Baseline shows BCR≈genuine (gap ~2 targets) and a **sampling ceiling ~28%**. Per OPS campaign methodology, this run **predates** `free_energy_strict` as the measured product default — **do not** cite 25.3% as proof the election fix worked; that remains **unmeasured** until a post-fix A/B. Route: E10 offline, wall oracle, then Wave 3 BCR. Softβ S1 is not the primary rate lever.

**Targets (fair claim methodology):**  
- Headline comparative: **S_top10** median over 10k bootstraps (3Dsig family).  
- Modern KPI: **S1** / **genuine** (rank-0, seed_echo=0, ≤2.0 Å).  
- Diagnostic: **BCR** (sampling ceiling).  

**Published anchors (do not mix):**  
- JCIM 2015 Table 2 Astex native FLRP: **top-1 = 45.2%**, **top-10 = 66.7%**.  
- 3Dsig 2017 presentation red medians (historical, not a current receipted campaign): FlexAID **~0.66** / FlexAIDdS **~0.69** (S_top10-style).  
- Production matrix: **`MC_st0r5.2_6.dat` MD5 `9dc93717dfed0698006d88dd6a9627bc`**.  
- Claim budget: **pop×gen = 1000×2000**, **R=10**, **seed-off**.

**Restore floors (from RESTORE_SUCCESS_SPEC):** genuine top-1 ≥ **45%** full-85 (or ≥3/6 discriminating subset as pre-check); historical v119 ~52% (caveats: seed-assist possible).

---

## 1. Election vs sampling (hard ceilings — do not re-confuse)

| Layer | What it can fix | What it cannot |
|-------|-----------------|----------------|
| **Sampling / search** | Raise **BCR** (near-native appears among heads) | — |
| **Election / ranking** | Convert BCR>0 into **S1 / S_top10** when good heads exist | Invent sub-2 Å if best head ~3 Å (frozen 1G9V often **~3.17 Å** ceiling) |
| **CF / wall / com** | Change which basins are CF-min and whether GA prefers burial | Alone, not fix gene-space niche starvation |

**Measured mechanisms (VERIFIED in source / on-disk):**

| Mechanism | Evidence | Implication |
|-----------|----------|-------------|
| **Search coverage dominant** | Clean probes: native often CF-min but blind BCR 5–7+ Å; near-native **freq~1** vs over-burial **freq~1666** | Primary for BCR=0 class |
| **ACF size bias (implicit population election)** | `cluster.cpp:178` uses `soft_beta::acf` → \(E_{\min}-T\ln Z\); header marks `acf` diagnostic; `free_energy_strict` collapses exact CF dups | At \(T>0\), size can dominate CF spread; **explicit** `ELECT_BY_POPULATION` OFF ≠ product free of size bias |
| **Wall saturates / com unbounded** | `soft_wall.h` soft-core; com floors deep negative | Enables burial free-lunch; wall fix before memetic |
| **Gene-space niche (`calc_rmsp`)** | `sig_share` from gene bounds; grid ordinal + angles | Basin starvation / BOOM wasted |
| **COM_BURIAL_CAP=-130 run** | OOM workers=6; 10/85; 1G9V elected 10.6 Å, BCR 2.34; totals still −1144 | **UNCITABLE**; per-optres/total not global bound; single-target-tuned |
| **Softβ DatasetRunner S1** | Default **OFF**; reorders heads only; BCR=0 ⇒ cannot create S1 | Not a sampling fix |

**Resolved slogan conflict:**  
- **REJECT** “elect by raw cluster frequency as the product default” for 1G9V-class (picks over-burial).  
- **KEEP** fixing **ACF multiplicity inflation** (shipped path at \(T>0\)) via `free_energy_strict` — different object than explicit population flag.

---

## 2. Full recommendation inventory

Tags: **KEEP** (do next / in sequence) · **DEFER** (after gates) · **REJECT** (do not chase) · **UNCITABLE** (incomplete/confounded evidence).

| ID | Recommendation | Tag | Grounding | Role |
|----|----------------|-----|-----------|------|
| **S1** | Search coverage: BOOM_INTERVAL + SIGMA_SCALE coupling; more coarse orientations; DoF budget with P1 | **KEEP** | ROOTCAUSE_CORRECTED; Claude Code P1; Opus M3 | Sampling — raise BCR |
| **S2** | Coarse-init CF-rank seeds (not CF&lt;0) | **KEEP** (landed on restore branches) | coarse_init fix | Sampling |
| **S3** | Hard-clash severity (not flat CF=10000) | **KEEP** (landed) | clashfix | Search gradient |
| **S4** | Water O.3 typing (not C.1) | **KEEP** (landed) | ca897577 lineage | Scoring hygiene |
| **S5** | Bigger pop / more gens as primary lever | **REJECT** as sole fix | Opus: basin visited once / not developed | Budget not main bottleneck |
| **S6** | CMA-ES as primary search (as-is) | **REJECT** until gene-0 discrete decode | Invalid CF=10000 probes | Search backend |
| **S7** | Gene-0 CMA decode + fair CMA vs BOOM A/B | **DEFER** | After E1/E2/search P1 | Sampling alt |
| **E1b** | Cluster emission: `acf` → `free_energy_strict` (flag-gated) | **KEEP** (first cheap election code) | cluster.cpp:178; SoftBetaFreeEnergy.h | Election |
| **E1a** | Explicit population-weighted election ON | **REJECT** for default / 1G9V | freq-1666 = 10.4 Å | Election |
| **E_softβ_S1** | DatasetRunner Softβ S1 default ON | **REJECT** as default; **KEEP** as opt-in A/B only | softbeta_election_policy | Election reorder |
| **E_softβ_default_ON_branch** | c297217c default ON soft-β | **DEFER** / do not merge blindly | main still default OFF | Election |
| **E2** | Un-cap / redesign steric wall | **KEEP** (after E1b design, before memetic) | soft_wall.h; Opus M2 | Scoring / enables refine |
| **E3** | Validity / admissibility gate on poses | **DEFER** | After E1b/E2 | Election hygiene |
| **E4** | Cartesian niche metric (replace gene rmsp) | **KEEP** (search wave, after E1b pilot) | gaboom sig_share / calc_rmsp | Sampling |
| **E5** | Memetic local refinement | **DEFER** — **only after E2** | Opus: refiner buries deeper if wall saturates | Sampling polish |
| **E7** | Wire `elec` to JSON | **KEEP** (independent, cheap) | Opus M5 | Scoring |
| **E10** | Offline independent rescore of heads (Vinardo/smina) | **KEEP** (do first — offline) | Opus E10 | Diagnostic split election vs score |
| **C1** | COM_BURIAL_CAP global product default | **REJECT** until redesign | per-optres ≠ total bound | Scoring |
| **C2** | COM_BURIAL_CAP=-130 full85 “validation” | **UNCITABLE** | OOM 10/85; CAP tuned to 1G9V native; empty git_commit | — |
| **C3** | COM_FLOOR / VCT_NORM canary (merged-but-OFF knobs) | **DEFER** clean serial A/B | DEEPDIVE armA | Scoring |
| **C4** | Polar desolv high weight alone | **REJECT** as sole lever | canaries 0/3 | Scoring |
| **M9** | Matrix production **9dc9** | **KEEP** | VCT audit; 72d7 rejected | Protocol |
| **M72** | Re-pin production to 72d7 | **REJECT** | packing fork; canary confound | Protocol |
| **H1** | Pocket HEM strip 1P2Y/1R9O only | **KEEP** (idempotent prep) | FINDINGS / probe_cf | Prep |
| **H2** | Strip distal 1G9V/1Q4G HEM as fix | **REJECT** | byte-identical CF | Prep |
| **P_proto** | Comparative A/B/C methodology (JCIM CF vs soft-β vs FlexAIDdS) | **KEEP** | COMPARATIVE_BENCHMARK_METHODOLOGY.md | Measurement |
| **P_metric** | Fail-closed genuine / S_top10 / no seed-echo | **KEEP** | aggregate_claim_metrics; bootstrap_3dsig | Measurement |
| **P_box** | Workers ≤2–4; serial science; no dual full85 | **KEEP** | OOM evidence | Ops |
| **P_receipt** | Always record binary SHA, commit, FLEXAIDDS_* scoring env | **KEEP** | com-cap gap | Ops |
| **V_sas** | sas_weight 1.0 as restore lever | **REJECT** as sole | red herring after cancel | Scoring |
| **V_vib** | Vib/Shannon physical ΔG as docking ranker default | **REJECT** for red-bar / S1 claim without labels | 3dsig contract | Thermodynamics |
| **R_restore_v119** | Target ≥45% genuine / aim ~52% | **KEEP** as **goal metric**, not a code patch | RESTORE_SUCCESS_SPEC | Goal |

---

## 3. Priority-ordered implementation sequence

### Box rules (every step)

- M3 Pro **~18 GB**: **one** science owner; **workers 2–4** max; `OMP_NUM_THREADS=1` per worker.  
- **No dual full85**. No `cmake --build` of mmap’d live binary.  
- Claim-style docks: matrix **9dc9**, **1000×2000**, **R=10** when citing rates, **seed-off**, label defined-cleft vs autonomous.

---

### Wave 0 — offline (no full-85 dock) — **start here**

| Step | Action | Acceptance test | Must not regress |
|------|--------|-----------------|------------------|
| **W0.1** | **E10:** Rescore frozen heads (1G9V com-cap BCR 2.34 vs elected 10.62; C0 partial successes) with independent score (smina/Vinardo or probe_cf panel) | Report: fraction of targets where independent scorer ranks near-native head above elected decoy | No product default change |
| **W0.2** | **E1b design:** env `FLEXAIDDS_ACF_STRICT=1` (name bikeshed OK) switches `cluster.cpp` emission from `acf` → `free_energy_strict` | Unit test: exact-duplicate CF members **do not** deepen G̃; default OFF = bit-identical | Default path parity |
| **W0.3** | Document receipt schema for scoring env (CAP, ACF_STRICT, TEMPER, matrix md5) | Template fields present in ops launcher or RUN_RECEIPT writer | — |

**First three KEEP actions implementable without full-85:** W0.1, W0.2, W0.3 (then W1 serial pilot).

---

### Wave 1 — flag-gated code + serial pilot (≤6–17 targets)

| Step | Action | Acceptance test | Must not regress |
|------|--------|-----------------|------------------|
| **W1.1** | Land **E1b** behind default-OFF flag; pilot on discriminating set **1G9V 1M2Z 1N1M 1J3J 1K3U 1L7F** + known goods **1HNN 1HP0 1HQ2** | On 1G9V-class: if BCR≤2.5 and size-dominated ACF was electing large wrong basin, **S1 improves or BCR-stable with better head CF rank**; log ACF vs strict G̃ | 1HNN/1HP0/1HQ2 genuine must not flip success→fail |
| **W1.2** | **E7** wire `elec` to JSON (default OFF) | Native CF oracle panel: elects when enabled does not mass-invert clean probes | Default OFF parity |
| **W1.3** | **Native-CF oracle A/B for E2 wall** (frozen poses, no GA) | Count CF inversions native vs decoy before/after wall change | No silent default ON |

---

### Wave 2 — scoring physics (wall before burial hacks)

| Step | Action | Acceptance test | Must not regress |
|------|--------|-----------------|------------------|
| **W2.1** | **E2** steric wall redesign / uncap (gated) | probe_cf: reduce com-over-burial invert rate on 1G9V **without** worsening 4/5 clean probes’ native-as-CF-min | Default OFF until pilot PASS |
| **W2.2** | **Do not** merge COM_BURIAL_CAP=-130 | N/A | — |
| **W2.3** | Optional COM_FLOOR / VCT_NORM **serial** canary **after** W2.1 design clarity | com magnitude no longer ~−3000 class on seed-able targets | Goods stay genuine |

---

### Wave 3 — sampling (BCR raisers)

| Step | Action | Acceptance test | Must not regress |
|------|--------|-----------------|------------------|
| **W3.1** | **S1** P1 anti-collapse: BOOM_INTERVAL × SIGMA_SCALE factorial (restarts=1 diversity bar) | Single-run near-native **freq>1** on 1J3J/1K3U; BCR drops | No SEC thrash; default OFF knobs |
| **W3.2** | **E4** Cartesian niche (or gene0-decoupled sharing) | Same targets: niche hits track Å proximity better than gene rmsp | — |
| **W3.3** | Coarse orientations env A/B (64 vs 256) on tight sites | BCR improvement on 1OF1/1J3J/1K3U | Cost linear — pilot only |
| **W3.4** | **E5** memetic **only if W2.1 PASS** | S1↑ on BCR-ready targets; native CF does not walk away by burial | Wall must hold |

---

### Wave 4 — measurement & claim

| Step | Action | Acceptance test |
|------|--------|-----------------|
| **W4.1** | Run comparative arms A/B/C per `COMPARATIVE_BENCHMARK_METHODOLOGY.md` | S_top10 bootstrap medians labeled; matrix 9dc9; R=10; seed-off |
| **W4.2** | Full-85 genuine claim only after pilot gates | genuine ≥45% or documented fail; **never** raw success column |
| **W4.3** | Report JCIM top-1 45.2% / top-10 66.7% / 3Dsig ~0.66–0.69 **separately** | No mixed bare % |

---

## 4. What “raising the fuck outta these numbers” means in practice

```
BCR floor (sampling)  ──W3──►  near-natives exist among heads
         │
         ▼
S1 / S_top10 (election) ──W0–W1──►  elect good heads (kill ACF size bias; honest Softβ)
         │
         ▼
Scoring landscape ──W2──►  stop rewarding over-burial (wall first; no CAP=-130)
         │
         ▼
Claim measurement ──W4──►  9dc9, R=10, seed-off, S_top10 + genuine labeled
```

**Rough ROI (qualitative, not a rate promise):**

| If bottleneck is… | Highest ROI next step |
|-------------------|------------------------|
| BCR=0 (never visit native) | W3.1 search P1 + E4 niche |
| BCR≤2, S1 fails (1G9V BCR 2.34, elect 10.6) | **W0.1 E10 + W1.1 E1b** then W2.1 wall |
| CF invert native vs decoy | W2.1 wall; not CAP=-130 |
| Reporting pollution | W4.3 metrics; kill seed-echo claims |

---

## 5. Explicit REJECT / UNCITABLE list (do not reopen without new evidence)

| Item | Why |
|------|-----|
| COM_BURIAL_CAP=-130 full85 tree | OOM, incomplete, single-target-tuned, receipt gaps → **UNCITABLE** |
| 72d7 as production | Rejected matrix fork |
| Explicit population election default ON | Picks over-burial basins |
| Softβ S1 as sampling fix | BCR=0 invariant |
| Polar-only high weight | 0/3 canaries |
| Distal HEM strip as 1G9V fix | No-op CF |
| Memetic before wall fix | Expected to worsen top-1 by burial |
| Dual full85 / workers≥6 on 18 GB | OOM proven |
| Comparing S1 single-run to JCIM top-10 66.7% or 3Dsig 0.66 as same number | Contract mismatch |

---

## 6. Immediate next 72 hours (concrete)

1. **E10 offline** on frozen 1G9V (+ 2–3 other election-gap targets) → write `workorders/E10_election_vs_scoring.md`.  
2. **Implement E1b** flag + unit test + default-OFF parity (no full dock).  
3. **Serial pilot** 6–9 targets, workers=2, 9dc9, defined-cleft if possible, ACF_STRICT on/off A/B.  
4. **Do not** re-fire com-cap full85 or softβ@72d7 autonomous as “entropy validation.”  
5. Keep **COM_BURIAL_CAP off main**.

---

## 7. References (on disk)

- `docs/implementation/COMPARATIVE_BENCHMARK_METHODOLOGY.md`  
- `docs/implementation/3dsig_red_pair_protocol.md`, `3dsig_shannon_ranking.md`, `softbeta_election_policy.md`  
- `METHODOLOGY.md`  
- `LIB/cluster.cpp` (~178), `LIB/SoftBetaFreeEnergy.h` (`acf` / `free_energy_strict`)  
- `LIB/soft_wall.h`  
- Workorders: `SYNTHESIS_opus5_audit_and_comcap_verdict_2026-07-24.md`, `CAUSAL_ANALYSIS_rate_regression_2026-07-24.md`, `ROOTCAUSE_CORRECTED_2026-07-24.md`, `RESTORE_SUCCESS_SPEC.md`, `FINDINGS_for_Grok_2026-07-24.md`
