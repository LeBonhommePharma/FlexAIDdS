# 3Dsig 2017 comparative methodology — reproduction contract

**Primary source (authoritative):**  
`/Users/lp.more/Downloads/Morency_LP_3Dsig_2017.pdf`  
L.-P. Morency, *The Impact of Conformational Entropy on the Accuracy of FlexAID in Binding Mode Prediction*, ISMB/ECCB 2017 — 3Dsig.

**Non-negotiable:** claim / comparative ranking and success statistics for “FlexAID vs FlexAID+entropy (FlexAIDdS)” must **reproduce this deck**, not invent a parallel protocol.

---

## 1. Physics / ranking (slides: “Introducing entropy in FlexAID’s scoring function”)

Density-based clustering groups poses into **binding modes**. Each mode is scored as a **soft-β free energy on the CF scoring proxy**:

\[
Z = \sum_{i\in\mathrm{mode}} e^{-\mathrm{CF}_i / T}
\quad
p_i = \frac{e^{-\mathrm{CF}_i / T}}{Z}
\]

\[
\tilde H = \sum_i p_i\,\mathrm{CF}_i
\quad
\tilde S = -\sum_i p_i\ln p_i
\quad
\tilde G = \tilde H - T\,\tilde S
\]

| Symbol | Meaning |
|--------|---------|
| \(\mathrm{CF}_i\) | Voronoi contact-function **scoring proxy** (a.u.), **not** experimental ΔG |
| \(T\) | Temperature in K; \(\beta=1/T\) (FlexAID soft-β, **not** \(1/(k_B T)\)) |
| Elect | **Lowest** \(\tilde G\) binding mode / head |

**Identity:** \(\tilde G = \tilde H - T\cdot\tilde S \equiv E_{\min}-T\ln Z_{\mathrm{local}}\) (cluster **ACF**).  
**Single implementation:** `LIB/SoftBetaFreeEnergy.h` (`flexaids::soft_beta::free_energy` / `acf`).  
**Local Z only:** each cluster/mode re-normalizes \(p_i\) over **its own members** — never the global `BindingPopulation::PartitionFunction` for ranking H/S/\(\tilde G\). Physical StatMech REMARK `free_energy` is separate; optional REMARK `soft_beta_G` logs the ranking objective.

### Code layers (must stay identical)

| Layer | Behavior | Log / API |
|--------|----------|-----------|
| `cluster.cpp` ACF | `soft_beta::acf(member CFs, T)`; emission order when \(T>0\) | `[ENTROPY_RANK]`; REMARK `soft_beta_G` |
| `BindingMode` classic H/S/F | Same \(\tilde G\) over **mode members** only (+ vib additive; not global Z) | BindingMode sort; REMARK `soft_beta_G` |
| `DatasetRunner` S1 | Same `soft_beta::free_energy` over heads + `.mcf` members | `[3DSIG-RANK]` |

| Env | Default | Meaning |
|-----|---------|---------|
| `FLEXAIDDS_ELECTION_SHANNON_F` | **0 (OFF)** | Elect by \(\tilde G=H-TS\) when set to `1` |
| `FLEXAIDDS_ELECTION_LEGACY_ZH` | 0 | Force legacy ZH / ≈ min-CF (not 3Dsig); same OFF path when Shannon unset |
| `FLEXAIDDS_ELECTION_SOFT_T` | 0 → resolve below | Soft-β \(T\) in K (env override) |
| `FLEXAIDDS_FORCE_CF_RANK_EMISSION` | 0 | Engine emits min-CF (rollback); classic SoftBeta path off |

**Default OFF until Astex pilot + SoftBeta identity.** Shannon S1 election stays **off** when both env knobs are unset. Enable with `FLEXAIDDS_ELECTION_SHANNON_F=1` after validation (`LIB/SoftBetaFreeEnergy.h`). Claim / 3Dsig launchers that need Shannon ON export that env (e.g. `scripts/run_C0_claim_clean.sh`).

**Soft-β \(T\) resolution** in `select_pose_freq_gated_pooled` (DatasetRunner S1):

1. `FLEXAIDDS_ELECTION_SOFT_T` if set and \(>0\) → log `source=env`
2. else **dock TEMPER** (`DockingConfig::temperature` / CONFIG `TEMPER`) if \(>0\) → log `source=dock`
3. else (legacy ZH only) `FLEXAIDDS_ELECTION_SCORE_TAU` if \(>0\) → log `source=env`
4. else **298 K** → log `source=fallback`

Log line: `[3DSIG-RANK] … T=… source=dock|env|fallback …`. Soft-β is \(\beta=1/T\) on the CF scoring proxy (not \(k_B\)).

**FlexAIDdS engine and DatasetRunner must not invent different ranking objectives.** Search still optimizes CF; when Shannon ranking is enabled, ranking/election uses \(\tilde G\) via `LIB/SoftBetaFreeEnergy.h`. Engine ACF and DatasetRunner election must share the same dock \(T\). Unit gates: `SoftBetaIdentity::*` and `BindingModeMatchesSoftBetaLocal` in `tests/test_classic_entropy_ranking.cpp`.

---

## 2. Benchmark methodology (slides: “Methodology used to benchmark FlexAID”)

Reuse **FlexAID JCIM 2015** comparative design (Gaudreault & Najmanovich, *J. Chem. Inf. Model.* 55:1323–1336, 2015), extended to entropy-on:

| Axis | 3Dsig 2017 PDF requirement |
|------|----------------------------|
| **Primary dataset** | Astex Diverse **N = 85** (native / cognate-pocket redock story in deck) |
| **Also reported in deck** | Astex Non-Native (N=1112), HAP2 flexibility set — secondary |
| **Per method, per case** | **10** independent simulations |
| **Budget each sim** | **2 000 000** energy evaluations |
| **Success (per case, bootstrap)** | RMSD **&lt; 2.0 Å** among the **top 10** predicted results |
| **Headline statistic** | **Median** success rate over **10 000** bootstrap resamples of the N cases (with replacement) |
| **Arms on barplots** | FlexAID · FlexAID+entropy (FlexAIDdS) · AutoDock Vina · FlexX · rDock |

### Mapping to current FlexAIDdS campaign knobs

| 3Dsig item | FlexAIDdS implementation |
|------------|---------------------------|
| 2 000 000 evals / sim | e.g. pop×gen = 2e6 per restart (e.g. 1000×2000) with **fixed gen**, pop×DoF only if documented |
| 10 sims / case | `FLEXAIDDS_RESTARTS=10` (or 10 independent jobs); **do not** silently use 5 without labeling |
| Top-10 success | Track success if **any of top-10 emitted modes** has RMSD &lt; 2 Å (**S_top10**); S1 = top-1 only is extra modern KPI |
| Bootstrap 10k median | `scripts/` analysis: resample cases, recompute success rate, report median + CI |
| Matrix | JCIM 2015 pin — **do not change mid-campaign** (see `docs/implementation/MATRIX_PIN_JCIM2015.md`) |
| Site / seed | Cognate pocket, **no native-pose seed** for claim-style fair compare |

### Success definitions (report all; label clearly)

| ID | Definition | 3Dsig deck |
|----|------------|------------|
| **S_top10** | Min RMSD among top-10 ranked modes ≤ 2 Å | **Primary in PDF bootstrap** |
| **S1** | Rank-0 / elected mode RMSD ≤ 2 Å | Modern claim KPI (stricter) |
| **S2** | S1 ∧ PoseBusters | Modern secondary (not in 2017 deck) |
| **BCR** | Best cluster-head RMSD ≤ 2 Å | Diagnostic sampling ceiling |

---

## 3. Density-based clustering (PDF)

- Binding modes = density-based clusters of similar poses (small movements, multi-conformer modes).  
- Production FO: **single** MinPts (literature Ankerst/Sander) — see `docs/implementation/fo_minpts_literature.md`.  
- **3Dsig red-pair arm B** must use that same single-pass FO (`CLUSTA FO` + engine `fo_choose_minpts`); **never** a triple MinPts ladder.  
- Engine emission names: FO dual-suffix `prefix_minPts_rank.pdb`; CF/DP `prefix_rank.pdb`. DatasetRunner must enumerate **both**.

---

## 4. Reproduction checklist (before claiming “3Dsig reproduced”)

- [ ] PDF open and methodology section re-read  
- [ ] Ranking: \(\tilde G=H-TS\) soft-β on CF; FlexAIDdS emission order ≡ DatasetRunner election objective  
- [ ] Log shows `[3DSIG-RANK]` / `[ENTROPY_RANK]` as appropriate  
- [ ] Matrix MD5 pinned and recorded in RUN_RECEIPT  
- [ ] Astex 85 (and secondary sets if claiming those bars)  
- [ ] 10 sims × 2e6 evals (or explicit deviation table)  
- [ ] S_top10 + 10k bootstrap median success rates for barplots  
- [ ] Arms labeled FlexAID / FlexAID+entropy / Vina / FlexX / rDock as applicable  

---

## 5. Red-pair execution (FlexAID vs FlexAIDdS bars only)

Operational freeze + serial launch for the **red** bars:

- **Protocol:** `docs/implementation/3dsig_red_pair_protocol.md`
- **Launcher:** `scripts/run_3dsig_red_pair_serial.sh` (A → B0 → B, `R=10`, `pop×gen=2e6`)
- **Archived figures:** `scripts/plot_3dsig_archived_bars.py`
- **Bootstrap metric:** `scripts/bootstrap_3dsig_s_top10.py`

C0 packaging campaigns are **not** the 2017 red-bar path until FO dual-suffix election is validated.

---

## 6. Related docs

- `docs/architecture/scoring_pipeline_schematic.md` — search vs score vs rank  
- `docs/classic_entropy_ranking.md` — engine emission  
- `docs/ensemble_pipeline.md` — 4-layer reproducibility  
- `LIB/SoftBetaFreeEnergy.h` — shared \(\tilde G\) / ACF math  
- `benchmarks/protocols/three_engine_entropy_comparison.md` — A/B0/B/C0 arms  
