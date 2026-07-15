# Literature synthesis — density clustering + binding-mode election

**Status:** swarm synthesis for FlexAID / FlexAIDdS production and 3Dsig red-pair.  
**Scope:** density-based binding modes, MinPts policy, soft-β \(\tilde G\) election, Density Peak (DP) role, and code-path mapping.  
**License hygiene:** clean-room Apache-2.0 implementation notes only; no GPL algorithm dumps or GPL-sourced snippets. FlexAID JCIM 2015 was historically GPLv3; FlexAIDdS is a separate Apache-2.0 modernization — cite the *science*, reimplement under Apache-2.0.  
**Primary sibling notes:**

| Doc | Role |
|-----|------|
| [`fo_minpts_literature.md`](fo_minpts_literature.md) | Single-pass FastOPTICS MinPts policy |
| [`3dsig_shannon_ranking.md`](3dsig_shannon_ranking.md) | Soft-β \(\tilde G\) ranking contract |
| [`3dsig_red_pair_protocol.md`](3dsig_red_pair_protocol.md) | Red FlexAID vs FlexAIDdS arms A/B0/B |
| [`docs/classic_entropy_ranking.md`](../classic_entropy_ranking.md) | Engine emission ACF / BindingMode F |
| [`docs/architecture/scoring_pipeline_schematic.md`](../architecture/scoring_pipeline_schematic.md) | Search vs cluster vs elect layers |

**Primary 3Dsig source:** L.-P. Morency, *The Impact of Conformational Entropy on the Accuracy of FlexAID in Binding Mode Prediction*, ISMB/ECCB 2017 — 3Dsig (deck PDF; local: `Morency_LP_3Dsig_2017.pdf`). Not a peer-reviewed article; methodology must be reproduced as stated in the deck + JCIM 2015 design.

---

## 1. Density-based binding modes vs CF centroid clustering

### 1.1 What a “binding mode” is (3Dsig + docking practice)

The 3Dsig deck defines binding modes as groups of **similar poses** that capture small molecular movements and multi-conformer basins, then scores **modes** (not single poses) with a soft-β free energy on the CF scoring proxy (slides: “Introducing entropy…” / “Outputs dynamic binding modes”). That is the conceptual bridge:

- **Search** finds low-CF poses (genetic algorithm).
- **Clustering** groups poses into basins (modes).
- **Election / ranking** chooses the preferred mode (and representative) under a free-energy-like objective when entropy is on.

Without clustering, entropy on the ensemble is ill-posed: you cannot define \(\tilde H\), \(\tilde S\) of a *mode* if every pose is its own singleton.

### 1.2 Density-based clustering (DBSCAN family → OPTICS / FastOPTICS)

Peer-reviewed foundation for density-based pose clustering:

| Paper | Venue | DOI | Idea used here |
|-------|--------|-----|----------------|
| Ester, Kriegel, Sander, Xu — **DBSCAN** | KDD 1996 | AAAI proceedings (KDD’96, pp. 226–231); commonly cited as *A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise* | Density-reachable clusters of arbitrary shape; MinPts + Eps; noise points; default **MinPts = 4** (2-D) |
| Sander, Ester, Kriegel, Xu — **GDBSCAN** | *Data Min. Knowl. Disc.* **2**:169–194 (1998) | [doi:10.1023/A:1009745219419](https://doi.org/10.1023/A:1009745219419) | Generalization; for **dim > 2**, rule of thumb **MinPts ≈ 2 · dim** |
| Ankerst, Breunig, Kriegel, Sander — **OPTICS** | *ACM SIGMOD Record* **28**(2):49–60 (1999) | [doi:10.1145/304181.304187](https://doi.org/10.1145/304181.304187) | Ordering + reachability plot; MinPts primary density param; experiments “good results” for **MinPts ∈ [10, 20]**; larger MinPts smooths reachability and reduces single-link chaining |

**Why density-based for docking ensembles:**

- Pose clouds in RMSD / gene space are **non-spherical** (elongated flexible-ligand corridors, multi-substate basins).
- Centroid / \(k\)-medoid schemes force a fixed \(k\) or a fixed RMSD radius and merge distinct basins or split dense wells.
- Density methods mark **noise** (sparse outliers) instead of forcing every pose into a mode — important when GA dumps many high-CF scrap poses.

**FlexAIDdS production entropy arm** uses **FastOPTICS** (`CLUSTA FO`) — an OPTICS-family implementation that builds a reachability ordering then extracts clusters at a single literature MinPts (see §2). Super-cluster helpers on energy projections are **not** a second pose-mode clustering (see `fo_minpts_literature.md`).

### 1.3 CF centroid / RMSD clustering (`CLUSTA CF`)

**CF path** (`LIB/cluster.cpp`) is the historical FlexAID baseline:

- Group poses by structural similarity with a **centroid / representative** logic driven by the CF landscape (and RMSD-style geometry among chromosomes).
- Modes still get member lists and frequencies; when \(T>0\), soft-β **ACF** can re-order emission (see §3).
- **3Dsig red arm A / B0** use CF clustering with **TEMPER 0** (no entropy ranking) — pure FlexAID-style min-CF story for the red “FlexAID” bar.

**Conceptual contrast:**

| Axis | CF centroid (`CLUSTA CF`) | Density FO (`CLUSTA FO`) |
|------|---------------------------|---------------------------|
| Cluster definition | Radius / centroid linkage around CF-competitive poses | Core-density + reachability (OPTICS) |
| Outliers | Often absorbed into nearest cluster | May remain noise / singleton |
| Multi-scale density | One structural threshold | Reachability encodes multi-density (single MinPts extracts one cut) |
| Historical red-bar role | FlexAID (entropy off) | FlexAIDdS (entropy on) |
| Emission names | `prefix_rank.pdb` | `prefix_minPts_rank.pdb` (dual suffix) |

**Important:** density clustering does **not** by itself add thermodynamics. It only defines the sets over which soft-β \(\tilde G\) is computed. Election is a separate layer (§3).

### 1.4 JCIM FlexAID 2015 context

Gaudreault & Najmanovich, *J. Chem. Inf. Model.* **55**:1323–1336 (2015), [doi:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078):

- Comparative docking methodology reused by 3Dsig (10 sims × 2 000 000 evals; bootstrap medians; Astex Diverse / non-native / HAP2-style flexibility stories).
- Establishes FlexAID as competitive when flexibility matters; **does not** introduce the 2017 entropy ranking — that is the 3Dsig deck extension.

Dataset anchors used in 3Dsig figures (cite when reporting bars):

- Hartshorn et al., *J. Med. Chem.* **50**:726–741 (2007), [doi:10.1021/jm061277y](https://doi.org/10.1021/jm061277y) — Astex Diverse (N = 85).
- Bootstrap medians follow Efron, *Biometrika* **68**:589–599 (1981), [doi:10.1093/biomet/68.3.589](https://doi.org/10.1093/biomet/68.3.589) (as cited on the 3Dsig deck).

---

## 2. Why **single** MinPts FastOPTICS (not a ladder) for production / red-pair B

### 2.1 Literature: one MinPts defines one density scale

OPTICS produces an ordering; **cluster extraction** still needs a density scale. Ankerst et al. (1999) treat **MinPts** as the primary density parameter and report a practical band **[10, 20]**. Sander et al. (1998) give **MinPts ≈ 2 · dim** for higher dimensions. Ester et al. (1996) give **MinPts = 4** as a 2-D default floor.

None of these papers mandate re-running the full algorithm at 3–5 MinPts values and unioning results as a *production* binding-mode definition. Multi-scale inspection of the **reachability plot** is for analysis; production needs one reproducible cut.

### 2.2 Engineering / science reasons to forbid the ladder

The legacy “triple MinPts ladder” (fixed 5/7/10 or ×1.5 scales) was a **testing-only** artifact in older trees. It is **forbidden** for production and for 3Dsig red-pair arm B because:

1. **Non-identifiability of modes** — union of clusters across MinPts double-counts basins and invents pseudo-modes that never coexist under one density definition.
2. **Inflated S_top10** — more emitted heads inflate “any of top-10 near native” without a fixed ranking budget.
3. **Irreproducible election** — DatasetRunner dual-suffix election becomes ambiguous when multiple MinPts suffixes represent the same physical basin.
4. **Log contract** — one `[FO-MINPTS]` line and one `Size of Population is K Binding Modes (minPts=N)` per clustering call (see `fo_minpts_literature.md`).
5. **Protocol freeze** — red-pair B CONFIG is `TEMPER 21` + `CLUSTA FO` only; no multi-FO CONFIG knobs.

### 2.3 FlexAIDdS single-MinPts rule (operational)

Implemented in `fo_choose_minpts()` (`LIB/FastOPTICS_cluster.cpp`) with constants in `LIB/ga_constants.h`:

1. \(\mathrm{dim}_{\mathrm{eff}} = \mathrm{clamp}(6 + f_{\mathrm{dih}}, 2, 20)\) (Sander-style effective dimension; ligand dihedrals from `FA->resligand->fdih` or IC gene count).
2. Prefer Ankerst **[10, 20]** when \(N\) is large enough; else Ester floor / feasibility caps (\(N/3\), hard max 50).
3. Sander \(2\cdot\dim_{\mathrm{eff}}\) enters the clamp; diversity softener only if CF diversity &lt; 5% (toward Ester 4).
4. **Exactly one** FastOPTICS + BindingPopulation pass. Super-cluster energy pre-filter / minibatch sampling (if enabled) are **not** extra FO pose clusterings.

**Red-pair B non-negotiable:** single literature MinPts only — `docs/implementation/3dsig_red_pair_protocol.md` §2.1.

---

## 3. Soft-β \(\tilde G = \tilde H - T\tilde S\) election (3Dsig) vs min-CF vs legacy Z+H

### 3.1 3Dsig soft-β free energy on CF (authoritative)

From the 3Dsig deck (“Introducing entropy in FlexAID’s scoring function”), for poses \(i\) in a mode:

\[
Z = \sum_{i \in \mathrm{mode}} e^{-\mathrm{CF}_i / T},
\quad
p_i = \frac{e^{-\mathrm{CF}_i / T}}{Z},
\quad
\tilde H = \sum_i p_i\,\mathrm{CF}_i,
\quad
\tilde S = -\sum_i p_i \ln p_i,
\quad
\tilde G = \tilde H - T\,\tilde S.
\]

| Symbol | Meaning |
|--------|---------|
| \(\mathrm{CF}_i\) | Voronoi contact-function **scoring proxy** (a.u.) — **not** experimental \(\Delta G\) |
| \(T\) | Temperature in K; \(\beta = 1/T\) (FlexAID soft-β), **not** \(1/(k_B T)\) |
| Elect | **Lowest** \(\tilde G\) mode / head |

**Analytic identity** (used in code for stability):

\[
\tilde G = E_{\min} - T \ln Z_{\mathrm{local}} \equiv \mathrm{ACF}
\]

with \(Z_{\mathrm{local}} = \sum_i \exp\bigl(-(E_i - E_{\min})/T\bigr)\). Shared implementation: `LIB/SoftBetaFreeEnergy.h`.

**Intuition (deck slide “Why considering entropy may be a good idea”):** a steep well with one very low CF can lose to a **wider** well with many moderately favorable poses — mode 2 preferred when entropy is counted.

### 3.2 Three election policies (do not mix in claims)

| Policy | Objective | When used | Logs / knobs |
|--------|-----------|-----------|--------------|
| **3Dsig soft-β** | \(\min \tilde G = \tilde H - T\tilde S\) over mode members | Production entropy ranking; red-pair B; DatasetRunner S1 default | `[ENTROPY_RANK]`, `[3DSIG-RANK]`; `FLEXAIDDS_ELECTION_SHANNON_F=1` |
| **min-CF** | Lowest head (or member) CF | FlexAID arm A / TEMPER 0; rollback | `FLEXAIDDS_FORCE_CF_RANK_EMISSION=1` or `force_cf_rank_emission` |
| **Legacy Z+H** | Heuristic \(Z\cdot e^{-\alpha H}\cdot\log1p(N)\)-style (≈ min-CF behavior in practice) | Explicit rollback only — **not** 3Dsig | `FLEXAIDDS_ELECTION_LEGACY_ZH=1` |

**Hard separation (AGENTS.md):**

- GA **search** always optimizes the CF proxy.
- **Ranking / election** may use \(\tilde G\) when entropy is on.
- Physical \(k_B\) StatMechEngine ledgers (F, H, −TS, Cv with true \(k_B\)) are **diagnostic** unless a validated full thermodynamic path is active — they must **not** silently elect S1.

### 3.3 Soft-\(T\) for red-pair B

Arm B freezes **TEMPER 21** (LP-optimized soft-\(T\) for the historical red bar). Override to 298 only with an explicit RUN_RECEIPT note. DatasetRunner soft-\(T\) defaults to dock \(T\) (`FLEXAIDDS_ELECTION_SOFT_T=0` → dock \(T\), else 298).

### 3.4 Vibrational add-on (FlexAIDdS extension)

Classic soft-β configurational \(\tilde G\) is the 3Dsig core. FlexAIDdS BindingMode ranking may add \((-T\cdot S_{\mathrm{vib}})\) from ENCoM/tENCoM **on top** of configurational F when enabled (`docs/classic_entropy_ranking.md`). Document vib ON/OFF in receipts; vib is **not** in the 2017 deck formula and must not be silently conflated with the red-bar \(\tilde G\).

---

## 4. Density Peak (Rodriguez–Laio): what it is; why **not** the historical FlexAID benchmark arm

### 4.1 Algorithm (peer-reviewed)

Rodriguez & Laio, *Science* **344**(6191):1492–1496 (2014), [doi:10.1126/science.1242072](https://doi.org/10.1126/science.1242072):

- Local density \(\rho_i\) from a cutoff kernel (Chi function / Gaussian variants).
- Cluster centers = points with high \(\rho\) **and** large distance \(\delta\) to any higher-density point.
- Remaining points assigned by following the nearest higher-density neighbor.

### 4.2 FlexAIDdS implementation

`LIB/DensityPeak_Cluster.cpp` (`CLUSTA DP`):

- Explicitly cites Science 344:1492–1496, Eq. (1) Chi kernel.
- Neighbor-rate band for cutoff distance \(d_c\) (classic 1–2% neighbor-rate heuristic from the paper’s community usage).
- Can elect the **density peak** as representative (`OUTPUT_CLUSTER_CENTER`), not lowest-CF member — different philosophy from CF and from FO heads after soft-β sort.

### 4.3 Why DP is **not** the historical FlexAID / 3Dsig red-pair entropy arm

| Reason | Detail |
|--------|--------|
| **3Dsig protocol** | Entropy arm = density-based FO modes + soft-β \(\tilde G\); red-pair freezes `CLUSTA FO`, forbids DP as entropy arm (`3dsig_red_pair_protocol.md`) |
| **JCIM 2015 / FlexAID bar** | Baseline FlexAID uses CF-style clustering without the 2017 entropy stack |
| **Different representative** | Peak density ≠ min-CF and ≠ soft-β mode optimum; changes S1/S_top10 semantics |
| **\(d_c\) sensitivity** | Cutoff choice is a free parameter; multiplies MinPts-like knobs and breaks frozen methodology |
| **Role in tree** | Optional algorithm for research / plumbing (`DPFO` pilots); **not** claim red bars |

**Allowed use of DP today:** diagnostic ablations, plumbing tests, optional third engine arm **only if** labeled and never substituted for FO in “FlexAID+entropy” claims.

---

## 5. Actionable mapping to FlexAIDdS code paths

### 5.1 Layer diagram (search → cluster → elect)

```text
GA (gaboom)  --optimizes-->  CF proxy (Vcontacts / vcfunction)
        |
        v
CLUSTA  { CF | FO | DP }   -->  cluster heads + Frequency + .mcf members
        |
        +-- T==0 or force_cf --> min-CF emission (FlexAID bar)
        |
        +-- T>0 classic soft-β --> ACF / BindingMode F  [ENTROPY_RANK]
        |
        v
DatasetRunner pool restarts --> elect min G̃  [3DSIG-RANK]  --> S1 / S_top10
```

See also `docs/architecture/scoring_pipeline_schematic.md`.

### 5.2 File / API map

| Concern | Path | Notes |
|---------|------|-------|
| CF centroid clustering | `LIB/cluster.cpp` | `void cluster(...)`; ACF re-sort when \(T>0\) and not `force_cf` |
| Soft-β math (shared) | `LIB/SoftBetaFreeEnergy.h` | \(\tilde G \equiv\) ACF; used by cluster, BindingMode, DatasetRunner |
| FastOPTICS entry | `LIB/FastOPTICS_cluster.cpp` | `fo_choose_minpts`, single pass, dual-suffix emission |
| MinPts constants | `LIB/ga_constants.h` | `GA_FOPTICS_*` Ankerst/Sander/Ester |
| OPTICS core / super-cluster | `LIB/fast_optics.cpp`, `LIB/FOPTICS.*`, `LIB/cpu_fast_optics.cpp` | Ordering + extraction |
| Density Peak | `LIB/DensityPeak_Cluster.cpp` | Rodriguez–Laio; **not** red-pair B |
| Dispatch by algorithm | `LIB/top.cpp` | `FO` / `DP` / default CF branch on `FA->clustering_algorithm` |
| BindingMode F | `LIB/BindingMode.cpp` | Mode free energy; optional vib |
| Classic emission contract | `docs/classic_entropy_ranking.md` | Rank-0 policy |
| DatasetRunner election | `LIB/DatasetRunner.cpp` | `select_pose_freq_gated_pooled`; dual-suffix FO enumeration |
| Config knobs | `LIB/config_defaults.h`, `LIB/config_parser.cpp`, `LIB/ProtocolConfig.*` | `clustering_algorithm`, `classic_entropy_ranking`, TEMPER |
| Env election flags | (process env) | `FLEXAIDDS_ELECTION_SHANNON_F`, `FLEXAIDDS_ELECTION_LEGACY_ZH`, `FLEXAIDDS_ELECTION_SOFT_T`, `FLEXAIDDS_FORCE_CF_RANK_EMISSION` |
| Red-pair launcher | `scripts/run_3dsig_red_pair_serial.sh` | A → B0 → B |
| CONFIG generation | `scripts/generate_flexaid_inp.py` | Emits `TEMPER` / `CLUSTA` |
| Bootstrap S_top10 | `scripts/bootstrap_3dsig_s_top10.py` | 10k median success |

### 5.3 CONFIG freeze for science arms

| Arm | Deck label | TEMPER | CLUSTA | Ranking |
|-----|------------|--------|--------|---------|
| **A** | FlexAID | 0 | CF | min-CF |
| **B0** | Master CF control | 0 | CF | min-CF |
| **B** | FlexAIDdS | 21 | **FO (single MinPts)** | soft-β \(\tilde G\) |

Matrix pin for red-pair: `MC_st0r5.2_6.dat` MD5 `72d7c7396702331d96ff12d18f831796` (see red-pair protocol).

### 5.4 Emission naming (election plumbing)

| Algorithm | Head PDB pattern | DatasetRunner must |
|-----------|------------------|--------------------|
| CF / DP | `prefix_rank.pdb` | Enumerate single-suffix ranks |
| FO | `prefix_minPts_rank.pdb` | Enumerate **dual-suffix**; one MinPts per run |

Until dual-suffix FO election is validated end-to-end, C0 packaging campaigns are **out of band** for 2017 red-bar claims (prefer classic A/B binaries per protocol).

---

## 6. Open questions / tests still needed

### 6.1 Protocol & packaging

- [ ] End-to-end FO dual-suffix election in DatasetRunner verified on pilot8 (logs: one `[FO-MINPTS]`, dual-suffix heads elect correctly).
- [ ] Arm B logs on live red-pair: exactly one `Size of Population … (minPts=N)` per clustering call.
- [ ] S_top10 + 10 000 bootstrap medians on Astex 85 vs archived targets (FlexAID **0.66**, FlexAIDdS **0.69** from deck labels) — or documented deviation table.
- [ ] Confirm vib OFF for strict 2017 \(\tilde G\) reproduction; if vib ON, label as “FlexAIDdS+vib” not raw 3Dsig.

### 6.2 Clustering science

- [ ] Sensitivity of S_top10 / S1 to `fo_choose_minpts` near Ankerst band edges (10 vs 20) at fixed ensemble — **analysis only**, not multi-MinPts production ladder.
- [ ] Compare FO vs CF clustering **with identical soft-β election** (isolate clustering vs ranking).
- [ ] DP ablation labeled as research: never mixed into red-pair B.
- [ ] Noise fraction under FO: are near-native poses discarded as noise on sparse high-DoF cases?

### 6.3 Soft-β / temperature

- [ ] TEMPER 21 vs 298 ablation on pilot set with fixed FO MinPts (receipt-noted).
- [ ] Identity check \(\tilde G = E_{\min}-T\ln Z\) vs explicit \(H-TS\) in unit tests (`SoftBetaFreeEnergy.h` already encodes both).
- [ ] Frequency gate (`freq > 1` preferred) impact on election_gap (BCR hit, S1 miss).

### 6.4 Correctness tests (existing + gaps)

| Test / tool | Status |
|-------------|--------|
| `tests/test_classic_entropy_ranking.cpp` | Classic ACF emission |
| `scripts/acf_vs_cf_ablation.py` + `tests/test_acf_vs_cf_ablation.py` | Election flip CF vs ACF |
| Soft-β numerical identity tests | Keep green when touching election |
| FO single-pass assertion in CI (grep one `[FO-MINPTS]` per run) | Recommended gate for arm B |
| Dual-suffix DatasetRunner integration test | **Still needed** for C0 claim packaging |

### 6.5 Licensing / provenance

- [ ] Keep all new clustering utilities Apache-2.0; no GPL FlexAID sources as copy-paste inspiration (`docs/licensing/clean-room-policy.md`).
- [ ] Cite JCIM / 3Dsig / OPTICS literature in papers; do not relicense deck content.

---

## 7. Canonical bibliography (verified DOIs / venues)

1. Ester M., Kriegel H.-P., Sander J., Xu X. A density-based algorithm for discovering clusters in large spatial databases with noise. *Proc. KDD’96* (AAAI), 1996, pp. 226–231.  
2. Sander J., Ester M., Kriegel H.-P., Xu X. Density-based clustering in spatial databases: the algorithm GDBSCAN and its applications. *Data Mining and Knowledge Discovery* **2**, 169–194 (1998). [doi:10.1023/A:1009745219419](https://doi.org/10.1023/A:1009745219419).  
3. Ankerst M., Breunig M.M., Kriegel H.-P., Sander J. OPTICS: ordering points to identify the clustering structure. *ACM SIGMOD Record* **28**(2), 49–60 (1999). [doi:10.1145/304181.304187](https://doi.org/10.1145/304181.304187).  
4. Rodriguez A., Laio A. Clustering by fast search and find of density peaks. *Science* **344**(6191), 1492–1496 (2014). [doi:10.1126/science.1242072](https://doi.org/10.1126/science.1242072).  
5. Gaudreault F., Najmanovich R.J. FlexAID: revisiting docking on non-native-complex structures. *J. Chem. Inf. Model.* **55**, 1323–1336 (2015). [doi:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078).  
6. Morency L.-P. The impact of conformational entropy on the accuracy of FlexAID in binding mode prediction. ISMB/ECCB 2017 — 3Dsig (conference deck).  
7. Hartshorn M.J. et al. Diverse, high-quality test set for the validation of protein–ligand docking performance. *J. Med. Chem.* **50**, 726–741 (2007). [doi:10.1021/jm061277y](https://doi.org/10.1021/jm061277y).  
8. Efron B. Nonparametric estimates of standard error: the jackknife, the bootstrap and other methods. *Biometrika* **68**, 589–599 (1981). [doi:10.1093/biomet/68.3.589](https://doi.org/10.1093/biomet/68.3.589).

*Optional secondary datasets referenced on the deck (cite if bars claimed):* Verdonk et al. non-native Astex set (*J. Chem. Inf. Model.* 2008); Gaudreault et al. HAP2 flexibility (*Bioinformatics* 2012).

---

## 8. One-page decision table (for agents)

| Question | Answer |
|----------|--------|
| What groups poses into modes for entropy? | Density FO for FlexAIDdS; CF centroid for classic FlexAID |
| Production FO MinPts? | **Single** literature MinPts (Ankerst + Sander + Ester floor) |
| Triple MinPts ladder? | **Forbidden** in production / red-pair B |
| Rank modes how for entropy claims? | Soft-β \(\tilde G=H-TS\) on CF; lowest wins |
| min-CF election? | FlexAID TEMPER 0 / force_cf rollback only |
| DP in red bars? | **No** |
| CF units thermodynamic ΔG? | **No** — scoring proxy unless full validated ledger |
| Proof logs for arm B | One `[FO-MINPTS]`; one `Size of Population … (minPts=N)`; `[ENTROPY_RANK]` / `[3DSIG-RANK]` |

---

*Document generated for FlexAIDdS literature-swarm synthesis. Do not invent citations beyond §7. Update sibling implementation notes if production MinPts or election defaults change — keep `AGENTS.md` terminology (CF proxy vs true ΔG) intact.*
