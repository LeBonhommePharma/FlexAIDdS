# Deep Science Audit — FlexAIDdS (26h window + current `main`)

**Date:** 2026-07-15  
**Tip audited:** `3b4a53ab1` (`main`)  
**Scope:** Ranking, clustering, election, search budget, 3Dsig/JCIM fidelity, thermodynamic language, success metrics.  
**Not in scope:** Homebrew, iCloud ops, Dependabot (except where they contaminate science provenance).

---

## 1. Executive scientific verdict

| Question | Answer |
|----------|--------|
| Does the **math** of soft-β free energy on CF make sense? | **YES** |
| Is that math **implemented as one identity** across engine + DatasetRunner? | **NO** |
| Can you claim “3Dsig Shannon free energy election” on current defaults? | **NO** |
| Can you claim physical binding free energy ΔG? | **NO** (correctly not, if language is disciplined) |
| Is FO production policy (single literature MinPts) scientifically sounder than the triple ladder? | **YES** |
| Is DoF budget (scale pop, fix gen) the right claim contract? | **YES** |
| Is Astex / 3Dsig red-pair **reproducible** from current pilots as written? | **NOT YET** |

**One-line verdict:**  
The **scientific idea** (search on CF proxy; rank modes by soft-β \(\tilde G=H-TS\) on CF; density-based modes; fixed-gen pop×DoF) is coherent and aligned with classic FlexAID / 3Dsig 2017. The **live packaging on `main`** breaks layer identity (T, members, FO sidecars, dual-suffix), so published success rates from this tip would not be a faithful realization of that idea.

---

## 2. Scientific layer model (must not be conflated)

| Layer | What it optimizes / reports | Units | Elects poses? |
|-------|----------------------------|-------|----------------|
| **L1 Search** | GA fitness = Voronoi **CF** (contact-function scoring proxy) | a.u. | No (explores) |
| **L2 Cluster** | Binding modes (CF / FO / DP) | geometry + CF | Groups poses |
| **L3 Rank (classic / 3Dsig)** | Soft-β \(\tilde G=\tilde H-T\tilde S\) on **CF** | CF a.u., \(T\) in K with \(\beta=1/T\) | **Yes** (emission + S1) |
| **L4 Thermo ledger** | StatMech / BindingMode physical F, H, S, Cv (+ vib) | kcal when \(k_B\) used | **No** (diagnostic unless explicitly wired) |

**AGENTS.md contract:** never sell L1/L3 as experimental ΔG. L3 is a **ranking objective on a scoring proxy**, not a validated free-energy of binding.

---

## 3. Soft-β ranking — theory (PASS)

### 3.1 Formula

For members \(\{E_i\}\) of one mode (CF values):

\[
Z=\sum_i e^{-(E_i-E_{\min})/T},\quad
p_i=\frac{e^{-(E_i-E_{\min})/T}}{Z}
\]

\[
\tilde H=\sum_i p_i E_i,\quad
\tilde S=-\sum_i p_i\ln p_i,\quad
\tilde G=\tilde H-T\tilde S
\]

**Identity (local partition):**

\[
\tilde G = E_{\min}-T\ln Z
\]

This is exactly the **ACF** form in `cluster.cpp` when \(T>0\):

```text
ACF = local_origin − ln(local_z) / β    with β = 1/T = FA->beta
```

**Soft-β:** \(\beta=1/T\) with \(T\) in kelvin as a **score temperature**, **not** \(1/(k_B T)\). CF is not kcal/mol. This is classic FlexAID / 3Dsig poster convention and is **scientifically legitimate as a ranking kernel** if applied consistently.

### 3.2 What the paper / deck requires

From `docs/implementation/3dsig_shannon_ranking.md` (3Dsig 2017 Morency):

- Modes = density-based clusters  
- Elect **lowest** \(\tilde G\)  
- Dataset: Astex Diverse N=85  
- 10 sims × 2e6 evals  
- Primary success: **S_top10** (any of top-10 RMSD &lt; 2 Å) + 10k bootstrap median  
- Arms: FlexAID (CF) vs FlexAID+entropy  

That is a coherent comparative methodology. The repo’s science problem is **implementation fidelity**, not the existence of the formula.

---

## 4. Soft-β ranking — implementation (FAIL identity)

Three code paths must agree. On current `main` they **do not**.

### 4.1 Layer map (current `main`)

| Layer | File | Objective | \(T\) source | Members | Local or global \(Z\)? |
|-------|------|-----------|--------------|---------|-------------------------|
| CF cluster ACF | `LIB/cluster.cpp` | \(E_{\min}-T\ln Z_{\mathrm{local}}\) | `FA->temperature` / `FA->beta` | Cluster members in chrom array | **Local** |
| BindingMode F | `LIB/BindingMode.cpp` | \(H-TS\) (+ vib + NATURaL) | `Population->Temperature` | Mode poses | **Global** \(Z=\) `PartitionFunction` |
| DatasetRunner S1 | `LIB/DatasetRunner.cpp` | \(\tilde G=H-TS\) min | **Hardcoded 298** unless env override | `.mcf` if present else head CF | **Local** over `.mcf` |

### 4.2 Finding S1 — Temperature mismatch (CRITICAL)

```cpp
// DatasetRunner.cpp ~924–929
double soft_T = 298.0;
if (proto.election_soft_T > 0.0)
    soft_T = proto.election_soft_T;
else if (!use_shannon_G && proto.election_score_tau > 0.0)
    soft_T = proto.election_score_tau;
// NEVER reads config.temperature / TEMPER / FA->temperature
```

Docs say: *“dock \(T\) (else 298)”*. Code does: **always 298** unless `FLEXAIDDS_ELECTION_SOFT_T` set.

**Arm B freeze uses TEMPER 21** (LP-optimized soft-β). Engine ranks modes with \(T=21\); DatasetRunner elects S1 with \(T=298\).

At \(T=298\), \(\exp(-\Delta\mathrm{CF}/T)\) is nearly flat for typical CF gaps of a few units → Shannon term is weak → election collapses toward **min-CF**.  
At \(T=21\), the same gaps are sharp → dense basins can beat deep sparse false minima.

**Scientific conclusion:** default S1 on main is **not** the same ranking objective as engine TEMPER-21 entropy ranking. Comparative “entropy helps” claims using DatasetRunner S1 at T=298 are **methodologically invalid** against arm B.

### 4.3 Finding S2 — FO has no member CF sidecar (CRITICAL)

DatasetRunner only gets multi-member \(\tilde S\) from `.mcf`:

```text
Member CFs: .mcf sidecar from cluster.cpp; else head CF only (S̃=0, G̃=CF).
```

`.mcf` is written in **`cluster.cpp`** (CF clustering path).  
**FastOPTICS / BindingMode emission path does not write `.mcf`.**

For FO (the entropy arm):

- Election sees **head CF only**  
- \(\tilde S = 0\), \(\tilde G = \mathrm{CF}\)  
- **Shannon free energy election is inactive** despite `[3DSIG-RANK]` logs  

**Scientific conclusion:** turning on `FLEXAIDDS_ELECTION_SHANNON_F` (default ON) **does not implement 3Dsig ranking for FO modes**. It rebrands min-CF with optional singleton inclusion.

### 4.4 Finding S3 — Global vs local \(Z\) in BindingMode (HIGH)

`BindingMode::compute_enthalpy/entropy` use:

```cpp
p_i = pose.boltzmann_weight / Population->PartitionFunction  // GLOBAL
```

`cluster.cpp` ACF and DatasetRunner soft_free_energy use **local** re-normalization inside the mode/cluster.

With global \(Z\):

- \(H\) and \(S\) for one mode depend on **other modes’** weights  
- \(\tilde G\) is **not** identical to local ACF  
- Dense modes and sparse modes interact through the global partition  

Unmerged branch `fix/softbeta-ranking-identity` correctly moves BindingMode to **local** SoftBeta (identity with ACF). That is the right science fix; it is **not on main**.

### 4.5 Finding S4 — Docs overclaim a shared header (HIGH)

`docs/implementation/3dsig_shannon_ranking.md` states:

> Shared implementation: `LIB/SoftBetaFreeEnergy.h`

**File does not exist on `main`.** Three inline implementations can drift. This is a reproducibility and review hazard.

### 4.6 Finding S5 — Default ON without validation (HIGH)

`FLEXAIDDS_ELECTION_SHANNON_F` defaults **true**. AGENTS.md: ranking changes require tests + explicit intent. Knobs-only ProtocolConfig tests do **not** validate Astex success-rate impact.

Also when Shannon path ON: `include_singletons` is forced ON → candidate pool changes vs freq-gated CF election.

### 4.7 Finding S6 — Numerical form of \(\tilde G\) (PASS with note)

DatasetRunner uses \(H-TS\) with log-sum-exp shift. Algebraically ≡ \(E_{\min}-T\ln Z\). Prefer ACF form for stability (as in SoftBeta header on the feature branch). Not a science bug at normal CF scales.

---

## 5. Clustering science (FO)

### 5.1 Single literature MinPts (PASS direction)

`fo_choose_minpts()` (`LIB/FastOPTICS_cluster.cpp`) implements one pass:

| Source | Rule | Use in code |
|--------|------|-------------|
| Ester 1996 | MinPts floor 4 | `GA_FOPTICS_MIN_POINTS` |
| Sander 1998 | \(\approx 2\cdot\dim\) | `dim = 6 + fdih` (SE(3)+torsions), cap 20 |
| Ankerst 1999 | MinPts ∈ [10,20] | clamp into band when \(N\) large enough |
| Ensemble | MinPts ≤ \(N/3\), &lt; \(N\) | feasibility |

**Why this is better than the triple ladder:** multi-MinPts emission invents **multiple artificial “populations”** of heads for one ensemble, contaminating S_top10 / BCR and making entropy arm non-comparable to a single density cut. Single pass is the scientifically correct production rule.

### 5.2 Caveats (MEDIUM)

1. **\(\dim_{\mathrm{eff}}=6+fdih\) is conformational DoF**, not FO ambient dimension (RMSD embedding uses ligand atom coords). Literature rule-of-thumb is applied to the **intrinsic** problem, not the ambient metric dim — defendable, but not a pure transcription of Sander.  
2. Ankerst band only fully active when \(N/3 \ge 10\) (roughly \(N\gtrsim 30\)). Small pilots (pop=200) may sit below full Ankerst band.  
3. **Diversity softener** (CF distinct ratio &lt; 0.05) is extra-literature; `diversity_ratio == 0` fails to soften (edge bug).  
4. High-\(N\) + high-dim climb can push MinPts above 20 (up to max constant) — still “Ankerst spirit,” not strict [10,20].  
5. Emission remains dual-suffix `prefix_minPts_rank.pdb` — packaging must enumerate it (incomplete on main for BCR).

### 5.3 FO vs DP pilot science (HIGH process)

DPFO small pilot (pop=200, gen=50, ~20k evals) is **not** a scientific ranking of FO vs DP; it is plumbing. Live FO packaging null RMSDs (−1 sentinels) made FO look empty while DP elected — packaging, not density algorithm superiority.

---

## 6. Search budget / DoF science (PASS with documentation debt)

### 6.1 Contract (correct)

Claim path (`FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1`, default):

\[
\mathrm{pop}_{\mathrm{eff}} = \mathrm{pop}_{\mathrm{base}} \times \max(1, n_{\mathrm{flex}}/4),\quad
n_{\mathrm{gen}} = n_{\mathrm{gen,base}}\ \mathrm{(fixed)}
\]

**Rationale:** high-DoF ligands need higher initial diversity \(H(X)\) over chromosomes, not longer trajectories. Premature collapse is a population problem.

Optional `BUDGET_SCALE` multiplies **pop** further for \(n_{\mathrm{genes}}\ge 14\).

Mode `0` = legacy gen-scale (forbidden for claim).  
Mode `-1` = fixed pop+gen (oracle-ceiling only).

### 6.2 Anti-pattern called out correctly

Treating CLI `1000×6000` plus `EVAL_SCALE_DIHEDRAL=-1` as “the claim budget” freezes search and **disables** DoF adaptation. That was a real agent mistake; docs fix (`68063cc9d`) is scientifically right.

### 6.3 Caveats

- Total evals **grow** with DoF (not iso-budget). Cross-engine A/B pilots that freeze 1000×2000 while C0 scales pop are **not effort-matched** on flexible ligands.  
- Receipts often store **base** pop/gen; effective budget lives in `[EVAL-BUDGET]` logs / per-target GA JSON — easy to misreport methods.  
- 3Dsig deck: **2e6 evals / sim**. Base \(1000\times 2000=2\mathrm{e}6\) matches; pop-scale **exceeds** 2e6 for flexible ligands — acceptable if documented, not if labeled “exact 3Dsig budget.”

---

## 7. TEMPER / soft-β temperature as a science parameter

### 7.1 Semantics

`TEMPER N` → `FA->temperature = N`, \(\beta=1/N\). Used for:

- Soft-β sampling (SMFREE) when enabled  
- ACF emission order  
- BindingMode classic ranking  

It is **not** “simulate at physiological 310 K with \(k_B\)”.

### 7.2 Freeze values

| Protocol era | Arm B TEMPER | Effect |
|--------------|--------------|--------|
| Early three-engine v1.0 | **298** | Softmax nearly flat → entropy ranking ≈ CF |
| Current 3Dsig red-pair / C0 B | **21** | Sharp soft-β; LP-optimized ranking temperature |

**Science note:** TEMPER 21 is a **hyperparameter of the ranking kernel**, not a claim of true thermodynamics at 21 K. Must be labeled that way in papers. Switching 298→21 mid-campaign without a protocol version is a **methods change**, not a “bugfix.”

---

## 8. 3Dsig / JCIM comparative fidelity

### 8.1 What the deck requires vs what the tree does

| Requirement | Status on main / pilots |
|-------------|-------------------------|
| Soft-β \(\tilde G\) ranking | Formula present; **T + FO members broken** for DatasetRunner S1 |
| Density modes (FO), not DP for entropy arm | Protocol correct; packaging incomplete |
| Astex Diverse 85 | Canonical dataset docs + apo strip strong |
| 10 sims × 2e6 evals | C0 gen=2000 base OK; R often 5 not 10; pop-scale alters evals |
| S_top10 + 10k bootstrap | **Not closed** in pilot parsers |
| Matrix pin JCIM-era | Protocol pins MD5 `72d7…` — good if enforced |
| No native seed for claim | Cognate no-seed stack improved; engine defaults still seed-capable |

### 8.2 Success metric science

| ID | Definition | Role |
|----|------------|------|
| **S_top10** | Min RMSD among top-10 ranked modes ≤/&lt; 2 Å | **3Dsig primary** |
| **S1** | Elected rank-0 RMSD ≤ 2 Å | Modern claim KPI (stricter) |
| **S2** | S1 ∧ PoseBusters | Modern secondary (AGENTS success) |
| **S3 / BCR** | Best cluster-head RMSD ≤ 2 Å | Sampling ceiling diagnostic only |

**Tensions:**

1. AGENTS: success = RMSD **and** PoseBusters → headline **S2**, not S1 alone.  
2. Deck: **S_top10**, not S1.  
3. Threshold **&lt;2.0** vs **≤2.0** not unified.  
4. Reporting S1/BCR as S_top10 is scientific fraud by metric substitution (even if accidental).

Admission aggregator after `9dbbd9fa9` is **sound for S1/S2/S3 claim gates** (fail-closed seeds, finite RMSD over flags). It does **not** implement S_top10 bootstrap.

---

## 9. Seed / oracle / fairness (HIGH for claim purity)

### 9.1 Oracle native seed (`26cb99276`)

Restoring gene-0 seed when `reference_ligand.file` is set is **correct for oracle-ceiling** diagnostics (was broken).  

**Risk:** engine defaults remain seed-friendly (`pose_seed_enabled`, non-zero seed fraction). Claim safety depends on **orchestration** setting seed_fraction=0 and empty reference for flood — fail-open if a script forgets.

### 9.2 No-seed cognate stack (`dea70ea88`)

Scientifically good direction: no crystal IC flood, scramble orientation, site-centered. Residual: **bound torsions** may still come from cognate ligand construction (pose-blind ≠ conformation-blind). Label methods as **cognate redock with random orientation**, not de novo pose generation.

### 9.3 Seed elitism / `_INI.pdb`

Seed-anchored elitism can put crystal `_INI.pdb` into the election pool. For **oracle** that is intentional; for **claim** it is contamination if enabled. Receipts have dual-truth on seed_elitism vs runtime override — provenance risk.

---

## 10. CF naming vs free energy (PASS with residual debt)

Positive: `best_score` documented as CF proxy; README overclaim “accurate ΔG” removed.  

Residual: CSV columns still named `predicted_dG`; affinity paths may map CF through −dG/RT → pKd theater; Python examples may still say “Best ΔG.”

**Science rule:** REMARK `free_energy` from StatMech is a **ledger quantity** under the model assumptions; REMARK CF / soft_beta G are ranking/search quantities. Never merge them in prose.

---

## 11. Apo integrity (PASS for Astex cognate strip)

Astex apo strip validation: **0 residual cognate ligand atoms** on 85 targets; apo≡deposit for 83/85 (ligand already outside apo for most). Strict gate still soft on missing files (ops). Science of receptor preparation for redock is **acceptably clean** for Astex Diverse cognate ligands.

Ligand-centered site PDBs (20 worst GetCleft outliers) are **scientifically correct** (6 Å shell of cognate SDF; crystal serial fidelity). They fix site definition, not ranking math.

---

## 12. What is *not* scientific free energy (do not overclaim)

Even when L3 works perfectly:

1. CF is a **Voronoi contact scoring proxy**, not MM/GBSA or experimental ΔG.  
2. Soft-β \(T\) is a **kernel temperature**, not thermodynamic temperature with \(k_B\).  
3. Shannon \(\tilde S\) over discrete GA poses is **not** full configurational entropy of continuous phase space.  
4. StatMech ledger F/S/Cv over the same ensemble is a **model free energy of the discrete sample**, not validated ITC ΔG.  
5. tENCoM / ENCoM vib terms are **corrections** with their own model assumptions.  
6. PoseBusters is a **pose quality** filter, not a free-energy proof.

Correct claim language:

> “Entropy-aware ranking of CF-scored poses via soft-β free energy of binding modes; success = RMSD (and PoseBusters) against crystal.”

Incorrect:

> “Computed binding free energy / true ΔG improved by entropy.”

---

## 13. Cross-layer failure mode (the load-bearing diagram)

```text
                    SEARCH (CF) ──────────────────────────────┐
                         │                                    │
              ┌──────────┴──────────┐                         │
              ▼                     ▼                         │
         CLUSTA CF              CLUSTA FO                     │
              │                     │                         │
         ACF local T=TEMPER    BindingMode F                  │
         + .mcf written        global Z + TEMPER              │
              │                     │ no .mcf                 │
              └──────────┬──────────┘                         │
                         ▼                                    │
              DatasetRunner S1 election                       │
              soft_T=298 default                              │
              G̃ from .mcf or CF only  ◄── FO: G̃=CF always  │
                         │                                    │
                         ▼                                    │
              S1 / S_top10 / claim tables                     │
                         │                                    │
              ≠ engine TEMPER-21 FO ranking                   │
              ≠ 3Dsig deck identity                           │
```

---

## 14. Prioritized science remediation

### P0 — Without these, no 3Dsig / entropy claim

| # | Fix | Why |
|---|-----|-----|
| 1 | Wire DatasetRunner `soft_T` ← dock `config.temperature` / TEMPER | T identity |
| 2 | Write `.mcf` (or equivalent member CF list) on FO / BindingMode emission | Shannon not vacuous |
| 3 | Merge shared SoftBeta (local \(Z\)) into cluster + BindingMode + DatasetRunner | One objective |
| 4 | Complete FO dual-suffix enumeration for election **and** BCR | Modes exist on disk |
| 5 | Default Shannon S1 **OFF** until Astex pilot at fixed protocol | Ranking change discipline |
| 6 | Emit top-10 mode RMSDs; bootstrap true S_top10 | Deck metric |

### P1 — Protocol honesty

| # | Fix |
|---|-----|
| 7 | Single freeze document: TEMPER B, gen, R, matrix MD5, eval_scale, soft_T |
| 8 | Receipts store **effective** pop/gen/T/MinPts/soft_T |
| 9 | Unify &lt;2.0 vs ≤2.0 |
| 10 | Fail-closed claim if FO log lacks single `[FO-MINPTS]` or shows ladder |
| 11 | Effort-match A/B vs C0 on high-DoF or label mismatch |

### P2 — Hardening

| # | Fix |
|---|-----|
| 12 | Engine default seed off for non-oracle |
| 13 | Unit tests: SoftBeta identity, T wiring, FO .mcf, dual-suffix BCR |
| 14 | Purge residual “ΔG” language on CF columns in user-facing outputs |

---

## 15. What you *can* defend scientifically today

If methods section is precise:

1. **Search** optimizes Voronoi CF (proxy).  
2. **Modes** from single-pass FO with literature-inspired MinPts (Ankerst/Sander/Ester composite heuristic).  
3. **Engine emission** at TEMPER&gt;0 uses local ACF soft-β (CF clustering path) or BindingMode classic F (with global-Z caveat).  
4. **DoF budget** scales population, fixed generations (claim path).  
5. **Apo receptors** for Astex Diverse are ligand-clean under residual-atom gate.  
6. **S1/S2/S3** definitions are enforceable in the claim aggregator (after seed fail-closed).  

You **cannot** yet defend:

- “DatasetRunner S1 = 3Dsig Shannon free energy on FO modes”  
- “Entropy arm and CF arm differ only by TEMPER/entropy ranking with matched election”  
- “Reproduced 3Dsig red-bar medians from this tree’s live packaging”  
- “Computed thermodynamic binding free energy”

---

## 16. Bottom line

The last 26 hours of science work **pointed the codebase in the right theoretical direction** (soft-β on CF, single FO MinPts, pop×DoF, 3Dsig protocol freeze language).  

The **implementation does not yet close the scientific loop**. The ranking identity across  
`cluster ACF` ↔ `BindingMode F` ↔ `DatasetRunner S1`  
is broken on temperature, partition locality, FO member energies, and FO packaging. Until P0 is closed and revalidated on Astex, treat entropy success claims from this tip as **unproven** and prefer classical FlexAID arm A/B0 binary paths with explicit TEMPER/CLUSTA receipts for any figure that must “make sense.”

---

## 17. Evidence anchors (code)

| Claim | Anchor |
|-------|--------|
| Soft-β formula + min G̃ | `LIB/DatasetRunner.cpp` ~904–1032 |
| soft_T = 298 default | `LIB/DatasetRunner.cpp` ~924–929 |
| `.mcf` only from cluster | `LIB/cluster.cpp` ~429+; DatasetRunner comment ~918 |
| ACF local LSE | `LIB/cluster.cpp` ~164–184 |
| BindingMode global Z | `LIB/BindingMode.cpp` ~296–350 |
| FO MinPts single pass | `LIB/FastOPTICS_cluster.cpp` `fo_choose_minpts` |
| Pop×DoF | `LIB/DatasetRunner.cpp` ~5582–5627 |
| SoftBeta header missing | no `LIB/SoftBetaFreeEnergy.h` on main |
| 3Dsig contract | `docs/implementation/3dsig_shannon_ranking.md` |

**Related swarm notes:** `docs/audit/26h-swarm/SUMMARY.md`, per-commit reports for `c82e6fc24`, `6ec671a92`, `0e39f3a0b`, `68063cc9d`, `dea70ea88`, `26cb99276`, `4e87c0b3c`.
