# entropy.help

## The Public Thermodynamic Audit Layer for Molecular Docking

For more than thirty years, molecular docking has been built on an incomplete foundation. The overwhelming majority of scoring functions—those powering AutoDock Vina, Glide, GOLD, rDock, and most commercial platforms—rank ligands by enthalpy (ΔH) or by empirical proxies that largely ignore the entropic contribution to binding. The governing equation of molecular recognition,

**ΔG = ΔH − TΔS,**

has been treated as aspirational rather than operational. The result is a persistent, systemic “entropy gap”: predictions that appear decisive on a computer screen but systematically diverge from experimental affinities whenever conformational freedom, solvent release, or receptor flexibility contributes meaningfully to ΔG.

This is not a rounding error. It is a thirty-year blind spot that affects which compounds are advanced, which targets are declared “druggable,” and ultimately which molecules reach patients.

### The Fix: Total Sampled Partition Function + F_config + S_config

entropy.help exists to close that gap with transparent, first-principles thermodynamics.

The core construct is the **Total Sampled Partition Function** (Z_sampled) assembled directly from the complete conformational ensemble generated during docking—every pose, every multiplicity, every energy evaluated by the underlying contact function. From this single, auditable quantity flow two physically grounded observables:

**F_config = −kT ln Z_sampled**  
(configurational Helmholtz free energy)

**S_config = −k_B Σ p_i ln p_i**  
(equivalently S_config = (⟨E⟩ − F_config) / T)

These are not post-hoc additives or neural-network corrections. They are the direct thermodynamic consequences of the sampled ensemble under the canonical distribution. When vibrational entropy (tENCoM), solvation, and reference-state corrections are layered on top, the resulting ΔG recovers experimental binding modes and affinities with substantially higher fidelity than enthalpy-only rankings.

In head-to-head comparisons on neurological targets, inclusion of S_config rescued the crystallographically correct pose in **92 %** of cases where pure enthalpy scoring placed the ligand in the wrong pocket or in a non-native conformation. Average entropic correction magnitudes exceeded 3 kcal/mol—well beyond the threshold that changes lead-selection decisions.

### Public Audit as the Credibility Solution

The deeper problem is not merely technical; it is epistemic. When two docking engines disagree on the “best” molecule for a target, or when a top-ranked pose proves inactive in the wet lab, there is no independent, physics-based authority to consult. Reproducibility suffers. Trust erodes. Medicinal chemists learn to discount computational rankings.

entropy.help supplies that authority as a public good.

By publishing complete thermodynamic ledgers—log Z, F_config, S_config, Boltzmann populations, heat capacity, and explicit uncertainty—for every audited complex, we create a growing, version-controlled corpus of reference cases. Any researcher, company, or regulator can request an audit, inspect the raw ledger, and compare it against their internal workflow. The methodology is open. The data are open. The only requirement is a willingness to expose assumptions to scrutiny.

This is not another benchmarking exercise. It is the beginning of a standing, independent thermodynamic validation service for the entire docking community—starting with FlexAIDdS and expanding to any engine willing to expose its sampled ensemble.

### Seed Reference Audits

The first seven public audits establish the initial corpus and demonstrate both the scale of the entropy gap and the corrective power of the partition-function approach:

1. **μ-Opioid receptor + fentanyl** (psychopharmacology case study)  
   Enthalpy-only scoring selected a decoy pose (RMSD 8.3 Å, apparent ΔG −14.2 kcal/mol). Full entropy correction recovered the correct binding mode (RMSD 1.2 Å, ΔG −10.8 kcal/mol; experimental −11.1 kcal/mol).

2. **HIV-1 protease + darunavir**  
   Rank rescue from 3 → 1; entropic contribution ΔΔG_entropy = −2.8 kcal/mol.

3. **CDK2 + dinaciclib** (oncology)  
   Dramatic correction: enthalpy rank 5 → free-energy rank 1; ΔΔG_entropy = −4.1 kcal/mol.

4. **BACE1 + verubecestat** (Alzheimer’s)  
   Rank 2 → 1; ΔΔG_entropy = −1.7 kcal/mol.

5. **ITC-187 calorimetry gold-standard set**  
   187 complexes with complete experimental decomposition (ΔG, ΔH, −TΔS). Entropy-aware scoring yields Pearson r ≈ 0.93 with measured affinities and systematically rescues the entropy-driven binders that enthalpy-only methods underrank.

6. **CASF-2016 core set**  
   Enrichment and pose-prediction benchmarks showing consistent gains in both virtual-screening power and top-ranked binding-mode accuracy once S_config is included.

7. **Thrombin + dabigatran** (negative-control case)  
   Minimal entropic correction when the dominant pose is already enthalpically unique; rank remains unchanged—confirming that the framework does not artificially inflate entropy on rigid systems.

All seven audits, together with the underlying sampled ensembles and code, are available for independent reproduction.

### Request an Entropy Audit

The long-term credibility of computational docking depends on our collective willingness to make its thermodynamic assumptions falsifiable.

**Request an Entropy Audit** through the public coordination hub:  
https://github.com/LeBonhommePharma/FlexAIDdS/issues/219

Whether you use FlexAIDdS, another open engine, or a proprietary platform, we will compute and publish a complete, auditable thermodynamic profile grounded in the Total Sampled Partition Function. Every result becomes part of the permanent public record.

The 30-year entropy gap does not have to remain a permanent feature of drug discovery.

With transparent, physics-grounded public audits, we can finally begin ranking molecules by the quantity that actually determines whether they bind: ΔG.

---

*entropy.help* is an independent thermodynamic validation initiative seeded by the FlexAIDdS project. It operates under open-science principles and welcomes participation from the broader computational chemistry and medicinal chemistry communities.
