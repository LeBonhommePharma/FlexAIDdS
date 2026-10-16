# BPD, Bupropion-Induced Rage & FlexAIDΔS Docking Design

**Research proposal draft — night session 2026-10-16 (04:00–06:00 EDT)**
**Author:** Louis-Philippe Morency (Le Bonhomme Pharma)
**Status:** Pre-experimental. Ligand SMILES curated; target PDB IDs to verify before pose work.
**Intended venue:** Nature Reviews Drug Discovery–style review + FlexAIDΔS benchmark proposal

---

## SLIDE 1 — Clinical Observation & Hypothesis

### Title
**Bupropion-Induced Rage in Undiagnosed BPD: A Natural Experiment in Catecholaminergic Push**

### Background
- Case 1: Mother, 66, **undiagnosed BPD** (suspected by ex-partner MD + self-observation over years). Naltrexone prescribed for alcohol dependence. Bupropion added → marked increase in rage episodes. No BPD-targeted medication in regimen.
- Case 2: 25-year-old medical student, top grades, no prior psychiatric history. Started bupropion 4–6 weeks prior → week-long rage episode aimed at mother and family **while physically separated from them**. Bupropion tapered gradually; escitalopram increased gradually → clinically flawless recovery.

### Hypothesis
BPD aggression = **serotonergic brake failure + catecholaminergic accelerator**.
- Bupropion (DAT/NET inhibitor, zero serotonergic action) pushes dopamine and noradrenaline without restoring 5-HT → rage amplification.
- Escitalopram (SERT inhibitor) restores the brake while the catecholaminergic driver is withdrawn → the combination is the intervention, not escitalopram alone.
- Naltrexone (non-selective opioid antagonist: mu, kappa, delta) may coincidentally modulate the **KOR dysphoria circuit** → dual benefit not planned by prescriber.

### Key Ligand Properties (RDKit-computed, 2026-10-16)

| Ligand | MW | LogP | TPSA | HBD | HBA | RotB |
|---|---:|---:|---:|---:|---:|---:|
| Bupropion | 239.7 | 3.30 | 29.1 | 1 | 2 | 3 |
| Naltrexone | 313.4 | 1.76 | 60.8 | 2 | 4 | 2 |
| Vafidemstat (scaffold approx) | 374.5 | 4.94 | 32.3 | 1 | 2 | 4 |
| Serotonin | 160.2 | 1.67 | 41.8 | 2 | 1 | 2 |
| Dopamine | 153.2 | 0.60 | 66.5 | 3 | 3 | 2 |
| Noradrenaline | 169.2 | 0.09 | 86.7 | 4 | 4 | 2 |
| 5-HIAA | 175.2 | 1.79 | 53.1 | 2 | 1 | 2 |
| Escitalopram | 324.4 | 3.81 | 36.3 | 0 | 3 | 5 |
| DynA(1–5) core | 513.6 | −0.20 | 196.9 | 7 | 7 | 15 |

**Notes:**
- Bupropion: lipophilic enough for easy BBB crossing; MW 239.7.
- Naltrexone: TPSA 60.8 — decent CNS penetration for an opioid antagonist.
- Vafidemstat: very lipophilic (LogP 4.94) — consistent with covalent FAD-pocket binding requiring deep access. **SMILES is scaffold approximation — replace with PubChem CID before pose analysis.**
- DynA(1–5) core: peptide; docking as-is is meaningless. Use active pentapeptide core or peptidomimetic.
- Escitalopram SMILES verified canonical: `CN(C)CCC[C@]1(c2ccc(F)cc2)OCc2cc(C#N)ccc21`.

---

## SLIDE 2 — FlexAIDΔS Docking Design

### Ligands (9)
Bupropion • Naltrexone • Vafidemstat • Serotonin • Dopamine • Noradrenaline • 5-HIAA • Escitalopram • DynA(1–5) core

### Targets (5 PDB structures + 1 off-target)

| Target | PDB | Key residues | Notes |
|---|---|---|---|
| 5-HT2A | 6A93 | Ser3.36, Ser5.46, Phe6.51, Phe6.52 | Psychedelic-bound; orthosteric site open |
| KOR | 4DJH | Asp3.32, Tyr3.33, Val6.55, Ile6.58 | JDTic antagonist-bound; deep pocket |
| DAT | 8Y2Z (cryo-EM) | Asp1.51, Ser1.52, Phe6.51 | Cocaine-site; extracellular vestibule |
| SERT | 5I71 | Asp1.51, Ser4.38, Phe6.51 | SSRI site; similar DAT scaffold |
| MAO-A | 2Z5X | FAD, Tyr4.05, Lys3.05 | Bupropion weak MAO-A inhibitor — off-target check |

### Method
- **FlexAIDΔS:** flexible side-chain rotamers + Shannon entropy collapse scoring.
- Scoring: ΔG = ΔH − TΔS. **Rank by entropy collapse magnitude, NOT raw binding energy alone.**
- **Negative control:** vafidemstat vs 5-HT2A / KOR / DAT / SERT — should NOT bind orthosterically (works upstream via KDM1A/LSD1 epigenetic mechanism). Negative result = positive evidence.
- **Positive controls:** serotonin@5-HT2A, dopamine@DAT, JDTic-like@KOR.
- **Dual-ligand scenario (new, from Case 2):** bupropion OFF DAT + escitalopram ON SERT, scored together for joint entropy collapse on 5-HT2A. Tests whether the clinical rescue (taper bupropion + ramp escitalopram) has a computable polypharmacology signature.

### Predictions
1. **Bupropion:** high DAT affinity, negligible SERT → catecholaminergic push without 5-HT brake. Explains rage in both cases.
2. **Naltrexone:** KOR blockade + residual mu affinity → dysphoria-circuit modulation. Explains possible dual benefit in Case 1.
3. **Vafidemstat:** no orthosteric binding at any target — negative across the board.
4. **Escitalopram:** high SERT affinity; entropy collapse on 5-HT2A may exceed single-ligand bupropion effect on DAT.
5. **Entropy ranking may reorder affinity-based rankings** — the drug that collapses receptor conformational entropy most stabilizes the circuit, regardless of Ki.

### Mechanistic bridge (CSF)
- Measure **5-HIAA in CSF** before/after bupropion in BPD patients; correlate with aggression scores.
- If 5-HIAA drops and aggression rises → serotonergic deficit is the bottleneck, not catecholaminergic excess per se.

---

## SLIDE 3 — Pipeline Context, The Gap & NRDD Hook

### Core claim
For a disorder with **ZERO approved drugs**, efficacy ranking by conformational entropy collapse identifies circuit-stabilizing polypharmacology that single-target affinity scoring has missed for 30+ years.

### Pipeline context (as of 2026-10)
- **Vafidemstat (Oryzon):** only Phase III BPD program. PORTICO Phase 2b complete (data at ECNP 2024); PORTICO-2 Phase 3 protocol submitted to FDA June 2025. Primary endpoint: **aggression** — the largest unmet-need symptom domain, the only domain with zero approved options.
- **KOR PET study (Biological Psychiatry, June 2026):** significantly lower KOR availability in BPD amygdala, ACC, hippocampus, insula, frontal cortex. Correlates with impulsivity and depressive symptoms. Hypothesis: mu/kappa imbalance — excess relative kappa signal.
- **Naltrexone retrospective (Timäus 2021):** OR 43.2 for BPD improvement — never confirmed in RCT. One German Phase 3 on naltrexone for dissociation (completed 2008) unpublished.
- **Network meta-analysis (Molecular Psychiatry, Feb 2026, Gerolymos et al.):** topiramate > aripiprazole > lamotrigine for hostility (high-certainty GRADE). **No drug reduces suicidal behavior.**
- No SSRI, antipsychotic, or mood stabilizer is BPD-approved. DBT remains the only robustly evidence-based treatment.

### Proposed trial (the gap)
- **RCT:** naltrexone vs placebo in **NON-drinking** BPD patients (isolates KOR effect from alcohol effect).
- Primary: self-injurious behavior frequency.
- Secondary: KOR PET before/after; CSF 5-HIAA as mechanistic bridge.
- Mechanistic sub-study: bupropion challenge → 5-HIAA trajectory → aggression correlation.

### NRDD angle
Entropy as a **ranking variable for off-label-to-approved translation** in disorders where symptom-targeted pharmacology has failed for three decades. FlexAIDΔS provides the computational engine: flexible receptors, Shannon entropy collapse, ΔG = ΔH − TΔS.

---

## BACKGROUND NOTES (full session, 2026-10-16)

### BPD primer
Borderline personality disorder: stable pattern of instability in relationships, self-image, emotions; intense fear of abandonment; rapid mood shifts; impulsivity (spending, substances, relationships). Not manipulation — a nervous system that regulates poorly. DBT is the most robust treatment.

### Projection mechanism
When self-image is unstable, there is no fixed reference point for observing one's own emotions → they are perceived in others. Splitting divides the world into all-good/all-bad; if you are not on the good side, the other must be. Every conversation becomes an involuntary diagnosis. 2019 *Personality Disorders* paper: BPD-trait individuals over-attribute negative intentions to ambiguous stimuli (measured, not impression).

### Epigenetics (order matters)
- **NOT:** trauma creates epigenetic marks.
- **YES:** pre-existing epigenetic vulnerability (inherited or installed in early years) makes trauma *possible as reprogramming* rather than a wound that heals. Trauma acts as amplifier that locks the circuit.
- Holocaust survivor studies: altered NR3C1 (glucocorticoid receptor) methylation transmitted to offspring.
- BPD: distinct methylation profiles on stress, serotonin, attachment genes.
- Changing epigenetics = reopening a plasticity window the disorder closed, not erasing trauma.
- **Heritability ~40–50%** (twin studies) — among the most heritable personality disorders. Probability, not sentence.

### Pharmacology landscape
- **No BPD-approved drug.** SSRIs most prescribed (impulsivity/irritability) — modest-to-null effects on core disorder. Atypicals (olanzapine, aripiprazole) target anger/aggression, effect sizes 0.3–0.5. Lamotrigine had a moment for emotional dysregulation; replication cooled it.
- **Endogenous opioid system:** self-injury releases endorphins → auto-reinforcing. Naltrexone tested with mixed-but-non-null results.
- **Vafidemstat:** covalent KDM1A/LSD1 inhibitor (FAD-dependent histone demethylase; removes methyl marks on H3K4me2 activation and H3K9me2 repression — bidirectional transcriptional switch). One daily dose due to near-irreversible FAD binding. Regulates serotonin, dopamine, fear-circuit genes (amygdala, PFC) upstream of receptors — regulation, not symptomatic targeting. Also implicated in long-term synaptic plasticity → may touch trauma epigenetics AND circuit flexibility simultaneously.
- **Kappa opioid receptor (KOR):** receptor of dysphoria. Endogenous ligand dynorphin. Activation → aversion, stress, sadness. Naltrexone blocks mu + kappa + delta. Automutilation as attempt to rebalance: beta-endorphin stimulates mu, counterbalancing kappa.
- **Bupropion warning:** DAT/NET inhibitor, zero serotonergic action. In BPD (already low 5-HT), removes brake without restoring steering. Feb 2026 network meta flagged methylphenidate similarly. **Two independent rage cases now documented (mother 66, med student 25), same drug, same 1–6 week window, rage aimed at closest attachment figures.** Case series candidate. Bupropion also weakly inhibits MAO-A (off-target check in docking design).
- **Bupropion + naltrexone combination:** both lower seizure threshold; bupropion worst in class. Vigilance if seizure history or alcohol use. Liver enzymes for naltrexone (transaminase elevation, dose-dependent).

### Naltrexone mechanism (alcohol)
Blocks mu-opioid receptors. Alcohol → beta-endorphin release → mu activation → removes GABAergic brake on VTA dopaminergic neurons → dopamine rises. Naltrexone cuts the relay. Works on alcohol because alcohol passes through the opioid-dopamine circuit; weak on dependencies that bypass it.

### Why nothing is approved (NRDD framing)
Thirty years of symptom-targeted pharmacology (SSRI → antipsychotic → mood stabilizer) has produced no BPD approval because no single target carries the disorder. The biology is distributed: 5-HT deficit, catecholamine excess, KOR/mu imbalance, HPA-axis epigenetic lock, attachment-circuit instability. **Entropy collapse as a cross-circuit ranking variable is the proposal.**

---

## FLAGS & CORRECTIONS BEFORE EXPERIMENTATION

1. **Vafidemstat SMILES is scaffold approximation** — replace with exact PubChem structure (CID to confirm) before any pose analysis.
2. **PDB IDs to verify as current:** 6A93 (5-HT2A), 4DJH (KOR), 8Y2Z (DAT cryo-EM), 5I71 (SERT), 2Z5X (MAO-A). Check for newer/higher-resolution entries.
3. **DynA(1–5) docking:** use peptidomimetic or active core, not full peptide (TPSA 196.9 — not drug-like).
4. **Case series ethics:** both rage cases are secondhand observation (mother's regimen; med student's course via conversation). Formal case report requires consent + clinical documentation access.
5. **No BPD diagnosis exists in either case** — "undiagnosed BPD" is clinical suspicion, not charted diagnosis. Frame accordingly in any manuscript.

---

## NEXT STEPS (for tomorrow's session)

- [ ] Verify vafidemstat exact structure (PubChem CID)
- [ ] Confirm/update PDB IDs; pull structures into FlexAIDΔS input pipeline
- [ ] Run single-ligand docking: bupropion, naltrexone, escitalopram, vafidemstat × 5 targets
- [ ] Run dual-ligand scenario: bupropion@DAT + escitalopram@SERT joint entropy collapse on 5-HT2A
- [ ] Rank all poses by Shannon entropy collapse magnitude
- [ ] Compare entropy ranking vs affinity ranking — document reordering
- [ ] Draft 5-HIAA/CSF mechanistic sub-study protocol
- [ ] Search published bupropion-aggression case reports (beyond these 2)
- [ ] NRDD-format outline: problem statement → target map → strategy recommendation

---

*Generated 2026-10-16 ~06:00 EDT. Session: Le Bonhomme Pharma × Eve/Grok. All RDKit descriptors computed in-session; literature citations (Timäus 2021, Gerolymos et al. 2026, KOR PET 2026) to be double-fact-checked before manuscript use.*
