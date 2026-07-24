# Scoring in FlexAID∆S: CF, ΔG_eff, and the Thermodynamic Engine

This document gives a complete account of the cost function, the NRGDock energy matrix, the thermodynamic ensemble scoring layer, and the physical filters applied to pose ranking in FlexAID∆S. The intended audience is a computational chemist or structural biologist who wants to understand exactly what the numbers mean.

---

## 1. The Voronoi Contact Function (CF)

### 1.1 What CF Measures

The core scoring proxy in FlexAID∆S is the **Voronoi contact function** (CF), derived from the Voronoi tessellation of atomic surfaces. Each heavy atom is assigned a Voronoi polyhedron; the contact area between two atoms is the area of the face shared by their polyhedra. CF is a sum over all intermolecular atom-pair contacts:

```
CF.com = Σ_{i ∈ ligand, j ∈ receptor} E(type_i, type_j) · A_{ij} / SAS_i
```

where `A_{ij}` is the Voronoi contact area between atoms `i` and `j`, `SAS_i` is the solvent-accessible surface of atom `i` (used as a normalization denominator), and `E(type_i, type_j)` is the entry in the 40×40 NRGDock energy matrix for the pair of atom types. Lower CF is better: a large negative CF.com means the ligand is making favorable, complementary contacts with the receptor across many well-matched atom-type pairs.

### 1.2 CF Component Breakdown

The full CF is a sum of terms:

```
CF = CF.com + CF.wal + CF.sas + CF.elec + CF.hbond + CF.pb_clash + CF.con
```

**CF.com** — Voronoi contact complementarity. The dominant term. Negative values indicate favorable packing; positive values indicate mismatched contacts (e.g., two donors forced into proximity without a bridging water).

**CF.wal** — Soft-wall steric repulsion. Atoms that overlap (interatomic distance < sum of vdW radii) accumulate a wall penalty. The default path uses a hard r^−12 potential, capped at **50 CF units per contact** (`WAL_CONTACT_CAP = 50.0` in `soft_wall.h`). The cap is critical: without it, a single deeply clashing atom pair can produce CF.wal values of 10^6 or more, which then overwhelm the partition function and cause the SIGSEGV that afflicted early versions on the 1M2Z Astex target. When `FLEXAIDDS_SOFTCORE_WAL` is set, a smooth parabolic plateau replaces the hard r^−12 below the contact floor, allowing the GA to escape clashes by gradient rather than by jumping discontinuously. When `FLEXAIDDS_WAL_COERCIVE` is set, the 50-unit cap is removed entirely, allowing deep clashes to dominate CF.com over-packing.

**CF.sas** — Solvent-accessible surface penalty. Ligand surface not in contact with the receptor contributes a desolvation cost (exposing non-polar surface to solvent is unfavorable relative to burial in a hydrophobic pocket).

**CF.elec** — Optional Coulomb electrostatic term, gated in the configuration. When active, uses the KCOULOMB = 332.0637 kcal·Å/(mol·e²) constant.

**CF.hbond** — Hydrogen-bond term. Evaluated by `HBondEvaluator`; default weight −2.5 (`FLEXAIDDS_HBOND_WEIGHT`). A negative value rewards correctly directed H-bonds; the weight is negative because H-bonds stabilize binding.

**CF.pb_clash** — PoseBusters intermolecular clash penalty. All-pairs check at a distance threshold of 0.75 × (vdW_i + vdW_j) (`pb_clash_ratio = 0.75`). Any atom pair below this threshold contributes `pb_clash_weight × severity^pb_clash_exponent` to CF. This term is **uncapped** by design — it must be able to overcome CF.com over-packing, unlike the WAL term which is capped. The receptor clash grid is built **once per dock session** (hoisted out of the per-eval loop) because the receptor is rigid; this hoist reduced wall-clock time for clash checking from O(N_evals × N_rec) to O(1 + N_evals × N_lig). Disabled by default (weight = 0.0); set `FLEXAIDDS_PB_CLASH_WEIGHT=1.0` to enable.

**CF.con** — Distance constraint penalty. Applied when constraints are defined in the configuration (e.g., pharmacophoric distance restraints). Each violated constraint adds KDIST to CF.

### 1.3 Implementation Notes

The Voronoi tessellation is computed by `Vcontacts()` in `Vcontacts.cpp`. The `vcfunction()` routine in `vcfunction.cpp` accumulates all CF terms per atom contact. The contacts buffer is cleared per-eval (default: O(MAX_ATOM_NUMBER) memset); `FLEXAIDDS_CONTACTS_EPOCH` replaces this with an O(1) epoch stamp to avoid the memset overhead on large receptors.

---

## 2. The NRGDock Energy Matrix

### 2.1 What the Matrix Encodes

The 40×40 matrix `MC_st0r5.2_6.dat` (loaded as `kEnergyMatrix` from `nrgrank_matrix.h`) is the statistical potential that converts Voronoi contact areas into binding scores. Each entry E(i, j) encodes the log-odds ratio of observing atom-type pair (i, j) in contact in a real PDB binding site versus the background frequency of that pair contact. Pairs that appear more frequently in real binding sites than expected at random receive negative entries (stabilizing); rare pairs receive positive entries (destabilizing).

The matrix is symmetric and was derived from a non-redundant set of PDB protein–ligand complexes, grouped by the 40-type SYBYL atom-type system. A "live" entry is one with a non-zero, non-trivial value fitted from sufficient PDB statistics. A "dead" entry (all zeros or with only 2–3 supporting contacts) provides effectively no signal.

### 2.2 How Contact Areas Are Weighted

The y-value lookup `get_yval(energy_matrix, area / SAS_i)` evaluates a piecewise function of the normalized contact area. As the fractional contact area `A_{ij}/SAS_i` grows (from a grazing contact to a fully buried face), the energy contribution changes continuously. Large contact areas between well-matched types (e.g., aromatic C–aromatic C) yield large stabilizing contributions; the same large areas between mismatched types (e.g., H-bond donor next to a non-polar atom) yield large destabilizing contributions.

This is what makes CF sensitive to shape complementarity *and* to electrostatic/chemical complementarity simultaneously: the Voronoi geometry provides the geometric weighting, and the energy matrix provides the chemical weighting.

### 2.3 Dead Rows and the Atom-Type Fixes

Several rows in `MC_st0r5.2_6.dat` are effectively dead — they were derived from too few PDB contacts to have meaningful fitted values, so interactions involving those atom types score near-zero regardless of the contact area. When a ligand atom is mis-typed into a dead row, that atom contributes nothing to CF.com, invisibly suppressing the discriminatory power of the score for that entire chemical group.

The key dead row is **type 8 (N.3)**: the aliphatic amine / sp3 nitrogen type. Almost no PDB entries were classified as N.3 contacts (the nitrogen was either protonated to N.4 or perceived as N.am), so this row accumulated fewer than a handful of live entries. Any sp3 nitrogen scored as type 8 would appear as a ghost — present geometrically but contributing zero energy. The fix maps N.3 → N.am (type 11), which has a well-populated row reflecting amine contacts as they actually appear in PDB structures.

The analogous problem with **type 26 (I, iodine)**: this row has only 3 live entries across 40 columns. Iodine is rare in PDB binding sites, so the matrix is nearly empty for type 26. Mapping I → BR (type 25) uses the bromine row, which is well-populated because bromine is far more common in PDB ligands and electronically similar to iodine (both are large polarizable halogens in the same group of the periodic table). This is a pragmatic approximation; a future version of the matrix would require a larger iodine-containing PDB training set to populate type 26 directly.

For **N.2 → N.ar** (type 10): the sp2 imine nitrogen (N.2 in SYBYL) is an H-bond *acceptor* — the lone pair is available for donation to a protein backbone NH or a charged residue. N.am (type 11), which was the previous mapping, is an H-bond *donor* (the NH of an amide). Mapping an acceptor nitrogen to a donor type reverses the chemical interpretation of the interaction and inverts the sign of the nearest-neighbor contacts in the matrix. The correct type is N.ar (type 10, aromatic/sp2 nitrogen), whose row was derived from pyridine-like nitrogens — the dominant chemical context for N.2 in drug-like molecules.

For **C.1 → C.2** (type 2): C.1 (sp carbon, type 1) has 10 live entries in the matrix but represents a rare functional group (alkynes, allenes, nitriles). C.2 (sp2 carbon, type 2) is broadly similar in polarity and packing geometry and has a richer training set. The fallback to type 2 is justified on both coverage and chemical grounds.

---

## 3. ΔG_eff: The Ensemble Free Energy

### 3.1 Motivation

The GA optimizes CF and converges a population toward the low-CF region of pose space. At the end of a run, that population is not a single pose but a distribution: many structurally similar low-CF members near the converged basin, plus some outliers from earlier generations. The question is: how should this distribution be collapsed to a single ranking score for the binding mode?

The classical answer is: take the best (lowest CF) member. This answer is wrong in two important ways. First, it discards information about the width of the basin — a narrow funnel and a broad flat plateau can produce the same minimum CF but have very different thermodynamic meanings. Second, it is noisy — the single-lowest member is sensitive to sampling fluctuations, particularly for flexible ligands with many near-degenerate torsion combinations.

ΔG_eff addresses both problems by treating the GA population as an empirical sample from a Boltzmann distribution and computing a free energy directly.

### 3.2 Formal Derivation

Define the Boltzmann weight of pose i in the converged population:

```
P_i = exp(−CF_i / T_eff) / Z
Z   = Σ_j exp(−CF_j / T_eff)
```

where T_eff is an effective temperature in CF units (default: 0.596). Then:

```
⟨CF⟩ = Σ_i P_i · CF_i       # Boltzmann-weighted mean CF (enthalpy proxy)
H     = −Σ_i P_i · ln P_i   # Shannon entropy of the pose distribution (nats)
ΔG_eff = ⟨CF⟩ − T_eff · H
```

This is formally analogous to the Helmholtz free energy F = ⟨E⟩ − TS, with T_eff playing the role of temperature and H (Shannon entropy of the pose distribution) playing the role of thermodynamic entropy S. The CF plays the role of energy.

The term −T_eff · H is negative (since H ≥ 0 and T_eff > 0), so it lowers ΔG_eff relative to ⟨CF⟩. A broader population (high H) gets a larger entropy bonus, but only if that breadth is real sampling diversity around a genuine basin. Critically, if the population is broad because it has converged to a false minimum — a flat region of CF space — then ⟨CF⟩ will also be high (poorly complementary), and the entropy bonus will not be large enough to overcome the large ⟨CF⟩. The two terms compete in exactly the right way to discriminate genuine narrow funnels from false flat plateaus.

### 3.3 The Two Calibrations

ΔG_eff is computed at two temperatures:

**T_eff = 0.596 (CF units)** — the scoring temperature. This is calibrated to the FlexAID CF scale: at T_eff = 0.596, the Boltzmann distribution over the GA population is genuinely dispersed (not collapsed to a near-delta function) while still being sensitive to the ~10 CF-unit range that separates good from poor poses on the Astex benchmark. This value was empirically tuned; at lower T the distribution collapses and ΔG_eff → min CF; at higher T the distribution flattens and ΔG_eff → ⟨CF⟩ uniformly.

**T = 21 (internal CF units, "kT_ISMB")** — the ISMB 2017 reporting calibration. This is a broader temperature at which the Boltzmann pose distribution is genuinely flat and the Shannon entropy H is close to log(N) for N population members. The T=21 pair (mean_CF_T21, dG_eff_T21) are *diagnostic reporting quantities only* — they do not enter G_bind or GA selection. They appear in output as `_T21` suffixed fields and in the ThermoWhiteboard I_ES, binding_regime, and CF_r2s diagnostics. Per the whiteboard, T=21 is baked into the *definition* of the left-hand quantities (ΔG₂₁, P_i(T=21)) — it is not a substitution; it names the specific observable.

### 3.4 G_bind: The Full Thermodynamic Free Energy

ΔG_eff is the *pose-population* free energy (configurational). The full binding free energy reported as `G_bind` additionally includes tENCoM vibrational entropy:

```
G_bind = T_eff · ⟨CF⟩_raw  −  TdS_shannon  +  TdS_vib
```

where:
- `T_eff · ⟨CF⟩_raw` — T_eff-weighted extensive enthalpy proxy (using the raw, un-normalized CF mean)
- `TdS_shannon` = T_eff · H — configurational entropy cost of the pose distribution (positive = unfavorable; a broad distribution costs entropy)
- `TdS_vib` = tencom_scale × (H_rep_bound − H_rep_ref) — vibrational entropy change upon binding, from the tENCoM elastic-network model. A negative value means the bound complex is *more rigid* than the separated receptor + ligand (entropy loss on binding); a positive value means binding softens the receptor (entropy gain). Physical values are typically in the range −5 to +5 nats; values outside this range are clamped.

The sign convention follows the standard thermodynamic convention: ΔG = ΔH − TΔS, where TΔS_shannon is the *cost* of configurational ordering (positive for a broad distribution that narrows upon binding) and TΔS_vib is the vibrational contribution (sign from actual receptor rigidification: negative if the receptor becomes more rigid, i.e., ΔS_vib < 0 means entropy loss, making it +TdS_vib in the G expression since G = H − TΔS).

This version of G_bind was restored from v88 (91.7% Astex-BCD) after a regression in v100 that inadvertently (1) attenuated the VCT signal by dividing by n_heavy (~12× attenuation), (2) flipped the sign of TdS_shannon, and (3) flipped the sign of TdS_vib, reducing accuracy to 9.4%.

### 3.5 Thermodynamic Impossibility Gate

When `FLEXAIDDS_THERMO_SCORE=1`, a physics filter is applied before election:

```
if ΔH > 0  AND  ΔS < 0:
    ΔG_eff ← +1000  (sentinel)
```

The reasoning is thermodynamic: from ΔG = ΔH − TΔS, if ΔH > 0 (endothermic) and ΔS < 0 (entropy loss), then −TΔS > 0 for every T > 0, making ΔG strictly positive at all physically realizable temperatures. Such a binding configuration cannot be spontaneous. Rather than propagating a nonsensical score through the clustering election, the gate assigns a large positive sentinel (+1000 CF units) so downstream clustering can never elect this configuration rank-0.

The ΔS source for this test is `TdS_vib` (the tENCoM vibrational entropy), because the population Shannon entropy H is always ≥ 0 by definition (each term −P_i ln P_i is non-negative for 0 < P_i ≤ 1), which would make the ΔS < 0 branch of the gate unreachable if H were used. TdS_vib is the only ΔS contribution in the engine that takes negative values (observed range: −1.86 to −1.95 nats on 1SG0/2GBP/1OF1 in test runs), so it is the physically meaningful discriminator for this gate.

The gate is a **zero-cost no-op when disabled**: the flag is read once (function-local static) and short-circuits before any per-pose iteration. In default benchmarking mode, the gate is off and the sentinel is never assigned.

---

## 4. T_eff = 0.596: Calibration

The value T_eff = 0.596 deserves explicit justification because it is not a physical temperature in kelvin — it is in the CF scoring unit system, which is an internal, dimensionless scale calibrated to PDB statistics.

The CF scale is such that a typical good pose (RMSD ≤ 2 Å, correct binding mode on the Astex set) scores approximately −30 to −50 CF units for CF.com, while the worst clashing decoys can score +100 or more. The energy *range* that needs to be meaningfully Boltzmann-weighted is roughly 50 CF units. For a Boltzmann factor to distinguish scores 10 CF units apart, T_eff must be on the order of 10 / ln(10) ≈ 4 CF units. At T_eff = 0.596, exp(−10/0.596) ≈ 10^{−7}, which would make the distribution nearly a delta function — except that the GA population is not drawn from a flat prior; it is already concentrated near the low-CF basin by the selection pressure. At the values actually seen in the converged population (spread of ~5–15 CF units around the basin), T_eff = 0.596 produces a distribution with H ≈ 2–5 nats, which is empirically the range where the entropy bonus discriminates true minima from false ones.

The complementary calibration at T = 21 corresponds to the ISMB 2017 whiteboard, where the reference temperature was set to produce ΔG values in the −10 to +10 range for the ITC-187 calorimetry benchmark, matching the kcal/mol scale of experimental ΔG measurements. These are not the same calibration — one is a scoring temperature for internal ranking; the other is a reporting temperature for comparing to calorimetry.

---

## 5. Cluster Election and Spread Guard

After the GA converges and the CF population is scored thermodynamically, `cluster.cpp` groups poses by RMSD and elects a representative from each cluster. The default election criterion is the cluster head's CF (or ΔG_eff when `FLEXAIDDS_THERMO_SCORE=1`).

The **two-gate spread guard** (`FLEXAIDDS_CLUSTER_SPREAD_MAX`) optionally demotes the rank-0 cluster head when *all three* of the following conditions hold simultaneously:

1. **Isolation**: The rank-0 head is more than `cluster_spread_max` Å from its top-4 peers (i.e., it sits in a separate region of pose space, not the consensus cluster).
2. **Minority population**: The rank-0 cluster holds less than `cluster_pop_min_fraction` (default 0.35) of the merged population across restarts. A dominant true minimum would attract most of the population.
3. **No restart consensus**: Fewer than `cluster_consensus_k` (default 3) independent restarts converge within `cluster_consensus_tau` (default 2.0 Å) of the rank-0 head. If rank-0 is a real minimum, independent restarts should find it reproducibly.

The three conditions are conjunctive: all must be true for demotion to occur. This was learned from experience: an earlier single-gate version (isolation alone) demoted valid shallow-binding-mode clusters and cost 11 targets on the Astex benchmark (64/85, a 24-point regression). The three-gate version restores all those targets while still catching genuine false minima.

---

## 6. Reporting Reference

The key fields in the output `[THERMO]` block and their derivation:

| Field | Derivation | Notes |
|:------|:-----------|:------|
| `G_bind` | T_eff · ⟨CF⟩_raw − TdS_shannon + TdS_vib | Primary G score; includes vibrational term |
| `H_vct` | ⟨CF⟩ / n_heavy | Intensive; per-heavy-atom enthalpy proxy (ITC-comparable) |
| `H_vct_raw` | ⟨CF⟩ = Σ P_i · CF_i | Extensive; enters G_bind directly |
| `TdS_shannon` | T_eff · (−Σ P_i ln P_i) | Configurational entropy cost |
| `TdS_vib` | tencom_scale · (H_rep_bound − H_rep_ref) | Vibrational ΔS from tENCoM |
| `dG_eff` | ⟨CF⟩ − T_eff · H at T_eff | Pose-population free energy at T_eff=0.596 |
| `dG_eff_T21` | ⟨CF⟩ − 21 · H at T=21 | Reporting-only; ISMB 2017 calibration |
| `I_ES` | (ΔH + TΔS) / (ΔH − TΔS) at T=21 | Enthalpy-entropy index ∈ [−1, +1] |
| `binding_regime` | classify(ΔH, ΔS, T=21) | `no_binding`, `enthalpy_driven`, `both_favorable`, `entropy_driven`, `borderline` |
| `thermo_impossible` | ΔH > 0 ∧ ΔS < 0 | Gate verdict; forces dG_eff → +1000 when true |
