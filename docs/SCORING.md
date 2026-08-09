# Scoring in FlexAID∆S: CF proxies and ensemble diagnostics

This document describes the cost function, NRGDock matrix, score-space ensemble
diagnostics, and their current integration boundaries. Names such as `dG_eff`
and `G_bind` are legacy wire fields, not evidence of physical free energy.

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

## 3. Score-space ensemble transform (legacy name: `dG_eff`)

The GA produces optimizer-selected CF records, not an enumerated equilibrium
ensemble. FlexAIDdS applies the following soft-min transform to some retained
records:

```
p_i       = exp(-CF_i / T_eff) / sum_j exp(-CF_j / T_eff)
H_sample  = -sum_i p_i ln(p_i)
G_tilde   = <CF>_p - T_eff H_sample
```

`T_eff` is expressed in the internal CF scale. `G_tilde` is useful as a
ranking/diagnostic proxy, but it is **not** a Helmholtz free energy or binding
free energy: CF is uncalibrated, the GA population has no canonical measure,
and exact K-fold duplication of the deposited records changes `G_tilde` by
`-T_eff ln(K)`. Sampling budget, deduplication, and cluster cardinality can
therefore change the value without changing the represented molecular states.

The historical field name `dG_eff` remains in some output for compatibility.
New consumers must inspect schema-v2 scientific provenance and report this
quantity as a CF/optimizer-sample proxy unless a separately validated energy
calibration and ensemble measure are supplied.

### 3.1 Recorded effective-temperature conventions

The code contains score-space constants such as `T_eff = 0.596` and the legacy
`T = 21` reporting convention. These are CF-scale parameters, not temperatures
in kelvin and not a conversion from CF to kcal/mol. Their empirical usefulness
must be established by a versioned benchmark receipt; this document does not
promote them to calorimetric calibration.

### 3.2 Shannon/tENCoM composition status

The legacy `G_bind` and `ShannonThermoStack` output mixes a CF proxy, Shannon
entropy of optimizer records, and a model-scale tENCoM diagnostic. The current
implementation does not provide a matched `Q_RL/(Q_R Q_L)` association cycle,
standard-state term, or validated unit conversion. It must therefore be
reported as `proxy_only`, never as physical `Delta G`, affinity, `Kd`, or `Ki`.

The legacy stack also receives a configurational free-energy value that already
contains `-T S_config` and subtracts a Shannon term again. Until the composition
path is repaired and independently tested, its numeric `deltaG` field is not an
additive thermodynamic ledger.

### 3.3 Impossibility predicate status

The predicate `dH > 0 && dS < 0` is mathematically valid only when both inputs
are commensurate state-function differences for the same physical cycle. The
current production inputs do not satisfy that contract. Moreover, the computed
sentinel is printed in the GA diagnostics but is not consumed by the later
exact-CF rescore, sort, or clustering path. `FLEXAIDDS_THERMO_SCORE=1` therefore
does **not** currently enforce a physics filter on elected poses.

---

## 4. Effective-temperature provenance

`T_eff = 0.596` and `T = 21` are historical CF-scale conventions. A numerical
constant that makes a score distribution convenient does not establish an
energy calibration. Any future statement that either convention predicts
calorimetric values must identify the calibration dataset, exact protocol and
code SHA, fixed denominator, held-out validation, uncertainty, and artifact
receipt. In the absence of that bundle, both values remain score parameters.

---

## 5. Cluster Election and Spread Guard

After the GA converges, the retained chromosomes are exactly rescored by CF,
sorted, and passed to clustering. Although thermodynamic-looking diagnostics
may be computed during the GA, the current top-level rescore boundary does not
promote their sentinel or `dG_eff` field into final clustering/election.

The **two-gate spread guard** (`FLEXAIDDS_CLUSTER_SPREAD_MAX`) optionally demotes the rank-0 cluster head when *all three* of the following conditions hold simultaneously:

1. **Isolation**: The rank-0 head is more than `cluster_spread_max` Å from its top-4 peers (i.e., it sits in a separate region of pose space, not the consensus cluster).
2. **Minority population**: The rank-0 cluster holds less than `cluster_pop_min_fraction` (default 0.35) of the merged population across restarts. A dominant true minimum would attract most of the population.
3. **No restart consensus**: Fewer than `cluster_consensus_k` (default 3) independent restarts converge within `cluster_consensus_tau` (default 2.0 Å) of the rank-0 head. If rank-0 is a real minimum, independent restarts should find it reproducibly.

The three conditions are conjunctive. Any accuracy statement about this guard
belongs in a fixed-denominator benchmark receipt, not in the scoring contract.

---

## 6. Reporting Reference

Legacy `[THERMO]` fields are retained for compatibility. Their current domain
and safe interpretation are:

| Field | Derivation | Notes |
|:------|:-----------|:------|
| `G_bind` | mixed legacy composition | `proxy_only`; not binding free energy |
| `H_vct` | ⟨CF⟩ / n_heavy | per-heavy-atom CF diagnostic; not calorimetric enthalpy |
| `H_vct_raw` | ⟨CF⟩ = Σ P_i · CF_i | CF-domain weighted mean |
| `TdS_shannon` | T_eff · (−Σ P_i ln P_i) | optimizer-sample Shannon diagnostic |
| `TdS_vib` | scaled difference of model-spectrum Shannon diagnostics | model scale; unmatched reference cycle |
| `dG_eff` | ⟨CF⟩ − T_eff · H | CF-domain soft-min proxy |
| `dG_eff_T21` | ⟨CF⟩ − 21 · H | legacy CF-domain reporting proxy |
| `I_ES` | ratio of proxy terms | diagnostic only; forbidden for ranking or physical claims |
| `binding_regime` | classification of proxy terms | label only; not a binding-state verdict |
| `thermo_impossible` | predicate over non-commensurate proxy terms | diagnostic only; not wired to final election |
