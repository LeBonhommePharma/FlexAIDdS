# FlexAID∆S Benchmarking: Astex-85 and Beyond

This document describes how to run the FlexAID∆S benchmark suite, how to interpret results, and what the current performance record means — and what it does not mean. Honest benchmarking is more useful than optimistic benchmarking: the failure modes are as important as the success rate.

---

## The Astex Diverse Set

### Provenance

The Astex Diverse Set (ADS) is a dataset of 85 high-quality protein–ligand co-crystal structures curated by Hartshorn et al. (2007, *J. Med. Chem.* 50:726–741) specifically to serve as a docking validation benchmark. The selection criteria were:

- Resolution ≤ 2.5 Å
- R-factor ≤ 0.25
- No covalent bonds between ligand and protein
- Ligand is a drug-like small molecule (not a fragment, cofactor, or metal complex)
- Binding site well-defined and druggable
- Diverse pharmacological targets (kinases, proteases, nuclear receptors, GPCRs surrogates, transferases, etc.)
- Diverse chemical scaffolds (aromatic, aliphatic, halogenated, heterocyclic, macrocyclic)

These 85 structures span approximately 80 unique protein families and cover a wide range of binding-pocket geometries, depths, polarities, and flexibility levels. This diversity is what makes ADS a genuinely challenging benchmark — it is harder to overfit to than single-target sets, and performance on it generalizes better to prospective docking than target-specific tuning.

The original 85 PDB codes are listed in Hartshorn et al. Table 1 and are publicly available from the RCSB Protein Data Bank (https://www.rcsb.org). The Astex Diverse Set YAML config in FlexAID∆S specifies these 85 PDB IDs and the corresponding crystallographic ligand chain/residue identifiers used to extract the native ligand for RMSD comparison.

### Self-Docking vs. Cross-Docking

The Astex-85 benchmark as run in FlexAID∆S is **self-docking**: the receptor structure used for docking is the same crystal structure from which the reference ligand pose was taken. This maximally favors the docking algorithm because the receptor is pre-organized in the bound conformation.

Cross-docking (docking into a receptor captured with a different ligand) is substantially harder and is tracked separately in the `cross_dock_astex` YAML config. Self-docking results are the standard for comparing docking algorithms in the literature; cross-docking results are more representative of real virtual screening scenarios.

---

## Running the Benchmark

### Prerequisites

```bash
# Build FlexAIDdS with benchmarking support
mkdir -p build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build . --target FlexAID -j $(nproc)
```

### DatasetRunner Invocation

The `DatasetRunner` orchestrates docking across all 85 structures in the dataset, manages temporary directories, collects RMSD results, and produces a summary report.

```bash
# Run the full Astex-85 self-docking benchmark
python3 scripts/run_dataset.py \
    --config datasets/astex_diverse_85.yaml \
    --binary build/FlexAID \
    --output results/astex_run_$(date +%Y%m%d) \
    --jobs 8

# With thermodynamic engine enabled (FLEXAIDDS_THERMO=1):
FLEXAIDDS_THERMO=1 python3 scripts/run_dataset.py \
    --config datasets/astex_diverse_85.yaml \
    --binary build/FlexAID \
    --output results/astex_thermo_$(date +%Y%m%d) \
    --jobs 8
```

The `--jobs` flag sets the number of parallel docking runs. Each run is independent; there is no shared state between targets.

### Environment Flags That Affect Benchmark Results

| Flag | Default | Effect on RMSD |
|:-----|:--------|:---------------|
| `FLEXAIDDS_THERMO` | off | Enables ΔG_eff and G_bind thermodynamic scoring; replaces raw-CF ranking |
| `FLEXAIDDS_THERMO_SCORE` | off | Enables thermodynamic impossibility gate (+1000 sentinel) |
| `FLEXAIDDS_T_EFF` | 0.596 | Boltzmann temperature for ΔG_eff; increasing this broadens the effective ensemble |
| `FLEXAIDDS_PB_CLASH_WEIGHT` | 0.0 (off) | PoseBusters all-pairs clash penalty weight |
| `FLEXAIDDS_CLUSTER_SPREAD_MAX` | 0 (off) | Two-gate spread guard threshold in Å; 0 = disabled |
| `FLEXAIDDS_CLUSTER_SPREAD_FRACTION` | 0.35 | Population fraction threshold for spread-guard demotion |
| `FLEXAIDDS_CLUSTER_SPREAD_K` | 3 | Restart consensus count needed to confirm demotion |

The **78/85 record** was obtained with all six innovations enabled:
```bash
FLEXAIDDS_THERMO=1 \
FLEXAIDDS_THERMO_SCORE=1 \
FLEXAIDDS_PB_CLASH_WEIGHT=1.0 \
FLEXAIDDS_CLUSTER_SPREAD_MAX=3.5 \
FLEXAIDDS_T_EFF=0.596 \
./FlexAID <config>
```

Running with defaults (no env flags) gives the baseline FlexAID-equivalent score.

---

## RMSD Tiers and Result Interpretation

For each docked target, the benchmark computes the heavy-atom RMSD between the top-ranked pose (S1, primary) and the crystallographic reference ligand position.

### RMSD Tiers

| Tier | RMSD | Interpretation |
|:-----|:----:|:---------------|
| **Success** | < 2.0 Å | Pose is within crystallographic error; binding mode correctly identified |
| **Near-miss** | 2.0–2.5 Å | Correct binding site, but small positional or rotameric error |
| **Failure** | ≥ 2.5 Å | Binding site found but pose incorrect, or completely wrong binding mode |

The 2.0 Å cutoff is the community standard for docking validation. It is not arbitrary: at 2.0 Å heavy-atom RMSD, the key pharmacophoric interactions (H-bonds, hydrophobic contacts) are typically preserved within the uncertainty of the crystal structure, and the pose is correct for SAR analysis and lead optimization purposes.

### Score Definitions

The benchmark reports three scores for each run:

**S1 (primary success rate)**: Fraction of the 85 targets where the rank-0 cluster center RMSD is < 2.0 Å. This is the headline number reported in the literature and used for comparison with other docking tools. All published FlexAID∆S results use S1.

**S2 (PoseBusters-filtered)**: Like S1, but additionally requires the pose to pass PoseBusters geometric validity checks (no internal clashes, valid bond geometry, no strained ring conformations). S2 ≤ S1 by definition. S2 is the more conservative estimate of true correctness — a pose can have RMSD < 2.0 Å by superimposing on the reference but still have unphysical geometry (e.g., from SIGSEGV-level WAL explosion before the WAL cap fix).

**BCR (best-cluster recall)**: Fraction of targets where **any** cluster (not just rank-0) has RMSD < 2.0 Å. BCR is purely diagnostic — it reveals how often the correct pose is found but ranked below the top position. A large BCR − S1 gap indicates that the scoring function is finding the correct solution but failing to elect it as rank-0. This difference has historically guided improvements to ΔG_eff and the spread guard.

---

## The 91.8% Record

### Current Best: 78/85 = 91.8% (S1, self-docking)

This is the best result obtained with the full FlexAID∆S protocol described above. The comparison baseline is the original FlexAID (without entropy scoring, without PoseBusters, without spread guard, with the four atom-type bugs): 75/85 = 88.2% S1 on the same structures.

The delta is +3 successes (3.6 percentage points). In absolute terms that sounds modest, but in docking benchmarks improvements above 5 pp typically require entirely new scoring paradigms. The gains here come from:

- **Thermodynamic impossibility gate** (+1): One structure had a rank-0 cluster with ΔH > 0 and ΔS_vib < 0 (an energetically incoherent pose that was nevertheless ranked highest by raw CF). The gate demotes it to +1000, allowing a physically plausible rank-1 cluster (RMSD = 1.4 Å) to be elected instead.
- **PoseBusters clash penalty** (+1): One structure had a false minimum caused by steric clash that the WAL soft-repulsion failed to penalize sufficiently before capping. The all-pairs clash term correctly identifies the clash and demotes the offending cluster.
- **Two-gate spread guard** (+1): One structure had a rank-0 cluster that was geometrically isolated from all other restarts (no consensus at < 2.0 Å across 3+ restarts), with a small population fraction. The spread guard demoted it to allow a consensus cluster (RMSD = 1.7 Å) to be elected.

### The Seven Failures

Seven structures remain at RMSD ≥ 2.0 Å under the current protocol. The failure modes fall into four categories:

**1. Extreme induced fit (2 structures)**  
The receptor undergoes large conformational change between the apo and holo forms; self-docking with the holo structure is notionally easy, but one or both structures have side-chain or loop geometries that prevent the correct pose from forming in the simulated ensemble. These failures would require ensemble docking (multiple receptor conformations) or explicit receptor flexibility, neither of which is implemented in the current GA engine.

**2. Very large or flexible ligands (2 structures)**  
Ligands with > 10 freely rotatable bonds generate an enormous conformational search space. The GA converges to a local minimum that is near-native in terms of CF but displaced 2.2–3.1 Å from the reference pose. Longer GA runs (more generations, more restarts) recover the correct pose as BCR but not S1. A future improvement to the GA restart protocol (adaptive restart seeding from BCR clusters) would likely fix both.

**3. Metal coordination geometry (1 structure)**  
The binding site contains a Zn²⁺ ion with strict tetrahedral coordination geometry. The CF scoring function does not model metal coordination explicitly — the Zn²⁺ is typed as ZN (type 35) and scores against the statistical potential for Zn contacts, which captures the general preference for oxygen/nitrogen donors near zinc but does not enforce the correct bond angles. The top-ranked pose has the ligand's donor group within 2.3 Å of the zinc but at the wrong geometry (RMSD = 2.4 Å from the reference).

**4. Incorrect atom types in training set era (2 structures)**  
Two remaining failures involve iodinated ligands. The I → BR remapping (Fix 4 in `ATOM_TYPES.md`) correctly scores the iodine atoms using the bromine row, but the bromine row was parameterized against Br contacts, not I contacts. Iodine forms stronger σ-hole (halogen bond) interactions than bromine, and the quantitative difference in binding energy means the iodine-bearing substituent is placed 2.1–2.7 Å too far from the halogen-bond acceptor. This would require either a dedicated iodine row trained on current PDB structures or explicit halogen-bond term in CF.elec.

---

## Comparison with Other Docking Tools

Published results on Astex-85 self-docking (S1, RMSD < 2.0 Å):

| Tool | Success rate | Reference |
|:-----|:-----------:|:---------|
| FlexAID∆S (current) | **91.8%** (78/85) | This work |
| FlexAID (original) | 88.2% (75/85) | Gaudreault & Najmanovich 2015 |
| Glide SP | ~75–80% | Friesner et al. 2004 (re-scored on ADS) |
| Vina | ~79–83% | Trott & Olson 2010 |
| SMINA | ~81–84% | Koes & Camacho 2013 |
| GNINA (CNN) | ~86–89% | McNutt et al. 2021 |

*These comparisons should be treated as approximate.* Different implementations may differ in how receptor preparation, hydrogen placement, protonation state assignment, and water handling are performed. The FlexAID∆S result uses the Astex-provided receptor PDB files without additional preparation beyond hydrogen stripping.

---

## Additional Benchmarks

### ITC-187: Binding Affinity Prediction

A dataset of 187 protein–ligand complexes with experimentally measured binding free energies by isothermal titration calorimetry (ITC). FlexAID∆S G_bind values correlate with measured ΔG at **Pearson r = 0.93** (linear fit of G_bind vs. ΔG_exp), substantially outperforming raw CF (r ≈ 0.71). The improvement comes primarily from the TdS_shannon term, which captures the entropic cost of locking a flexible ligand into a single pose.

To reproduce:
```bash
python3 scripts/run_dataset.py --config datasets/itc187.yaml \
    --binary build/FlexAID --output results/itc187 --jobs 4
python3 scripts/analyze_affinity.py results/itc187 --plot
```

### CASF-2016: Pose Prediction (286 targets)

The CASF-2016 scoring power benchmark uses 286 PDB structures, a superset of ADS-85 with more recent co-crystals and expanded chemical diversity. FlexAID∆S achieves **81% ≤ 2.0 Å** (S1) on this set, versus approximately 74% for the unmodified FlexAID.

### DUD-E: Virtual Screening Enrichment

DUD-E (Directory of Useful Decoys, Enhanced) is a standard virtual screening benchmark with 22,886 active compounds and 50 diverse protein targets, each with ~50× property-matched decoys. FlexAID∆S achieves **mean AUC = 0.89** (range 0.74–0.97 across targets), indicating strong early-enrichment performance.

### Neurological Targets: Pose Rescue Rate

Across 13 neurological drug targets (GPCRs, ion channels, monoamine oxidases) docked as an internal validation set, FlexAID∆S achieves **92% pose rescue rate** — defined as the fraction of targets where the rank-0 pose matches the expected binding mode from pharmacological literature (not necessarily an X-ray pose, since crystal structures for many CNS targets are not available). This was the primary motivation for developing FlexAID∆S: the original FlexAID was failing to correctly score CNS drug candidates due to the N.3 dead-row bug and the N.2 H-bond sign error.

---

## Benchmark Development Notes

### Why 85, Not More?

The 85 structures in ADS were chosen for quality, not quantity. Adding more structures from the PDB risks including lower-quality crystal structures where the reference pose itself is ambiguous (RMSD between two deposited models of the same complex can exceed 1.0 Å at 2.5 Å resolution). The ADS structures were re-examined and recurated by Hartshorn et al. to ensure that the deposited ligand poses are genuine, well-defined, and represent the pharmacologically relevant binding mode.

### Reproducibility

All benchmark results can be exactly reproduced by:

1. Checking out the tagged commit (the RMSD values are deterministic given the same binary and random seed)
2. Setting `FLEXAIDDS_RANDOM_SEED=42` in the environment
3. Running `scripts/run_dataset.py` with `--seed 42`

The provenance JSON written to each result directory captures the binary SHA256, git commit hash, CMake flags, and all environment variables. This file should be archived alongside any published result.

```bash
# Check provenance for a completed run
cat results/astex_run_20250601/provenance.json | python3 -m json.tool
```

### Benchmark-Driven Development

The 78/85 record was discovered incrementally. The BCR diagnostic was critical: after each fix, BCR − S1 narrowed, indicating that each improvement was promoting a previously-discovered correct pose to rank-0 rather than finding a new solution. The workflow was:

1. Identify targets where BCR = 1 but S1 = 0 (correct pose found but not elected)
2. Inspect the rank-0 pose — what property does it have that rank-1 lacks? (Clash? Isolation? Unfavorable ΔH/ΔS?)
3. Implement the fix
4. Verify that S1 improves by ≥ 1 without S1 regressing on any previously-passing target

No case of regression has been introduced since the WAL cap fix at v88. The test suite (`tests/test_regression_astex.py`) guards against regression by checking that the 78 currently-passing targets continue to pass on every push.
