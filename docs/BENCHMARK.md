# FlexAID∆S Benchmarking: Astex-85 and Beyond

This document describes how to run the FlexAID∆S benchmark suite and how to interpret
results. Honest benchmarking is more useful than optimistic benchmarking: the failure
modes are as important as the success rate.

**This repository currently publishes no Astex-85 docking-power rate.** See the
withdrawal notice below and `README.md`.

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
cmake --build . --target FlexAID benchmark_datasets -j $(nproc)
```

`scripts/run_dataset.py` and `scripts/analyze_affinity.py` **do not exist** in this tree.
Do not copy commands that name them.

### DatasetRunner invocation

The C++ `benchmark_datasets` binary (CMake target) orchestrates docking across the
dataset, manages temporary directories, collects RMSD results, and writes a summary.
The skill wrapper adds manifest capture and packaging. The one-command reviewer path
is `scripts/reproduce_astex85.sh` (default **blind**; see `REPRODUCIBILITY.md`).

```bash
# Canonical campaign binary (METHODOLOGY.md §0 / §3)
./build/benchmark_datasets \
    --benchmark "crossdock_json:benchmarks/datasets/benchmark_astex_native_85.json" \
    --output results/astex_run \
    --threads 4 \
    --omp-threads 2

# Skill wrapper (manifest + optional --package)
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
    --dataset astex_diverse --tier 1 --dry-run

# Blind reviewer reproduction (SEED_ELITISM=0; not an oracle ceiling)
bash scripts/reproduce_astex85.sh
```

Dataset YAML: `benchmarks/datasets/astex_diverse.yaml` (`docking_mode: self_docking`).
Each target run is independent; there is no shared state between targets.

Follow `METHODOLOGY.md` §0 (matrix pin, seed, budget, in-place RMSD) and §3
(autonomous / blind; no seed elitism). A rate is publishable only with a provenance
receipt. `FLEXAIDDS_NATIVE_SEED_FRAC` is a **dead knob** on today's DatasetRunner path
(the runner always emits `seed_fraction: 0.0`). The live oracle lever is
`FLEXAIDDS_SEED_ELITISM=1` (injects `_INI.pdb`). Do not export either on a docking-power run.

### Environment flags that affect ranking (method documentation)

These flags are listed so a receipt can record what was on. They are **not** a published
recipe for a success rate.

| Flag | Default | Effect |
|:-----|:--------|:-------|
| `FLEXAIDDS_THERMO` | off | Enables thermodynamic scoring overlay |
| `FLEXAIDDS_THERMO_SCORE` | off | Enables thermodynamic impossibility gate (+1000 sentinel) |
| `FLEXAIDDS_T_EFF` | 0.596 | Boltzmann temperature for ΔG_eff |
| `FLEXAIDDS_PB_CLASH_WEIGHT` | 0.0 (off) | PoseBusters all-pairs clash penalty weight |
| `FLEXAIDDS_CLUSTER_SPREAD_MAX` | 0 (off) | Two-gate spread guard threshold in Å; 0 = disabled |
| `FLEXAIDDS_CLUSTER_SPREAD_FRACTION` | 0.35 | Population fraction threshold for spread-guard demotion |
| `FLEXAIDDS_CLUSTER_SPREAD_K` | 3 | Restart consensus count needed to confirm demotion |
| `FLEXAIDDS_SEED_ELITISM` | off | If `1`, injects `_INI.pdb` — **oracle ceiling, not docking power** |

The withdrawn 78/85 figure was obtained with several of these innovations enabled
(retained here as method documentation only — not a claim; see the withdrawal notice):

```bash
FLEXAIDDS_THERMO=1 \
FLEXAIDDS_THERMO_SCORE=1 \
FLEXAIDDS_PB_CLASH_WEIGHT=1.0 \
FLEXAIDDS_CLUSTER_SPREAD_MAX=3.5 \
FLEXAIDDS_T_EFF=0.596 \
./FlexAID <config>
```

Running with defaults (no env flags) is the baseline FlexAID-equivalent path.

---

## RMSD tiers and result interpretation

For each docked target, the benchmark computes the heavy-atom RMSD between the top-ranked pose (S1, primary) and the crystallographic reference ligand position. Success is defined in `METHODOLOGY.md` §0 / §3: rank-0 in-place RMSD **`<= 2.0 Å`**.

### RMSD tiers

| Tier | RMSD | Interpretation |
|:-----|:----:|:---------------|
| **Success** | **<= 2.0 Å** | Pose is within the community docking-power cutoff; binding mode correctly identified |
| **Near-miss** | > 2.0–2.5 Å | Correct binding site, but small positional or rotameric error |
| **Failure** | ≥ 2.5 Å | Binding site found but pose incorrect, or completely wrong binding mode |

The 2.0 Å cutoff is the community standard for docking validation. It is not arbitrary: at 2.0 Å heavy-atom RMSD, the key pharmacophoric interactions (H-bonds, hydrophobic contacts) are typically preserved within the uncertainty of the crystal structure, and the pose is correct for SAR analysis and lead optimization purposes.

### Score definitions

The benchmark reports three scores for each run. **This repository currently publishes none of them as a FlexAID∆S rate.**

**S1 (primary success rate)**: Fraction of the 85-target denominator where the rank-0 cluster center RMSD is **<= 2.0 Å**. This is the headline metric used in the literature for comparison with other docking tools.

**S2 (PoseBusters-filtered)**: Like S1, but additionally requires the pose to pass PoseBusters geometric validity checks (no internal clashes, valid bond geometry, no strained ring conformations). S2 ≤ S1 by definition. S2 is the more conservative estimate of true correctness — a pose can have RMSD <= 2.0 Å by superimposing on the reference but still have unphysical geometry.

**BCR (best-cluster recall)**: Fraction of targets where **any** cluster (not just rank-0) has RMSD <= 2.0 Å. BCR is purely diagnostic — it reveals how often the correct pose is found but ranked below the top position. A large BCR − S1 gap indicates that the scoring function is finding the correct solution but failing to elect it as rank-0. That gap is expected under RMSD-blind CF ranking; it is not a licence to elect by RMSD.

---

## Astex-85 success rate — WITHDRAWN pending verification

> **No Astex-85 success rate is published by this repository at present.**
>
> The figures previously stated here (78/85 = 91.8%) are **withdrawn**. They could not be
> reproduced from a receipted, blind, unseeded run on the current engine, and the
> reproduction path they referenced (`scripts/run_dataset.py`, `scripts/analyze_affinity.py`)
> does not exist in the tree. This matches the status already stated in `README.md`:
> *unverified / pending receipt*.
>
> A rate will be republished only when it is produced by a run that satisfies
> `METHODOLOGY.md` §0 — blind, `native_pose_seeded=0`, no seed elitism, a fixed 85-target
> denominator, and a provenance receipt pinning engine hash, matrix hash, and input hashes.
>
> Do not cite 91.8%, 94.1%, or 88.2% as current FlexAID(∆S) docking power. The 88.2%
> (75/85) figure previously listed as FlexAID 2015 S1 is a **literature misquote**;
> Gaudreault & Najmanovich 2015 JCIM **Table 2** is **top-1 45.2% / top-10 66.7%**
> (`METHODOLOGY.md` §3). 94.1% (80/85) was an oracle ceiling (`REPRODUCIBILITY.md`).
>
> The methodology, diagnostics and failure-mode *categories* below remain valid as *method*
> documentation. Any number appearing in them is illustrative of the analysis, not a claim.

---

## Comparison with other docking tools

Published literature on Astex-85 self-docking. FlexAID∆S has **no verified rate** in this
repository. FlexAID 2015 numbers are JCIM Table 2 (do not swap top-1 / top-10 labels).

| Tool | Success rate | Reference |
|:-----|:-----------:|:---------|
| FlexAID∆S (current) | *withdrawn — no verified rate published* | `README.md`; this page |
| FlexAID (2015) | **top-1 45.2% / top-10 66.7%** | Gaudreault & Najmanovich 2015 JCIM Table 2; `METHODOLOGY.md` §3 |
| Glide SP | ~75–80% | Friesner et al. 2004 (re-scored on ADS) |
| Vina | ~79–83% | Trott & Olson 2010 |
| SMINA | ~81–84% | Koes & Camacho 2013 |
| GNINA (CNN) | ~86–89% | McNutt et al. 2021 |

The figure 88.2% (75/85) previously listed in this table as FlexAID 2015 S1 is **withdrawn
as a literature misquote** — not JCIM Table 2.

*Literature comparisons for Glide / Vina / SMINA / GNINA should be treated as approximate.*
Different implementations may differ in receptor preparation, hydrogen placement, protonation,
and water handling.

---

## Additional benchmarks — unverified / no receipt

> **Unverified / no receipt.** The Pearson **r = 0.93** (ITC-187), CASF-2016 **81%**,
> DUD-E **AUC 0.89**, and CNS **92%** figures previously stated here are **withdrawn**.
> No provenance receipt for those numbers exists in this repository. Do not cite them
> as FlexAID∆S performance. Dataset YAML files may exist; that is not a result.

### ITC-187: binding affinity prediction (intended campaign)

A dataset of 187 protein–ligand complexes with experimentally measured binding free energies
by isothermal titration calorimetry. YAML: `benchmarks/datasets/itc187.yaml`.
`scripts/analyze_affinity.py` does not exist. When this campaign is run, use `benchmark_datasets` or
`.grok/skills/flexaidds/scripts/dataset_runner.py` and archive a `METHODOLOGY.md` §0 receipt
before quoting any correlation.

### CASF-2016, DUD-E, neurological targets

These remain *intended* evaluation surfaces (pose prediction, virtual-screening enrichment,
internal CNS pose-rescue). Previously quoted percentages (CASF 81%, DUD-E AUC 0.89, CNS 92%)
are withdrawn until a receipted run exists.

---

## Benchmark development notes

### Why 85, not more?

The 85 structures in ADS were chosen for quality, not quantity. Adding more structures from the PDB risks including lower-quality crystal structures where the reference pose itself is ambiguous (RMSD between two deposited models of the same complex can exceed 1.0 Å at 2.5 Å resolution). The ADS structures were re-examined and recurated by Hartshorn et al. to ensure that the deposited ligand poses are genuine, well-defined, and represent the pharmacologically relevant binding mode.

### Reproducibility

A docking-power run is reproducible only when it is **blind** and receipted:

1. Follow `METHODOLOGY.md` §0 and §3 (autonomous, `native_pose_seeded=0`, no seed elitism).
2. Pin engine hash, energy-matrix hash (`MC_st0r5.2_6.dat`), and input hashes in provenance JSON.
3. Run `benchmark_datasets`, `.grok/skills/flexaidds/scripts/dataset_runner.py`, or
   `scripts/reproduce_astex85.sh` (default blind). Do not invoke `scripts/run_dataset.py`
   (it does not exist).

The provenance JSON written to each result directory captures the binary SHA256, git commit hash, CMake flags, and all environment variables. This file should be archived alongside any published result.

```bash
# Check provenance for a completed run
python3 -m json.tool results/astex_run/provenance.json
```

### Benchmark-driven development

The withdrawn 78/85 figure was arrived at incrementally; that history is **not** a current
pass list. BCR − S1 remains a valid *diagnostic* (correct pose found but not elected). Do
not elect the near-native cluster by RMSD in order to close that gap — that converts S1
into an oracle. There is no in-tree `tests/test_regression_astex.py` guarding
"currently-passing" targets.
