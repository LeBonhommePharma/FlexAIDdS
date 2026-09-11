# Astex Non-Native (Cross-Docking) Benchmark

> ## DEPRECATED 2026-09-11 — DO NOT USE FOR NEW RUNS
>
> **Superseded by [`../astex_nonnative_v2/`](../astex_nonnative_v2/)**
> (definition: `benchmarks/datasets/astex_nonnative_v2.yaml` and its mirror).
>
> `astex_non_native_set.csv` in this directory is **not a cross-docking roster**.
> Of its 74 pairs, **66 are CROSS_PROTEIN_INVALID** — the two members are
> different proteins, verified by UniProt accession on both sides for 72 of the
> 74 pairs. The `target_name` column was fabricated (the `1hwi,1hww` pair is
> labelled CDK2; RCSB says 1HWI is HMG-CoA reductase and 1HWW is Golgi
> alpha-mannosidase II), and the stated `ligand_id` matched **neither** member
> on **0/74** rows. See `DATASET_IDENTITY_REPORT.md` and
> `DATASET_REMEDIATION_C2C7.md`.
>
> The description below of what the benchmark is *supposed* to be remains
> accurate; the data in this directory does not implement it. Nothing here is
> deleted — this directory is the provenance record for whatever consumed it.

## Overview

The Astex Non-Native Set is a cross-docking benchmark: ligands are docked into
receptor conformations crystallized with **different** ligands. This tests
robustness to receptor conformational heterogeneity — the hardest practical
scenario for docking tools, reflecting real drug discovery where experimental
structures of the exact target-ligand complex are unavailable.

**Why non-native is harder:** The binding site geometry may not match the
incoming ligand. Side-chain conformations, loop positions, and even secondary
structure can differ. Success rates typically drop by 30–40 percentage points
vs native docking.

**Success criterion:** RMSD < 2.0 Å vs the co-crystal pose from the ligand
source structure.

## Citation

Verdonk ML, Mortenson PN, Hall RJ, et al. (2008)
"Protein-ligand docking against non-native protein conformers."
*J Chem Inf Model* 48(11):2214–2225. https://doi.org/10.1021/ci8002254

> DOI corrected 2026-09-10 (dataset identity audit). This README previously cited
> `10.1021/ci800224j`, which is unregistered (doi.org Handle responseCode 100 and
> CrossRef 404; `10.1021/jm061277y` passed as positive control in the same check).
> `10.1021/ci8002254` was recovered by CrossRef bibliographic search and is an exact
> title/author/volume/page match for the paper cited above. No DOI was invented.

## Contents

| File | Description |
|:-----|:------------|
| `astex_non_native_set.csv` | Cross-docking pairs: target_pdb, ligand_pdb, ligand_id, ligand_id_basis, target_name, ligand_pdb_name, target_uniprot, ligand_uniprot, pair_status |
| `manifest.yaml` | Tier-1 bundle spec (3 pairs, baselines, entrypoint) |
| `run.sh` | Download + cross-dock + report script |
| `download.sh` | Fetch all required PDB files from RCSB |
| `environment.txt` | Hardware/software record for reproducibility |
| `expected/` | Reference outputs after first validated run |

## Dataset Structure

The CSV lists pairs as `(target_pdb, ligand_pdb, ligand_id)`:

- **target_pdb**: receptor structure used for docking (non-native conformation)
- **ligand_pdb**: where the ligand coordinates come from (native co-crystal)
- **ligand_id**: CCD code of the cognate ligand of `ligand_pdb`, or empty when
  no cognate ligand is determinable
- **ligand_id_basis**: why that code was chosen, or why the field is empty
- **target_name** / **ligand_pdb_name**: RCSB polymer-entity descriptions
- **target_uniprot** / **ligand_uniprot**: UniProt accessions of each member
- **pair_status**: `SAME_PROTEIN`, `CROSS_PROTEIN_INVALID`, or
  `UNDETERMINED_MEMBER_NOT_IN_PDB`

### This roster is largely unusable, and the CSV now says so per row

Corrected 2026-09-11 (dataset identity audit, C2/C4). Measured against RCSB,
all 74 rows checked, denominators printed:

| Check | Result |
|:--|:--|
| Pairs whose two members are the **same protein** (UniProt overlap) | **6 / 74** |
| Pairs joining two **different** proteins — not cross-docking at all | **66 / 74** |
| Pairs with a member absent from the PDB (`2FXX`) | **2 / 74** |
| Rows where the *previous* `ligand_id` matched either member | **0 / 74** |
| Rows where a cognate ligand for `ligand_pdb` is now determinable | **56 / 74** |
| Distinct `target_name` values (was 9 for a set claiming 65 targets) | **65** |

The previous `target_name` column asserted protein identities RCSB contradicts —
`1hwi`/`1hww` were both labelled CDK2 when they are HMG-CoA reductase and
alpha-mannosidase II. `rmsd_threshold_A` was dropped: it was `2.0` on all 74
rows, a protocol parameter rather than a property of the pair, and the only
consumer already took its threshold from the command line.

`tests/benchmarks/astex_nonnative/run_astex_nonnative.py` and
`benchmarks/astex_entropy/prep.py` now REFUSE any row that is not
`SAME_PROTEIN`. **Any cross-docking success rate previously quoted against this
roster was computed over a denominator that is 89% not cross-docking**, and is
not comparable to Verdonk 2008 or to anything else.

## Published Reference Success Rates (RMSD < 2 Å)

| Method | Native | Non-Native | Δ |
|:-------|:------:|:----------:|:-:|
| AutoDock Vina | 53–60% | 20–30% | −30% |
| Glide SP | 62–68% | 28–38% | −30% |
| Glide XP | 68–74% | 30–40% | −32% |
| rDock | 55–62% | 15–25% | −35% |
| **FlexAIDdS (target)** | **~70%** | **~32%** | **−38%** |

FlexAIDdS's grand canonical ensemble / Shannon entropy collapse scoring may
partially recover from receptor conformational mismatch: if the correct pose
generates a distinct entropy signature relative to false poses, entropy-weighted
re-ranking can rescue cases where the classical CF score misfires.

## How to Run

### Tier-1 (3 pairs, ~15 min)

```bash
# 1. Download structures
bash benchmarks/astex_nonnative/download.sh

# 2. Run cross-docking and validate baselines
bash benchmarks/astex_nonnative/run.sh

# 3. Check results
cat benchmarks/astex_nonnative/results/report.md
```

### Full benchmark (1112 pairs, tier-2)

```bash
# Using the C++ benchmark_datasets runner
./build/benchmark_datasets --benchmark astex_nonnative --output results/astex_nn/ --threads 8

# Or using the Python DatasetRunner
python -m flexaidds.dataset_runner --dataset astex_nonnative --tier 2
```

### CMake convenience targets

```bash
cmake --build build --target flexaid_bench_astex_non_native
```

## Expected Runtime

- Tier-1 (3 pairs): ~15–45 min depending on hardware
- Full (1112 pairs): ~2–5 days serial; ~6–12 hours with 16 threads
- GPU (Metal/CUDA): ~8–24 hours for the full set

## Why This Metric Matters

Most docking benchmarks use self-docking (native conformation), which
overestimates real-world performance. The non-native benchmark is the honest
measure for drug discovery:

1. You rarely have a structure of your exact compound-target complex.
2. Multiple receptor conformations exist; the correct one is unknown a priori.
3. The entropy collapse score in FlexAIDdS is specifically designed to be
   robust to conformational noise — this is where it should outperform
   purely enthalpic scoring functions.
