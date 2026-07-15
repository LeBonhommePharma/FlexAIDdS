# Astex Diverse Set Benchmark

## Overview

The Astex Diverse Set is the gold-standard native-pose docking benchmark: 85
high-quality protein-ligand co-crystal structures selected for diversity across
drug-target families and medicinal-chemistry relevance. Each complex provides a
receptor structure and a co-crystallized ligand pose that serves as the
experimental ground truth.

**Success criterion:** RMSD < 2.0 Å vs the co-crystal pose.

## Citation

Hartshorn MJ, Verdonk ML, Chessari G, et al. (2007)
"Diverse, high-quality test set for the validation of protein-ligand docking performance."
*J Med Chem* 50(4):726–741. https://doi.org/10.1021/jm061277y

Original set: Nissink JWM, Murray C, Hartshorn M, et al. (2002)
*Proteins* 49(4):457–471. https://doi.org/10.1002/prot.10232

## Canonical data path

**Structures live in `astex_diverse/<PDB>/` (this directory’s nested tree).**
That path is the repository source of truth for the 85-complex set.

- Full map of duplicates / deprecations: [`../datasets/CANONICAL.md`](../datasets/CANONICAL.md)
- Checksums: [`../datasets/astex_diverse_sha256.csv`](../datasets/astex_diverse_sha256.csv)
- Nested `data/astex_diverse/` is **deprecated** (do not use for new work)

## Contents

| File | Description |
|:-----|:------------|
| `astex_diverse/` | **Canonical** per-PDB structures (apo, ligand, binding site, deposit) |
| `astex_diverse_set.csv` | 85 complexes: PDB ID, ligand ID, resolution, RMSD threshold, citation |
| `manifest.yaml` | Tier-1 bundle spec (5 targets, baselines, entrypoint) |
| `run.sh` | Download + dock + report script |
| `download.sh` | Fetch PDB files from RCSB; extract apo receptors + ligands via `benchmark_datasets` |
| `environment.txt` | Hardware/software record for reproducibility |
| `expected/` | Reference outputs after first validated run |
| `data/` | **Deprecated** nested copy + tier-1 loose PDBs — see `data/DEPRECATION.md` |
| `structures/` | Symlink alias layer only (not a second store) |

## Dataset

82 unique PDB entries (some ligand IDs repeat across structurally related complexes).
Major target families:

| Ligand | Protein family | # complexes |
|:-------|:---------------|:-----------:|
| MK1 | CDK2 | 14 |
| PLN | Thymidine kinase | 17 |
| 2AN | p38α MAP kinase | 8 |
| STU | CDK2 / thrombin (staurosporine) | 7 |
| BRL | CDK2 | 6 |
| GW5 | PPARγ | 6 |
| IDS | CDK2 | 7 |
| Others | Various | 20 |

## Published Reference Success Rates (top-1, RMSD < 2 Å)

| Method | Success rate |
|:-------|:-----------:|
| AutoDock Vina | 53–60% |
| Glide SP | 62–68% |
| Glide XP | 68–74% |
| GOLD ChemScore | 62–70% |
| FlexAID (no entropy) | ~55–62% |
| **FlexAIDdS (target)** | **~70%** |

## How to Run

### Tier-1 (5 targets, ~10 min)

```bash
# 1. Download structures and extract ligands
bash benchmarks/astex_diverse/download.sh

# 2. Run docking and validate baselines
bash benchmarks/astex_diverse/run.sh

# 3. Check results
cat benchmarks/astex_diverse/results/report.md
```

### Full benchmark (85 complexes, tier-2)

```bash
# Using the C++ benchmark_datasets runner
./build/benchmark_datasets --benchmark astex --output results/astex/ --threads 8

# Or using the Python DatasetRunner
python -m flexaidds.dataset_runner --dataset astex_diverse --tier 2
```

### CMake convenience targets

```bash
cmake --build build --target flexaid_bench_astex_diverse
```

## Expected Runtime

- Tier-1 (5 targets): ~10–30 min depending on hardware
- Full (85 complexes): ~4–12 hours (serial); ~1–2 hours with 8 threads
- GPU (Metal/CUDA): ~30–60 min for the full set

## Notes

- Structures are fetched directly from RCSB; no registration required.
- Ligand SDF extraction is automated via `benchmark_datasets --prepare-only`.
- The `astex_diverse_set.csv` lists all 85 entries with resolution and ligand ID.
- Entropy rescue is the fraction of complexes where ΔS correction moves a
  sub-2Å pose from rank > 1 to rank 1 (FlexAIDdS-specific metric).
