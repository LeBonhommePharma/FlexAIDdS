# Astex Non-Native (Cross-Docking) Benchmark

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
*J Chem Inf Model* 48(11):2214–2225. https://doi.org/10.1021/ci800224j

## Contents

| File | Description |
|:-----|:------------|
| `astex_non_native_set.csv` | Cross-docking pairs: target_pdb, ligand_pdb, ligand_id, threshold, target_name |
| `manifest.yaml` | Tier-1 bundle spec (3 pairs, baselines, entrypoint) |
| `run.sh` | Download + cross-dock + report script |
| `download.sh` | Fetch all required PDB files from RCSB |
| `environment.txt` | Hardware/software record for reproducibility |
| `expected/` | Reference outputs after first validated run |

## Dataset Structure

The CSV lists pairs as `(target_pdb, ligand_pdb, ligand_id)`:

- **target_pdb**: receptor structure used for docking (non-native conformation)
- **ligand_pdb**: where the ligand coordinates come from (native co-crystal)
- **ligand_id**: three-letter heteroatom code

The full set covers 65 protein families grouped by ligand scaffold. The CSV in
this bundle contains 73 representative pairs across 8 target families (CDK2,
ERK2, PPARγ, p38α, thymidine kinase, factor Xa, CDK2-staurosporine, CDK2-BRL,
CDK2-IDS).

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
