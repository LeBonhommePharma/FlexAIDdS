# FlexAIDdS Benchmark Suite

Automated distributed benchmarking system for FlexAIDdS across all major
molecular docking datasets.  Evaluates entropy-aware scoring against gold
standards for scoring power, docking power, virtual screening enrichment,
and the primary FlexAIDdS discriminator: **ΔS rescue rate** (Shannon Energy
Collapse metric).

---

## Directory Structure

```
benchmarks/
├── DatasetRunner.py      # Orchestrator: discovery, distribution, aggregation
├── metrics.py            # Pure-computation metric functions
├── ligand_prep.py        # SMILES → Mol2 → .inp preparation pipeline
├── run.py                # CLI entry point
├── __init__.py           # Package init
├── datasets/             # YAML dataset configurations + canonical docs
│   ├── CANONICAL.md      # Which Astex tree is authoritative (+ deprecations)
│   ├── astex_diverse.yaml # Astex Diverse Set (FlexAID 2015 validation #1)
│   ├── astex_diverse_sha256.csv  # Key-file checksums for canonical tree
│   ├── astex_diverse_manifest.json
│   ├── astex_nonnative.yaml # Astex Non-Native Set (FlexAID 2015 validation #2)
│   ├── hap2.yaml         # HAP2 holo/apo benchmark (FlexAID 2015 validation #3)
│   ├── casf2016.yaml     # CASF-2016 (285 targets, 4-power benchmark)
│   ├── itc187.yaml       # ITC-187 calorimetry gold standard
│   ├── dude37.yaml       # DUD-E 37-target cross-docking
│   ├── lsd_docking.yaml  # lsd.docking.org 11-target prospective
│   ├── psychopharm23.yaml # NRGlab 23 CNS targets
│   ├── naloxone_ss.yaml  # Naloxone / Serotonin Syndrome focus + entropy collapse
│   ├── muv.yaml          # MUV 17-target unbiased VS
│   └── erds_specificity.yaml  # ERDS Z-score specificity
├── astex_diverse/        # Astex harness + CANONICAL structures
│   └── astex_diverse/    # <PDB>/ apo, ligand, binding_site (source of truth)
└── smoke/                # Lightweight smoke validation (CI sanity)
    ├── run.sh
    ├── manifest.yaml
    └── README.md
```

**Astex structures:** use `benchmarks/astex_diverse/astex_diverse/` only.
Nested `benchmarks/astex_diverse/data/astex_diverse/` is deprecated
(`benchmarks/datasets/CANONICAL.md`). Large `astex_repro/` campaign outputs
are local/gitignored and are not the dataset.

---

## Quick Start

### Install dependencies

```bash
pip install numpy scipy pyyaml
# Optional for ligand preparation:
conda install -c conda-forge rdkit
# Optional for distributed execution:
pip install mpi4py
```

### Run a tier-1 (fast) benchmark

```bash
# 5 CASF targets — ~2 minutes
python -m benchmarks.run --dataset casf2016 --tier 1

# Specific metric only
python -m benchmarks.run --dataset itc187 --tier 1 --metric entropy_rescue_rate

# Dry run (no docking binary required — synthetic scores)
python -m benchmarks.run --all --tier 1 --dry-run
```

### Run a full tier-2 benchmark

```bash
# All datasets
python -m benchmarks.run --all --tier 2 --bootstrap

# Specific datasets
python -m benchmarks.run --dataset casf2016 --dataset psychopharm23 --tier 2
```

### Distributed MPI run

```bash
mpirun -n 8 python -m benchmarks.run --all --tier 2 --distributed --nodes 8
```

### Combined distributed test + benchmark pass

```bash
MPI_PROCS=4 CTEST_JOBS=8 ./scripts/run_distributed_validation.sh
```

### Use as a library

```python
from benchmarks.DatasetRunner import DatasetRunner

runner = DatasetRunner(
    datasets_dir="benchmarks/datasets",
    results_dir="results/my_run",
    binary="/path/to/FlexAID",
    bootstrap_ci=True,
)

# Run one dataset
result = runner.run_single("itc187", tier=2)
print(f"Pearson r = {result.metrics['scoring_power_pearson_r']:.3f}")
print(f"Rescue rate = {result.metrics['entropy_rescue_rate']:.3f}")

# Run all and save report
report = runner.run_all(tier=2)
report.save("results/my_run/final_report")
```

---

## Benchmark Tiers

| Tier | Trigger | Targets | Timeout | Purpose |
|------|---------|---------|---------|---------|
| **1** | PR to `master` | 5 CASF targets | 10 min | Sanity gate — catches obvious regressions |
| **2** | Merge to `master` + weekly cron | All datasets | 6 h | Full statistical validation |

Tier-1 runs gate merges: a regression (measured metric < baseline × 0.95)
marks the PR check as failed.

---

## Datasets

### FlexAID 2015 Validation Sequence

These are the first-class docking benchmarks for FlexAID lineage comparisons,
in the order used by Gaudreault & Najmanovich (2015):

1. **Astex Diverse Set** (`astex_diverse`) — 85 native holo complexes for self-docking.
2. **Astex Non-Native Set** (`astex_nonnative`) — 65 targets / 1112 structures for cross-docking into non-native receptor conformations.
3. **HAP2** (`hap2`) — FlexAID-native holo/apo protein-pair benchmark from the 2015 JCIM paper.

Run all discovered datasets with `--all`; the packaged dataset runner uses
`benchmark_order` metadata so these three run and display first.

### CASF-2016
- **285 complexes** from PDBbind Core Set 2016
- Four evaluation powers: scoring, ranking, docking, screening
- Download: [pdbbind.org.cn/casf.asp](http://www.pdbbind.org.cn/casf.asp)
- FlexAIDdS targets: r = 0.88, docking power = 82%, EF1% = 12.5

### ITC-187
- **187 complexes** with full ITC thermodynamic decomposition (ΔG, ΔH, −TΔS)
- Primary validation for entropy-aware scoring
- FlexAIDdS targets: r = 0.93, RMSE = 1.4 kcal/mol
- **Entropy rescue rate target: 72%**

### DUD-E 37
- **37 targets** × 3 receptor states (holo / apo / AlphaFold2)
- Cross-docking VS enrichment
- Download: [dude.docking.org](https://dude.docking.org/)

### lsd.docking.org 11-Target Subset
- **11 prospective targets** from ultra-large library docking campaigns
- Confirmed experimental hits provide ground truth
- Access: [lsd.docking.org](https://lsd.docking.org/) (academic registration)

### NRGlab Psychopharmacology 23-Target
- **23 CNS targets**: GPCRs, monoamine transporters, ion channels
- Internal NRGlab dataset — set `FLEXAIDDS_BENCHMARK_DATA` to data root
- FlexAIDdS target: **92% entropy rescue rate**

### MUV
- **17 targets** from Maximum Unbiased Validation
- Actives and decoys share identical physicochemical distributions
- Download: [pharmvip.com](https://www.pharmvip.com/muv)

### ERDS Specificity Z-Score
- **15 targets** × 50,000 Enamine REAL background compounds per target
- Validates that ΔS correction improves signal/noise rather than shifting scores
- Z-score target: < −2.5

---

## Metrics

| Metric | Function | Description |
|--------|----------|-------------|
| `entropy_rescue_rate` | `entropy_rescue_rate()` | **Primary ΔS discriminator** — % complexes where entropy correction rescues the crystal pose from rank > 3 to rank ≤ 3 |
| `ef_1pct` | `enrichment_factor(fraction=0.01)` | Enrichment factor at 1% of ranked list |
| `ef_5pct` | `enrichment_factor(fraction=0.05)` | Enrichment factor at 5% |
| `log_auc` | `log_auc()` | Logarithmic AUC (early enrichment emphasis) |
| `scoring_power_pearson_r` | `scoring_power()` | Pearson r vs experimental ΔG |
| `scoring_power_rmse` | `scoring_power()` | RMSE vs experimental ΔG (kcal/mol) |
| `docking_power_top1` | `docking_power(top_n=1)` | % targets with top-1 pose RMSD < 2.0 Å |
| `docking_power_top3` | `docking_power(top_n=3)` | % targets with any top-3 pose RMSD < 2.0 Å |
| `target_specificity_zscore` | `target_specificity_zscore()` | Z-score of binder scores vs random background |
| `hit_rate_top10` | `hit_rate_top_n(n=10)` | Fraction of actives in top-10 predictions |

All scalar metrics support **95% bootstrap CIs** via `bootstrap_ci()`.

---

## CI/CD Integration

### Tier-1 (`.github/workflows/benchmark-tier1.yml`)

- Triggers on every PR to `master`
- 10-minute hard timeout
- Posts results as a sticky PR comment
- Build is cached (re-used across PRs with identical source hash)
- Binary auto-detection: falls back to `--dry-run` if binary is absent

### Tier-2 (`.github/workflows/benchmark-tier2.yml`)

- Triggers on merge to `master` + weekly Sunday cron
- Full 6-hour window, all datasets
- Dataset files cached from Zenodo/object storage
- Regression check: exits non-zero if any metric drops > 5% below baseline
- Report committed as a commit comment on `master`
- Optional MPI multi-node job (disabled until HPC runner is configured)

---

## Adding a New Dataset

1. Create `benchmarks/datasets/<slug>.yaml` following the existing templates.
2. Populate `targets`, `metrics`, and `expected_baselines`.
3. Set `tier1_subset_size` (≤ 10 for fast CI runs).
4. Register a Zenodo DOI if the dataset is publicly distributable.
5. Update `benchmark-tier2.yml` to add a cache step for the new dataset.

---

## Environment Variables

| Variable | Description |
|----------|-------------|
| `FLEXAIDDS_BINARY` | Path to `FlexAID` executable |
| `FLEXAIDDS_BENCHMARK_DATA` | Root directory for downloaded dataset files |

---

## Reproducibility Policy

Each benchmark bundle in this directory must reach at minimum
`reproducibility_level: replayable-from-repository-artifacts` before
scientific claims can be based on its outputs.  See `smoke/manifest.yaml`
for the template structure.
