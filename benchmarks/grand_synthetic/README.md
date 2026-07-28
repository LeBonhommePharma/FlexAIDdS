# Grand Canonical Synthetic Exact-Solution Fixtures (P4)

**Location**: `benchmarks/grand_synthetic/`

**Purpose**: Provide machine-verifiable ground-truth cases for GrandPartitionFunction
(Ξ, p_bind, occupancy, selectivity) with **known canonical log_Z**, concentrations,
and **analytically computed** expected outputs. These are independent of docking
GA sampling; used to validate harness, Python fallback, C++ GPF, and HW parity.

**Scientific requirements (AGENTS.md + GPF plan)**:
- Analytical exact (use log-sum-exp identical to implementation).
- Reproducibility: pin T=298.0 K, c°=1 M, kB consistent with code (0.001987206 kcal/mol/K).
- Objectivity: report within documented fp tolerance (1e-9 relative for log quantities).
- Later integration: loadable by `test_grand_partition.cpp` extensions, `grand_calibrate.py`,
  `validate_benchmark_results.py` extensions, and Python `flexaidds` grand module.
- Do not claim real ΔG; these test the partition math only.

**How generated**:
- Exact values computed offline with Python math.log / logsumexp equivalent.
- See `compute_grand_exact.py` (to be added in P4 run step) or inline formulas in competition_example.yaml.
- All values use natural log, dimensionless partition functions.

**Usage**:
- Python: load json, feed log_Z + c to pure-Py GPF impl, assert match.
- C++: embed or load in future `TEST(GrandPartition, SyntheticExact)` cases.
- Harness: `python benchmarks/grand_synthetic/grand_calibrate.py --synthetic ...` (or from package; see also top-level scripts/ for related)

**Files**:
- `dual_ligand_exact.json`: basic 2-ligand equal conc + conc-varied cases.
- `multi_ligand_exact.json`: 3+ ligands + extreme ratio for stability. (Used for exact 3L verification: p_bind 0.900/0.090/0.009 at 1uM equal, log_Xi matches hand-calc.)
- `grand_calibrate.py`: harness for synthetic + runner integration.
- `grand_exact_cases.csv`: reference values.

**P3+ outputs in campaigns**: When using --conc or competition YAML, runner emits *_grand_summary.csv (and grand_summary in JSON) with log_Z/conc_M/log_Xi/p_bind per ligand. See DatasetRunner and README top-level for usage.
- `*.csv`: tabular for easy loading in validate/calibrate scripts.
- `README.md`: this file.

**Verification gate (P4)**:
- Fixtures must be loadable and math self-consistent (pA+pB+p_empty == 1.0 within 1e-12).
- Later: run harness dry-run and match to printed expected.

**Temperature & constants**:
T = 298.0 K
RT ≈ 0.592187 kcal mol⁻¹ (for ΔG = -RT lnZ comments only)
c° = 1.0 M

See GPF_IMPLEMENTATION_PLAN.md §P4 and docs/GrandPartitionFunction_Report.md §2 for formulas.
Apache-2.0. Aligns with reproducibility and zero-overclaim rules.