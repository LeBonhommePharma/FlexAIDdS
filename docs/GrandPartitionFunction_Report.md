# Grand Partition Function for Competitive Drug Binding

## Technical Report — FlexAIDdS v2.0

**Author**: Louis-Philippe Morency, PhD Candidate
**Affiliation**: Université de Montréal, NRGlab
**Date**: 2026-04-14
**Document scope**: Analysis of features, novelties, and drug discovery applications of the `GrandPartitionFunction` module

---

## 1. Executive Summary

The `GrandPartitionFunction` module implements a **grand canonical ensemble** framework for competitive ligand binding within the FlexAIDdS molecular docking engine. Unlike conventional docking tools that score ligands independently, this module treats the receptor binding site as a thermodynamic system where multiple ligands compete for occupation. The approach yields concentration-dependent binding probabilities, pairwise selectivity ratios, and thermodynamically rigorous free energy rankings — all derived from partition functions computed over conformational ensembles sampled by the genetic algorithm.

Key results:

- **Competitive selectivity** computed as a first-class thermodynamic quantity (not post-hoc score subtraction)
- **Concentration-dependent occupancy** via fugacity/activity terms (`z_i = c_i / c°`)
- **Numerically stable** log-space arithmetic handling partition functions spanning 10^500+
- **Cross-ligand knowledge transfer** through Bayesian conformer priors and binding-site hotspots
- **Thread-safe concurrent accumulation** supporting parallel multi-ligand docking campaigns

---

## 2. Theoretical Foundation

### 2.1 Grand Canonical Ensemble

The binding site is modeled as a system that can be in one of N+1 states: empty (apo) or occupied by one of N ligands. The grand partition function is:

```
Ξ = 1 + Σᵢ zᵢ·Zᵢ
```

where:

| Symbol | Definition | Units |
|--------|-----------|-------|
| Ξ | Grand partition function | dimensionless |
| 1 | Empty site contribution | dimensionless |
| zᵢ | Fugacity = cᵢ / c° | dimensionless |
| cᵢ | Ligand i concentration | M |
| c° | Standard state concentration | 1 M |
| Zᵢ | Canonical partition function of ligand i | dimensionless |

The canonical partition function Zᵢ is computed by `StatMechEngine` from the Boltzmann-weighted conformational ensemble:

```
Zᵢ = Σ_j exp(−β·E_j) · w_j
```

where the sum runs over all sampled poses j with energies E_j and weights w_j, and β = 1/(k_B·T).

### 2.2 Thermodynamic Observables

From Ξ, all thermodynamic observables follow:

**Binding probability** (concentration-dependent):

```
p(ligand_i) = zᵢ·Zᵢ / Ξ = exp(ln zᵢ + ln Zᵢ − ln Ξ)
```

**Empty probability** (apo site):

```
p(empty) = 1 / Ξ = exp(−ln Ξ)
```

**Intrinsic Helmholtz free energy** (concentration-independent):

```
F_i = −kT · ln Zᵢ    [kcal/mol]
```

**Binding free energy** (vs. reference state):

```
ΔG_bind = F_bound − F_ref
```

**Selectivity ratio** (overflow-safe):

```
S(A/B) = (z_A·Z_A) / (z_B·Z_B) = exp(ln Z_A − ln Z_B + ln(c_A/c_B))
```

**Log-selectivity** (apparent, concentration-weighted):

```
ln S(A/B) = ln Z_A − ln Z_B + ln(c_A/c_B)
```

**Intrinsic log-selectivity** (concentration-independent):

```
ln S_int(A/B) = ln Z_A − ln Z_B
```

This quantity is useful for SAR/potency comparisons where concentration effects are irrelevant.

### 2.3 Log-Space Arithmetic

All partition function values are stored as `ln Z`. The log-sum-exp identity computes `ln Ξ` without evaluating any `exp()` argument larger than `max_val`:

```
ln Ξ = max_val + ln(Σᵢ exp(xᵢ − max_val))
```

where the set {xᵢ} includes 0 (empty site) and all `ln(zᵢ·Zᵢ)` values. This guarantees numerical stability for partition functions spanning hundreds of orders of magnitude.

---

## 3. Implementation Architecture

### 3.1 Class Hierarchy

```
TargetServer
├── GrandPartitionFunction    (owns GPF state + mutex)
├── TargetKnowledgeBase       (conformer priors, hotspots)
├── TargetValidation          (structure validation)
└── Session management        (atomic session IDs)
```

`GrandPartitionFunction` owns its own `std::mutex` — all public methods are thread-safe. `TargetServer` holds a separate mutex only for knowledge-base operations. This avoids nested locking and deadlock risk.

### 3.2 Data Storage

Each ligand is stored as a `LigandEntry`:

```cpp
struct LigandEntry {
    double log_Z;    // intrinsic ln(Z_i)
    double log_c;    // ln(c_i/c°) — stored directly to avoid subtraction drift
    double log_zZ;   // ln(z_i · Z_i) = log_c + log_Z (cached for performance)
};
```

The triple-field design separates:
- **Intrinsic quantities** (`log_Z`): used for free energy, ranking — concentration-independent
- **Concentration** (`log_c`): stored directly to prevent floating-point drift during repeated `overwrite_ligand` calls
- **Grand canonical quantities** (`log_zZ`): used for probabilities, selectivity, Ξ — concentration-dependent

### 3.3 Ligand Update Semantics

| Operation | Effect on log_Z | Effect on log_zZ | Use case |
|-----------|----------------|-------------------|----------|
| `add_ligand` | Set | Set (= log_c + log_Z) | New ligand |
| `overwrite_ligand` | Replace | Replace (preserves log_c) | Re-dock with better protocol |
| `merge_ligand` | log_sum_exp | log_sum_exp | Combine independent ensembles |
| `remove_ligand` | Removed | Removed | Discard ligand |

**Merge** computes `ln(Z_old + Z_new)` via log-sum-exp, correctly combining independent conformational ensembles without double-counting. Concentration is preserved through both overwrite and merge.

### 3.4 Numerical Stability Guarantees

| Scenario | Protection |
|----------|-----------|
| Partition functions spanning 10^500 | Log-space storage, log-sum-exp |
| Selectivity ratios > 10^300 | `log_selectivity()` returns raw Δ; `selectivity()` returns ±Inf/0.0 beyond ±700 |
| Empty ensemble | Returns ln Ξ = 0, p(empty) = 1.0 |
| Extreme concentration ratios | Absorbed into log_zZ, no overflow |
| No heap allocation for Ξ | Single-pass, no std::vector |
| Physically impossible concentrations | Values > 1000 M rejected with `std::invalid_argument` |

### 3.5 Public API Summary

| Method | Returns | Description |
|--------|---------|-------------|
| `add_ligand(name, log_Z, c_M)` | void | Register ligand with partition function and concentration |
| `add_ligand(name, engine, c_M)` | void | Register ligand from StatMechEngine |
| `add_or_overwrite(name, log_Z, c_M)` | void | Atomic insert-or-replace (TOCTOU-safe) |
| `overwrite_ligand(name, new_log_Z)` | void | Replace Z, preserve concentration |
| `merge_ligand(name, new_log_Z)` | void | Combine ensembles via log-sum-exp |
| `remove_ligand(name)` | void | Remove ligand from ensemble |
| `log_Xi()` | `double` | ln(Ξ) — grand partition function |
| `binding_probability(name)` | `double` | p(ligand bound), concentration-weighted |
| `empty_probability()` | `double` | p(apo site) = 1/Ξ |
| `F_bound(name)` | `double` | Helmholtz free energy of bound ensemble (kcal/mol) |
| `delta_G_bind(name, F_ref)` | `double` | ΔG_bind = F_bound − F_ref |
| `selectivity(a, b)` | `double` | Apparent selectivity ratio (concentration-weighted) |
| `log_selectivity(a, b)` | `double` | ln S_apparent(a/b), overflow-safe |
| `log_intrinsic_selectivity(a, b)` | `double` | ln(Z_A/Z_B), concentration-independent |
| `rank()` | `vector<LigandRank>` | Ligands sorted by ΔG ascending |
| `all_log_Z()` | `vector<pair>` | All ligand names with intrinsic ln Z |
| `all_log_zZ()` | `vector<pair>` | All ligand names with concentration-weighted ln(z·Z) |

All query methods are `[[nodiscard]]`. Copy and move constructors are `= delete` (mutex ownership).

---

## 4. Drug Discovery Applications

### 4.1 Selectivity Screening

Traditional workflow: dock ligand A to target T, dock ligand A to off-target O, compare scores post-hoc.

GPF workflow: dock all ligands to target T in one campaign, compute pairwise selectivity directly.

**Advantage**: Selectivity accounts for the full conformational ensemble of each ligand-target pair, not just the best pose. A ligand that adopts a single favorable pose has different selectivity than one that populates ten moderately favorable poses — the entropy contribution is captured automatically.

### 4.2 Concentration-Dependent Occupancy

Setting ligand concentrations to physiological values enables:

- **Dose-response prediction**: Sweep c_i to generate occupancy curves
- **Therapeutic window estimation**: Compare occupancy at therapeutic vs toxic concentrations
- **Combination therapy modeling**: Multiple drugs at different concentrations competing for the same site
- **Tissue-specific predictions**: Different drug concentrations in plasma vs CSF vs intracellular

Example: Drug A (1 nM affinity, 10 μM plasma concentration) vs Drug B (100 pM affinity, 100 nM plasma concentration) — the GPF correctly predicts which drug actually occupies the binding site in vivo.

### 4.3 Lead Optimization with Ensemble Feedback

The `TargetServer` accumulates knowledge across ligand docking sessions:

1. **Conformer priors**: Bayesian posterior over receptor conformations weighted by binding probability. Well-binding ligands "vote" on the receptor conformation, improving priors for subsequent ligands.

2. **Binding site hotspots**: Energy-weighted grid points accumulated across all sessions, building a consensus pharmacophore map.

3. **Iterative refinement**: `overwrite_ligand` for re-docking with better protocols; `merge_ligand` for combining expanded ensembles.

This creates a positive feedback loop: better priors → better docking → better priors.

### 4.4 Virtual Screening Campaigns

For large-scale virtual screening:

1. Dock all candidates (thread-safe concurrent sessions)
2. Rank by intrinsic free energy (concentration-independent)
3. Apply physiological concentrations to compute actual occupancy
4. Compute selectivity matrices for lead candidates
5. Re-dock top hits with finer grids (overwrite or merge)

The thread-safe `create_session()` / `register_result()` pattern supports parallel docking of 50+ ligands with automatic GPF accumulation.

---

## 5. Comparison with Existing Tools

| Feature | FlexAIDdS GPF | AutoDock Vina | Glide | GOLD | CADD (general) |
|---------|:---:|:---:|:---:|:---:|:---:|
| Grand canonical ensemble | Yes | No | No | No | No |
| Concentration-dependent occupancy | Yes | No | No | No | Rare |
| Selectivity as first-class output | Yes | No | No | No | No |
| Partition function from ensemble | Yes | No | Partial | No | Rare |
| Entropy from conformational diversity | Yes | No | Partial | No | Rare |
| Cross-ligand knowledge transfer | Yes | No | No | No | No |
| Log-space numerical stability | Yes | N/A | N/A | N/A | N/A |
| Thread-safe concurrent docking | Yes | Limited | No | No | Rare |
| Open source | Yes (Apache-2.0) | Yes (GPL) | No (commercial) | No (commercial) | Varies |

**Key differentiator**: Most docking tools optimize for the single best pose of a single ligand. FlexAIDdS's GPF treats the full conformational ensemble of all competing ligands simultaneously as a grand canonical thermodynamic system.

---

## 6. Limitations and Future Directions

### 6.1 Current Limitations

| Limitation | Impact | Potential solution |
|-----------|--------|-------------------|
| Single binding site | Cannot model multi-site targets (allosteric + orthosteric) | Product of grand partition functions |
| No cooperativity | Assumes independent single-site binding | Multi-site Hamiltonian |
| Fixed temperature | Cannot sweep temperature for melting curves | Re-instantiate GPF per temperature |
| Implicit solvent in "1" | Empty site does not model explicit water competition | WaterMap-style solvent thermodynamics |
| No ligand-ligand interactions | Cannot model self-competition or dimerization | Higher-order terms in Ξ |
| Standard state 1 M | Dilute nM–pM concentrations dominated by log(c) term | Physically correct; user must account for it |

### 6.2 Proposed Extensions

1. **Multi-site grand partition function**: `Ξ = Π_sites Ξ_site` for targets with multiple independent binding sites, enabling allosteric modulator screening.

2. **Temperature sweeps**: Integrate with existing parallel tempering framework to generate full van't Hoff plots (ln K vs 1/T) for entropy/enthalpy decomposition.

3. **Water displacement scoring**: Replace the "1" (empty site) with `Z_water = Π_sites exp(−β·G_water)` to model the thermodynamic cost of displacing ordered water molecules.

4. **Time-dependent occupancy**: Combine with pharmacokinetic concentration-time profiles [c(t)] to predict dynamic target engagement.

5. **Multi-target selectivity matrices**: Extend to compute selectivity across multiple targets simultaneously (one GPF per target, cross-target comparison).

---

## 7. Software Engineering Summary

### 7.1 Dependencies

- **Internal**: `statmech.h` (StatMechEngine, kB_kcal constant), `TargetKnowledgeBase.h`, `TargetValidation.h`, `flexaid.h` (FA_Global, atom, resid)
- **External**: C++26 standard library only (`<cmath>`, `<mutex>`, `<unordered_map>`, `<vector>`, `<string>`, `<algorithm>`, `<limits>`, `<stdexcept>`)
- **No external dependencies**: No Eigen, Boost, or third-party libraries required for GPF functionality

### 7.2 Thread Safety

- All public methods of `GrandPartitionFunction` are guarded by `std::mutex`
- `TargetServer` uses separate mutex for knowledge-base operations (avoids nested locks)
- Session IDs via `std::atomic<int>`
- No deadlock risk (single-lock ownership pattern)

### 7.3 Test Coverage

| Test file | Tests | Coverage |
|-----------|-------|----------|
| `test_grand_partition.cpp` | 29 | Construction, single/multi ligand, competitive binding, equal ligands, ranking, overwrite/merge, remove, error handling, extreme values, StatMech integration, concentration, log-selectivity, intrinsic selectivity, all_log_zZ, impossible concentration guard, ΔG_bind |
| `test_target_server.cpp` | 12 | Construction, validation, session management, GPF integration, re-docking, knowledge base, concurrent registration |
| `test_multi_site_gpf.cpp` | 11+ | Independent/coupled sites, cooperativity (pos/neg), occupancy, marginal probs (unit tests only; full multi-cleft wiring is P5 follow-up) |

### P3+ Emission: CSV/JSON + --conc + Competition YAML

- Python `DatasetRunner` (and CLI `--conc` / `default_conc_M`) + per-ligand `ligand_concs` from YAML emit:
  - `*_grand_summary.csv` with columns: ligand, log_Z, conc_M, log_Xi, p_bind
  - grand_summary embedded in per-dataset JSON reports + markdown summary tables
- C++ binary: `--conc` / `--concentration` sets default for TargetServer sessions; YAML competition_sets / ligands with `conc_M` parsed in runner.
- `benchmarks/datasets/competition_example.yaml` provides per-ligand `conc_M` examples (nM-uM) + expected for validation harness.
- `compute_grand_partition` (Python) + C++ GPF produce matching Ξ, p_bind etc.
- Sidecar/REMARK: [GRAND] lines and grand_*.json when competitive context active (from prior chunks).
- Repro: use `grand_calibrate.py --synthetic` (in benchmarks/grand_synthetic/) and exact fixtures (multi_ligand_exact.json verified 3L equal/varied conc cases).

MultiSiteGPF (product over sites + coupling matrix) remains available for future multi-cleft/allosteric (non-breaking; not yet auto-wired from cleft detection into campaigns — documented follow-up).

See also: GPF_IMPLEMENTATION_PLAN.md (P3/P5), runner.py _save + grand emission paths, competition_example.yaml.

---

## 8. Conclusions

The `GrandPartitionFunction` module provides a thermodynamically rigorous framework for competitive ligand binding that goes beyond conventional single-ligand, single-pose docking. Its key contributions to the drug discovery workflow are:

1. **Principled selectivity prediction** grounded in statistical mechanics rather than score subtraction
2. **Physiologically relevant predictions** through concentration-dependent binding probabilities
3. **Ensemble-level thermodynamics** capturing entropic contributions from conformational diversity
4. **Scalable concurrent screening** with thread-safe accumulation across parallel docking sessions
5. **Iterative refinement** through knowledge accumulation across ligand campaigns

The implementation maintains numerical stability across 1000+ orders of magnitude through log-space arithmetic, requires no external dependencies beyond the C++ standard library, and integrates seamlessly with the existing FlexAIDdS pipeline through the `TargetServer` interface.

---

*Generated for FlexAIDdS v2.0 — Le Bonhomme Pharma*
