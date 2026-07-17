# Theory: Canonical (NVT) vs Grand-Canonical (μVT) in FlexAIDdS

```
// ============================================================================
//  MAKE COMPETITIVE BINDING GREAT AGAIN
// ============================================================================
```

> **Code map**
> - NVT \(Z\): `LIB/statmech.*` (`StatMechEngine`)
> - Abstract surface: `LIB/PartitionFunctionBase.h`
> - Multi-N + competitive μVT façade: `LIB/GrandCanonicalEngine.*`
> - Single-site competitive Ξ (production): `LIB/GrandPartitionFunction.*`
> - Multi-site: `LIB/MultiSiteGPF.*`
> - Audit: `docs/audit/grand_canonical_muVT_phase0.md`

## One-paragraph summary

Docking samples poses at fixed particle number (one ligand geometry in the
site). That is a **canonical** ensemble problem: compute \(Z(N{=}1,V,T)\) from
the GA pose energies. Competitive pharmacology needs **concentrations** and
**multiple species**, which is a **grand-canonical** problem: form
\(\Xi(\mu,V,T)\) from per-occupancy or per-ligand \(Z\) values and fugacities
\(\lambda = e^{\beta\mu}\) (or \(z_i = c_i/c^\circ\)).

## Canonical NVT (what StatMechEngine does)

\[
Z(N,V,T) = \sum_i n_i\, e^{-\beta E_i}, \qquad
F = -k_BT\ln Z, \qquad
S = (\langle E\rangle - F)/T
\]

- Implementation: `statmech::StatMechEngine` (`LIB/statmech.h`).
- Numerics: log-sum-exp via `flexaids::log_sum_exp_dispatch`
  (Shannon-stabilized kernel; AVX-512 / Metal / scalar).
- Adapter: `flexaids::CanonicalPartitionAdapter` exposes the same
  `PartitionFunctionBase` surface without changing `StatMechEngine` ABI.
- Used by `BindingMode` for per-mode Helmholtz free energy and the
  thermodynamic ledger (`ThermodynamicBreakdown`).
- **GA ranking** uses the CF/contact-function scoring proxy, not \(F\).

## Grand-canonical μVT

### (A) Multi-N occupancy — `GrandCanonicalEngine`

\[
\Xi(\mu,V,T) = \sum_{N=0}^{N_{\max}} \lambda^{N}\, Z(N,V,T), \qquad
\lambda = e^{\beta\mu}
\]

Observables:

| Quantity | Formula | API |
|----------|---------|-----|
| \(\ln\Xi\) | \(\mathrm{lse}_N(N\ln\lambda + \ln Z_N)\) | `log_Xi_multiN()` |
| \(\langle N\rangle\) | \((1/\beta)\,\partial\ln\Xi/\partial\mu\) | `mean_N_multiN()` |
| \(\mathrm{Var}(N)\) | \(\langle N^2\rangle - \langle N\rangle^2\) | `var_N_multiN()` |
| \(p(N)\) | \(\lambda^N Z_N / \Xi\) | `occupancy_probability(N)` |

Outer \(N\) summation is OpenMP-parallel for large \(N_{\max}\); the log-sum-exp
kernel is the same Shannon dispatch used for NVT \(Z\).

### (B) Multi-species competitive binding — `GrandPartitionFunction`

For a single site that is empty or occupied by one of \(M\) ligands:

\[
\Xi = 1 + \sum_{i=1}^{M} z_i Z_i, \qquad
z_i = \frac{c_i}{c^\circ},\quad c^\circ = 1\,\mathrm{M}
\]

Equivalently (\(\mu_i = k_BT\ln(c_i/c^\circ)\), \(\mu^\circ=0\)):

\[
\Xi(\boldsymbol{\mu},V,T) = 1 + \sum_i e^{\beta\mu_i} Z_i.
\]

| Quantity | Formula | API |
|----------|---------|-----|
| \(\ln\Xi\) | log-sum-exp\((0,\ln(z_i Z_i))\) | `log_Xi()` |
| \(p_i\) | \(z_i Z_i / \Xi\) | `binding_probability(name)` |
| \(p_\mathrm{empty}\) | \(1/\Xi\) | `empty_probability()` |
| \(\langle N\rangle\) | \(1 - p_\mathrm{empty}\) | `mean_occupancy()` / `mean_N()` |
| Selectivity (apparent) | \((z_A Z_A)/(z_B Z_B)\) | `selectivity(A,B)` |
| Selectivity (intrinsic) | \(Z_A/Z_B\) | `log_intrinsic_selectivity` |
| Mixing entropy | \(-k_B\sum_\alpha p_\alpha\ln p_\alpha\) | `mixing_entropy()` |
| Ligand entropy collapse | \(1 - S_\mathrm{lig}/\ln M\) | `ligand_entropy_collapse()` |

`GrandCanonicalEngine` **composes** `GrandPartitionFunction` for channel (B)
and does not reimplement that Ξ math.

### Analogy to Shannon entropy collapse

When one ligand dominates the bound population, the Shannon mixing entropy
over species crashes toward zero — the same mathematical signature as:

- GA diversity collapse (population entropy → 0 under premature convergence)
- Binding-mode configurational entropy collapse when poses concentrate

`ligand_entropy_collapse()` is **diagnostic only** and never used for ranking.

## Architecture (v1)

```
PartitionFunctionBase
├── CanonicalPartitionAdapter   (NVT: wraps log_Z / StatMechEngine)
└── GrandCanonicalEngine        (μVT)
      ├── multi-N channel       Σ_N λ^N Z_N   [Shannon LSE + OpenMP]
      └── competitive channel   GrandPartitionFunction
            └── MultiSiteGPF / TargetServer (site registry)
```

## NRGsuite / multi-ligand concentration API

After ligands are registered (docked → \(\ln Z_i\) known), retitrate without
re-docking:

```text
set_concentration([L1, L2, ...])   # conceptual NRGsuite command
```

C++:

```cpp
flexaids::GrandCanonicalEngine gce(310.0);
gce.add_competitive_ligand("fentanyl", log_Z_f, 1e-9);
gce.add_competitive_ligand("naloxone", log_Z_n, 1e-6);
gce.set_competitive_concentrations({"fentanyl", "naloxone"}, {1e-9, 1e-4});
// or TargetServer:
server.set_concentration("naloxone", 1e-4);
```

Python (pure helper always; C++ `_core` when built with μVT sources):

```python
from flexaidds.grand_canonical import set_concentration, CompetitiveSite, plot_occupancy_curve
site = CompetitiveSite(temperature_K=310.0)
site.add("fentanyl", log_Z=..., c_M=1e-9)
site.add("naloxone", log_Z=..., c_M=1e-6)
set_concentration(site, {"naloxone": 1e-4})
print(site.binding_probability("naloxone"), site.mean_N())
curve = site.occupancy_vs_concentration("naloxone", [1e-9, 1e-7, 1e-5, 1e-3])
plot_occupancy_curve(curve, title="Naloxone titration (μVT)")
```

## Scoring-loop boundary (Phase 2)

- `VoronoiCFBatch` / GA: **one ligand geometry per chromosome** (CF proxy).
- Multi-ligand GA individuals may carry a `LigandVector` (names + concentrations)
  as **metadata** for post-hoc Ξ; fugacity is **not** folded into CF scores
  without an explicit feature flag + tests (`AGENTS.md`).
- `BindingMode` can tag `ligand_concentration_M` / species name for ledger feed.

## CMake

```bash
cmake -B build -DBUILD_TESTING=ON -DFLEXAIDS_ENABLE_MUVT=ON
cmake --build build -j
ctest --test-dir build -R 'Grand(Partition|Canonical)' --output-on-failure
```

## Ranking guardrail

- **Search / pose order:** CF Voronoi proxy (`VoronoiCFBatch`, GA).
- **Post-hoc ledger:** NVT `StatMechEngine` + μVT `GrandCanonicalEngine` / GPF.
- Changing concentrations updates \(\Xi\) and \(p_i\), **not** the GA
  chromosome ranking, unless an experimental reweight path is explicitly
  enabled and tested.

## Toy validation systems (synthetic log_Z)

| System | Ligands | Test |
|--------|---------|------|
| MOR | fentanyl + naloxone | `GrandPartition.MORnaloxone_NVT_vs_muVT`, `GrandCanonicalEngine.Benchmark_MOR_*` |
| 5-HT2A | 5-MeO-DMT + 5-HT | `GrandCanonicalEngine.Benchmark_5HT2A_*` |

These use synthetic partition functions — **not** experimental \(K_i\) claims.
Compare predicted \(\langle N\rangle\) / apparent \(K_i\) trends only after
real docked \(\ln Z_i\) and calibration (`AffinityCalibration`).

## References in-repo

- `LIB/PartitionFunctionBase.h` — abstract NVT/μVT surface
- `LIB/GrandCanonicalEngine.h` — multi-N + competitive façade
- `LIB/GrandPartitionFunction.h` — competitive single-site Ξ
- `tests/test_grand_partition.cpp`, `tests/test_grand_canonical_engine.cpp`
- `docs/dev/thermo_source_map.md`, `docs/audit/grand_canonical_muVT_phase0.md`
