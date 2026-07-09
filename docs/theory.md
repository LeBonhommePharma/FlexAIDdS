# Theory: Canonical (NVT) vs Grand-Canonical (μVT) in FlexAIDdS

> **Source of truth for code:** `LIB/GrandPartitionFunction.*` **is** the
> single-site grand-canonical (μVT) competitive-binding engine. There is no
> separate parallel `GrandCanonicalEngine` class. Canonical NVT lives in
> `LIB/statmech.*` (`StatMechEngine`). Multi-site extensions use
> `LIB/MultiSiteGPF.*`.

## One-paragraph summary

Docking samples poses at fixed particle number (one ligand molecule in the
site). That is a **canonical** ensemble problem: compute \(Z(N{=}1,V,T)\) from
the GA pose energies. Competitive pharmacology needs **concentrations** and
**multiple species**, which is a **grand-canonical** problem: form
\(\Xi(\mu,V,T)\) from the per-ligand \(Z_i\) and fugacities
\(z_i = c_i/c^\circ\).

## Canonical NVT (what StatMechEngine does)

\[
Z(N,V,T) = \sum_i n_i\, e^{-\beta E_i}, \qquad
F = -k_BT\ln Z, \qquad
S = (\langle E\rangle - F)/T
\]

- Implementation: `statmech::StatMechEngine` (`LIB/statmech.h`).
- Numerics: log-sum-exp via `flexaids::log_sum_exp_dispatch`
  (Shannon-stabilized kernel; AVX-512 / Metal / scalar).
- Used by `BindingMode` for per-mode Helmholtz free energy and the
  thermodynamic ledger (`ThermodynamicBreakdown`).
- **GA ranking** uses the CF/contact-function scoring proxy, not \(F\).

## Grand-canonical μVT (what GrandPartitionFunction does)

For a single site that is empty or occupied by one of \(M\) ligands:

\[
\Xi = 1 + \sum_{i=1}^{M} z_i Z_i, \qquad
z_i = \frac{c_i}{c^\circ},\quad c^\circ = 1\,\mathrm{M}
\]

Equivalently with chemical potentials (ideal solution,
\(\mu_i = \mu_i^\circ + k_BT\ln(c_i/c^\circ)\)):

\[
\Xi(\boldsymbol{\mu},V,T) = 1 + \sum_i e^{\beta\mu_i} Z_i
\quad(\mu^\circ=0\text{ convention in code}).
\]

Observables:

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

### Analogy to Shannon entropy collapse

When one ligand dominates the bound population, the Shannon mixing entropy
over species crashes toward zero — the same mathematical signature as:

- GA diversity collapse (population entropy → 0 under premature convergence)
- Binding-mode configurational entropy collapse when poses concentrate

See also README “entropy dualism” notes. The metric
`ligand_entropy_collapse()` is **diagnostic only** and never used for ranking.

## What is *not* reimplemented

| Component | Role | Do not duplicate |
|-----------|------|------------------|
| `StatMechEngine` | Canonical \(Z\), \(F\), \(S\), \(C_v\) | — |
| `GrandPartitionFunction` | Single-site \(\Xi\), concentrations, selectivity | No second GC engine |
| `MultiSiteGPF` | Multiple sites + cooperativity \(\omega\) | — |
| `TargetServer` | Owns receptor + live \(\Xi\) across ligands | — |
| UnifiedHardwareDispatch LSE | Numerics for both NVT and μVT | — |

## NRGsuite / multi-ligand concentration API

After ligands are registered (docked → \(\ln Z_i\) known), retitrate without
re-docking:

```text
set_concentration([L1, L2, ...])   # conceptual NRGsuite command
```

C++:

```cpp
target::GrandPartitionFunction& xi = server.grand_partition();
xi.set_concentrations({"fentanyl", "naloxone"}, {1e-9, 1e-6});  // M
// or:
server.set_concentration("naloxone", 1e-4);
```

Python (pure helper; works without `_core` GC bindings):

```python
from flexaidds.grand_canonical import set_concentration, CompetitiveSite
site = CompetitiveSite(T=310.0)
site.add("fentanyl", log_Z=..., c_M=1e-9)
site.add("naloxone", log_Z=..., c_M=1e-6)
set_concentration(site, {"naloxone": 1e-4})
print(site.binding_probability("naloxone"), site.mean_N())
```

## Ranking guardrail

- **Search / pose order:** CF Voronoi proxy (`VoronoiCFBatch`, GA).
- **Post-hoc ledger:** NVT `StatMechEngine` + μVT `GrandPartitionFunction`.
- Changing concentrations updates \(\Xi\) and \(p_i\), **not** the GA
  chromosome ranking, unless an experimental reweight path is explicitly
  enabled and tested.

## References in-repo

- `LIB/GrandPartitionFunction.h` — API contract and thermodynamic comments
- `tests/test_grand_partition.cpp` — analytic + MOR/naloxone toy NVT vs μVT
- `docs/dev/thermo_source_map.md` — file map
- `docs/dev/thermo_invariants.md` — numerical invariants
