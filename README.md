<div align="center">

# FlexAID∆S

**Thermodynamically-aware molecular docking engine**

[![CI](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml/badge.svg)](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![C++26](https://img.shields.io/badge/C%2B%2B-26-blue.svg)](https://en.cppreference.com/w/cpp/26)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9-3776AB.svg)](https://www.python.org/)
[![Astex-85](https://img.shields.io/badge/Astex--85-78%2F85%20%3D%2091.8%25-brightgreen.svg)](#benchmark-astex-85)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5b00078-blue)](https://doi.org/10.1021/acs.jcim.5b00078)

**[Installation](docs/INSTALLATION.md)** ·
**[User Guide](docs/USERGUIDE.md)** ·
**[Scoring](docs/SCORING.md)** ·
**[Atom Types](docs/ATOM_TYPES.md)** ·
**[Benchmarks](docs/BENCHMARK.md)** ·
**[Changelog](VERSION.md)**

</div>

---

FlexAID∆S is a thermodynamically-aware molecular docking engine for structure-based drug design. Where conventional docking programs report the single lowest-energy pose, FlexAID∆S treats the entire GA-sampled pose population as a statistical ensemble and extracts a free energy — balancing enthalpy against conformational entropy — to identify binding modes that are both energetically favorable *and* physically accessible at finite temperature. The result is a docking engine that is particularly effective on targets where the correct binding mode is not the lowest-enthalpy pose, a systematic failure mode of classical scoring functions.

FlexAID∆S descends from [FlexAID](https://doi.org/10.1021/acs.jcim.5b00078) (Gaudreault & Najmanovich, *J. Chem. Inf. Model.* 2015) and extends it with a thermodynamic scoring layer, corrected atom-type assignments, a physical-realism clash penalty, and a two-gate spread guard that prevents false minima from dominating the ranked output. On the Astex-85 benchmark (85 diverse protein–ligand co-crystal structures), FlexAID∆S achieves **78/85 = 91.8%** at RMSD ≤ 2.0 Å, compared to 75/85 = 88.2% for the original FlexAID.

---

## The Physics

### The Contact Function

The genetic algorithm searches ligand pose space — six rigid-body degrees of freedom plus all rotatable torsions — and evaluates each pose with the **Voronoi contact function** (CF). CF measures shape complementarity by decomposing the molecular surface into Voronoi polyhedra and integrating the contact area between atom pairs, weighted by a 40×40 energy matrix (`MC_st0r5.2_6.dat`) trained on PDB-derived contact statistics across 40 SYBYL atom types. Lower CF is better — a perfect complementary fit between a ligand and its receptor pocket approaches the global minimum.

The total CF for a pose is:

```
CF = CF.com  +  CF.wal  +  CF.sas  +  CF.elec  +  CF.hbond  +  CF.pb_clash  +  CF.con
```

where `CF.com` is the Voronoi contact complementarity, `CF.wal` is the soft-wall steric repulsion (capped at 50 CF units per contact to prevent numerical blow-up), `CF.sas` is an accessible-surface area term, `CF.elec` is an optional electrostatic term, `CF.hbond` a hydrogen-bond term, `CF.pb_clash` the PoseBusters intermolecular clash penalty, and `CF.con` a distance-constraint term.

The energy matrix maps atom-type pairs to statistical potentials derived from contact frequencies in the PDB. A pair that appears more often in real binding sites than in a random background gets a negative entry (stabilizing); one that appears less often gets a positive entry (destabilizing). With 40 atom types and full symmetry, there are 820 unique interaction parameters. See [docs/SCORING.md](docs/SCORING.md) for the full derivation.

### Why min(CF) Is Not Enough

Ranking by the single lowest-CF pose conflates two distinct phenomena. A narrow, deep funnel — one correct binding mode with many similar low-energy neighbors — is genuinely good binding. A broad flat landscape — many structurally diverse poses all scoring near the same CF minimum — is a false minimum: the ligand is weakly complementary to the receptor across a large region of pose space, and the *apparent* lowest CF is a sampling artifact rather than a true thermodynamic minimum. Classical min(CF) cannot distinguish these cases.

The GA amplifies the problem. Because selection pressure drives chromosomes toward the lowest-CF region, a false minimum can accumulate a large sub-population before the algorithm converges, even though the physiologically correct pose exists elsewhere in the landscape and would score higher under min(CF) than the false minimum's best member.

### ΔG_eff: Ensemble Free Energy

FlexAID∆S replaces the single-pose ranking criterion with a **Boltzmann ensemble free energy**, ΔG_eff, computed over the converged GA population:

```
P_i    =  exp(−CF_i / T_eff) / Z        [Boltzmann weight of pose i]
Z      =  Σ_i exp(−CF_i / T_eff)        [partition function]
⟨CF⟩  =  Σ_i P_i · CF_i                [Boltzmann-weighted mean enthalpy proxy]
H      = −Σ_i P_i · ln P_i             [Shannon entropy of the pose distribution, nats]
ΔG_eff =  ⟨CF⟩ − T_eff · H
```

The Shannon entropy term H penalizes a wide, diffuse pose distribution. A ligand with many structurally distinct low-CF poses (broad landscape) has high H, which *increases* ΔG_eff relative to the enthalpy alone, correctly demoting it. A ligand with a narrow, well-converged population centered on one pose (deep funnel) has low H, leaving ΔG_eff close to ⟨CF⟩, correctly promoting it.

The effective temperature T_eff = 0.596 (in CF scoring units) is calibrated so that at room temperature the Boltzmann distribution is neither collapsed to a delta function nor uniformly flat, matching the ISMB 2017 calibration for the FlexAID scoring scale. A second, broader calibration at T = 21 (internal CF units, corresponding to the ISMB 2017 whiteboard convention) is computed as a reporting diagnostic.

In addition to the Shannon term, a **vibrational entropy correction** (TdS_vib) from the tENCoM elastic-network model measures how much the receptor's torsional flexibility changes upon binding, contributing to G_bind:

```
G_bind = T_eff · ⟨CF⟩_raw  −  TdS_shannon  +  TdS_vib
```

Here TdS_shannon is the configurational entropy cost (positive = unfavorable; a more disordered bound-state distribution costs entropy) and TdS_vib is the vibrational entropy change upon binding (typically negative for a ligand that rigidifies the receptor, i.e., stabilizing).

### Thermodynamic Impossibility Gate

When `FLEXAIDDS_THERMO_SCORE=1` is active, FlexAID∆S applies a physics filter: any pose for which ΔH > 0 and ΔS < 0 *simultaneously* is flagged as **thermodynamically impossible**. From ΔG = ΔH − TΔS, if ΔH > 0 and ΔS < 0 then −TΔS > 0 for all T > 0, making ΔG strictly positive at every temperature. Such a configuration cannot bind spontaneously under any physically realizable condition. Poses failing this gate receive a sentinel ΔG_eff = +1000, ensuring they can never be elected rank-0. The ΔS source for this test is the vibrational entropy term TdS_vib (the Shannon population entropy is always ≥ 0 by construction and would make the ΔS < 0 arm unreachable).

---

## What's New vs. FlexAID

| Feature | FlexAID (2015) | FlexAID∆S |
|:--------|:--------------:|:---------:|
| Pose ranking criterion | min(CF) | ΔG_eff = ⟨CF⟩ − T·H |
| Ensemble thermodynamics | ✗ | ✓ Helmholtz F, entropy S, Cv |
| Thermodynamic impossibility gate | ✗ | ✓ ΔH > 0 ∧ ΔS < 0 → sentinel +1000 |
| Intermolecular clash detection | Approximated (23× undercounting) | ✓ Full all-pairs PoseBusters penalty |
| Receptor clash grid | Rebuilt every CF eval | ✓ Loop-invariant hoist (once per dock) |
| Spread guard (false-minima demotion) | ✗ | ✓ Two-gate: distance + frequency + consensus |
| Atom type: N.2 (sp2 imine) | → N.am (donor, wrong sign) | → N.ar (acceptor, correct) |
| Atom type: N.3 (amine) | → N.3/type-8 (dead matrix row) | → N.am/type-11 (live) |
| Atom type: C.1 (sp carbon) | → C.1/type-1 (sparse) | → C.2/type-2 (better sampled) |
| Atom type: I (iodine) | → type-26 (3 live entries) | → BR/type-25 (full halogen row) |
| WAL repulsion | Unbounded (SIGSEGV on extreme clashes) | ✓ Capped at 50 CF units per contact |
| Vibrational entropy (tENCoM) | ✗ | ✓ Torsional elastic-network model |
| Astex-85 success rate | 75/85 = 88.2% | **78/85 = 91.8%** |

---

## Quick Start

### Build

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The build produces:
- `build/FlexAIDdS` — the main docking executable (LTO-optimized)
- `build/FlexAID` — legacy-compatible interface
- `build/tENCoM` — standalone vibrational entropy tool

### Dock a Ligand

```bash
# Flexible docking with default settings
./build/FlexAIDdS receptor.pdb ligand.mol2

# With a JSON configuration file
./build/FlexAIDdS receptor.pdb ligand.mol2 -c config.json

# Rigid-body screening mode (faster, no torsion sampling)
./build/FlexAIDdS receptor.pdb ligand.mol2 --rigid

# Legacy two-file interface (FlexAID compatible)
./build/FlexAID config.inp ga.inp output_prefix
```

Supported ligand formats: MOL2 (SYBYL), SDF/MOL (V2000 and V3000). Supported receptor formats: PDB, mmCIF.

### Example JSON Configuration

```json
{
  "ga": {
    "num_chromosomes": 1000,
    "num_generations": 6000,
    "num_restarts": 5,
    "temperature": 0.596
  },
  "flexibility": {
    "ligand_torsions": true,
    "ring_conformers": true,
    "chiral_centers": true
  },
  "scoring": {
    "pb_clash_weight": 1.0,
    "hbond_weight": -2.5
  },
  "thermodynamics": {
    "temperature": 298,
    "clustering_algorithm": "DP"
  }
}
```

### Interpret Output

FlexAID∆S writes ranked PDB files and a REMARK-annotated summary for each binding mode. The key fields in the `[THERMO]` block:

```
REMARK  G_bind       = -8.41   # ΔG_bind: primary ranking criterion (kcal/mol)
REMARK  H_vct_raw    = -6.23   # Boltzmann-weighted mean CF (enthalpy proxy)
REMARK  TdS_shannon  =  1.94   # Configurational entropy cost (T·H, positive = unfavorable)
REMARK  TdS_vib      = -0.24   # Vibrational entropy gain (negative = stabilizing)
REMARK  dG_eff       = -4.31   # ΔG_eff = <CF> − T_eff·H (at T_eff=0.596)
REMARK  dG_eff_T21   = -5.88   # ΔG_eff at ISMB 2017 calibration T=21
REMARK  binding_regime = enthalpy_driven
```

### Python Analysis

```bash
pip install -e ./python
```

```python
import flexaidds as fd

# Load results from a completed docking run
docking = fd.load_results("path/to/results/")

# Inspect top binding modes
for mode in docking.binding_modes[:3]:
    print(f"Rank {mode.rank}: ΔG = {mode.free_energy:.2f}  RMSD = {mode.rmsd:.2f} Å")

# Compute thermodynamics directly on an energy array
from flexaidds import StatMechEngine
engine = StatMechEngine(temperature=298.0)
engine.add_energies(pose_energies)
th = engine.compute()
print(f"F = {th.free_energy:.3f}  S = {th.entropy:.4f}  Cv = {th.heat_capacity:.4f}")
```

### PyMOL Visualization

Install via Plugin Manager → `pymol_plugin/`. Then:

```
flexaids_load path/to/results/
flexaids_color_boltzmann       # color poses by Boltzmann weight
flexaids_thermo                # show thermodynamic breakdown panel
```

---

## Benchmark: Astex-85

The **Astex Diverse Set** (85 diverse protein–ligand co-crystal complexes, covering 55 therapeutic targets and 12 protein families) is the standard benchmark for evaluating flexible docking accuracy. The success criterion is RMSD ≤ 2.0 Å between the predicted top-ranked pose and the crystallographic reference.

| Engine | Astex-85 Success Rate | Condition |
|:-------|:---------------------:|:----------|
| FlexAID (2015, published) | 75/85 = 88.2% | Cognate redock, no native seed |
| FlexAID∆S (this work) | **78/85 = 91.8%** | Cognate redock, no native seed |

The three additional successes recovered by FlexAID∆S relative to the baseline are targets where min(CF) elected an incorrect binding mode that ΔG_eff correctly demoted. The atom-type fixes (N.2 → N.ar, N.3 → N.am) account for at least one of the three; the PoseBusters clash penalty accounts for another; the spread guard accounts for the third.

The benchmark is fully reproducible. See [docs/BENCHMARK.md](docs/BENCHMARK.md) for dataset provenance, exact commands, and per-target results.

---

## Protocol Flags

FlexAID∆S exposes its scientific innovations as environment-variable flags so each mechanism can be enabled, disabled, or tuned independently for ablation studies and benchmarking.

| Flag | Default | Effect |
|:-----|:-------:|:-------|
| `FLEXAIDDS_THERMO_SCORE` | OFF | Promote ΔG_eff = ⟨CF⟩ − T·H to primary ranking criterion (replaces min CF) |
| `FLEXAIDDS_PB_CLASH_WEIGHT` | `0.0` | PoseBusters intermolecular clash penalty weight; set to `1.0` to enable |
| `FLEXAIDDS_PB_CLASH_RATIO` | `0.75` | Clash threshold as fraction of summed vdW radii |
| `FLEXAIDDS_WAL_COERCIVE` | OFF | Remove WAL_CONTACT_CAP=50 ceiling; deep clashes override CF.com |
| `FLEXAIDDS_WAL_STIFF` | `0` | Override soft-wall stiffness k (default 50) for sweep experiments |
| `FLEXAIDDS_T_EFF` | `0.596` | Effective temperature in CF units for Boltzmann pose weights |
| `FLEXAIDDS_REPORT_T` | `21.0` | Reporting temperature for ISMB 2017 whiteboard diagnostics |
| `FLEXAIDDS_CLUSTER_SPREAD_MAX` | `0.0` | Activate spread guard; set to pocket radius in Å (e.g. `8.0`) |
| `FLEXAIDDS_CLUSTER_POP_MIN_FRACTION` | `0.35` | Minimum rank-0 population fraction below which demotion is eligible |
| `FLEXAIDDS_CLUSTER_CONSENSUS_K` | `3` | Minimum restarts that must agree with rank-0 to veto demotion |
| `FLEXAIDDS_CLUSTER_CONSENSUS_TAU` | `2.0` | RMSD radius (Å) within which a restart head counts as agreeing |
| `FLEXAIDDS_SOFTCORE_WAL` | OFF | Enable soft-core (parabolic) wall instead of r^-12 hard wall |
| `FLEXAIDDS_CONTACTS_EPOCH` | OFF | O(1) contacts buffer clear (epoch counter, vs. O(N) memset) |
| `FLEXAIDDS_EVAL_SCALE_DIHEDRAL` | `1` | GA budget scaling: 1=pop-scale (default), 0=gen-scale, -1=fixed |
| `FLEXAIDDS_RESTARTS` | `5` | Number of independent GA restarts per target |
| `FLEXAIDDS_PARALLEL_RESTARTS` | ON | Launch restart workers concurrently |
| `FLEXAIDDS_USE_SHANNON` | OFF | Enable ShannonThermoStack configurational entropy computation |
| `FLEXAIDDS_RING_FLEX` | OFF | Enable non-aromatic ring pucker sampling (LigandRingFlex) |
| `FLEXAIDDS_HBOND_WEIGHT` | `-2.5` | Hydrogen-bond term coefficient |

---

## Architecture

```
  ┌─────────────────────────────────────────────────────────────────┐
  │                        Input Layer                              │
  │  PDB receptor  │  MOL2 / SDF ligand  │  JSON config            │
  │  (PDB, mmCIF)  │  (V2000, V3000)     │  (or legacy .inp)       │
  └────────────────┬────────────────────────────────────────────────┘
                   │  atom typing (40-type NRGDock)
                   ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │                   Genetic Algorithm (gaboom)                    │
  │  • Population: rigid-body + torsion chromosomes                 │
  │  • Selection pressure: CF cost function (lower = better)        │
  │  • Restarts × parallel workers, configurable budget scaling     │
  └────────────────────────────┬────────────────────────────────────┘
                               │  per-eval scoring
                               ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │              Voronoi Contact Function (vcfunction)              │
  │  Voronoi tessellation → contact areas → energy matrix lookup    │
  │  + soft-wall repulsion  (capped at 50 CF units)                 │
  │  + PoseBusters clash penalty  (all-pairs, loop-invariant grid)  │
  │  + H-bond / GIST / metal-coordination optional terms            │
  └────────────────────────────┬────────────────────────────────────┘
                               │  converged population
                               ▼
  ┌────────────────────────────────────────┐  ┌──────────────────────────┐
  │      ThermodynamicEngine               │  │    Pose Validation        │
  │  ⟨CF⟩ = Σ P_i · CF_i                  │  │  RMSD ≤ 2.0 Å: S1        │
  │  H    = −Σ P_i · ln P_i               │  │  PoseBusters bust: S2    │
  │  ΔG_eff = ⟨CF⟩ − T_eff · H           │  │  BCR (sampling ceiling)  │
  │  G_bind = T·⟨CF⟩_raw − TdS + TdS_vib │  └──────────────────────────┘
  │  Impossibility gate (ΔH>0 ∧ ΔS<0)    │
  └────────────────────────────────────────┘
                               │  ranked binding modes
                               ▼
  ┌────────────────────────────────────────┐  ┌──────────────────────────┐
  │       StatMechEngine (statmech.cpp)    │  │    tENCoM (tENCoM/)       │
  │  Canonical Z, Helmholtz F, Cv         │  │  Torsional ENM: ΔS_vib   │
  │  WHAM, parallel tempering, TI         │  │  Hessian → normal modes   │
  │  AVX-512 / Eigen / OpenMP dispatch    │  │  B-factor vibrational     │
  └────────────────────────────────────────┘  └──────────────────────────┘
                               │
                               ▼
  ┌─────────────────────────────────────────────────────────────────┐
  │                      Output & Analysis                          │
  │  Ranked PDB files  │  REMARK-annotated thermodynamics           │
  │  Python API (flexaidds)  │  PyMOL plugin  │  CLI inspector      │
  └─────────────────────────────────────────────────────────────────┘
```

---

## Building with Optional Features

```bash
# Standard release build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# With GoogleTest unit tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure

# With Python bindings
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --build build --parallel

# Apple Silicon with Metal GPU acceleration
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_METAL=ON

# NVIDIA GPU with CUDA
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_CUDA=ON

# HPC cluster with AVX-512 + OpenMP + MPI
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_AVX512=ON -DFLEXAIDS_USE_OPENMP=ON -DFLEXAIDS_USE_MPI=ON
```

**Compiler requirements:** C++26 (GCC ≥ 14, Clang ≥ 18, Apple Clang ≥ 16/Xcode 16, MSVC ≥ 19.40). CMake ≥ 3.28. Eigen3 recommended.

**Runtime hardware dispatch** (where built): CUDA → Metal → AVX-512 → AVX2 → OpenMP → scalar. No configuration required; the engine selects the fastest available backend at startup.

---

## Testing

```bash
# C++ tests (GoogleTest)
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure

# Python tests (pytest)
cd python && pip install -e . && pytest tests/ -q
```

Tests marked `@requires_core` skip gracefully when the C++ `_core` extension is not built. The CI matrix covers Linux (GCC 14, Clang 18), macOS (Apple Clang), and a Python bindings smoke test.

---

## Publications

If you use FlexAID or FlexAID∆S, please cite:

> Gaudreault F & Najmanovich RJ (2015). FlexAID: Revisiting Docking on Non-Native-Complex Structures. *J. Chem. Inf. Model.* **55**(7):1323–1336. [DOI:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)

Related work:

- Gaudreault F, Morency LP & Najmanovich RJ (2015). NRGsuite. *Bioinformatics* **31**(23):3856–3858. [DOI:10.1093/bioinformatics/btv458](https://doi.org/10.1093/bioinformatics/btv458)
- Frappier V et al. (2015). ENCoM. *Proteins* **83**(11):2073–2082. [DOI:10.1002/prot.24922](https://doi.org/10.1002/prot.24922)
- Morency LP & Najmanovich RJ (2026). FlexAID∆S — *methods manuscript in preparation*

---

## License

Apache License 2.0 — free for academic and commercial use.

Copyright © 2026 Le Bonhomme Pharma · Louis-Philippe Morency

No GPL dependencies. See [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) for the full dependency license matrix.
