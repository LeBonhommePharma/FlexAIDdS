<div align="center">

# FlexAID∆S

**Thermodynamically-aware molecular docking engine**

[![CI](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml/badge.svg)](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![C++26](https://img.shields.io/badge/C%2B%2B-26-blue.svg)](https://en.cppreference.com/w/cpp/26)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9-3776AB.svg)](https://www.python.org/)
[![Astex-85](https://img.shields.io/badge/Astex--85-unverified%20%7C%20pending%20receipt-lightgrey.svg)](#benchmark-astex-85)
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

FlexAID∆S descends from [FlexAID](https://doi.org/10.1021/acs.jcim.5b00078) (Gaudreault & Najmanovich, *J. Chem. Inf. Model.* 2015) and extends it with a score-space ensemble diagnostic layer, corrected atom-type assignments, a physical-realism clash penalty, and a two-gate spread guard that demotes false minima in the ranked output. Its accuracy on the Astex-85 benchmark is **unverified/pending**: this repository publishes no receipted success rate. See [Benchmark: Astex-85](#benchmark-astex-85).

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

### Score-space ensemble diagnostic

FlexAID∆S can summarize retained GA records with a soft-min/Shannon
transform (legacy field name `dG_eff`):

```
p_i      = exp(−CF_i / T_eff) / Z
H_sample = −Σ_i p_i ln(p_i)
G_tilde  = ⟨CF⟩_p − T_eff · H_sample
```

This is an optimizer-sample, arbitrary-CF-unit diagnostic. It is not a
canonical partition function: the GA records do not define an equilibrium
measure, exact duplication changes `G_tilde`, and CF has no validated kcal/mol
conversion. `T_eff = 0.596` and the legacy `T = 21` are internal score
parameters, not temperatures in kelvin.

tENCoM and Shannon diagnostics are available separately, but the legacy
combined `G_bind`/`deltaG` composition is not a matched association cycle and
currently double-counts configurational entropy. It must not be interpreted as
physical binding free energy, affinity, `Kd`, or `Ki`.

### Diagnostic impossibility predicate

`FLEXAIDDS_THERMO_SCORE=1` computes and prints a legacy `dH > 0 && dS < 0`
sentinel. The inputs are not commensurate physical state differences, and the
sentinel is not consumed by the later exact-CF rescore/clustering boundary. The
flag therefore does not currently enforce a physics filter on the elected pose.

---

## What's New vs. FlexAID

| Feature | FlexAID (2015) | FlexAID∆S |
|:--------|:--------------:|:---------:|
| Pose ranking criterion | min(CF) | Cluster-local soft-β `G̃ = H̃ − T·S̃` over the same CF samples (still CF-bound; no physical energy enters) |
| Ensemble analysis | ✗ | ✓ fail-closed proxy/canonical provenance ledger |
| Impossibility predicate | ✗ | diagnostic only; not wired to final election |
| Intermolecular clash detection | Approximated (23× undercounting) | ✓ Full all-pairs PoseBusters penalty |
| Receptor clash grid | Rebuilt every CF eval | ✓ Loop-invariant hoist (once per dock) |
| Spread guard (false-minima demotion) | ✗ | ✓ Two-gate: distance + frequency + consensus |
| Atom type: N.2 (sp2 imine) | → N.am (donor, wrong sign) | → N.ar (acceptor, correct) |
| Atom type: N.3 (amine) | → N.3/type-8 (dead matrix row) | → N.am/type-11 (live) |
| Atom type: C.1 (sp carbon) | → C.1/type-1 (sparse) | → C.2/type-2 (better sampled) |
| Atom type: I (iodine) | → type-26 (3 live entries) | → BR/type-25 (full halogen row) |
| WAL repulsion | Unbounded (SIGSEGV on extreme clashes) | ✓ Capped at 50 CF units per contact |
| Vibrational diagnostic (tENCoM) | ✗ | ✓ torsional elastic-network model scale |
| Astex-85 success rate | See claim receipt | See fixed-denominator claim receipt |

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

Every key below is read by `LIB/config_parser.cpp`; unknown keys are silently ignored, so
mistyped names fail open rather than erroring.

```json
{
  "ga": {
    "num_chromosomes": 1000,
    "num_generations": 2000,
    "seed": 0
  },
  "flexibility": {
    "ligand_torsions": true,
    "intramolecular": true
  },
  "scoring": {
    "pb_clash_weight": 1.0,
    "hbond_weight": -2.5
  },
  "thermodynamics": {
    "temperature": 300,
    "clustering_algorithm": "CF"
  },
  "output": {
    "max_results": 10
  }
}
```

`thermodynamics.temperature` is an integer score-scale parameter, not kelvin. Restart count and
ring-pucker sampling are not JSON config keys — they are set through `FLEXAIDDS_RESTARTS` and
`FLEXAIDDS_RING_FLEX` (see [Protocol Flags](#protocol-flags)).

### Interpret Output

FlexAID∆S writes ranked PDB files and a REMARK-annotated summary for each binding mode. Every numeric field below is in arbitrary CF units over optimizer samples — the emitted schema says so explicitly, and none of it is a physical ΔG, affinity, `Kd`, or `Ki`:

```
REMARK thermo_schema_version = 2
REMARK thermo_claim_validity = proxy_only
REMARK thermo_energy_domain = cf_arbitrary_units
REMARK thermo_ensemble_measure = optimizer_samples
REMARK thermo_reference_state = bound_only
REMARK proxy_free_energy = <float>   # ensemble transform over the mode's CF samples
REMARK free_energy = <float>         # deprecated alias of proxy_free_energy
REMARK soft_beta_G = <float>         # ranking objective G̃ (see below)
REMARK enthalpy = <float>            # weighted mean CF of the mode
REMARK entropy = <float>             # −Σ p ln p, per unit of the score-scale parameter
REMARK heat_capacity = <float>
REMARK temperature = <float>         # score-scale parameter, not a measured temperature
REMARK binding_mode = <int>
REMARK pose_rank = <int>
```

**What actually ranks binding modes.** `BindingMode::compute_energy()` (`LIB/BindingMode.cpp`)
is the ranking objective, and it is CF-bound: in the default classic path it returns the
mode-local soft-β free energy `G̃ = H̃ − T·S̃` computed over the member poses' own `CF` values
(optionally re-weighted by `CF.pb_clash` at election), plus a pose-independent vibrational
diagnostic and an optional NATURaL constant. `LIB/cluster.cpp` sorts clusters by the same
soft-β quantity (`ACF`, over `app_evalue`), so emission order and the S1 election agree with it.
Because both corrections are constants for a given receptor, they cancel in every difference
and never change the order. There is no `G_bind`, no `dG_eff`, and no `binding_regime` REMARK:
those legacy names are not emitted by the engine.

### Python Analysis

```bash
pip install -e ./python
```

```python
import flexaidds as fd

# Load results from a completed docking run
docking = fd.load_results("path/to/results/")

# Inspect top binding modes (proxy_free_energy is in CF units, not kcal/mol)
for mode in docking.binding_modes[:3]:
    print(f"Rank {mode.rank}: proxy F = {mode.proxy_free_energy}  soft_beta_G = {mode.soft_beta_G}")

# Compute ensemble statistics directly on a score array.
# Keyword is temperature_K; samples are added with add_sample/add_samples.
from flexaidds import StatMechEngine
engine = StatMechEngine(temperature_K=298.0)
engine.add_samples(pose_scores)
th = engine.compute()
print(f"F = {th.free_energy:.3f}  S = {th.entropy:.4f}  Cv = {th.heat_capacity:.4f}")
# Interpretation gate: proxy_only unless calibrated-energy provenance is attached.
print(th.claim_validity.value)
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

The **Astex Diverse Set** (85 diverse protein–ligand co-crystal complexes, covering 55 therapeutic
targets and 12 protein families) is the standard benchmark for evaluating flexible docking accuracy.

**Status: unverified / pending receipt. This repository publishes no validated Astex-85 success
rate for FlexAID∆S.** The previously advertised headline was withdrawn because the repository's own
audit notes record two defects that make it unusable as evidence:

- Runs were native-pose seeded, so a large share of "successes" are seed echoes — poses that
  reproduce the seeded input rather than a blind prediction.
- The harness `success` column means "the docking job ran", not "RMSD below the acceptance
  threshold", so the printed rate overstates accuracy.

A replacement number is deliberately not derived or estimated here. It becomes publishable only
after a blind (unseeded) campaign is deposited with a provenance receipt — binary digest, raw
ensemble digest, git object ID, runner command — and registered in
[`docs/entropy-help/audits/audits.json`](docs/entropy-help/audits/audits.json), at which point
`scripts/validate_thermo_claims.py` will allow the claim to appear on this page.

See [docs/BENCHMARK.md](docs/BENCHMARK.md) for dataset provenance and the exact commands.

---

## Protocol Flags

FlexAID∆S exposes its scientific innovations as environment-variable flags so each mechanism can be enabled, disabled, or tuned independently for ablation studies and benchmarking.

| Flag | Default | Effect |
|:-----|:-------:|:-------|
| `FLEXAIDDS_THERMO_SCORE` | OFF | Computes and prints the legacy `dH > 0 ∧ dS < 0` sentinel. Diagnostic only: the exact-CF rescore/clustering boundary does not read it, so it does not change the elected pose or the ranking criterion |
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
| `FLEXAIDDS_USE_SHANNON` | OFF | Presence-gated GA monitor that pools ligand ANM mode eigenvalues across cluster representatives and reports an ω-space Shannon diagnostic. The eigenvalues never enter CF, fitness, or election |
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
  │   Mode election (BindingMode/cluster)  │  │    Pose Validation        │
  │  ALL QUANTITIES IN CF UNITS            │  │  RMSD acceptance: S1     │
  │  H̃  = weighted mean CF of the mode    │  │  PoseBusters bust: S2    │
  │  S̃  = −Σ p_i ln p_i over mode members │  │  BCR (sampling ceiling)  │
  │  G̃  = H̃ − T·S̃   ← RANKING OBJECTIVE  │  └──────────────────────────┘
  │  T is a score-scale parameter, not K   │
  │  (+ pose-independent vib / NATURaL     │
  │   constants: cancel in every diff)     │
  └────────────────────────────────────────┘
                               │  ranked binding modes
                               ▼
  ┌────────────────────────────────────────┐  ┌──────────────────────────┐
  │  StatMechEngine (statmech.cpp)         │  │    tENCoM (tENCoM/)       │
  │  DIAGNOSTIC ONLY on this path          │  │  DIAGNOSTIC ONLY          │
  │  Canonical-form Z, F, Cv over CF       │  │  Torsional ENM model      │
  │  samples; emitted as proxy_only        │  │  scale; no kcal/mol       │
  │  WHAM / tempering / TI available for   │  │  calibration, not summed  │
  │  calibrated inputs supplied elsewhere  │  │  into a binding cycle     │
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
