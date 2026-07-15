<div align="center">

# FlexAID∆S

**Entropy-aware molecular docking engine**  
Genetic algorithm search · Voronoi contact-function scoring · statistical-mechanics ensemble analysis

[![CI](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml/badge.svg)](https://github.com/LeBonhommePharma/FlexAIDdS/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![C++26](https://img.shields.io/badge/C%2B%2B-26-blue.svg)](https://en.cppreference.com/w/cpp/26)
[![Python](https://img.shields.io/badge/python-%E2%89%A5%203.9-3776AB.svg)](https://www.python.org/)
[![Version](https://img.shields.io/badge/version-2.0.3-brightgreen.svg)](VERSION.md)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jcim.5b00078-blue)](https://doi.org/10.1021/acs.jcim.5b00078)

**[Installation](docs/INSTALLATION.md)** ·
**[User guide](docs/USERGUIDE.md)** ·
**[Support matrix](docs/SUPPORT_MATRIX.md)** ·
**[Reproducibility](docs/REPRODUCIBILITY.md)** ·
**[Benchmarks](docs/BENCHMARKS.md)** ·
**[Changelog](VERSION.md)** ·
**[Website](https://lebonhommepharma.github.io/FlexAIDdS/)**

</div>

---

## Intended use

**FlexAID∆S** (FlexAID with ∆S) is a production-oriented docking and ensemble-analysis stack for structure-based design, computational chemistry, and industrial R&D workflows. It extends the published [FlexAID](https://doi.org/10.1021/acs.jcim.5b00078) flexible-docking lineage with:

| Layer | Role |
|:------|:-----|
| **Search** | Genetic algorithm (GA) exploration of ligand pose and conformation |
| **Scoring proxy** | Voronoi **contact function (CF)** for ranking during search |
| **Ensemble analysis** | Partition-function / Shannon / vibrational terms over the sampled ensemble |
| **Validation** | RMSD admission, optional PoseBusters (`bust`), tENCoM / Eigen diagnostics |

**Primary deliverables**

| Artifact | Description |
|:---------|:------------|
| `FlexAIDdS` | Release docking executable (LTO-oriented build) |
| `FlexAID` | Legacy-compatible docking executable |
| `tENCoM` | Torsional elastic-network vibrational-entropy tool |
| `flexaidds` (Python) | Results I/O, analysis API, CLI inspector |
| `benchmark_datasets` | Campaign runner for Astex / protocolized benchmarks |

License: **Apache-2.0** (academic and commercial use). See [License](#license--compliance).

---

## Scientific integrity (read before citing numbers)

Docking scores and ensemble free-energy estimates are **not** automatically experimental binding free energies.

| Term | Meaning in this codebase |
|:-----|:-------------------------|
| **CF / contact-function scoring proxy** | Geometry-based score used to drive and rank the GA search (Voronoi CF / Vcontacts). |
| **Ensemble-derived free energy estimate** | Helmholtz-style *F*, entropy *S*, heat capacity *C<sub>v</sub>* from the sampled ensemble (StatMech / BindingMode). Requires the thermodynamic path to be enabled and validated for the claim. |
| **Thermodynamic ledger** | Structured breakdown (*F*, *H*, −*TS*, *C<sub>v</sub>*, Boltzmann weights) — reporting layer, not a guarantee of wet-lab Δ*G*. |

**Benchmark admission (TIER-1 claim path)** — see [`benchmarks/protocols/admission_metrics_contract.md`](benchmarks/protocols/admission_metrics_contract.md):

| Metric | Definition | Role |
|:-------|:-----------|:-----|
| **S1** | Elected pose RMSD ≤ 2.0 Å (Hungarian / protocol RMSD) | Primary success |
| **S2** | S1 ∧ official PoseBusters pass (`bust`) | Modern secondary |
| **S3 / BCR** | Any-cluster RMSD ≤ 2.0 Å | Sampling ceiling — **diagnostic only** |

Native pose seeding is **forbidden** on claim runs (`seed_echo` / `native_pose_seeded` must be zero). DoF search budget modulates **population (chromosomes)**, not generations — base CLI is often `pop=1000`, `gen=6000`; effective pop is reported in `[EVAL-BUDGET]` logs (`FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1`). See [AGENTS.md](AGENTS.md) and the three-engine protocol.

Full thermodynamics reference: [`docs/thermodynamics.md`](docs/thermodynamics.md).

---

## Support boundary

This repository is a full research platform. **Not every module is a supported product surface.**

### Supported (production path)

- `FlexAIDdS` / `FlexAID` / `tENCoM` CLI binaries  
- Documented JSON / legacy config docking workflows  
- `flexaidds` Python package (pure-Python + optional `_core` acceleration)  
- GoogleTest + pytest suites gated by CI  
- Reproducibility and benchmark contracts under [`benchmarks/`](benchmarks/)  
- Apache-2.0 first-party code; permissive third-party only ([`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md))

### Experimental / non-contract

- Swift / HealthKit-style device bridges  
- TypeScript / PWA dashboards  
- Fleet / multi-tenant distributed orchestration  
- NATURaL co-translational workflows  
- Accelerator paths not listed in the [support matrix](docs/SUPPORT_MATRIX.md)

Authoritative lists:  
[`docs/VALIDATED_CAPABILITIES.md`](docs/VALIDATED_CAPABILITIES.md) ·
[`docs/EXPERIMENTAL_CAPABILITIES.md`](docs/EXPERIMENTAL_CAPABILITIES.md) ·
[`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md)

---

## Quick start

### Build (C++26)

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

```bash
# Flexible dock (defaults include ligand flexibility; temperature via config)
./build/FlexAIDdS receptor.pdb ligand.mol2
# or legacy two-file interface:
# ./build/FlexAID config.inp ga.inp output_prefix
```

### Python analysis package

```bash
pip install -e ./python
```

```python
import flexaidds as fd

run = fd.load_results("path/to/results")
print(fd.__version__)  # 2.0.3
```

### macOS Homebrew

```bash
brew tap lebonhommepharma/flexaidds https://github.com/LeBonhommePharma/FlexAIDdS
# Homebrew 6+ tap trust (formula-scoped; required when HOMEBREW_REQUIRE_TAP_TRUST is set):
brew trust --formula lebonhommepharma/flexaidds/flexaidds
brew install lebonhommepharma/flexaidds/flexaidds
# Metal (stable v2.0.3+ on main/tag — no feature branch):
#   brew install --build-from-source --with-metal lebonhommepharma/flexaidds/flexaidds
```

Native tools and the Python package are separate installs. Full platform notes (including recovery if a stale HEAD branch breaks reinstall): [`docs/INSTALLATION.md`](docs/INSTALLATION.md).

---

## Capabilities

### Docking engine

- Genetic algorithm search with configurable population, generations, restarts, and diversity controls  
- **Voronoi contact-function (CF)** scoring for shape complementarity  
- Dead-end elimination (DEE) pruning of ligand conformational space  
- Batch CF evaluation (`VoronoiCFBatch`) with OpenMP / SIMD paths  
- Clustering: centroid, FastOPTICS, density-peak  
- Multi-format ligands: MOL2, SDF/MOL, SMILES (where enabled); receptors PDB / CIF  
- Full ligand flexibility by default: torsions, ring conformers, R/S centers  
- Optional GIST / H-bond / metal-ion terms where configured  

### Ensemble thermodynamics

- Canonical ensemble partition function *Z*, Helmholtz *F*, Shannon configurational entropy *S*, *C<sub>v</sub>*  
- Binding-mode clustering and Boltzmann reweighting  
- Grand-canonical paths for competitive / concentration-aware analysis (research)  
- tENCoM vibrational corrections (relative unless calibrated — see [`docs/TENCOM_ENTROPY_CALIBRATION.md`](docs/TENCOM_ENTROPY_CALIBRATION.md))  
- Log-sum-exp numerical stability throughout  

### Pose validation

- In-tree **NativePoseQC** (`LIB/PoseBust`) — clean-room diagnostic suite (Apache-2.0)  
- Optional subprocess bridge to official **PoseBusters** `bust` (BSD, not vendored) for claim-ready **S2**  
- Standalone product: **[PoseBust](https://github.com/LeBonhommePharma/PoseBust)** (C++26, independent of this monorepo)  

### Hardware

Runtime preference order (where built): **CUDA → Metal → AVX-512 → AVX2 → OpenMP → scalar**.  
Platforms: Linux (GCC/Clang), macOS (Clang, Apple Silicon Metal), Windows (MSVC). See [support matrix](docs/SUPPORT_MATRIX.md).

---

## Architecture

```text
  Receptor / ligand I/O          Genetic algorithm           CF scoring
  (PDB, MOL2, SDF, …)    --->    (gaboom / restarts)   --->  (Voronoi / Vcontacts)
                                        |
                                        v
                              Pose ensemble + clustering
                              (BindingMode / election)
                                        |
                    +-------------------+-------------------+
                    |                                       |
                    v                                       v
           Thermodynamic ledger                    Pose validation
           (StatMech, Shannon, tENCoM)             (RMSD · PoseBust · bust)
```

**Ranking during search** uses the CF scoring proxy. Thermodynamic quantities are derived from the **sampled ensemble** after (or alongside) search, depending on configuration. Election policies for claim campaigns are protocol-defined (e.g. CF-only arms vs entropy-ranked arms).

---

## Build system

### Requirements

| Component | Minimum |
|:----------|:--------|
| Compiler | C++26 (GCC ≥ 14, Clang ≥ 18, Apple Clang ≥ 16 / Xcode 16, MSVC ≥ 19.40) |
| CMake | ≥ 3.28 |
| Recommended | Eigen3 |
| Optional | OpenMP, CUDA, Metal (macOS), pybind11, MPI |

### Common configurations

```bash
# Release + tests
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure

# Python bindings
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_PYTHON_BINDINGS=ON
cmake --build build --parallel

# Apple Silicon Metal
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFLEXAIDS_USE_METAL=ON
```

### Selected CMake options

| Option | Default | Purpose |
|:-------|:--------|:--------|
| `BUILD_TESTING` | OFF | GoogleTest suite |
| `BUILD_PYTHON_BINDINGS` | OFF | pybind11 `_core` |
| `FLEXAIDS_USE_OPENMP` | ON | Threading |
| `FLEXAIDS_USE_AVX2` | ON | SIMD |
| `FLEXAIDS_USE_AVX512` | OFF | SIMD (HPC) |
| `FLEXAIDS_USE_CUDA` | OFF | NVIDIA GPU |
| `FLEXAIDS_USE_METAL` | OFF | Apple GPU |
| `FLEXAIDS_USE_MPI` | OFF | Distributed domain decomposition |
| `ENABLE_TENCOM_TOOL` | ON | `tENCoM` binary |

After changing sources or `CMakeLists.txt`, always reconfigure and rebuild; do not assume linking still succeeds.

---

## Usage

### Command line

```bash
./build/FlexAIDdS receptor.pdb ligand.mol2
./build/FlexAIDdS receptor.pdb ligand.mol2 -c config.json
./build/FlexAIDdS receptor.pdb ligand.mol2 --rigid   # screening mode
```

JSON configuration (illustrative):

```json
{
  "thermodynamics": { "temperature": 298, "clustering_algorithm": "DP" },
  "ga": { "num_chromosomes": 1000, "num_generations": 6000 },
  "flexibility": { "ligand_torsions": true, "ring_conformers": true }
}
```

Legacy:

```bash
./build/FlexAID config.inp ga.inp output_prefix
```

### Python

```python
import flexaidds as fd

# High-level docking (requires bindings / engine as documented)
results = fd.dock(
    receptor="receptor.pdb",
    ligand="ligand.mol2",
    compute_entropy=True,
)

docking = fd.load_results("output_prefix")
for mode in docking.binding_modes:
    print(mode.rank, mode.free_energy, mode.entropy)
```

```python
from flexaidds import StatMechEngine

engine = StatMechEngine(temperature=298.0)
engine.add_energies(pose_energies)
thermo = engine.compute()
print(thermo.free_energy, thermo.entropy)
```

CLI inspector:

```bash
python -m flexaidds /path/to/results/
python -m flexaidds /path/to/results/ --top 5 --json
```

### PyMOL

Install via Plugin Manager → `pymol_plugin/`. Commands include `flexaids_load`, `flexaids_show_ensemble`, `flexaids_color_boltzmann`, `flexaids_thermo`. Requires `pip install -e python/`.

### tENCoM

```bash
tENCoM reference.pdb target.pdb [-T 298] [-o prefix]
```

Vibrational entropy is a **relative** heuristic unless a validated eigenvalue-to-frequency calibration is supplied.

---

## Benchmarking & campaigns

Protocolized campaigns (Astex Diverse, three-engine comparison, admission metrics) live under [`benchmarks/`](benchmarks/) and agent contracts in [`AGENTS.md`](AGENTS.md).

**Operational norms for production campaigns**

- Pin interaction matrix MD5 (e.g. `MC_st0r5.2_6.dat`) in every arm’s receipt  
- Cognate / defined-cleft redock: **no native pose seed** on claim paths  
- Report **S1** and **S2**; do not sell BCR/S3 as abstract success  
- **All campaign results → iCloud Drive** (`$FLEXAIDDS_RESULTS` under CloudDocs; ~2 TB quota) — see [`docs/ICLOUD_BENCHMARK_STORAGE.md`](docs/ICLOUD_BENCHMARK_STORAGE.md)  
- Stage Mach-O binaries on **local** disk only; never new claim OUT under `~/flexaidds_results`  
- On memory-constrained hosts: **one heavy GA process at a time**  
- Aggregate from on-disk CSV/JSON — never from chat memory  

Three-engine protocol: [`benchmarks/protocols/three_engine_entropy_comparison.md`](benchmarks/protocols/three_engine_entropy_comparison.md).  
Benchmark skill (agents): [`.agents/skills/flexaidds-benchmarking/SKILL.md`](.agents/skills/flexaidds-benchmarking/SKILL.md).

---

## Reproducibility

A number in documentation is **not** automatically repository-reproducible. A claim is reproducible when a bundle under [`benchmarks/`](benchmarks/) provides dataset provenance, pinned binaries/matrices, exact commands, and metric definitions.

Policy: [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md) · [`benchmarks/README.md`](benchmarks/README.md)

---

## Quality assurance

### C++

```bash
cmake -S . -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

### Python

```bash
pip install -e ./python
pytest tests/ -q
```

Tests marked `@requires_core` skip when the C++ extension is absent.

### CI

GitHub Actions: multi-compiler C++ builds, pure-Python tests, Python binding smoke, license hygiene. See [`.github/workflows/`](.github/workflows/) and [`docs/SUPPORT_MATRIX.md`](docs/SUPPORT_MATRIX.md).

---

## Documentation map

| Document | Audience |
|:---------|:---------|
| [Installation](docs/INSTALLATION.md) | Deployers, IT, scientists |
| [User guide](docs/USERGUIDE.md) | End users, API consumers |
| [Support matrix](docs/SUPPORT_MATRIX.md) | Platform owners |
| [Thermodynamics](docs/thermodynamics.md) | Methods / modelers |
| [Benchmarks](docs/BENCHMARKS.md) | Validation leads |
| [Admission metrics](benchmarks/protocols/admission_metrics_contract.md) | Claim authors |
| [Clean-room policy](docs/licensing/clean-room-policy.md) | Legal / contributors |
| [VERSION.md](VERSION.md) | Release managers |
| [AGENTS.md](AGENTS.md) | Automation & AI coding agents |

Website: [lebonhommepharma.github.io/FlexAIDdS](https://lebonhommepharma.github.io/FlexAIDdS/)

---

## Related software

| Project | Relation |
|:--------|:---------|
| [PoseBust](https://github.com/LeBonhommePharma/PoseBust) | Standalone C++26 pose validation (NativePoseQC + optional `bust`) |
| [NRGsuite](https://doi.org/10.1093/bioinformatics/btv458) | PyMOL docking UI lineage |
| [Shannon](https://github.com/LeBonhommePharma/Shannon) | Shared Shannon-entropy methodology outside docking (separate product) |

---

## Publications

If you use FlexAID or FlexAID∆S, please cite:

> Gaudreault F & Najmanovich RJ (2015). FlexAID: Revisiting Docking on Non-Native-Complex Structures.  
> *J. Chem. Inf. Model.* 55(7):1323–1336.  
> [DOI:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)

Related:

- Gaudreault F, Morency LP & Najmanovich RJ (2015). NRGsuite. *Bioinformatics* 31(23):3856–3858. [DOI:10.1093/bioinformatics/btv458](https://doi.org/10.1093/bioinformatics/btv458)  
- Frappier V et al. (2015). ENCoM. *Proteins* 83(11):2073–2082. [DOI:10.1002/prot.24922](https://doi.org/10.1002/prot.24922)  
- Morency LP & Najmanovich RJ (2026). FlexAID∆S methods manuscript — *in preparation*

---

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) and [AGENTS.md](AGENTS.md) before opening PRs.

| Rule | Requirement |
|:-----|:------------|
| License of contributions | Apache-2.0 |
| Allowed dependencies | Apache-2.0, BSD, MIT, MPL-2.0, PSF |
| Forbidden | GPL / AGPL (including as implementation inspiration) |
| Verification | Build and tests green before merge; no silent ranking changes without tests + flag |

AI agent skill (optional): [`.grok/skills/flexaidds/SKILL.md`](.grok/skills/flexaidds/SKILL.md)

---

## License & compliance

**Apache License 2.0** — free for academic and commercial use.

| File | Purpose |
|:-----|:--------|
| [LICENSE](LICENSE) | Full Apache-2.0 text |
| [NOTICE](NOTICE) | Copyright and attribution |
| [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md) | Dependency and optional-tool licenses |
| [docs/licensing/LICENSE_MATRIX.md](docs/licensing/LICENSE_MATRIX.md) | Compatibility matrix |
| [docs/licensing/clean-room-policy.md](docs/licensing/clean-room-policy.md) | GPL isolation policy |

Copyright © 2026 Le Bonhomme Pharma · Louis-Philippe Morency

---

## Maintainers

**Le Bonhomme Pharma** · [GitHub](https://github.com/LeBonhommePharma)  
**Louis-Philippe Morency** — project lead  

Issues: [github.com/LeBonhommePharma/FlexAIDdS/issues](https://github.com/LeBonhommePharma/FlexAIDdS/issues)

---

*FlexAID∆S is provided as-is under Apache-2.0. Validate all numerical claims against pinned benchmarks and experimental data before regulatory or clinical decision use.*
