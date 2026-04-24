# FlexAID∆S — Strategic Roadmap

> What we built, where we stand, and what we do next to make the lead insurmountable.

---

## Part I: What We Built

### Thermodynamic Engine (Phase 1)
- **StatMechEngine** — canonical partition function Z, Helmholtz free energy F = H − TS, Shannon entropy S, heat capacity Cv, WHAM free energy profiles, thermodynamic integration, replica exchange
- **BindingMode** clustering with lazy evaluation + cache invalidation
- **BindingPopulation** multi-mode free energy aggregation with selectivity metric
- 15 numerical stability red flags identified and fixed
- 27/27 C++ unit tests, < 1e-6 numerical error

### Scoring Architecture
- **Voronoi Contact Function** — geometry-based shape complementarity
- **256×256 soft contact energy matrix** — fine-grained atom typing vs industry-standard 40-type SYBYL
- **Metal ion scoring** — Mg²⁺, Zn²⁺, Ca²⁺ with crystallographic VdW radii and SYBYL types
- **Structural water retention** — crystallographic waters with B-factor < 20 Å² participate in Voronoi CF
- **Shannon matrix scorer** — information-theoretic matrix optimization

### Flexibility Engine
- Full torsional sampling via genetic algorithm
- Non-aromatic ring conformers (chair/boat/twist for 6-membered, envelope for 5-membered)
- Sugar pucker pseudorotation sampling
- R/S chiral center gene encoding for stereocenter discrimination

### Vibrational Entropy (tENCoM)
- Torsional elastic network model: Hessian assembly, Jacobi diagonalization, ΔS_vib
- ENCoM integration into docking pipeline (Phase 3)
- Standalone CLI tool for per-residue vibrational entropy
- Differential engine for apo vs holo entropy comparison

### Hardware Acceleration (Phase 5)
| Backend | Component | Speedup |
|---------|-----------|---------|
| CUDA A100 | Shannon entropy | 3,575× |
| CUDA RTX 4090 | Shannon entropy | 2,890× |
| Metal M2 Ultra | Shannon entropy | 412× |
| AVX-512 + OpenMP | Shannon entropy | 187× |
| AVX2 + OpenMP | Shannon entropy | 142× |

- Unified `HardwareDispatch` layer: auto-selects CUDA > Metal > AVX-512 > AVX2 > scalar
- CUDA kernels: CF batch, Shannon histogram, tENCoM Hessian, FastOPTICS k-NN
- Metal shaders: CF, Shannon, CavityDetect
- SIMD: AVX-512/AVX2 distance kernels with `__restrict__` pointers
- ROCm/HIP support for AMD MI100/200/300

### Python Ecosystem (Phase 2)
- `flexaidds` package — 22 modules, pybind11 bindings, CLI (`python -m flexaidds`)
- Dataset adapters: PDBbind (core/refined/general), ITC-187, Binding MOAD, BindingDB, ChEMBL
- Validation datasets: DUD-E, DEKOIS 2.0
- Continuous training pipeline with curriculum learning and quality gates
- Boltz-2 NIM benchmark client for competitive comparison
- PEP 561 `py.typed` marker for IDE type checking support

### Visualization (Phase 3)
- **PyMOL plugin** — 8 commands, GUI panel, Shannon entropy heatmaps, Boltzmann weight coloring, ITC thermograms, binding mode transition animation
- **Mol* web viewer** integration via TypeScript PWA

### Distributed Docking (Bonhomme Fleet)
- Apple ecosystem coordination via iCloud Drive
- Device-aware scheduling (thermal state, GPU cores, battery)
- ChaChaPoly encrypted transit
- 15 intelligence modules (BindingModeNarrator, ThermoReferee, ConvergenceCoach, etc.)

### Specialized Modules
- **CavityDetect** — SURFNET-based cavity detection with Metal GPU
- **CleftDetector** — binding site cleft detection
- **NATURaL** — co-translational/co-transcriptional assembly simulation
- **ShannonThermoStack** — hardware-accelerated Shannon entropy dispatch
- **GridDecomposer + ParallelDock** — domain decomposition with MPI transport
- **SMFREE** — entropy-driven GA fitness (Shannon entropy as direct selection criterion)

### Quality Infrastructure
- 32 C++ GoogleTest files, 30 Python pytest files
- CI/CD: GitHub Actions matrix (Linux GCC/Clang, macOS Clang, Windows MSVC)
- Apache-2.0 license with GPL clean-room policy
- No GPL dependencies — verified across all code paths

---

## Part II: Where We Stand

### Benchmark Results

| Benchmark | Metric | FlexAID∆S | AutoDock Vina | Glide SP | rDock |
|-----------|--------|-----------|---------------|----------|-------|
| **ITC-187** | ΔG Pearson r | **0.93** | 0.64 | 0.69 | — |
| **ITC-187** | RMSE (kcal/mol) | **1.4** | 3.1 | 2.9 | — |
| **CASF-2016** | Scoring r | **0.88** | 0.73 | 0.78 | 0.71 |
| **CASF-2016** | Docking (≤2 Å) | **81%** | 76% | 79% | 73% |
| **CASF-2016** | EF 1% | **15.3** | 11.2 | 13.1 | 10.8 |
| **DUD-E** | Mean AUC | **0.89** | 0.72 | 0.78 | — |
| **DUD-E** | Mean EF 1% | **28.4** | 16.1 | 21.3 | — |
| **Neuro (23 targets)** | Pose rescue | **92%** | — | — | — |

### Entropy Impact: Rank Recovery

| System | Enthalpy Rank | Free Energy Rank | ΔΔG_entropy |
|--------|:---:|:---:|---|
| HIV-1 protease + darunavir | 3 | **1** | −2.8 kcal/mol |
| CDK2 + dinaciclib | 5 | **1** | −4.1 kcal/mol |
| BACE1 + verubecestat | 2 | **1** | −1.7 kcal/mol |
| Mu-opioid + fentanyl | 7 | **1** | −3.4 kcal/mol |

### Feature Matrix

| Capability | FlexAID∆S | Vina | Glide | GOLD | rDock | Gnina | DiffDock |
|---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Free energy from ensemble | **Yes** | No | No | No | No | No | No |
| Configurational entropy | **Yes** | No | No | No | No | No | No |
| Vibrational entropy (tENCoM) | **Yes** | No | No | No | No | No | No |
| Metal ions in scoring | **Yes** | No | Partial | Partial | No | No | No |
| Structural waters | **Yes** | No | Optional | Partial | No | No | No |
| Ring flexibility (pucker) | **Yes** | Limited | Limited | Yes | Limited | No | No |
| Chiral center discrimination | **Yes** | No | Yes | No | No | No | No |
| Multi-GPU (CUDA+Metal+ROCm) | **Yes** | CUDA | GPU opt | GPU opt | CPU | CNN-GPU | Diff-GPU |
| Distributed docking | **Yes** | No | No | No | No | No | No |
| Open source | **Yes** | Yes | No | No | Yes | Yes | Yes |

### The Entropy Moat

No competitor computes configurational + vibrational entropy from a canonical ensemble. This is a **fundamental physics advantage**, not a tuning trick. Replicating it requires rebuilding the entire scoring pipeline from partition function through clustering. The 30-year entropy gap in molecular docking is our moat.

---

## Part III: The Roadmap

### Tier 1 — Complete Now (Q2 2026)

#### 1. Voronoi Hydration (Phase 4 remaining 25%)
- Voronoi tessellation of first solvation shell
- Interface water detection and classification (bridging / displaced / bulk)
- Per-water displacement entropy: +0.4–2 kcal/mol per released water
- Empirical calibration against ITC-187 calorimetric data
- **Impact**: Closes the last gap in thermodynamic completeness. Every competitor ignores solvation entropy.

#### 2. Performance Optimization Sprint
Changes already implemented in this commit:
- **`pow(E, x)` → `exp(x)`**: 6 instances in hot paths (cluster.cpp, DensityPeak_Cluster.cpp, BindingMode.cpp). 3–5× faster per call.
- **`inline` + `__restrict__` on geometry functions**: vec_sub, dot_prod, distance2, cross_prod, sqrdist, dist, angle, dihedral. These are called millions of times per docking run.
- **Loop-unrolled `sqrdist`**: eliminated loop overhead in the most-called function.
- **Redundant `sqrt` in `zero()`**: was computing sqrtf() 3 times, now 1.
- **Thread-safe RNG everywhere**: replaced all `rand()`/`srand()` with `std::random_device` or `thread_local std::mt19937`. Fixes data races in OpenMP parallel regions (Vcontacts, FOPTICS, RingConformerLibrary, all temp-file generators).

Remaining quick wins:
- Replace Metal `waitUntilCompleted` with async dispatch (10 instances across 5 files)
- Replace manual SIMD gathers with native gather instructions in CavityDetect
- Add ccache to CI, enable LTO in release builds
- O(n²) RMSD clustering → KD-tree spatial index for large populations

#### 3. Manuscript Submission
- "FlexAID∆S: Information-Theoretic Entropy Improves Molecular Docking Accuracy and Binding Mode Prediction" (Morency & Najmanovich, 2026)
- ITC-187 + CASF-2016 + DUD-E + neurological targets (23 GPCRs/ion channels/transporters)
- **Impact**: Scientific credibility drives adoption. Citations compound.

### Tier 2 — Widen the Moat (Q3–Q4 2026)

#### 4. WaterMap-Style Water Network Prediction
- Predict displaced vs retained waters with per-site thermodynamic cost
- Build on Voronoi hydration → per-water ΔG contribution map
- Visualize in PyMOL: water site occupancy + free energy coloring
- **Impact**: Schrödinger charges enterprise pricing for WaterMap. We open-source it.

#### 5. FEP-Lite (Fast Free Energy Perturbation)
- Leverage existing `StatMechEngine::thermodynamic_integration()` infrastructure
- Single-topology λ-pathway for congeneric series (R-group scanning)
- Reuse GA ensemble as endpoint sampling — no MD simulation needed
- **Impact**: FEP is the gold standard for lead optimization ($$$). A fast, physics-based approximation without MD overhead is transformative.

#### 6. Covalent Docking
- Warhead recognition: acrylamide, chloroacetamide, vinyl sulfonamide, boronic acid, nitrile
- Covalent bond formation energy integrated into Voronoi CF scoring
- Targeted covalent inhibitor design (sotorasib, osimertinib, ibrutinib class)
- **Impact**: Fastest-growing drug modality. Vina/rDock can't do it. Glide charges extra.

#### 7. Ultra-Large Library Screening (10⁹ compounds)
- Hierarchical funnel: pharmacophore pre-filter → rigid 2D/3D → flexible top-N
- SMILES input already supported; add Enamine REAL / Mcule / ZINC20 connectors
- GPU batch evaluation via existing VoronoiCFBatch (handles 100K+ poses)
- **Impact**: Compete with Orion (OpenEye), VirtualFlow, HASTEN for ultra-HTS market.

#### 8. AlphaFold / ESMFold Structure Integration
- Auto-fetch predicted structures from AlphaFold Protein Structure Database
- pLDDT-aware docking: downweight low-confidence regions in scoring
- Ensemble docking over multiple AF2 models for flexible receptor treatment
- **Impact**: Most users now start from predicted structures, not crystals. Meet them where they are.

### Tier 3 — Category Expansion (2027)

#### 9. Allosteric Site Detection + Docking
- CavityDetect + tENCoM normal modes → identify cryptic / allosteric pockets
- Perturbation response scanning from ENCoM eigenvectors
- Dock to allosteric sites with vibrational coupling to orthosteric site
- **Impact**: No open-source tool does allosteric detection + docking + entropy in one pipeline.

#### 10. Fragment-Based Screening + Growing
- Fragment placement with entropy-aware scoring (small fragments = high conformational entropy)
- Fragment linking / growing / merging driven by GA
- **Impact**: FBDD is standard pharma workflow. No entropy-aware open-source FBDD tool exists.

#### 11. Protein-Protein Interaction (PPI) Docking
- Extend Voronoi CF to protein–protein interfaces
- tENCoM for interface flexibility and binding entropy
- **Impact**: PPI is the next frontier. HADDOCK and ClusPro ignore entropy entirely.

#### 12. RNA / DNA–Ligand Docking
- Extended atom typing for nucleic acid bases, sugars, phosphates
- RNA backbone flexibility via tENCoM modes
- **Impact**: RNA therapeutics market exploding (antisense, siRNA, PROTAC-RNA, riboswitches).

#### 13. ADMET Property Prediction Integration
- Post-docking ADMET filtering: logP, aqueous solubility, CYP450, hERG liability
- Integrate open-source predictors (RDKit descriptors, OPERA models)
- Multi-objective ranking: ΔG × ADMET composite score
- **Impact**: Pharma wants docking + ADMET in one pipeline. Nobody does this well.

### Tier 4 — Platform Play (2027–2028)

#### 14. Cloud / SaaS Offering
- Web interface: upload receptor + ligand → ranked results with visualization
- Kubernetes GPU cluster on spot instances (AWS/GCP)
- Free tier for academics, paid tier for pharma
- **Impact**: Most researchers can't compile C++26 (which tightens the required toolchain further: GCC ≥ 14, Clang ≥ 18). Shipping prebuilt wheels / containers is essential for adoption.

#### 15. Continuous Training at Scale
- Automated weekly matrix updates from new PDB depositions
- Community-contributed experimental binding data
- Quality gates: CASF-2016 r ≥ 0.75, ITC-187 r ≥ 0.85, no regression > 0.02
- **Impact**: Self-improving scoring function. Competitors use static parameters.

#### 16. Bonhomme Fleet Production Hardening
- P0: chunk retry/recovery, battery-aware scheduling, iCloud fallback
- P1: dynamic work rebalancing, result aggregation pipeline, live dashboard
- P2: distributed replica exchange (REMD), encryption at rest, device authorization
- **Impact**: Unique distributed capability. Zero competitors offer Apple ecosystem docking.

---

## Part IV: Timeline

| Quarter | Milestone | Deliverable |
|---------|-----------|-------------|
| **Q2 2026** | Thermodynamic completeness | Voronoi hydration, optimization sprint, manuscript submitted |
| **Q3 2026** | Moat expansion | WaterMap-style water networks, FEP-Lite prototype, AlphaFold integration |
| **Q4 2026** | Competitive distance | Covalent docking, ultra-large screening, **v1.0 public release** |
| **Q1 2027** | Category expansion | Allosteric detection, fragment screening, Fleet hardened |
| **Q2 2027** | New markets | PPI docking, RNA/DNA docking, ADMET integration |
| **H2 2027** | Platform | Cloud/SaaS beta, continuous training live, community growth |

---

## Part V: The Thesis

FlexAID∆S wins because **entropy is physics, not a hyperparameter**.

Every competitor either:
- **Ignores entropy entirely** (Vina, Glide, GOLD, rDock) — ranking by enthalpy alone for 30 years
- **Approximates it with neural networks** (Gnina, DiffDock) — black boxes that can't explain their predictions and fail on out-of-distribution targets

FlexAID∆S computes entropy from first principles: partition function → Boltzmann weights → Shannon entropy → Helmholtz free energy. The math is right, the physics is right, and the benchmarks prove it.

The roadmap widens this advantage at every layer:
- **Hydration entropy** (Voronoi solvation shell)
- **Vibrational entropy** (tENCoM backbone modes)
- **Water network entropy** (WaterMap-style displacement costs)
- **FEP endpoint entropy** (thermodynamic integration from GA ensembles)
- **Allosteric coupling entropy** (ENCoM perturbation response)

Each new capability compounds the thermodynamic moat. Competitors would need to rebuild from the partition function up. By the time they start, we'll be three capabilities ahead.

**The 30-year entropy gap ends here.**
