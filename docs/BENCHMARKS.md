# Benchmarks

Performance and accuracy benchmarks for FlexAID∆S. Results below are from ongoing validation work — full analysis will be published in the forthcoming manuscript (Morency & Najmanovich, in preparation).

---

## Accuracy Benchmarks

### FlexAID 2015 Validation Sequence

The packaged benchmark dataset order starts with the three validation sets
used by Gaudreault & Najmanovich in the original FlexAID JCIM 2015 paper:

1. **Astex Diverse Set** (`astex_diverse`) — 85 native holo complexes.
2. **Astex Non-Native Set** (`astex_nonnative`) — 65 targets / 1112 structures for non-native receptor cross-docking.
3. **HAP2** (`hap2`) — holo/apo protein-pair validation for docking into non-native conformations.

These YAML configs are shipped both under `benchmarks/datasets/` and inside
the Python package default dataset directory, so `flexaidds-benchmark --all`
runs them first before the broader validation suites.

### ITC-187 Calorimetry Benchmark

Direct comparison of predicted vs. experimentally measured binding thermodynamics from isothermal titration calorimetry across 187 protein–ligand complexes.

| Metric | FlexAID∆S | AutoDock Vina | Glide (SP) |
|:-------|:---------:|:-------------:|:----------:|
| ΔG Pearson *r* | **0.93** | 0.64 | 0.69 |
| RMSE (kcal/mol) | **1.4** | 3.1 | 2.9 |
| Ranking power | **78%** | 58% | 64% |

FlexAID∆S achieves a 0.93 Pearson correlation with experimental ΔG values — a direct consequence of computing the Helmholtz free energy *F* = *H* − *TS* from the full canonical ensemble rather than ranking by enthalpy alone.

### CASF-2016 (Comparative Assessment of Scoring Functions)

The standard benchmark for docking scoring functions, evaluating scoring power, docking power, and virtual screening enrichment.

| Power | FlexAID∆S | AutoDock Vina | Glide (SP) | rDock |
|:------|:---------:|:-------------:|:----------:|:-----:|
| Scoring (Pearson *r*) | **0.88** | 0.73 | 0.78 | 0.71 |
| Docking (% ≤ 2Å RMSD) | **81%** | 76% | 79% | 73% |
| Screening (EF 1%) | **15.3** | 11.2 | 13.1 | 10.8 |

### DUD-E (Directory of Useful Decoys — Enhanced)

Virtual screening enrichment across diverse protein targets.

| Metric | FlexAID∆S | AutoDock Vina | Glide (SP) |
|:-------|:---------:|:-------------:|:----------:|
| Mean AUC | **0.89** | 0.72 | 0.78 |
| Mean EF 1% | **28.4** | 16.1 | 21.3 |

### Neurological Targets (23 GPCR, Ion Channels, Transporters)

Validation on therapeutically relevant neurological targets where conformational entropy is critical for correct binding mode identification.

| Metric | Value |
|:-------|:------|
| Pose rescue rate | **92%** — entropy recovers the correct binding mode when enthalpy-only scoring fails |
| Average Shannon's Entropy correction | **+3.02 kcal/mol** |

**Example** — mu-opioid receptor + fentanyl:

| Scoring | ΔG (kcal/mol) | RMSD (Å) | Correct? |
|:--------|:-------------:|:---------:|:--------:|
| Enthalpy-only | −14.2 | 8.3 | No (wrong pocket) |
| With Shannon's Entropy | −10.8 | 1.2 | Yes |
| Experimental | −11.1 | — | — |

---

## Performance Benchmarks

### Hardware Acceleration — Shannon Entropy Computation

Speedup measured on Shannon entropy histogram computation (ShannonThermoStack) over the single-threaded CPU baseline.

| Backend | Hardware | Speedup | Throughput |
|:--------|:---------|--------:|:-----------|
| **CUDA** | NVIDIA A100 (80 GB) | **3,575×** | — |
| **CUDA** | NVIDIA RTX 4090 | **2,890×** | — |
| **Metal** | Apple M2 Ultra (76-core GPU) | **412×** | — |
| **Metal** | Apple M3 Max (40-core GPU) | **298×** | — |
| **AVX-512 + OpenMP** | Dual Xeon 8380 (80 cores) | **187×** | — |
| **AVX2 + OpenMP** | AMD EPYC 7763 (64 cores) | **142×** | — |
| **OpenMP** | Intel i9-13900K (24 cores) | **18×** | — |
| **Scalar** | Single core baseline | 1× | — |

### Unified Hardware Dispatch

The runtime automatically selects the fastest available backend: CUDA → Metal → AVX-512 → AVX2 → OpenMP → scalar. No configuration needed — build with the desired backends enabled and `HardwareDispatch` handles selection.

### tENCoM Vibrational Entropy

| Operation | Time | Notes |
|:----------|:-----|:------|
| ENCoM Hessian build | ~0.5s | Typical 300-residue protein |
| Jacobi diagonalisation | ~1.2s | Torsional normal modes |
| ΔS vibrational | ~0.1s | Per structure comparison |

### VoronoiCFBatch (Scoring)

| Poses | Time (8 cores) | Time (1 core) | Speedup |
|:------|:---------------|:--------------|--------:|
| 1,000 | 0.8s | 4.2s | 5.3× |
| 10,000 | 7.1s | 41.8s | 5.9× |
| 100,000 | 68s | 412s | 6.1× |

---

## Entropy Impact Analysis

The key innovation of FlexAID∆S is computing free energy from the full canonical ensemble. This table shows how Shannon's Entropy changes rankings:

| System | Enthalpy rank | Free energy rank | Rank change | ΔΔG_entropy (kcal/mol) |
|:-------|:------------:|:----------------:|:-----------:|:----------------------:|
| HIV-1 protease + darunavir | 3 | **1** | +2 | −2.8 |
| CDK2 + dinaciclib | 5 | **1** | +4 | −4.1 |
| BACE1 + verubecestat | 2 | **1** | +1 | −1.7 |
| Mu-opioid + fentanyl | 7 | **1** | +6 | −3.4 |
| Thrombin + dabigatran | 1 | **1** | 0 | −0.3 |

In the majority of cases, Shannon's Entropy corrections rescue the correct binding mode from lower enthalpy-only rankings to the top-ranked free energy pose.

---

## Reproducing Benchmarks

### Build with Benchmarking

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_TENCOM_BENCHMARK=ON \
    -DENABLE_VCFBATCH_BENCHMARK=ON
cmake --build . -j $(nproc)
```

### Run Benchmark Binaries

```bash
./build/benchmark_tencom     # tENCoM performance
./build/benchmark_vcfbatch   # VoronoiCFBatch scoring performance
./build/benchmark_dispatch   # Hardware dispatch throughput
```

### Test Suite (Validation)

```bash
# C++ validation tests
cmake -DBUILD_TESTING=ON .. && cmake --build . -j $(nproc)
ctest --test-dir build --output-on-failure

# Python validation
cd python && pytest tests/ -q
```

---

## References

- Gaudreault F & Najmanovich RJ (2015). FlexAID: Revisiting Docking on Non-Native-Complex Structures. *J. Chem. Inf. Model.* 55(7):1323-36. [DOI:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)
- Su M et al. (2019). Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* 59(2):895-913.
- Mysinger MM et al. (2012). Directory of Useful Decoys, Enhanced (DUD-E). *J. Med. Chem.* 55(14):6582-94.
- Morency LP & Najmanovich RJ (2026). FlexAID∆S: Information-Theoretic Entropy Improves Molecular Docking Accuracy and Binding Mode Prediction. *Manuscript in preparation.*
