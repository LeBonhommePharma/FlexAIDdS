# Benchmarks

This page tracks benchmark intent, benchmark entry points, and claim maturity for
FlexAIDdS. It does not present accuracy, enrichment, or hardware-speedup numbers
as repository-reproducible unless the corresponding benchmark bundle contains
dataset provenance, immutable inputs, exact commands, seeds, metric scripts, and
expected outputs.

The claim gate follows [REPRODUCIBILITY.md](REPRODUCIBILITY.md). Manuscript or
exploratory values may exist outside the repository, but they must not be mixed
with replayable repository evidence.

---

## Current validation status

The current repository artifacts do not yet support final public benchmark
tables for ITC-187, CASF-2016, DUD-E, neurological target pose rescue, hardware
dispatch speedups, or VoronoiCFBatch throughput.

Known local audit status:

| Check | Status | Notes |
|:------|:-------|:------|
| Existing C++ tests from the build tree | Passing | 52 tests passed locally during audit. |
| Fresh benchmark-enabled local build | Failing | `benchmark_vcfbatch` currently needs its scoring dependency linked before this target can be treated as replayable. |
| Python tests | Failing | 3 tests failed during audit: one parallel benchmark mock test and two result-loading error-path tests. |
| Astex committed report | Failing placeholder | `benchmark_results/astex_diverse_report.md` records 0/85 successful systems. |
| Expected benchmark outputs | Missing | Several `benchmarks/*/expected/` directories contain placeholders rather than validated reference outputs. |

Until these items are corrected, any numeric performance or accuracy result
should be labeled **preliminary** or **external/manuscript-only**, not
repository-reproducible.

---

## Claim maturity by benchmark family

| Benchmark family | Intended metric family | Repository maturity |
|:-----------------|:-----------------------|:--------------------|
| ITC-187 calorimetry | Delta-G correlation, RMSE, ranking power | Preliminary only until replayable inputs, expected outputs, and metric scripts are committed. |
| CASF-2016 | Scoring, docking, and screening power | Preliminary only until CASF acquisition/preprocessing and expected results are committed. |
| DUD-E | AUC and enrichment factors | Preliminary only until target/decoy provenance and metric scripts are committed. |
| Neurological targets | Pose rescue and entropy contribution analysis | Preliminary only until target list, structures, ligands, configs, and reference results are committed. |
| Hardware dispatch | Shannon entropy throughput by backend | Preliminary only until benchmark binaries build cleanly and machine metadata is recorded. |
| tENCoM vibrational entropy | Runtime and differential entropy behavior | Preliminary only until calibration status and reference systems are documented. |
| VoronoiCFBatch | Batch scoring throughput | Blocked until the benchmark target links and runs cleanly. |

---

## Reproducible benchmark bundle requirements

Every benchmark promoted to repository-reproducible must include:

- dataset provenance or an acquisition script
- immutable identifiers, checksums, or archived input snapshots
- preprocessing commands
- exact command lines
- fixed seeds where stochastic sampling is used
- metric calculation scripts
- expected outputs or metric snapshots
- git SHA and machine/environment details

A benchmark should not publish final table values in this page until its bundle
satisfies those requirements.

---

## Build with benchmarking

```bash
cmake .. -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_TENCOM_BENCHMARK=ON \
    -DENABLE_VCFBATCH_BENCHMARK=ON
cmake --build . -j $(nproc)
```

At the time this page was updated, the `benchmark_vcfbatch` target still needed
a linkage fix before the benchmark-enabled build could be used as validation
evidence.

---

## Run benchmark binaries

```bash
./build/benchmark_tencom
./build/benchmark_vcfbatch
./build/benchmark_dispatch
```

Record the executable path, git SHA, compiler, CPU/GPU model, enabled backends,
thread count, seed values, input bundle checksums, and output artifact hashes
for every run intended to support a public claim.

---

## Test suite

```bash
# C++ validation tests
cmake -DBUILD_TESTING=ON .. && cmake --build . -j $(nproc)
ctest --test-dir build --output-on-failure

# Python validation
cd python && pytest tests/ -q
```

Treat green tests as necessary but not sufficient for scientific benchmark
claims. The benchmark bundle still has to prove provenance, replayability, and
metric calculation integrity.

---

## References

- Gaudreault F & Najmanovich RJ (2015). FlexAID: Revisiting Docking on Non-Native-Complex Structures. *J. Chem. Inf. Model.* 55(7):1323-36. [DOI:10.1021/acs.jcim.5b00078](https://doi.org/10.1021/acs.jcim.5b00078)
- Su M et al. (2019). Comparative Assessment of Scoring Functions: The CASF-2016 Update. *J. Chem. Inf. Model.* 59(2):895-913.
- Mysinger MM et al. (2012). Directory of Useful Decoys, Enhanced (DUD-E). *J. Med. Chem.* 55(14):6582-94.
