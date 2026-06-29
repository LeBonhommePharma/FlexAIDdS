# Linux CPU/CUDA Performance Baseline Runbook

**Audience:** Operator on Ubuntu x86-64 (local machine, HPC node, or CI)  
**macOS note:** Dock timing harvest is read-only safe on macOS; **build + microbench capture require Linux.**

Cross-reference: `docs/PERFORMANCE_SWARM_AUDIT_2026.md`, `results/perf_swarm/README.md`

---

## Overview

| Artifact | Platform | How populated |
|----------|----------|---------------|
| `baseline_linux_cpu.json` | AVX2 + OpenMP | Microbench + optional dock harvest |
| `baseline_linux_cuda.json` | CUDA + AVX2 | GPU machine only |
| `baseline_macos_metal.json` | Metal (done) | 182 dock timings harvested on macOS |

Two metric families share schema v1.0.0:

1. **Microbenchmarks** — `benchmark_tencom`, `benchmark_vcfbatch` → `benchmarks[]`
2. **Dock timings** — `stderr.log` TIMING SUMMARY → `dock_timings_harvested`

---

## 1. Linux CPU — dependencies

```bash
sudo apt-get update
sudo apt-get install -y \
  cmake ninja-build libeigen3-dev libomp-dev gcc-14 g++-14 git
export CC=gcc-14 CXX=g++-14
```

---

## 2. Linux CPU — exact CMake flags

```bash
cd /path/to/FlexAIDdS

cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=gcc-14 \
  -DCMAKE_CXX_COMPILER=g++-14 \
  -DFLEXAIDS_USE_CUDA=OFF \
  -DFLEXAIDS_USE_METAL=OFF \
  -DFLEXAIDS_USE_AVX2=ON \
  -DFLEXAIDS_USE_OPENMP=ON \
  -DBUILD_TESTING=ON \
  -DENABLE_TENCOM_BENCHMARK=ON \
  -DENABLE_VCFBATCH_BENCHMARK=ON

cmake --build build -j "$(nproc)" \
  --target benchmark_tencom benchmark_vcfbatch
```

Shortcut (same flags):

```bash
chmod +x benchmarks/linux/build_linux_cpu.sh
./benchmarks/linux/build_linux_cpu.sh
```

---

## 3. Run microbenchmarks

```bash
cd build
./benchmark_tencom | tee tencom_bench.txt
./benchmark_vcfbatch 200 20 | tee vcfbatch_bench.txt
```

**Reference row:** TeNCoM uses `N_res=200` for `build_ms_full` / `sample_ms_full`.  
**VCF batch:** pop=200, genes=20; metric is OpenMP serial/parallel speedup.

---

## 4. Capture microbench baseline

```bash
python3 scripts/capture_microbench_baseline.py \
  --label linux_cpu \
  --merge results/perf_swarm/baseline_linux_cpu.json \
  --tencom build/tencom_bench.txt \
  --vcfbatch build/vcfbatch_bench.txt
```

Commit the updated `results/perf_swarm/baseline_linux_cpu.json` when metrics are real numbers (not `null`).

---

## 5. Harvest dock timings (optional, campaign results)

After a `benchmark_datasets` campaign on Linux:

```bash
export FLEXAIDDS_RESULTS_ROOT=/path/to/campaign/results

python3 scripts/harvest_perf_baselines.py \
  --label linux_cpu \
  --results-root "${FLEXAIDDS_RESULTS_ROOT}" \
  --out results/perf_swarm
```

This writes `baseline_linux_cpu.json` **or** overwrites if same label — use `--merge` flow:

1. Capture microbench first (step 4).
2. Harvest dock timings to a temp file, then merge `dock_timings_harvested` manually, **or** re-run harvest with a distinct label (`linux_cpu_dock`) and merge in a follow-up PR.

Compare pre/post dock campaigns:

```bash
python3 scripts/compare_dock_timings.py \
  --baseline results/perf_swarm/baseline_linux_cpu.json \
  --current results/perf_swarm/baseline_linux_cpu_post_p1.json
```

---

## 6. Linux CUDA — manual (no GitHub-hosted GPU)

```bash
# NVIDIA driver + CUDA toolkit (nvcc in PATH)
./benchmarks/linux/build_linux_cpu.sh --cuda

cd build
./benchmark_tencom | tee tencom_bench.txt
./benchmark_vcfbatch 200 20 | tee vcfbatch_bench.txt

python3 scripts/capture_microbench_baseline.py \
  --label linux_cuda \
  --merge results/perf_swarm/baseline_linux_cuda.json \
  --tencom build/tencom_bench.txt \
  --vcfbatch build/vcfbatch_bench.txt
```

`.github/workflows/perf.yml` includes a disabled `benchmark-cuda` job template for `runs-on: [self-hosted, linux, cuda]`. Set `if: false` → `if: true` when a GPU runner is registered.

---

## 7. CI workflow (GitHub Actions)

**Workflow:** `.github/workflows/perf.yml` (manual `workflow_dispatch`)

| Input | Default | Purpose |
|-------|---------|---------|
| `baseline_sha` | `main` | Commit with populated `baseline_linux_cpu.json` |
| `platform` | `linux_cpu` | Baseline file name suffix |

**Dispatch from GitHub UI:** Actions → Performance Regression → Run workflow.

**Compare locally against committed baseline:**

```bash
python3 scripts/compare_perf_baseline.py \
  --baseline results/perf_swarm/baseline_linux_cpu.json \
  --tencom build/tencom_bench.txt \
  --vcfbatch build/vcfbatch_bench.txt
```

Regression threshold: **±5%** per metric (`scripts/compare_perf_baseline.py`).

---

## 8. Tier-1 paired smoke (accuracy + timing gate)

For SoA / P0 validation on Linux (longer run):

```bash
# Build with benchmark_datasets enabled
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DFLEXAIDS_USE_AVX2=ON -DFLEXAIDS_USE_OPENMP=ON \
  -DENABLE_BENCHMARK_DATASETS=ON \
  -DFLEXAIDS_USE_CUDA=OFF -DFLEXAIDS_USE_METAL=OFF
cmake --build build -j "$(nproc)" --target benchmark_datasets FlexAIDdS

python3 scripts/run_perf_validation_smoke.py \
  --results-root /path/to/perf_swarm_validation \
  --build-scalar build \
  --nice 19
```

---

## 9. Checklist before marking Linux baselines “done”

- [ ] `baseline_linux_cpu.json` — `benchmarks[].metrics` populated (no `null`)
- [ ] `baseline_linux_cpu.json` — `git.commit` matches baseline capture commit
- [ ] `perf.yml` compare passes at `baseline_sha` with ≤5% drift
- [ ] (Optional) `dock_timings_harvested.count` > 0 for end-to-end dock regression
- [ ] (Optional) `baseline_linux_cuda.json` on GPU host

---

## 10. macOS operator (read-only harvest only)

```bash
# Harvest existing campaign logs without launching docks
python3 scripts/harvest_perf_baselines.py \
  --label macos_metal \
  --results-root "${FLEXAIDDS_RESULTS_ROOT:-$HOME/Documents/PhD/Programs/FlexAIDdS/results}"

python3 scripts/compare_dock_timings.py \
  --baseline results/perf_swarm/baseline_macos_metal.json \
  --current results/perf_swarm/baseline_macos_metal_post_p0.json
```

Do **not** expect Linux microbench or CUDA numbers from macOS.