# FlexAIDdS Full Benchmarking Plan
## MacBook Pro 14" M3 Pro · 18 GB Unified RAM · macOS 14+

> **Status:** Ready to execute on iCloud+ 2TB renewal
> **Branch:** `master` (post PR #190 — all C-1–C-5 thermodynamic fixes merged)
> **Primary convergence metric:** Shannon Energy Collapse H(X) < 2 bits
>
> Core thermodynamic APIs report Shannon entropy in **nats**. Convert only at
> reporting or convergence-boundary code with `H_bits = H_nats / ln(2)`.
> **Thermodynamic engine:** Grand Canonical Ensemble (log Ξ, F_bound, binding selectivity)

---

## 1. Hardware Profile & Optimal Configuration

### M3 Pro Topology
| Resource | Spec | FlexAIDdS Mapping |
|---|---|---|
| P-cores | 6 × Firestorm (3.7 GHz) | OpenMP heavy threads — physics, energy eval |
| E-cores | 6 × Icestorm (2.8 GHz) | I/O, PDB parsing, logging |
| GPU | 18-core Apple GPU (Metal 3) | Future Metal compute shaders (currently CPU-side) |
| Unified RAM | 18 GB LPDDR5 @ 150 GB/s | All energy matrices in-core, no swap |
| SIMD | Neon 128-bit (≡ SSE4.2) | `_NEON_` flag replaces AVX2 path |
| L2 cache | 12 MB shared | Fits ~150k energy grid points per tile |

### Recommended Environment Variables
```bash
export OMP_NUM_THREADS=6
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_WAIT_POLICY=passive   # avoid busy-spin on E-core stalls
export FLEXAID_MAX_MEM_MB=16384  # reserve 2 GB for OS
export SHANNON_TRACE_LEVEL=2     # 0=off 1=final 2=per-step 3=debug
export FLEXAID_SIMD=NEON         # disable AVX512 (not on M3), enable NEON
export FLEXAID_SEED=42
```

### ulimits
```bash
ulimit -n 65536   # file descriptors (85 complexes × multiple outputs)
ulimit -s unlimited   # stack (tENCoM deep recursion)
```

---

## 2. Benchmark Suite Inventory

### 2.1 Astex Diverse Set — Native Pose Prediction
- **File:** `benchmarks/astex_diverse/astex_diverse_set.csv`
- **Complexes:** 85 (Hartshorn 2007, J Med Chem 50:726)
- **Success criterion:** top-ranked pose RMSD < 2.0 Å vs. crystal structure
- **Targets:** FlexAIDdS ≥ **65%** (Vina baseline = 58%); hard floor ≥ 58%
- **Key metric:** Per-complex H(X) at convergence + RMSD correlation

### 2.2 Astex Non-Native Set — Cross-Docking
- **File:** `benchmarks/astex_nonnative/astex_non_native_set.csv`
- **Pairs:** 73 cross-docking pairs across 6 target families:
  - CDK2: 1hwi↔1hww, 1oq5↔1sq5
  - ERK2: 1n2v↔2bsm
  - PPAR: 1r1h↔2cl7↔2vd0↔2vd1
  - p38 MAP kinase: 1t40↔1t46↔1y6b↔2cs2
  - Thymidine kinase (HSV-1 pairs)
  - Factor Xa (serine protease pairs)
- **Target:** FlexAIDdS ≥ **35%** (Vina cross-docking baseline ~23%)
- **Key metric:** log Ξ(correct pocket) − log Ξ(alternatives) as selectivity signal

### 2.3 Recommended Additional Datasets (Phase 2)

| Dataset | Complexes | Purpose | Priority |
|---|---|---|---|
| CASF-2016 core set | 285 | Scoring/ranking power | HIGH |
| PDBbind v2020 refined | 5316 | Affinity correlation (Pearson R) | HIGH |
| DUD-E (40 targets) | ~100k decoys | Enrichment / ROC-AUC | MEDIUM |
| LIT-PCBA (15 targets) | ~1.6M | Realistic VS benchmark | MEDIUM |
| Merck FEP benchmark | 8 series | ΔΔG accuracy vs. FEP+ | LOW (future) |

> **Immediate action:** CASF-2016 and PDBbind v2020 core set are publication-critical. Add before thesis submission.

### 2.4 Unit Tests (Gate)
- **Status:** 52/52 passing (2026-05-08, post PR #190)
- **Role:** Must pass before any benchmark run. Abort if any fail.
- **Runtime:** ~48 s on M3 Pro (DatasetRunnerTests dominates)

---

## 3. Shannon Energy Collapse Methodology

### 3.1 Convergence Definition (Rigorous)

```text
H(X_t)      = -sum_i p(x_i,t) * log2(p(x_i,t))    [bits]
H_nats(X_t) = -sum_i p(x_i,t) * ln(p(x_i,t))      [nats]
H_bits      = H_nats / ln(2)
```

where `p(xᵢ,t)` is the Boltzmann-weighted probability of pose cluster *i* at step *t*.
Core scoring and thermodynamic code should keep the natural-log form and convert
to bits only at convergence/reporting boundaries.

- **Convergence threshold:** H(X) < **2.0 bits** — >75% Boltzmann weight in ≤ 2 clusters
- **Hard convergence (thesis-quality):** H(X) < **1.0 bit** — >50% weight in single dominant pose

Physical analogy: entropy collapse = supersaturated solution precipitating one crystal. The search has found its energy funnel — analogous to radar pre-compensation locking onto a single target bearing.

### 3.2 Per-Step Trace Instrumentation

```bash
export SHANNON_TRACE_LEVEL=2   # per-step CSV output
```

Output format per complex:
```
step, H_bits, n_clusters, top_pose_weight, top_pose_rmsd_vs_crystal
0,    6.32,   64,         0.016,           N/A
100,  4.71,   38,         0.031,           3.2
500,  2.88,   14,         0.083,           1.8
1000, 1.43,    5,         0.241,           1.1
1500, 0.87,    2,         0.551,           0.9    ← CONVERGED
```

### 3.3 H(X) vs. RMSD Correlation Analysis

**Hypothesis:** Spearman ρ(H_final, RMSD_top1) > 0 — lower entropy at convergence → better pose.

Physical mechanism: early entropy collapse = narrow energy funnel = correct binding mode. High H_final = rugged landscape or wrong binding site.

**Test:** Spearman rank correlation, 95% bootstrap CI, n=85, 10,000 resamples.

### 3.4 H(X) Profiles by Target Family

| Family | Expected H_final | Rationale |
|---|---|---|
| Kinases (ATP pocket) | < 1.0 bit | Narrow hydrophobic, fast collapse |
| Proteases (S1 pocket) | < 1.5 bits | Well-defined geometry |
| GPCRs | 1.5–3.0 bits | Flexible TM bundle, broader landscape |
| Nuclear receptors | 1.0–2.0 bits | Buried but large cavity |

Plot: H(X) vs. normalized step (0→1), mean ± SD per family.

### 3.5 Shannon-Weighted Success Rate (SWSR)

```
SR_H = Σᵢ (1/H_final,i · successᵢ) / Σᵢ (1/H_final,i)
```

- SR_H > SR: correctly calibrated (high confidence = high accuracy) ✓
- SR_H < SR: overconfident on hard cases — investigate scoring function

```python
def swsr(successes, H_finals):
    w = 1.0 / np.maximum(H_finals, 0.1)   # floor at 0.1 to avoid inf
    w /= w.sum()
    return float(np.dot(w, successes.astype(float)))
```

---

## 4. Thermodynamic Validation Suite

**Gate:** Run before any benchmark. Any failure = abort the run.

### C-1: Boltzmann Weight Integrity
```python
assert all(w > 0.0 for w in weights),   "C-1 FAIL: zero/negative weights"
assert abs(sum(weights) - 1.0) < 1e-9,  "C-1 FAIL: weights don't sum to 1"
assert min(weights) > 1e-300,            "C-1 FAIL: underflow (int multiplicity cast)"
```

### C-2: β = 1/(kB·T) at 300K = 1.6774 mol/kcal
```python
EXPECTED_BETA = 1.0 / (0.001987204 * 300.0)   # 1.6774 mol/kcal
rel_err = abs(beta_reported - EXPECTED_BETA) / EXPECTED_BETA
assert rel_err < 1e-6, f"C-2 FAIL: β={beta_reported:.6f} (expected {EXPECTED_BETA:.6f})"
```

### C-3: Vibrational Entropy — No 2π Bias
```python
# S_vib = kB·(1 - ln(ħω/kBT)), classical limit
# C-3 bug: added ln(2π) ≈ 1.838/mode ≈ 109 kcal/mol per 100 modes
x = (1.4388 * omega_cm) / T    # ħω/kBT (hc/kB = 1.4388 cm·K)
S_correct = kB_kcal * (1.0 - math.log(x))
assert abs(S_reported - S_correct) < 1e-4, "C-3 FAIL: 2π bias present"
```

### C-4: Heat Capacity Cv = Var(E)/(kB·T²)
```python
# Correct: kB·T²  NOT  (kBT)² — off by 503× if bug present
Cv_computed = var_E / (kB * T * T)
Cv_expected = N_particles * 1.5 * kB   # monatomic ideal gas reference
assert abs(Cv_computed - Cv_expected) / Cv_expected < 0.05, \
    f"C-4 FAIL: Cv={Cv_computed:.4f}, expected {Cv_expected:.4f}"
```

### C-5: Torsional Entropy Linear in 1/ω
```python
# S_tor = kB·T/ω  →  S(ω=50)/S(ω=100) = 2.0, NOT 4.0 (quadratic bug)
ratio = S_omega50 / S_omega100
assert abs(ratio - 2.0) < 0.01, f"C-5 FAIL: ratio={ratio:.3f} (expected 2.0)"
```

### Grand Partition Function Sanity
```python
assert abs(log_Xi(c=0.0)) < 1e-12,  "GPF FAIL: log Ξ(c=0) ≠ 0"
assert log_Xi(1e-6) < log_Xi(1e-3) < log_Xi(1.0) < log_Xi(1e3), \
    "GPF FAIL: not monotonically increasing"
```

---

## 5. Storage Budget

### Per-Complex Output (one run)
| File type | Size/complex | ×85 complexes |
|---|---|---|
| Top-10 poses (PDB) | ~250 KB | 21 MB |
| Energy matrix (binary) | ~2 MB | 170 MB |
| Shannon trace (LEVEL=2, CSV) | ~500 KB | 43 MB |
| Trajectory dump (if enabled) | ~50 MB | 4.3 GB |
| Grid/NRGMAP files | ~10 MB | 850 MB |
| **Total (no trajectory)** | | **~1.1 GB** |
| **Total (with trajectory)** | | **~5.4 GB** |

### Full Suite Estimate
| Dataset | Runs | No traj | With traj |
|---|---|---|---|
| Astex Diverse (85) | 85 | 1.1 GB | 5.4 GB |
| Astex Non-Native (73 pairs) | 73 | 0.9 GB | 4.6 GB |
| CASF-2016 (285) | 285 | 3.6 GB | 18 GB |
| PDBbind v2020 (285 × 5 seeds) | 1425 | 18 GB | 90 GB |
| DUD-E (~4k actives) | ~4000 | 50 GB | 250 GB |
| LIT-PCBA (~15k actives) | ~15000 | 185 GB | 925 GB |
| Raw PDB downloads | — | ~200 GB | — |
| **Phase 1 (Astex + CASF)** | | **~6 GB** | **~28 GB** |
| **Full production** | | **~460 GB** | **~1.3 TB** |

**Verdict:** 2 TB iCloud+ is correct. Phase 1 (thesis minimum) fits in ~30 GB. Use `SHANNON_TRACE_LEVEL=1` (final H only) for DUD-E/LIT-PCBA; `LEVEL=2` (per-step) for Astex + CASF-2016 only. Compress PDB source files with `zstd -T0` (~4:1 ratio).

---

## 6. One-Shot Execution Script

```bash
#!/usr/bin/env bash
# run_full_benchmark.sh — FlexAIDdS full benchmark, MacBook Pro M3 Pro
# Usage: caffeinate -i ./run_full_benchmark.sh [PHASE=1|2] [SEED=42]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD_DIR="${REPO_ROOT}/build"
RESULTS_DIR="${REPO_ROOT}/benchmark_results/$(date +%Y%m%d_%H%M%S)"
SEED="${SEED:-42}"
PHASE="${PHASE:-1}"

# M3 Pro optimal
export OMP_NUM_THREADS=6
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export OMP_WAIT_POLICY=passive
export FLEXAID_MAX_MEM_MB=16384
export FLEXAID_SIMD=NEON
export FLEXAID_SEED="${SEED}"
ulimit -n 65536; ulimit -s unlimited
mkdir -p "${RESULTS_DIR}"
echo "Results: ${RESULTS_DIR} | Phase: ${PHASE} | Seed: ${SEED}"

# ── Step 0: Unit tests (GATE) ──────────────────────────────────────────────────
cd "${BUILD_DIR}"
ctest --output-on-failure -j6 --timeout 120 \
    || { echo "UNIT TEST FAILURE — aborting"; exit 1; }

# ── Step 1: Thermodynamic validation ──────────────────────────────────────────
export SHANNON_TRACE_LEVEL=2
python3 "${REPO_ROOT}/tests/thermodynamic_validation.py" \
    --build-dir "${BUILD_DIR}" \
    --output "${RESULTS_DIR}/thermo_validation.json"

# ── Step 2: Astex Diverse (85 complexes) ──────────────────────────────────────
time "${REPO_ROOT}/benchmarks/astex_diverse/run_astex_diverse.sh" \
    --build-dir "${BUILD_DIR}" --output-dir "${RESULTS_DIR}/astex_diverse" \
    --seed "${SEED}" --nthreads 6
python3 "${REPO_ROOT}/tests/benchmarks/astex_diverse/evaluate.py" \
    "${RESULTS_DIR}/astex_diverse" \
    --rmsd-threshold 2.0 --bootstrap-n 10000 --shannon-weighted \
    --output "${RESULTS_DIR}/astex_diverse_report.json"

# ── Step 3: Astex Non-Native (73 pairs) ───────────────────────────────────────
export SHANNON_TRACE_LEVEL=1   # final H only for cross-docking
time "${REPO_ROOT}/benchmarks/astex_nonnative/run_astex_non_native.sh" \
    --build-dir "${BUILD_DIR}" --output-dir "${RESULTS_DIR}/astex_nonnative" \
    --seed "${SEED}" --nthreads 6
python3 "${REPO_ROOT}/tests/benchmarks/astex_nonnative/evaluate.py" \
    "${RESULTS_DIR}/astex_nonnative" \
    --rmsd-threshold 2.0 --bootstrap-n 10000 \
    --output "${RESULTS_DIR}/astex_nonnative_report.json"

# ── Step 4: CASF-2016 (if downloaded) ─────────────────────────────────────────
if [[ -d "${REPO_ROOT}/benchmarks/casf2016" ]]; then
    export SHANNON_TRACE_LEVEL=2
    time "${REPO_ROOT}/benchmarks/casf2016/run_casf2016.sh" \
        --build-dir "${BUILD_DIR}" --output-dir "${RESULTS_DIR}/casf2016" \
        --seed "${SEED}" --nthreads 6
fi

# ── Step 5: Unified report ─────────────────────────────────────────────────────
python3 "${REPO_ROOT}/scripts/generate_benchmark_report.py" \
    --results-dir "${RESULTS_DIR}" \
    --output-html "${RESULTS_DIR}/REPORT.html" \
    --output-md   "${RESULTS_DIR}/REPORT.md"
echo "══════ BENCHMARK COMPLETE → ${RESULTS_DIR}/REPORT.md ══════"
cat "${RESULTS_DIR}/REPORT.md"
```

### 6.1 Computational Timing Profile — Source Analysis

> **Derived from `gaboom.cpp`, `CoarseScreen.cpp`, `MIFGrid.h`, `generate_grid.cpp`, `partition_grid.cpp`**
> **All estimates for M3 Pro Firestorm 3.7 GHz, OMP_NUM_THREADS=6, no GPU offload**

#### Total eval_chromosome calls per complex

| Event | Chromosome count |
|---|---|
| Initial population | 1,000 |
| GA loop (500 gen × 1000 chrom, BOOM fraction=1.0) | 500,000 |
| Grid repartitions (10 events × 900 repop chrom, popszpartition=100) | 9,000 |
| **Total per complex** | **510,000** |

All via `#pragma omp parallel for schedule(dynamic)` across 6 P-cores.
Dirty-tracking fires when `dirty_atm.size() < natm/2` (normal modes off) → only <10% of atoms restored between consecutive evals.

#### Per-step wall-clock (single-thread estimates)

| Step | Cost | Notes |
|---|---|---|
| **vcfunction (Vcontacts scoring)** | **~1–5 ms** | Voronoi cell construction + LJ contact areas; dominant hot path |
| ic2cf (internal → Cartesian) | ~1–5 µs | Rotation matrix chain for torsion genes |
| Thread workspace restore (dirty) | ~2–5 µs | ~300 atom struct memcopy |
| optres ptr redirect + cf clear | ~1–2 µs | natm loop + n_optres × 9 field zero |
| Shannon entropy check (per gen) | ~0.1 ms | O(1000 chrom) H(X) computation |
| InStreamClustering (per event) | ~3–10 ms | O(N_elite²) pairwise RMSD, every 100 gen |
| partition_grid (top-100 elites) | ~5 ms | std::map over popszpartition=100 |
| slice_grid (midpoint insertion) | ~10 ms | O(N_pts) midpoint map iteration |
| compute_mif (OMP, adapted grid) | ~8–15 ms | LJ probe scan 6-thread, cutoff=6.0 Å |
| generate_grid (0.375 Å, ~79k pts) | ~20 ms | Triple while-loop sphere masking |
| CoarseScreener setup | ~5 ms | 6.56 Å cells × 39 atom types × 27 neighbors, OMP |
| build_clash_grid (0.25 Å) | ~2–5 ms | OMP over bounding-box voxels |
| PDB parse + data structure init | ~2–5 s | File I/O + atom/residue array setup |
| FastOPTICS/DBSCAN clustering | ~100–300 ms | O(N_pop²) final ensemble reduction |

#### Per-phase breakdown (6 OMP P-cores, ~80% parallel efficiency)

| Phase | Events | Single-thread cost | Effective w/ 6 OMP | Phase total |
|---|---|---|---|---|
| **P0 — Setup** | 1× | — | — | **~5–7 s** |
| PDB parse + init | 1 | 2–5 s | serial | 3 s |
| generate_grid + MIF + CDF | 1 | 45 ms | serial | 45 ms |
| CoarseScreen (grid + clash) | 1 | 10 ms | serial | 10 ms |
| Initial pop eval (1000 chrom OMP) | 1,000 | ~3 ms | 0.6 ms eff. | **0.6 s** |
| **P1 — GA loop (500 gen)** | 500,000 | ~3 ms | 0.6 ms eff. | **~5 min** |
| eval_chromosome dominates | | | | |
| Entropy checks (every 10 gen, ×50) | 50 | 0.1 ms | serial | 5 ms |
| InStreamClustering (every 100 gen, ×5) | 5 | ~5 ms | serial | 25 ms |
| Stagnation check (every gen, ×500) | 500 | 0.01 ms | serial | 5 ms |
| **P2 — Grid repartitions (every 50 gen, ×10)** | 10× | — | — | **~15–20 s** |
| QuickSort + partition_grid + slice | 10 | 25 ms | serial | 250 ms |
| compute_mif (adapted grid, OMP) | 10 | 12 ms | 2.5 ms | 25 ms |
| Repopulate 900 chrom (OMP) | 9,000 | ~3 ms | 0.6 ms eff. | **5.4 s** |
| **P3 — Post-processing** | 1× | — | — | **~400 ms** |
| FastOPTICS/DBSCAN | 1 | 200 ms | — | 200 ms |
| Grand canonical log Ξ | 1 | 50 ms | — | 50 ms |
| Write output files | 1 | 100 ms | — | 100 ms |
| **TOTAL PER COMPLEX** | | | | **~5–8 min** |

#### Cost attribution

```
vcfunction (Vcontacts)    ≈ 92–95% of wall time   ← the ONLY lever that matters
Grid repartition repops   ≈  3–5%
Setup/parse               ≈  1–2%
Everything else           ≈ <1%
```

> **Vectorization target:** Vcontacts on Neon SIMD is the singular optimization opportunity.
> `partition_grid`, `compute_mif`, and `slice_grid` combined < 5% — do not optimize those.

---

### 6.2 Wall-Clock Estimates (M3 Pro, 6 threads)

**Source-analysis range:** 5–8 min/complex. Calibrate with the pilot run in Section 9 before committing to the full suite.

| Step | 5 min/complex | 8 min/complex |
|---|---|---|
| Unit tests + thermo validation | ~53 s | ~53 s |
| Pilot calibration (5 Astex complexes) | ~25 min | ~40 min |
| Astex Diverse (85 complexes) | ~7 h | ~11.3 h |
| Astex Non-Native (73 pairs) | ~6 h | ~9.7 h |
| CASF-2016 (285 complexes) | ~23.8 h | ~38 h |
| **Phase 1 minimum (Astex Diverse only)** | **~7 h** | **~11.3 h** |
| **Phase 1 full (+ CASF-2016)** | **~31 h** | **~49 h** |

> Run with `caffeinate -i ./run_full_benchmark.sh` to prevent sleep.
> After the pilot, replace the range above with the empirical `wall_time_s` mean × N_complexes.

---

## 7. Statistical Methodology

### 7.1 Bootstrap CI (always — n=85 too small for normal approximation)
```python
import numpy as np

def bootstrap_sr(successes, n_boot=10000, ci=0.95):
    n = len(successes)
    boot = [np.random.choice(successes, n, replace=True).mean()
            for _ in range(n_boot)]
    a = (1 - ci) / 2
    return successes.mean(), np.quantile(boot, a), np.quantile(boot, 1-a)

# Example: 57/85 → SR = 67.1% (95% CI: 57.6%–76.5%)
```

### 7.2 Effect Size vs. Vina Baseline (Fisher's Exact)
```python
from scipy.stats import fisher_exact

def vs_baseline(n_success, n_total, n_base_success, n_base_total):
    table = [[n_success,      n_total      - n_success],
             [n_base_success, n_base_total - n_base_success]]
    OR, p = fisher_exact(table, alternative='greater')
    delta = n_success/n_total - n_base_success/n_base_total
    return delta, OR, p

# Astex: FlexAIDdS vs. Vina (Vina: 49/85 = 57.6%)
delta, OR, p = vs_baseline(57, 85, 49, 85)
# → ΔRMSD-SR = +9.4pp, OR = 1.54, p = 0.041
```

### 7.3 Shannon-Weighted Success Rate
```python
def swsr(successes, H_finals):
    w = 1.0 / np.maximum(H_finals, 0.1)   # floor at 0.1 to avoid inf
    w /= w.sum()
    return float(np.dot(w, successes.astype(float)))
# SWSR > SR: calibrated. SWSR < SR: overconfident failures → investigate.
```

### 7.4 Spearman ρ(H_final, RMSD)
```python
from scipy.stats import spearmanr
rho, p = spearmanr(H_finals, RMSDs)
# Expected: ρ > 0 (higher entropy at convergence → worse pose)
# |ρ| < 0.2 AND p > 0.05 → H(X) not predictive → investigate scoring function
```

### 7.5 Report Full RMSD Distribution
Do not report binary SR only. Include:
- Fraction with RMSD < 1.0 Å (near-crystal)
- Fraction with RMSD < 2.0 Å (standard success — primary)
- Fraction with RMSD < 3.0 Å (near-success)
- Fraction with RMSD > 5.0 Å (catastrophic failures)
- Median RMSD across all complexes

---

## 8. Publication-Readiness Checklist

### Code Integrity
- [x] C-1–C-5 thermodynamic bugs fixed (PR #190, merged 2026-05-08)
- [x] 52/52 unit tests passing on master (2026-05-08)
- [ ] Numerical validation suite passes (Section 4 — C-1 through C-5 + GPF)
- [ ] Git tag: `git tag v1.0.0-benchmark` before benchmark run

### Reproducibility
- [ ] Fixed seed=42 documented in results header
- [ ] Exact commit hash recorded in results directory
- [ ] Hardware + macOS version documented
- [ ] Input PDBs from Astex archive (DOI: 10.1021/jm060066m) — no re-processed structures
- [ ] Zenodo or OSF deposit with DOI (required for thesis citation)

### Statistical Rigor
- [ ] Bootstrap CI on all success rates (not Wilson/normal approximation)
- [ ] Fisher's exact test vs. Vina baseline (p < 0.05 threshold)
- [ ] SWSR calibration check (SR_H vs. SR direction)
- [ ] Spearman ρ(H_final, RMSD) with bootstrap CI
- [ ] Full RMSD distribution reported (not just binary)

### Baseline Comparisons
- [ ] Vina 1.2.7 results on same 85 complexes (reproduce or cite Hartshorn 2007 Table 3)
- [ ] Explicit statement: "Vina baseline from Hartshorn 2007, reproduced with Vina 1.2.7"
- [ ] Optional: Smina, GNINA-CNN for extended comparison figure

### Limitations Section (Required in Thesis)
- [ ] Astex Diverse (2007) structures are widely known — not a blind challenge
- [ ] Cross-docking 73 pairs across 6 families — not a random sample
- [ ] No explicit solvation (PBSA/GBSA) — Grand Canonical implicit treatment described
- [ ] M3 Pro timings not transferable to HPC cluster benchmarks

### Data Availability
- [x] Benchmark CSVs committed in `benchmarks/` (PR #192)
- [x] Evaluation scripts in `tests/benchmarks/` (PR #192)
- [ ] Raw results (poses, H traces, energy matrices) deposited on Zenodo (DOI)

---

## 9. Pre-Flight Checklist (Day-of Execution)

```bash
#!/usr/bin/env bash
cd ~/Projects/FlexAIDdS

# 1. Clean master, correct commit
git status
git log --oneline -3

# 2. Build is current
ls -la build/FlexAIDdS
cmake --build build --parallel 6 --target FlexAIDdS

# 3. 52/52 unit tests pass
cd build && ctest --output-on-failure -j6 && cd ..

# 4. Storage ≥ 50 GB free
df -h ~/

# 5. Pilot calibration — 5 Astex complexes, measure wall_time_s
# Recalibrates the 5–8 min/complex estimate from source analysis (Section 6.1)
PILOT_COMPLEXES=$(head -6 benchmarks/astex_diverse/astex_diverse_set.csv | tail -5 | cut -d',' -f1)
for pdb in $PILOT_COMPLEXES; do
    time ./build/FlexAIDdS         --input  benchmarks/astex_diverse/structures/${pdb}/${pdb}_protein.pdb         --ligand benchmarks/astex_diverse/structures/${pdb}/${pdb}_ligand.mol2         --output /tmp/pilot_${pdb}         --seed 42 --nthreads 6 2>&1 | tee /tmp/pilot_${pdb}.log
done
# Extract wall_time_s from DockingResult output and compute mean:
grep -h "wall_time_s" /tmp/pilot_*.log | awk '{sum+=$2; n++} END {
    mean=sum/n;
    printf "Pilot mean: %.1f s/complex  →  %.1f h for 85  →  %.1f h for 85+73+285
",
    mean, 85*mean/3600, (85+73+285)*mean/3600}'

# 6. Tag before run
git tag "benchmark-$(date +%Y%m%d)" -m "Pre-benchmark snapshot"

# 7. Run (prevent sleep)
caffeinate -i ./run_full_benchmark.sh PHASE=1 SEED=42
```

---

## 10. Phase Roadmap

| Phase | Datasets | Goal | When |
|---|---|---|---|
| **0** | Unit tests + thermo validation | Gate — always first | Always |
| **1** | Astex Diverse + Non-Native | Thesis minimum | iCloud+ renewal |
| **1b** | + CASF-2016 (285) | Scoring/ranking power | After CASF download |
| **2** | + PDBbind v2020 refined | Affinity correlation (Pearson R) | Post-thesis prep |
| **3** | + DUD-E + LIT-PCBA | Enrichment / virtual screening | Paper submission |

---

*Generated: 2026-05-08 | Updated: 2026-05-08 (timing profile from source analysis — gaboom.cpp, CoarseScreen.cpp, MIFGrid.h)*
*FlexAIDdS master post-PR #190 + PR #192 | M3 Pro benchmarking plan v1.1*
*Shannon Energy Collapse H(X) < 2 bits: primary convergence metric throughout all phases*
*vcfunction (Vcontacts) = 92–95% of wall time — sole vectorization target for performance*
