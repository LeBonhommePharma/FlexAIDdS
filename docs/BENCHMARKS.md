# Benchmarks

Performance and accuracy benchmarks for FlexAID∆S.

**Claim maturity** follows `docs/REPRODUCIBILITY.md`:

| Label | Meaning |
|-------|---------|
| **Repository-reproducible** | Replayable bundle + commands + artifacts in-repo |
| **Published external** | Peer-reviewed literature (e.g. FlexAID JCIM 2015) |
| **Preliminary / target** | Appears in docs or design notes but **not** backed by a current claim receipt in this checkout |

Agents and skills **must not** present preliminary numbers as measured campaign success. CF/`best_score` is a scoring proxy; ensemble F is not experimental ΔG unless the full ledger path is active and validated.

---

## Accuracy Benchmarks

### FlexAID 2015 Validation Sequence

The packaged benchmark dataset order starts with the three validation sets
used by Gaudreault & Najmanovich in the original FlexAID JCIM 2015 paper:

1. **Astex Diverse Set** (`astex_diverse`, `docking_mode: self_docking`) — 85 native holo complexes.
2. **Astex Non-Native Set** (`astex_nonnative`, `docking_mode: cross_docking`) — 65 targets / ~1112 structures for non-native receptor cross-docking.
3. **HAP2** (`hap2`, `docking_mode: cross_docking`) — holo/apo protein-pair validation for docking into non-native conformations.

These YAML configs ship under `benchmarks/datasets/` (and the Python package mirror). Every YAML declares an explicit `docking_mode`; `validate_dataset_semantics.py` fails closed on contradictions.

**Published external (FlexAID 2015 JCIM, FLRP ideal subset, top-1)** — cite the paper, not as FlexAIDdS live rates:

| Set | Published top-1 (FLRP) | Source field in YAML |
|-----|------------------------|----------------------|
| Astex native | ~0.45 | `astex_diverse.yaml` → `published_baselines` |
| Astex non-native | ~0.39 | `astex_nonnative.yaml` |
| HAP2 | ~0.22 | `hap2.yaml` |

Live FlexAIDdS / classic three-engine **claim** rates require on-disk `result.csv` / `RUN_RECEIPT` + RMSD/PoseBusters gates (see skill deception-proof contract). Do not overwrite published FlexAID numbers with unreceipted memory.

### ITC-187 Calorimetry Benchmark — **PRELIMINARY / TARGET**

Design target for entropy-aware affinity correlation (see `benchmarks/datasets/itc187.yaml` comments). **Not** a published FlexAIDdS Pearson *r*, RMSE, or ranking-power rate — unverified / pending receipt.

| Metric | FlexAID∆S (target / preliminary) | Notes |
|:-------|:--------------------------------:|:------|
| ΔG Pearson *r* | **0.93 (target)** | Training/gate threshold language in tooling; not a live claim without ITC receipt package |
| RMSE (kcal/mol) | **1.4 (target)** | Same |
| Ranking power | **78% (target)** | Same |

Comparator columns previously listed for Vina/Glide are **literature context only**, not a head-to-head receipt in this repository.

Do **not** say “FlexAID∆S achieves 0.93 Pearson” in claim language until a `REPRODUCIBILITY_MANIFEST` + affinity table for ITC-187 is produced under the admission contract.

### CASF-2016 — **PRELIMINARY / TARGET**

Scoring / docking / screening powers in older drafts (r≈0.88, docking ~81%, EF1% ~15) are **not published FlexAIDdS rates** — unverified / pending receipt, not repository-reproducible claim packages. CASF YAML is `docking_mode: affinity_scoring` for scoring-power framing; pose success still needs RMSD+PoseBusters when pose claims are made.

### DUD-E / neurological vignettes — **PRELIMINARY**

Mean AUC / EF and mu-opioid narrative tables in prior docs are illustrative or design-stage. Treat as experimental until a benchmark bundle exists under `benchmarks/` with expected outputs.

---

## Performance Benchmarks

### Hardware Acceleration — Shannon Entropy Computation

Speedup measured on Shannon entropy histogram computation (ShannonThermoStack) over the single-threaded CPU baseline. Treat as **hardware microbench** (not docking success rates).

| Backend | Hardware | Speedup (reported) |
|:--------|:---------|--------:|
| **CUDA** | NVIDIA A100 (80 GB) | **3,575×** |
| **CUDA** | NVIDIA RTX 4090 | **2,890×** |
| **Metal** | Apple M2 Ultra (76-core GPU) | **412×** |
| **Metal** | Apple M3 Max (40-core GPU) | **298×** |
| **AVX-512 + OpenMP** | Dual Xeon 8380 (80 cores) | **187×** |
| **AVX2 + OpenMP** | AMD EPYC 7763 (64 cores) | **142×** |
| **OpenMP** | Intel i9-13900K (24 cores) | **18×** |
| **Scalar** | Single core baseline | 1× |

Re-run `./build/benchmark_dispatch` (or Shannon unit benches) on the current machine before quoting new hardware numbers.

### Unified Hardware Dispatch

The runtime selects among built backends: CUDA → Metal → AVX-512 → AVX2 → OpenMP → scalar when enabled at compile time.

### tENCoM Vibrational Entropy (order-of-magnitude)

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

## Entropy Impact Analysis — **ILLUSTRATIVE**

Historical ranking-change vignettes (HIV-PR, CDK2, etc.) are **not** substitute for Astex claim packages. Prefer ensemble ledger language (F, H, −TS) over “true ΔG”.

---

## Reproducing Benchmarks

### Build with Benchmarking

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_TENCOM_BENCHMARK=ON \
    -DENABLE_VCFBATCH_BENCHMARK=ON
cmake --build build -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
python3 .grok/skills/flexaidds/scripts/resolve_build.py --sync-env --write-pin
export FLEXAIDDS_REQUIRE_BUILD=1
```

### Dataset semantics + skill preflight

```bash
python3 .grok/skills/flexaidds/scripts/validate_dataset_semantics.py
python3 .grok/skills/flexaidds/scripts/validate_skill.py
python3 .grok/skills/flexaidds/scripts/ensure_docking_data.py --check
```

### Any target / any ligand (fast path)

```bash
# Local files
python3 .grok/skills/flexaidds/scripts/dock_any.py \
  --receptor receptor.pdb --ligand ligand.mol2 --temperature 298.15

# Self-docking from RCSB (HET residue code)
python3 .grok/skills/flexaidds/scripts/dock_any.py \
  --pdb 1STP --ligand-res BTN --dry-run
```

### Campaigns

```bash
# Dry-run first
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
  --dataset astex_diverse --tier 1 --dry-run --resume --package

# Local OUT (never claim from iCloud-only live GA trees)
python3 .grok/skills/flexaidds/scripts/dataset_runner.py \
  --dataset astex_diverse --tier 2 --resume --package \
  --results-dir "${FLEXAIDDS_LOCAL_ROOT:-$HOME/flexaidds_results}/benchmarks"
```

### Test Suite (Validation)

```bash
cmake -B build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build -j "$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
ctest --test-dir build --output-on-failure
python3 -m pytest tests/test_flexaid_skill.py tests/test_dataset_semantics.py -q
```

See also: `METHODOLOGY.md`, `AGENTS.md` § Benchmark storage, `.grok/skills/flexaidds/SKILL.md` § Deception-proof claim contract.
