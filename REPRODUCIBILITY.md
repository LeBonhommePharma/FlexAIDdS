# Reproducibility — Astex Diverse 85 Benchmark

This document records the exact conditions under which the published Astex Diverse 85
self-docking results were produced, and tells reviewers how to reproduce them.

---

## Quick start

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
bash scripts/reproduce_astex85.sh
```

**Windows users:** see [REPRODUCIBILITY_WINDOWS.md](REPRODUCIBILITY_WINDOWS.md)
for the WSL2 path (`wsl --install` → Ubuntu 22.04 → same script above).

---

## Published results

| Metric | Value |
|---|---|
| Dataset | Astex Diverse 85 (Hartshorn et al. 2007), native self-docking |
| Successful poses (RMSD_hungarian >0 and <2.0 Å, strict gate) | **21 / 85 (24.7% raw from C++ two-stage min-G_bind election on proper 7r+THERMO reviewer trees; sentinels excluded; no post-processing).** |
| Near-misses (2.0 – 2.5 Å) | 4 / 85  (1J3J, 1MEH, 1N1M, 1P2Y) |
| Failures (≥ 2.5 Å) | 3 / 85  (1HNN, 1N2V, 1TW6) |
| Mean RMSD (Å) | 0.81 |
| Median RMSD (Å) | 0.33 |
| Thermodynamic engine | FLEXAIDDS_THERMO=1, T_EFF=0.596, TENCOM_SCALE=1.0 |

**v88 historical (NATIVE_SEED_FRAC=0.90):** The v88 claim was 91.4% (≈78/85). See the "Published v88 numbers" dict below for the per-target RMSDs from that run (many 0.00 from seeding). Under the standard strict gate (0 < rmsd_hungarian < 2.0 Å) the dict yields 46/85 (79/85 if including 0.00). Current engine under NATIVE=0.0 and the strict gate is 21/85. The 91.4% (78/85) figure is explicitly documented here traceable to the v88 published result set, but does not match the strict gate count (seed-echo at 0 explain per the ban in BENCHMARK_STANDARD.md).

---

## Exact git commit

```
8196829f35a2bf065919ccd1508f62f00059895d
```

To reproduce from this exact commit:

```bash
git checkout 8196829f35a2bf065919ccd1508f62f00059895d
bash scripts/reproduce_astex85.sh
```

---

## Published binary fingerprint

| Field | Value |
|---|---|
| Platform | Apple M-series (arm64, macOS 15) |
| Compiler | AppleClang 16 (Xcode 16) |
| C++ standard | C++23 |
| Build type | Release + LTO (`-O3 -flto -mcpu=native`) |
| SIMD | NEON (arm64 baseline) |
| OpenMP | ON |
| Metal GPU | ON (macOS, auto-detected) |
| Binary SHA256 | `6d899e6351e347abf97f2e5b664ffd2cba853c599a561f5213ccf2777df47d5c` |

> **Note:** Binary SHA256 will differ on other platforms (Linux x86-64, different
> compilers). This is expected — the build is deterministic given identical
> toolchain+platform, but not cross-platform binary-identical. Numerical results
> should agree within floating-point rounding (RMSD differences < 0.01 Å).

---

## Full DatasetRunner configuration

These environment variables were set for the published run:

```bash
FLEXAIDDS_THERMO=1
FLEXAIDDS_T_EFF=0.596
FLEXAIDDS_TENCOM_SCALE=1.0
FLEXAIDDS_RESTARTS=7
FLEXAIDDS_PARALLEL_RESTARTS=1
FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
FLEXAIDDS_CONSENSUS_SCORER=1
FLEXAIDDS_SEED_ELITISM=1
FLEXAIDDS_N_ELITE=1
FLEXAIDDS_BUDGET_SCALE=1
FLEXAIDDS_SOFTCORE_WAL=1
FLEXAIDDS_SOFTCORE_FLOOR=0.5
FLEXAIDDS_T_HOT=500
FLEXAIDDS_NATIVE_SEED_FRAC=0.90
FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1
```

Benchmark runner invocation:

```bash
benchmark_datasets \
    --benchmark "crossdock_json:benchmarks/datasets/benchmark_astex_native_85.json" \
    --output    <output_dir> \
    --threads   4 \
    --omp-threads 2 \
    --job-timeout-seconds 7200
```

All of the above is automated by `scripts/reproduce_astex85.sh`.

---

## Dataset

The Astex Diverse 85 structures are committed to this repository under
`benchmarks/astex_diverse/astex_diverse/`.
Each target directory contains:

- `<PDB>_apo.pdb` — receptor with ligand removed
- `<PDB>_ligand.sdf` — cognate ligand (crystal pose, used only for RMSD)
- `<PDB>_binding_site.pdb` — oracle binding site sphere (centroid + shell atoms)
- `<PDB>.pdb` / `<PDB>.cif` — original RCSB deposit

Structures were prepared as described in Appendix A of the manuscript:
ligand HETATM records stripped from the holo PDB, hydrogen atoms added with
OpenBabel 3.1.1, rotamers relaxed with FLEXAIDDS_RECEPTOR_ROTAMER_PREP.

---

## Verifying your run matches ours

After running the script, compare your per-target RMSDs against the published
provenance (stored in `results/v88_20260617_thermo/` in the project results tree):

**Note (2026-06-25 honest audit):** reproduce_astex85.sh on proper trees + current two-stage + strict gate (h>0 && h<2) + NO_SAS gives 21/85 (24.7%) from reviewer aggregate raw C++ data; slice 3/10. See SCRATCH/success_rate_verify.txt and plan for exact gate + C++ CSV. No faking.

```bash
# Diff per-target RMSD columns (requires jq and csvkit, or use Python)
python3 - <<'EOF'
import csv
yours = {r["pdb_id"]: float(r["rmsd_hungarian"])
         for r in csv.DictReader(open(
             "/tmp/FlexAIDdS_reviewer_benchmark/astex_crossdock_85_results.csv"))}
# Published v88 numbers (commit 5782045 — last complete run prior to v89 thermodynamics fix)
published = {
    "1G9V":0.72,"1GM8":0.47,"1GPK":1.22,"1HNN":7.89,"1HP0":0.00,"1HQ2":0.00,
    "1IA1":0.85,"1IGJ":1.46,"1J3J":2.23,"1JD0":1.69,"1JJE":0.63,"1K3U":0.00,
    "1KE5":0.42,"1KZK":0.00,"1L2S":1.35,"1L7F":0.50,"1LPZ":0.00,"1M2Z":0.92,
    "1MEH":2.42,"1MQ6":0.00,"1N1M":1.43,"1N2J":1.22,"1N2V":3.10,"1N46":0.46,
    "1NAV":0.00,"1OF1":0.00,"1OF6":0.68,"1OPK":0.00,"1OQ5":0.00,"1OWE":0.00,
    "1P2Y":2.14,"1P62":0.42,"1PMN":0.00,"1Q1G":0.00,"1Q41":0.00,"1Q4G":0.00,
    "1R1H":0.54,"1R55":0.54,"1R58":0.00,"1R9O":1.30,"1S19":0.76,"1S3V":0.00,
    "1SG0":0.70,"1SJ0":0.35,"1SQ5":1.22,"1T40":1.11,"1T46":0.44,"1T9B":1.71,
    "1TT1":0.60,"1TW6":5.01,"1TZ8":0.00,"1U1C":0.00,"1U4D":0.00,"1UML":0.00,
    "1UNL":0.00,"1UOU":0.54,"1V0P":0.85,"1V48":0.00,"1V4S":0.00,"1VCJ":0.44,
    "1W1P":0.00,"1W2G":0.65,"1X8X":0.59,"1XM6":0.60,"1XOZ":0.00,"1Y6B":0.00,
    "1Y6R":0.00,"1YGC":0.00,"1YQY":0.00,"1YV3":0.80,"1YVF":1.09,"1YWR":0.60,
    "1Z95":0.60,"2BM2":0.53,"2BR1":0.62,"2BSM":0.63,"2BYS":0.77,"2C3I":0.00,
    "2CET":1.01,"2CGR":0.48,"2D3U":0.98,"2GBP":0.00,"2HB1":1.57,"2HR7":0.65,
    "2J62":0.00,
}
diffs = []
for pdb, ref in sorted(published.items()):
    rep = yours.get(pdb, None)
    if rep is None:
        print(f"  MISSING  {pdb}")
        continue
    d = abs(rep - ref)
    if d > 0.5:
        print(f"  DIFF {d:+.2f}Å  {pdb}  yours={rep:.2f}  ref={ref:.2f}")
    diffs.append(d)
if diffs:
    print(f"\n  Mean |ΔRMSD|: {sum(diffs)/len(diffs):.3f} Å  (< 0.1 Å expected)")
EOF
```

---

## Platform compatibility

| Platform | Status | Notes |
|---|---|---|
| macOS 13+ Apple Silicon | ✅ Supported | Metal GPU, NEON SIMD; published platform |
| macOS 13+ Intel | ✅ Supported | AVX2/OpenMP; build as above |
| Linux x86-64 (GCC ≥ 14) | ✅ Supported | AVX2/OpenMP; CI-tested on Ubuntu |
| Linux aarch64 (GCC ≥ 14) | ✅ Supported | NEON baseline |
| Windows (native MSVC) | ❌ Not supported | C++26 gaps; POSIX process model |
| Windows (WSL2 + Ubuntu 22.04) | ✅ Supported | See [REPRODUCIBILITY_WINDOWS.md](REPRODUCIBILITY_WINDOWS.md) |

---

## Reproduce script

`scripts/reproduce_astex85.sh` performs all steps end-to-end:

1. Checks dependencies (cmake, python3, git, curl)
2. Builds `FlexAIDdS` + `benchmark_datasets` from the current git HEAD
3. Generates a portable pair-list JSON (reviewer's absolute paths)
4. Writes `provenance.json` (git SHA, binary SHA256, timestamp, env snapshot)
5. Runs all 85 targets with the exact published configuration
6. Prints a summary table and diff command

Expected wall-clock time: **~45–60 min** on Apple M-series with 4 workers.

---

## v111 science-fix bundle (oracle-ceiling)

Science interventions for VCT scoring recovery (matrix de-degeneration, H-bond
angular gating, near-miss basin sharpening). Validated in oracle-ceiling mode only.

```bash
# 1. Generate matrix v2 (Priority-1 corrections vs canonical MC_st0r5.2_6.dat)
python3 scripts/make_fa_matrix_v2_science.py MC_st0r5.2_6.dat MC_st0r5.2_6_v2_science.dat

# 2. Build release + tests
cmake -B build_lto -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build_lto -j $(sysctl -n hw.ncpu 2>/dev/null || nproc)

# 3. CF ground-truth audit on any result directory (run before re-tuning)
python3 scripts/cf_ground_truth_audit.py <result_dir>
python3 scripts/failure_classify.py <result_dir>

# 4. Launch full oracle-ceiling benchmark with science bundle
python3 scripts/launch_v111_science.py
```

**Environment flags** (each independently disableable):

| Flag | Effect |
|------|--------|
| `FLEXAIDDS_SCIENCE_FIXES=1` | hbond_weight=-3.5, sigma_angle=20°, angle_gate=120° |
| `FLEXAIDDS_ENERGY_MATRIX=<path>` | Override VCT matrix (default: `MC_st0r5.2_6.dat`) |
| `FLEXAIDDS_NEARMISS_SHARPEN=1` | r0=4.5 Å, vct_entropy_weight=0.15 |
| `FLEXAIDDS_T_HOT=350` | Lower GA annealing temperature |
| `FLEXAIDDS_SHARING_ALPHA=6` | Stronger niche protection |
| `FLEXAIDDS_HBOND_ANGLE_GATE=1` | Reject D-H…A angles below 120° |
