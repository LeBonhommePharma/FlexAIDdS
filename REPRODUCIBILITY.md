# Reproducibility — Astex Diverse 85 Benchmark

This document tells reviewers how to run a **blind** Astex Diverse 85 self-docking
campaign and how to read historical provenance from a **withdrawn** oracle run.

**This repository publishes no Astex-85 success rate.** Measurement rules live in
`METHODOLOGY.md` §0 (environment invariants, in-place RMSD) and §3 (autonomous /
blind accuracy gate). Success ⇔ rank-0 in-place RMSD **`<= 2.0 Å`**.

Do not treat this file as “how to reproduce the published 80/85”. That figure is
withdrawn (oracle ceiling, not docking power).

---

## Quick start

```bash
git clone https://github.com/LeBonhommePharma/FlexAIDdS.git
cd FlexAIDdS
bash scripts/reproduce_astex85.sh
```

Default path is **blind**: `FLEXAIDDS_SEED_ELITISM=0` and `FLEXAIDDS_NATIVE_SEED_FRAC=0`.
`FLEXAIDDS_NATIVE_SEED_FRAC` is a dead knob on today's DatasetRunner path (always
emits `seed_fraction: 0.0`); the live oracle lever is `FLEXAIDDS_SEED_ELITISM=1`
(injects `_INI.pdb`). The default script sets both to off.

Optional `bash scripts/reproduce_astex85.sh --oracle-ceiling` enables native-pose
elitism and prints **`ORACLE CEILING — not docking power`**. Do not cite that arm
as S1.

**Windows users:** see [REPRODUCIBILITY_WINDOWS.md](REPRODUCIBILITY_WINDOWS.md)
for the WSL2 path (`wsl --install` → Ubuntu 22.04 → same script above).

---

## Published results — WITHDRAWN

**This repository publishes no Astex-85 success rate at present.**

The table previously here reported 80/85 (94.1 %). It is **withdrawn** for a specific and
disqualifying reason: it was produced with `FLEXAIDDS_SEED_ELITISM=1` and
`FLEXAIDDS_NATIVE_SEED_FRAC=0.90` (90 % of the genetic-algorithm population notionally
seeded with the crystal pose; on today's DatasetRunner the fraction knob is dead and
elitism injects `_INI.pdb`). That is an **oracle ceiling**, not docking power, and
`METHODOLOGY.md` §0 forbids reporting it as a success rate:

> *NEVER report seed-elitism numbers as docking power.*

`scripts/reproduce_astex85.sh` defaults to blind. It does not compare output against
94.1 % or 80/85. Do not document `SEED_ELITISM=1` / `NATIVE_SEED_FRAC=0.90` as the
current recipe — those assignments appear only on `--oracle-ceiling` and in the
historical block below.

A rate will be republished only from a run that is blind (`native_pose_seeded=0`,
`seed_echo=0`), uses a fixed 85-target denominator, and carries a provenance receipt
pinning engine hash, energy-matrix hash, and input hashes.

---

## Historical provenance of the withdrawn oracle run

The pin fields below describe the **withdrawn** 80/85 (94.1 %) campaign only. They
are not the current reproduce recipe and not a docking-power target.

### Exact git commit

```
8196829f35a2bf065919ccd1508f62f00059895d
```

### Binary fingerprint (withdrawn run)

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
> toolchain+platform, but not cross-platform binary-identical.

### Environment of the withdrawn run (not current)

These variables were set for the withdrawn oracle campaign. **Do not export
`FLEXAIDDS_SEED_ELITISM=1` or `FLEXAIDDS_NATIVE_SEED_FRAC=0.90` on a docking-power
run.**

```bash
FLEXAIDDS_THERMO=1
FLEXAIDDS_T_EFF=0.596
FLEXAIDDS_TENCOM_SCALE=1.0
FLEXAIDDS_RESTARTS=7
FLEXAIDDS_PARALLEL_RESTARTS=1
FLEXAIDDS_EVAL_SCALE_DIHEDRAL=1
FLEXAIDDS_CONSENSUS_SCORER=1
FLEXAIDDS_SEED_ELITISM=1          # oracle — withdrawn; not the default script
FLEXAIDDS_N_ELITE=1
FLEXAIDDS_BUDGET_SCALE=1
FLEXAIDDS_SOFTCORE_WAL=1
FLEXAIDDS_SOFTCORE_FLOOR=0.5
FLEXAIDDS_T_HOT=500
FLEXAIDDS_NATIVE_SEED_FRAC=0.90   # dead knob today; withdrawn recipe only
FLEXAIDDS_RECEPTOR_ROTAMER_PREP=1
```

Benchmark runner invocation (same binary; default script now blinds the two seed knobs):

```bash
benchmark_datasets \
    --benchmark "crossdock_json:benchmarks/datasets/benchmark_astex_native_85.json" \
    --output    <output_dir> \
    --threads   4 \
    --omp-threads 2 \
    --job-timeout-seconds 7200
```

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

## The 0.00 Å per-target table is a seed-echo artefact

The table below is **not a reproduction target**. Many entries are **0.00 Å**
because the withdrawn campaign injected the crystal pose (`SEED_ELITISM=1` /
native-seed elitism). That is a **seed-echo artefact**, not docking power.
`METHODOLOGY.md` §0 forbids reporting seed-elitism / `_INI.pdb` RMSD as the result.

Do not diff a blind run against these numbers. They are retained only as provenance
of the withdrawn oracle campaign (commit `5782045` / v88 tree, last complete run
prior to the v89 thermodynamics fix).

```text
1G9V 0.72  1GM8 0.47  1GPK 1.22  1HNN 7.89  1HP0 0.00  1HQ2 0.00
1IA1 0.85  1IGJ 1.46  1J3J 2.23  1JD0 1.69  1JJE 0.63  1K3U 0.00
1KE5 0.42  1KZK 0.00  1L2S 1.35  1L7F 0.50  1LPZ 0.00  1M2Z 0.92
1MEH 2.42  1MQ6 0.00  1N1M 1.43  1N2J 1.22  1N2V 3.10  1N46 0.46
1NAV 0.00  1OF1 0.00  1OF6 0.68  1OPK 0.00  1OQ5 0.00  1OWE 0.00
1P2Y 2.14  1P62 0.42  1PMN 0.00  1Q1G 0.00  1Q41 0.00  1Q4G 0.00
1R1H 0.54  1R55 0.54  1R58 0.00  1R9O 1.30  1S19 0.76  1S3V 0.00
1SG0 0.70  1SJ0 0.35  1SQ5 1.22  1T40 1.11  1T46 0.44  1T9B 1.71
1TT1 0.60  1TW6 5.01  1TZ8 0.00  1U1C 0.00  1U4D 0.00  1UML 0.00
1UNL 0.00  1UOU 0.54  1V0P 0.85  1V48 0.00  1V4S 0.00  1VCJ 0.44
1W1P 0.00  1W2G 0.65  1X8X 0.59  1XM6 0.60  1XOZ 0.00  1Y6B 0.00
1Y6R 0.00  1YGC 0.00  1YQY 0.00  1YV3 0.80  1YVF 1.09  1YWR 0.60
1Z95 0.60  2BM2 0.53  2BR1 0.62  2BSM 0.63  2BYS 0.77  2C3I 0.00
2CET 1.01  2CGR 0.48  2D3U 0.98  2GBP 0.00  2HB1 1.57  2HR7 0.65
2J62 0.00
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

`scripts/reproduce_astex85.sh` performs all steps end-to-end on the **blind** default:

1. Checks dependencies (cmake, python3, git, curl)
2. Builds `FlexAIDdS` + `benchmark_datasets` from the current git HEAD
3. Generates a portable pair-list JSON (reviewer's absolute paths)
4. Writes `provenance.json` (git SHA, binary SHA256, timestamp, env snapshot, arm label)
5. Runs all 85 targets with `SEED_ELITISM=0` and `NATIVE_SEED_FRAC=0`
6. Prints observed S1 with the RMSD instrument name (`rmsd_hungarian` vs
   `rmsd_to_crystal`) and `N_denominator=85` — not a comparison to the withdrawn 80/85

Expected wall-clock time: **~45–60 min** on Apple M-series with 4 workers.

`--oracle-ceiling` is labelled **not docking power** and must not be used to
republish a rate.
