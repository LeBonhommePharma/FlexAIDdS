# Reproduction Gap Audit — Why 81.2% / 83.5% Is Hard to Reproduce (2026-06-29)

> Investigates v50b → v124 → v126 regression chain. All rates: oracle-ceiling Astex 85, ≤2.0 Å.

## Executive Summary

**The 81.2% figure is the wrong anchor.** Measured oracle self-dock rates:

| Run | Success | Notes |
|-----|---------|-------|
| v50b (`efc4f5d`) | **71/85 (83.5%)** | Pinned binary + matrix MD5 `72d7c739` |
| v121 SMFREE smoke | 61/85 (71.8%) | Post-PSHARE revert |
| v122b | 63/85 (74.1%) | v50b cfg partial restore |
| **v124 full85** | **78/85 (91.8%)** | 85/85 dirs; **not** 50/85 (stale CSV) |
| v126 partial | trending ~58% mid-run | honest selector + `ba5364d3` binary |

**Cannot reproduce 81.2%** if targeting TIER-2 cross-dock while running oracle self-dock JSON — different benchmarks.

## v50b Protocol Snapshot

From `scripts/launch_v50b.py` + `v50b_20260614_consensus5r`:

- `FLEXAIDDS_CONSENSUS_SCORER=1`, `SEED_ELITISM=1`, 5-restart
- SMFREE fitness (overflow → degenerate GA)
- `RECEPTOR_ROTAMER_PREP=0`, `NATIVE_SEED_FRAC=0.90`
- 298 K, **no** `--mode oracle-ceiling` CLI flag
- Matrix MD5 `72d7c739` (HEAD: `9dc93717`)
- Binary SHA `dbfaca09…` @ `efc4f5d`

## Knob Diff vs HEAD (v124/v126)

| Knob | v50b | v124/v126 HEAD |
|------|------|----------------|
| hbond_search in JSON | OFF | ON (v123+) |
| sas_weight | 1.0 (implicit) | 1.0 explicit |
| normalize_area | — | true |
| Temperature | 298 K | 300 K (v124) / 298 K (v126) |
| boltzmann_composite | exp (overflow) | logsumexp (v126) |
| H-bond/VCT patch | absent | `ba5364d3` (v126 binary) |
| consensus-guard | absent | v124+ |

## Ranked Root Causes

1. **Wrong benchmark tier cited** (81.2% ≠ v50b self-dock 83.5%)
2. **Scoring / dock_config drift** (hbond_search, matrix MD5, temperature)
3. **Selector regime change** (logsumexp + consensus-guard vs freq+consensus accident)
4. **PSHARE regression** (v117; fixed v121 but path-dependent)
5. **Stale aggregate CSV** (v124 looked like 50/85 until per-target recount → 78/85)

## Target-level Flips (v50b pass → v124 fail)

Notable: **1M2Z, 1XM6, 1GPK, 1R58, 1T9B** — investigate as Level-2 false-minima / scoring sensitivity targets.

Many v124 "wins" are `ini_elitism` with RMSD≈0 but high `best_cluster_rmsd` (seed echo, not GA basin).

## Minimal v128 Reproduction Recipe

1. `git checkout efc4f5d` → build `build_lto` → stamp `/tmp/FlexAIDdS_v128`
2. Verify matrix MD5 `72d7c739…`
3. Copy `launch_v50b.py` → `launch_v128_v50b_repro.py`
4. Env: identical to v50b (`CONSENSUS_SCORER=1`, no `RECEPTOR_ROTAMER_PREP`, etc.)
5. CLI: `--temperature 298`, no `--mode oracle-ceiling`, timeout 5400s
6. **Gate**: ≥71/85 on per-target `result.csv`
7. Then ablate one knob per run: hbond (v123) → logsumexp (a4056163) → H-bond/VCT (ba5364d3)

## Recommended v128 Experiment

**v128_v50b_repro** — pinned `efc4f5d` binary + matrix, v50b env exactly. This isolates whether HEAD drift alone explains gap vs v50b. Do **not** use as thesis headline number — use as ablation baseline.

## Thesis Integration

- v50b 83.5% under degenerate SMFREE = seed + consensus accident (mechanism now documented)
- v124 91.8% = oracle ceiling with consensus-guard + mixed binary resume
- v126 full-85 = first **honest** selector measurement — the only number that resolves the thesis