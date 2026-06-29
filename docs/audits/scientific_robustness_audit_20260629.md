# Scientific Robustness Audit — v124–v127 (2026-06-29)

> Cross-validated against repo files, benchmark results, and published docking literature.
> Oracle-ceiling Astex 85 self-dock; ≤2.0 Å success criterion unless noted.

## Executive Summary

1. **Protocol**: v124–v127 use `ORACLE_CEILING` (crystal IC, seed elitism, oracle site) — a **ceiling** measurement, not blind/autonomous docking. Do not compare to blind cross-dock literature (~20–35%) without explicit labeling.

2. **v124 final (verified)**: **78/85 (91.8%)** from per-target `result.csv` (85/85 dirs). Aggregate `astex_diverse_results.csv` is **stale** (59 rows) — do not cite for thesis.

3. **v50b reference (same oracle self-dock family)**: **71/85 (83.5%)** per `~/flexaidds_results/v50b_20260614_consensus5r/` at commit `efc4f5d`. The **81.2%** figure in `BENCHMARK_STANDARD.md` refers to a **different tier** (TIER-2 oracle cross-dock) — not this launch path.

4. **v126 (Option B)**: logsumexp-stable `boltzmann_composite` is **numerically correct**; log-space ranking preserves composite order. First run where freq-gated selector is not masked by overflow + consensus fallback.

5. **Consensus**: Code default `FLEXAIDDS_CONSENSUS_SCORER=0` since `ce8f3368`; v126/v127 launch scripts set `=1` globally. Runtime verified: 46/46 completed v126 targets had `[CONSENSUS]` in `stderr.log` (2026-06-29).

6. **Thermodynamic claims**: `predicted_dG`, Shannon fields are **CF-derived reporting proxies** — not calibrated ΔG (see `docs/dev/thermo_invariants.md`, `AGENTS.md`).

7. **v127 / `ba5364d3`**: H-bond cosine gate + C.ar stacking changes scoring kernel — **validation run**, not neutral packaging; Astex overfitting risk if interpreted as confirmation.

8. **Thesis framing**: v50b's high oracle rate under SMFREE overflow is a **documented degenerate-accident** (random GA + seed + consensus). v126 measures the **honest free-energy selector** regime. Either outcome vs v50b/v124 band is publishable.

## Literature Benchmarks (orientation)

| Source | Setting | Typical success |
|--------|---------|-----------------|
| Hartshorn et al. 2007 (Astex/GOLD) | Native redock, unbiased | ~80% order-of-magnitude |
| Blind cross-dock benchmarks | No crystal seed | ~20–35% |
| Oracle-ceiling self-dock (this repo) | 90% native seed + elitism | **83–92%** band (v50b/v124) |

## Protocol Audit (v124–v127)

| Version | Key change | Scientific role |
|---------|------------|-----------------|
| v124 | Consensus-guard for INI | CF heuristic; protects seed from worse cluster |
| v125 | Consensus default OFF | Diagnostic: overflow breaks selector |
| v126 | logsumexp composite | Principled relative free-energy ranking |
| v127 | HEAD + H-bond/VCT | Scoring-kernel ablation baseline (not selector-only) |

## Risk Register

| Priority | Risk |
|----------|------|
| P0 | Citing stale aggregate CSV instead of 85 `result.csv` files |
| P0 | Conflating TIER-2 81.2% with oracle self-dock 83.5%/91.8% |
| P0 | Overclaiming `predicted_dG` as true thermodynamics |
| P1 | Mixed-binary v124 resume (59 @ old SHA + 26 @ HEAD) |
| P1 | `ba5364d3` scoring interpreted as validation without ablation |
| P2 | Many wins are `ini_elitism` RMSD≈0 — disclose seed-echo fraction |

## Recommendations

- **v126 full-85** is the control experiment for the honest selector thesis.
- Regenerate v124 `astex_diverse_results.csv` from 85 per-target files before publication tables.
- Run **hbond-off ablation** (v127b @ `a4056163`) if v126 < 78%.
- Run **TIER-1** (`NATIVE_SEED_FRAC=0`) before claiming autonomous docking performance.

## References (repo)

- `benchmarks/BENCHMARK_STANDARD.md`, `docs/dev/thermo_invariants.md`, `AGENTS.md`
- `LIB/DatasetRunner.cpp` (~833 logsumexp, ~6454 consensus, ~6586 guard)
- `scripts/launch_v126_smoke.py`, `scripts/launch_v127_full85.py` (env_snapshot)