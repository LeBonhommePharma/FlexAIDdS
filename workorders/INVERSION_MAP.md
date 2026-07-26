# Native–Elected CF inversion map — PASS

**Written:** 2026-07-26T01:43:56.472558+00:00  
**Pilot poses:** `/Users/lp.more/flexaidds_results/pilot_w1_boom_interval_20260725_134740`  
**One variable:** pose role under fixed production LOCCLF (native vs elected)  
**ε:** 0.5 CF units  
**OUT:** `/Users/lp.more/flexaidds_results/workorders/inversion_map_20260725_213932`  

**Verdict:** **PASS** — classified 8/8; LOCKED=3 MISS=5 TIED=0

**Next lever (guided):** SEARCH-MISS dominant → FLEXAIDDS_COARSE_ORIENTATIONS=256 W1 pilot (matrix 9dc9)

## Counts

| Class | N |
|-------|--:|
| SCORING-LOCKED | 3 |
| SEARCH-MISS | 5 |
| TIED | 0 |
| of MISS: sampling_ceiling (BCR>2) | 5 |
| of MISS: election_pool (BCR≤2) | 0 |

## Per-target

| PDB | class | subclass | CF native | CF elected | dCF(n−e) | elect RMSD | BCR |
|-----|-------|----------|----------:|-----------:|---------:|-----------:|----:|
| 1J3J | SEARCH-MISS | sampling_ceiling | -51.244 | -48.892 | -2.352 | 62.2229 | 22.9559 |
| 1K3U | SEARCH-MISS | sampling_ceiling | -138.811 | -32.443 | -106.369 | 11.4681 | 11.7827 |
| 1L7F | SEARCH-MISS | sampling_ceiling | -118.989 | -42.490 | -76.499 | 3.9227 | 3.9633 |
| 1N1M | SEARCH-MISS | sampling_ceiling | -57.585 | -37.111 | -20.474 | 5.6602 | 4.0427 |
| 1M2Z | SEARCH-MISS | sampling_ceiling | -117.745 | 118.466 | -236.211 | 13.7872 | 13.0393 |
| 1OQ5 | SCORING-LOCKED | false_min_attractor | -107.870 | -125.757 | +17.886 | 3.9457 | 1.6481 |
| 1SQ5 | SCORING-LOCKED | false_min_attractor | -135.211 | -164.592 | +29.381 | 5.0963 | 1.6467 |
| 1YGC | SCORING-LOCKED | false_min_attractor | -146.591 | -216.077 | +69.486 | 1.7516 | 1.0103 |

## Cadence

- Phase: 1.5 diagnosis (score-only) / roadmap Phase 1 extension  
- Genuine/BCR: **not claimed** this step  
- **PASS**  
- No full-85; no memetic; no WAL re-panel  

Methodology: `NEXT_CAMPAIGN_STEP.md` · BENCHMARKING_ROADMAP Phase 1 then Phase 4 sampling.

## Regime interpretation (campaign-critical)

| Bucket | Targets | Meaning |
|--------|---------|---------|
| **SEARCH-MISS** (5) | 1J3J 1K3U 1L7F 1N1M 1M2Z | Crystal scores **better** CF than elected; GA never found a competitive near-native basin. **Sampling problem.** |
| **SCORING-LOCKED** (3) | 1OQ5 1SQ5 1YGC | Elected scores **better** CF than crystal SDF while BCR is ≤2 Å (near-native exists in pool). **Scoring prefers non-crystal / decoy attractor over crystal coordinates.** |

This **splits** the E10 “sampling primary” statement: clean probes are sampling-limited; **election-gap near-misses are scoring-limited** (native not CF-min even when near-native poses exist).

### Implications for Phase 4
1. **Do not** run BOOM/coarse on SCORING-LOCKED codes expecting genuine wins without a scoring change.  
2. **Do** run sampling levers (**COARSE_ORIENTATIONS**, niche, etc.) on **SEARCH-MISS** codes only.  
3. pb_clash / burial work targets the SCORING-LOCKED class (strong decoys), not the clean-probe SEARCH-MISS class.

### Liveness / tool notes
- ops/gates configs missing for 1OQ5/1SQ5/1YGC — used pilot `dock_config.json` (DatasetRunner production emit).  
- CF gate smoke: n=5 scored, 0 skip (`cf_gate_probe_cf.sh` --config+--ligand) → Phase 0 tool PASS.  
- 1M2Z native CF = −117.745 (matches methodology production LOCCLF).

### Provenance
- Pilot poses: `/Users/lp.more/flexaidds_results/pilot_w1_boom_interval_20260725_134740`  
- ε = 0.5  
- Live OUT: `/Users/lp.more/flexaidds_results/workorders/inversion_map_20260725_213932`  
- Script: `scripts/native_elected_cf_inversion_map.py`  

**Blocks:** full-85, memetic, WAL re-panel, interval-only BOOM still hold.
