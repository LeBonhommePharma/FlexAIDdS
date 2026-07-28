# Pilot8 dual-zero — Native vs Elected CF inversion map

**Date:** 2026-07-28  
**Panel:** pilot8 `{1G9V,1GPK,1MEH,1P62,1Q4G,1R9O,1T40,2BYS}` from arm A  
`comparative_pilot8_20260728`  
**Status:** **COMPLETE — 8/8 classified SEARCH-MISS (sampling_ceiling)**  
**Does not claim docking success rates; does not authorize full-85.**

## Purpose

After P3 **SCIENCE HOLD** (BCR=0/8 and S_top10=0/8 on reconstruction arm A), answer per target:

> Under a **fixed production LOCCLF** surface, is the **crystal native** better, worse, or tied vs the **elected** pose on the CF proxy?

That splits **SCORING-LOCKED** (landscape prefers decoy) from **SEARCH-MISS** (native is CF-better but was not found/elected).

## Method

| Item | Value |
|------|--------|
| One variable | Pose role: crystal ligand SDF vs pilot elected PDB |
| Fixed surface | Production dock_config template (S4 leaf axes: VCT r0=7, sas_weight=0.4, hbond on) |
| Scorer | `probe_cf --mode direct` + FlexAIDdS C binary + data dir matrix/defs |
| Env | `FLEXAIDDS_WAL_COERCIVE=0`; no BOOM/PB_CLASH coercive |
| ε | 0.5 CF units |
| dCF | CF_native − CF_elected (negative ⇒ native better) |
| Script | `scripts/native_elected_cf_inversion_map.py --pilot8` |

**OUT:** `~/flexaidds_results/campaigns/three_engine/analysis/pilot8_inversion_map_20260728/`  
(`inversion_map.json`, `inversion_map.csv`, `inversion_map.md`)

## Counts

| Class | N |
|-------|--:|
| SCORING-LOCKED | **0** |
| SEARCH-MISS | **8** |
| TIED | **0** |
| of MISS: sampling_ceiling (BCR>2) | **8** |
| of MISS: election_pool (BCR≤2) | **0** |

## Per-target (modern LOCCLF re-score)

| PDB | class | subclass | CF native | CF elected (probe) | dCF(n−e) | elect RMSD | BCR | GA REMARK CF |
|-----|-------|----------|----------:|-------------------:|---------:|-----------:|----:|-------------:|
| 1G9V | SEARCH-MISS | sampling_ceiling | −49.6 | ~1.80e5 | ≪0 | 11.03 | 5.95 | (see receipts) |
| 1GPK | SEARCH-MISS | sampling_ceiling | −107.3 | ~1.59e5 | ≪0 | 8.60 | 4.46 | |
| 1MEH | SEARCH-MISS | sampling_ceiling | −79.8 | ~2.30e5 | ≪0 | 13.28 | 6.00 | |
| 1P62 | SEARCH-MISS | sampling_ceiling | −144.7 | ~8.87e4 | ≪0 | 4.90 | 4.69 | **−1236.7** |
| 1Q4G | SEARCH-MISS | sampling_ceiling | −129.3 | ~1.48e5 | ≪0 | 7.57 | 4.43 | |
| 1R9O | SEARCH-MISS | sampling_ceiling | −117.1 | ~1.18e5 | ≪0 | 7.60 | 5.91 | |
| 1T40 | SEARCH-MISS | sampling_ceiling | −212.4 | ~2.59e5 | ≪0 | 7.71 | 5.13 | |
| 2BYS | SEARCH-MISS | sampling_ceiling | −157.4 | ~9.09e4 | ≪0 | 10.79 | 4.25 | |

## Caveats (mandatory)

1. **Arm A poses are reconstruction FlexAID**, not historical JCIM SHA. Scores here are **modern FlexAIDdS LOCCLF** on those poses.
2. **Cross-engine re-score pathology:** GA REMARK CF on elected 1P62 is **−1236** (favourable under A binary), but `probe_cf` re-score yields huge **CF.wal** (~8e4). Absolute CF_elected is **not** comparable to GA REMARK; ordinal ranking (native ≪ elected under modern LOCCLF) still holds.
3. **BCR always >2** on this panel → subclass is **sampling_ceiling**, not election-pool (near-native was not in the scanned pool).
4. **No claim** of genuine/BCR rates; **full85 still forbidden** while P3=hold.

## Interpretation

- **0 SCORING-LOCKED** under modern LOCCLF: the dual-zero pilot is **not** explained by “crystal is CF-worse than the elected decoy” on this surface.
- **8/8 SEARCH-MISS + sampling_ceiling:** crystal is CF-better; GA/election never produced a ≤2 Å pose. That matches P3 hold as **sampling/prep ceiling**, not “landscape hates native.”
- The wall blow-up on re-score of complex PDBs flags a **pose I/O / typing** risk when bridging reconstruction FlexAID outputs into modern `probe_cf` — treat absolute elected CF as diagnostic only.

## Next lever (guided, not launched here)

1. **Primary:** sampling/prep for native-like placement (coarse orient / cleft / budget) **only if** modern C-arm or leaf-config docks are the claim path — not BOOM thrash.  
2. **Do not** full-85 until P3 leaves hold with multi-arm interpretable success.  
3. Optional: score-only map on **arm C** pilot poses when available (same engine as probe_cf) for apples-to-apples CF.  
4. Prep/cleft audit remains open if sampling still fails after geometry-faithful re-score.

## Related

- `P3_SCIENCE_HOLD.md` (dual-zero N=8)  
- `scripts/native_elected_cf_inversion_map.py`  
- Prior panel map: `workorders/INVERSION_MAP.md` (different codes)
