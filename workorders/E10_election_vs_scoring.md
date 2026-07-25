# E10 — election vs scoring / sampling (STEP 1)

**Date:** 2026-07-25  
**Phase:** W0.1 / E10 offline (no re-dock)  
**Script:** `scripts/e10_election_vs_scoring.py`  
**Campaign:** `~/flexaidds_results/v_autonomous_20260724_160919` (named rsync from archive batch 20260725T095624Z; materialize **in progress** during first cut)

## One variable

None (offline diagnostic on frozen heads).

## Results (first cut — partial materialize)

| Metric | Value |
|--------|------:|
| Targets analyzed | **8** (partial local tree; full 79–85 still rsyncing) |
| Election-gap (BCR≤2.5 Å, elected>2.0, seed_echo=0) | **1 / 8** |
| Size-bias suspects (soft_β_G ≪ CF, high freq) | **8 / 8** |

Machine outputs: see implementer scratch `e10/e10.json`, `e10/e10.csv`, `e10/e10.md`.

### Per-target (partial)

| PDB | rmsd | BCR | gap? | size_bias? |
|-----|-----:|----:|:----:|:----------:|
| 1G9V | 4.50 | 2.06 | Y | Y |
| 1GM8 | 3.58 | 3.41 | n | Y |
| 1GPK | 3.21 | 3.45 | n | Y |
| 1HNN | 1.58 | 1.42 | n | Y |
| 1HP0 | 3.84 | 3.80 | n | Y |
| 1HQ2 | 1.51 | 1.94 | n | Y |
| 1IA1 | 2.65 | 2.76 | n | Y |
| 1IGJ | 74.08 | 29.67 | n | Y |

## Interpretation (gate)

- **Independent CF-better-than-elected near-native heads:** rare in this slice (at most 1 head note on 1IA1); **election_gap fraction = 1/8 = 12.5%** among analyzed — **not** a large election-dominated failure mode.
- **Size-bias_suspect = 8/8** confirms this frozen run used **multiplicity-sensitive ACF** (`soft_beta_G` hundreds of units below pose CF with high `freq`) — consistent with OPS note that the 25.3% baseline **predates** measured `free_energy_strict` product default. **Do not** cite 25.3% as proof election fix worked.
- **Sampling still dominant** when BCR≫2 (1IGJ, 1GM8, 1GPK, 1HP0).

## ACCEPT vs methodology STEP 1

| Criterion | Result |
|-----------|--------|
| E10 script on science branch | PASS |
| Run on local baseline leaf (no CloudDocs find) | PASS (partial N=8 while rsync continues) |
| Fraction independent scorer prefers near-native over elected | **Small** (election not primary wall) |
| STOP before sampling? | **NO** — proceed to STEP 2 wall oracle; re-run E10 when full materialize completes |

## Reporting cadence

- **Phase:** STEP 1 E10  
- **One variable:** n/a (offline)  
- **Genuine / BCR / gap:** deferred to full materialize; partial election_gap **1/8**  
- **PASS/FAIL:** **PASS** to continue (small election fraction; size-bias documents legacy ACF)
