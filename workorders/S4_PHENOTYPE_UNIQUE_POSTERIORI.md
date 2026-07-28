# S4 PHENOTYPE_UNIQUE near-miss — a posteriori (FINAL)

**OUT:** `/Users/lp.more/flexaidds_results/s4_pheno_unique_near_miss_20260727_211213`  
**Completed:** 2026-07-28T02:59:45Z (`ALL_ARMS_DONE`) · eval 2026-07-28T03:00:57Z  
**One variable:** `FLEXAIDDS_PHENOTYPE_UNIQUE=1` vs control unset  
**Binary SHA256 (both arms):** `afd5cf42d8cb726de5b92fb66431095360f6161a088945eba110299cd09e4f57`  
**Git tip:** `9569df16` · matrix **9dc9** · R=2 · NO_SEC=1 · seed-echo 0  

## Protocol

| Item | Value |
|------|--------|
| Panel | NEAR_MISS `1N1M`, `1L7F` |
| Control | all FLEXAIDDS_NEW_SEARCH_* unset |
| Treatment | `FLEXAIDDS_PHENOTYPE_UNIQUE=1` only |
| Workers | 2 · OMP=1 |
| Sol #9 | acquired then released after eval |

## Metrics

| arm | code | elect | BCR | elected_cf | wall_s |
|-----|------|------:|----:|-----------:|-------:|
| control | 1L7F | 3.9907 | 3.9907 | −157.729 | 3161 |
| control | 1N1M | 6.3999 | 4.1954 | −99.314 | 2477 |
| pheno_unique | 1L7F | **4.3053** | **3.8964** | −152.151 | 3249 |
| pheno_unique | 1N1M | 6.3999 | **4.1753** | −99.314 | 2610 |

| Δ | 1L7F | 1N1M | mean |
|---|-----:|-----:|-----:|
| dBCR (tx−ctrl) | **−0.094** | **−0.020** | **−0.057** |
| d_elect | **+0.315** | 0.000 | — |

## L4

| arm | n `[NEW-SEARCH-ARCH]` |
|-----|----------------------:|
| control | **0** |
| pheno_unique | **4** |

L4 **PASS** (treatment live, control zero).

## Judgment

| Layer | Result |
|-------|--------|
| L4 | **PASS** |
| Magnitude floor (mean_dBCR ≤ −0.5 or BCR&lt;2 or elect≤2.5) | **FAIL** (mean_dBCR = −0.057) |
| Elect regression | 1L7F elect worsened by +0.315 Å (&gt;0.25) |
| status | **PASS_LIVENESS** |
| `ACCEPT_S4_PHENO` | **False** |

### Interpretation

Phenotype-unique classic mutate **fires and is slightly BCR-helpful** on both targets (small negative dBCR) but far short of the −0.5 Å floor. 1N1M elect remains locked at **6.3999** (same attractor as G4.1/election/G4.3). 1L7F elect **regressed** while BCR improved slightly — ranking still elects a worse pose than control.

Does **not** unlock full-85. Does **not** promote PHENOTYPE_UNIQUE into claim recipe.

## Pins

`validate-pins --out $OUT` → **PINS_OK** (shared binary SHA both arms).  
Sol #9 lock **released** after eval.

## Flip residual

Phase-4 architecture residual still open for optional **B-only** pilot (`FLEXAIDDS_BASIN_REINJECT=1`) as a separate one-var run. Publication path: **P2 oracle + A/B/C binary staging** before any comparative N=85 — **not** full-85 from this pilot.

## Non-claims

- Not genuine top-1 / PoseBusters / claim_ready success.  
- CF/contact-function scoring proxy only.  
- Not permission to enable PHENOTYPE_UNIQUE by default.
