# R2 PhD 2016 — Exhaustive Original Results Catalog

**Source (read-only):** `/Users/lp.more/Documents/PhD/Seminaire/R2_PhD_2016/Files`  
**Generated:** 2026-07-15  
**Machine JSON:** [`docs/implementation/data/R2_PhD_2016_results_catalog.json`](data/R2_PhD_2016_results_catalog.json)  
**Slim job-success JSON:** [`docs/figures/3dsig_2017/R2_PhD_2016_job_success.json`](../figures/3dsig_2017/R2_PhD_2016_job_success.json)

Original R2 files were **not modified**. Large `training_st0r5.2.surfaces` (~886 MB) and `R_packages/summary/` (~1400 RMSD lists) were inventoried shallowly (no content deep-walk).

---

## 1. Metric definitions

| Metric | Definition |
|--------|------------|
| **S_topK (std)** | Fraction of cases where **any** pose with rank ∈ [0, K) has RMSD ≤ 2.0 Å |
| **S_oracle** | Lowest RMSD among all poses/ranks for the case ≤ 2.0 Å |
| **Bootstrap** | `ASTEX_bootstrap.out` (`concatenate_results.pl`): resamples of **oracle** success over 84 cases |
| **Pooled** | All minPts + repeats merged per case before scoring |
| **Row-median** | Median of per-(minPts, repeat) success rates |

Astex Diverse expected **N = 84**. Success threshold **RMSD ≤ 2.0 Å**.

---

## 2. Primary table — job → success

Source preference: `ASTEX_final` pooled → `ALL_astex` → `FlexAID_CF+FO` → `BEST`.

| Job | Algo | TEMP | CLUS | minPts | Notes | Source | n | S_top1 | S_top10 | S_oracle | Boot median | Boot 95% CI |
|----:|:----:|-----:|-----:|:-------|:------|:-------|--:|-------:|--------:|---------:|------------:|:------------|
| 2180 | FO | 13 | 2.5 | 15 | — | ALL_astex | 77 | 0.4545 | 0.5065 | 0.5065 | n/a | — |
| 3771 | FO | 13 | 2.0 | 15 | — | ALL_astex | 80 | 0.3250 | 0.3875 | 0.3875 | n/a | — |
| 8572 | CF | 13 | 2.0 | REF | — | FlexAID_CF+FO | 84 | 0.0119 | 0.2024 | 0.2024 | n/a | — |
| 9156 | FO | 13 | 2.0 | — | — | ALL_astex | 26 | 0.4615 | 0.5769 | 0.5769 | n/a | — |
| 9763 | CF | 0 | 2.0 | REF | — | FlexAID_CF+FO | 84 | 0.1190 | 0.3214 | 0.3214 | n/a | — |
| 18852 | FO | 13 | 2.0 | 10,15,22 | 10 repeats | ASTEX_final | 83 | 0.8313 | 0.8675 | 0.8675 | 0.8929 | [0.821, 0.952] |
| 21046 | FO | 13 | 2.5 | 12 | CF < 0 filter, 10 repeats | ASTEX_final | 84 | 0.8571 | 0.8571 | 0.8571 | 0.8571 | [0.774, 0.929] |
| 21822 | CF | 0 | 2.0 | REF | 10 repeats | ASTEX_final | 84 | 0.6548 | 0.9643 | 0.9643 | 0.9643 | [0.917, 1.000] |
| 28081 | CF | 13 | 2.0 | REF | MC_st0r5.2 matrix, 10 repeats | ASTEX_final | 84 | 0.5595 | 0.9643 | 0.9643 | 0.9643 | [0.923, 1.000] |
| 29332 | CF | 13 | 2.0 | REF | NRG_mat_BEST_13 matrix, 10 repeats | ASTEX_final | 84 | 0.3452 | 0.7857 | 0.7857 | 0.7857 | [0.697, 0.857] |
| 30281 | FO | 13 | 2.5 | — | — | FlexAID_CF+FO | 34 | 0.5882 | 0.6471 | 0.6471 | n/a | — |
| 31578 | FO | 13 | 3.0 | — | — | ALL_astex | 24 | 0.4167 | 0.5417 | 0.5417 | n/a | — |
| 31955 | FO | 13 | 2.0 | 12 | CF < 0 filter, 10 repeats | ASTEX_final | 84 | 0.8333 | 0.8690 | 0.8690 | 0.8690 | [0.792, 0.929] |

### Compact handoff table (job → median success)

| Job | Algo | S_top1 | S_top10 | S_oracle | Bootstrap median |
|----:|:----:|-------:|--------:|---------:|-----------------:|
| 2180 | FO | 0.4545 | 0.5065 | 0.5065 | n/a |
| 3771 | FO | 0.3250 | 0.3875 | 0.3875 | n/a |
| 8572 | CF | 0.0119 | 0.2024 | 0.2024 | n/a |
| 9156 | FO | 0.4615 | 0.5769 | 0.5769 | n/a |
| 9763 | CF | 0.1190 | 0.3214 | 0.3214 | n/a |
| 18852 | FO | 0.8313 | 0.8675 | 0.8675 | 0.8929 |
| 21046 | FO | 0.8571 | 0.8571 | 0.8571 | 0.8571 |
| 21822 | CF | 0.6548 | 0.9643 | 0.9643 | 0.9643 |
| 28081 | CF | 0.5595 | 0.9643 | 0.9643 | 0.9643 |
| 29332 | CF | 0.3452 | 0.7857 | 0.7857 | 0.7857 |
| 30281 | FO | 0.5882 | 0.6471 | 0.6471 | n/a |
| 31578 | FO | 0.4167 | 0.5417 | 0.5417 | n/a |
| 31955 | FO | 0.8333 | 0.8690 | 0.8690 | 0.8690 |

---

## 3. Full `jobsID.info` mapping

| Job | Algo | TEMP | CLUS RMSD | minPoints | nRepeats | Matrix / notes | Description |
|----:|:----:|-----:|----------:|:----------|---------:|:---------------|:------------|
| 2180 | FO | 13 | 2.5 | 15 | — | — | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.5, TEMP = 13, minPoints init = 15 |
| 3771 | FO | 13 | 2.0 | 15 | — | — | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.0, TEMP = 13, minPoints init = 15 |
| 8572 | CF | 13 | 2.0 | REF | — | — | Astex Diverse set with CF cluster, CLUS RMSD = 2.0, TEMP = 13 |
| 9156 | FO | 13 | 2.0 | None | — | — | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.0, TEMP = 13 |
| 9763 | CF | 0 | 2.0 | REF | — | — | Astex Diverse set with CF cluster, CLUS RMSD = 2.0, TEMP = 0 |
| 18852 | FO | 13 | 2.0 | 10,15,22 | 10 | 10 repeats | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.0, TEMP = 13, minPoints = 10,15,22, nRepeats = 10 |
| 21046 | FO | 13 | 2.5 | 12 | 10 | CF < 0 filter, 10 repeats | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.5, TEMP = 13, minPoints = 12, nRepeats = 10; # CF < 0 |
| 21822 | CF | 0 | 2.0 | REF | 10 | 10 repeats | Astex Diverse set with CF cluster, CLUS RMSD = 2.0, TEMP = 0, nRepeats = 10 |
| 28081 | CF | 13 | 2.0 | REF | 10 | MC_st0r5.2 matrix, 10 repeats | Astex Diverse set with CF cluster, CLUS RMSD = 2.0, TEMP = 13; MC_st0r5.2 (10 rep) |
| 29332 | CF | 13 | 2.0 | REF | 10 | NRG_mat_BEST_13 matrix, 10 repeats | Astex Diverse set with CF cluster, CLUS RMSD = 2.0, TEMP = 13; NRG_mat_BEST_13 (10 rep) |
| 30281 | FO | 13 | 2.5 | None | — | — | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.5, TEMP = 13 |
| 31578 | FO | 13 | 3.0 | None | — | — | Astex Diverse set with FastOPTICS, CLUS RMSD = 3.0, TEMP = 13 |
| 31955 | FO | 13 | 2.0 | 12 | 10 | CF < 0 filter, 10 repeats | Astex Diverse set with FastOPTICS, CLUS RMSD = 2.0, TEMP = 13, minPoints = 12, nRepeats = 10; # CF < 0 |

**FO vs CF / TEMP / matrix**

- **FastOPTICS (FO):** 3771, 2180, 9156, 30281, 31578, 18852, 21046, 31955
- **CF cluster (CF):** 8572, 9763, 21822, 28081, 29332
- **TEMP=0:** 9763, 21822 — **TEMP=13:** remaining jobsID entries
- **Matrices:** 28081 = `MC_st0r5.2`; 29332 = `NRG_mat_BEST_13`; 21046/31955 = CF < 0 filter

---

## 4. `ASTEX_bootstrap.out` — all columns

- Planned iterations (script): 10000
- Raw data rows: **463** (incomplete file; last row truncated)
- Rows used (full 6 columns): **462**
- Header prints each job ID twice; statistics use unique first 6 columns.
- Values = **oracle success rate** on bootstrap resamples of 84 cases.

| Job | n | Median | Mean | Std | P2.5 | P25 | P75 | P97.5 | Min | Max |
|----:|--:|-------:|-----:|----:|-----:|----:|----:|------:|----:|----:|
| 18852 | 462 | 0.892857 | 0.891105 | 0.035050 | 0.8214 | 0.8690 | 0.9167 | 0.9524 | 0.7500 | 0.9643 |
| 21046 | 462 | 0.857143 | 0.855674 | 0.040274 | 0.7738 | 0.8333 | 0.8810 | 0.9286 | 0.7143 | 0.9643 |
| 21822 | 462 | 0.964286 | 0.963848 | 0.019762 | 0.9167 | 0.9524 | 0.9762 | 1.0000 | 0.8929 | 1.0000 |
| 28081 | 462 | 0.964286 | 0.964415 | 0.019798 | 0.9229 | 0.9524 | 0.9762 | 1.0000 | 0.8810 | 1.0000 |
| 29332 | 462 | 0.785714 | 0.783215 | 0.043498 | 0.6967 | 0.7500 | 0.8095 | 0.8571 | 0.6548 | 0.9048 |
| 31955 | 462 | 0.869048 | 0.866342 | 0.037594 | 0.7920 | 0.8452 | 0.8929 | 0.9286 | 0.7262 | 0.9524 |

**Bootstrap ranking (median oracle):** 28081 ≈ 21822 (0.9643) ≫ 18852 (0.8929) > 31955 (0.8690) > 21046 (0.8571) ≫ 29332 (0.7857).

---

## 5. Success rates from compiled pose tables

### 5.1 `ASTEX_final.dat` (primary multi-repeat set)

- Poses: **54091**; jobs: 18852, 21046, 21822, 28081, 29332, 31955

#### Per (job, minPts) — median over repeats

| Job | minPts | n_rows | n_cases | S_top1 | S_top10 | S_oracle |
|----:|:-------|-------:|--------:|-------:|--------:|---------:|
| 18852 | 10 | 10 | 68 | 0.4511 | 0.5946 | 0.5946 |
| 18852 | 15 | 10 | 68 | 0.4786 | 0.5846 | 0.5846 |
| 18852 | 22 | 10 | 68 | 0.4063 | 0.5901 | 0.5901 |
| 21046 | 12 | 9 | 80 | 0.5600 | 0.6267 | 0.6267 |
| 21822 | REF | 10 | 84 | 0.2679 | 0.6845 | 0.6845 |
| 28081 | REF | 10 | 84 | 0.1786 | 0.7083 | 0.7083 |
| 29332 | REF | 10 | 84 | 0.1250 | 0.5714 | 0.5714 |
| 31955 | 12 | 10 | 83 | 0.5223 | 0.5891 | 0.5891 |

#### Per job — pooled (all minPts/repeats)

| Job | n_cases | S_top1 | S_top3 | S_top5 | S_top10 | S_oracle | median best RMSD |
|----:|--------:|-------:|-------:|-------:|--------:|---------:|-----------------:|
| 18852 | 83 | 0.8313 | 0.8554 | 0.8675 | 0.8675 | 0.8675 | 0.634 |
| 21046 | 84 | 0.8571 | 0.8571 | 0.8571 | 0.8571 | 0.8571 | 0.666 |
| 21822 | 84 | 0.6548 | 0.8571 | 0.9048 | 0.9643 | 0.9643 | 0.632 |
| 28081 | 84 | 0.5595 | 0.9048 | 0.9286 | 0.9643 | 0.9643 | 0.646 |
| 29332 | 84 | 0.3452 | 0.6905 | 0.7619 | 0.7857 | 0.7857 | 0.814 |
| 31955 | 84 | 0.8333 | 0.8571 | 0.8571 | 0.8690 | 0.8690 | 0.574 |

### 5.2 `ALL_astex.dat`

- Poses: **63100**; jobs: 18852, 21046, 2180, 21822, 28081, 29332, 31578, 31955, 3771, 9156

| Job | minPts | n_rows | n_cases | S_top1 | S_top10 | S_oracle |
|----:|:-------|-------:|--------:|-------:|--------:|---------:|
| 2180 | 15 | 1 | 77 | 0.2987 | 0.4156 | 0.4156 |
| 2180 | 22 | 1 | 77 | 0.3506 | 0.4416 | 0.4416 |
| 2180 | 33 | 1 | 77 | 0.2597 | 0.4156 | 0.4156 |
| 2180 | 49 | 1 | 77 | 0.3117 | 0.4416 | 0.4416 |
| 2180 | 73 | 1 | 77 | 0.2727 | 0.4156 | 0.4156 |
| 3771 | 15 | 1 | 80 | 0.2500 | 0.3375 | 0.3375 |
| 3771 | 22 | 1 | 80 | 0.2625 | 0.3000 | 0.3000 |
| 3771 | 33 | 1 | 80 | 0.2000 | 0.3375 | 0.3375 |
| 3771 | 49 | 1 | 80 | 0.1875 | 0.3375 | 0.3375 |
| 3771 | 73 | 1 | 79 | 0.2025 | 0.3165 | 0.3165 |
| 9156 | 25 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 9156 | 28 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 9156 | 31 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 9156 | 32 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 36 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 37 | 1 | 3 | 0.6667 | 0.6667 | 0.6667 |
| 9156 | 40 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 9156 | 41 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 9156 | 42 | 1 | 6 | 0.6667 | 0.8333 | 0.8333 |
| 9156 | 46 | 1 | 3 | 0.3333 | 1.0000 | 1.0000 |
| 9156 | 48 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 54 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 9156 | 55 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 9156 | 60 | 1 | 3 | 0.0000 | 0.3333 | 0.3333 |
| 9156 | 61 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 9156 | 63 | 1 | 6 | 0.5000 | 0.8333 | 0.8333 |
| 9156 | 65 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 69 | 1 | 3 | 0.6667 | 1.0000 | 1.0000 |
| 9156 | 72 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 74 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 9156 | 79 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 81 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 9156 | 82 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 9156 | 83 | 1 | 1 | 0.0000 | 1.0000 | 1.0000 |
| 9156 | 88 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 90 | 1 | 3 | 0.0000 | 0.3333 | 0.3333 |
| 9156 | 91 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 9156 | 94 | 1 | 6 | 0.5000 | 0.6667 | 0.6667 |
| 9156 | 96 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 9156 | 97 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 18852 | 10 | 10 | 68 | 0.4511 | 0.5946 | 0.5946 |
| 18852 | 15 | 10 | 68 | 0.4786 | 0.5846 | 0.5846 |
| 18852 | 22 | 10 | 68 | 0.4063 | 0.5901 | 0.5901 |
| 21046 | 12 | 9 | 80 | 0.5600 | 0.6267 | 0.6267 |
| 21822 | REF | 10 | 84 | 0.2679 | 0.6845 | 0.6845 |
| 28081 | REF | 10 | 84 | 0.1786 | 0.7083 | 0.7083 |
| 29332 | REF | 10 | 84 | 0.1250 | 0.5714 | 0.5714 |
| 31578 | 22 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 31578 | 33 | 1 | 2 | 0.5000 | 1.0000 | 1.0000 |
| 31578 | 35 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 31578 | 36 | 1 | 2 | 1.0000 | 1.0000 | 1.0000 |
| 31578 | 38 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 42 | 1 | 3 | 0.6667 | 0.6667 | 0.6667 |
| 31578 | 43 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 31578 | 45 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 31578 | 49 | 1 | 2 | 0.5000 | 1.0000 | 1.0000 |
| 31578 | 50 | 1 | 2 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 52 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 31578 | 53 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 31578 | 54 | 1 | 4 | 0.7500 | 0.7500 | 0.7500 |
| 31578 | 57 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 62 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 63 | 1 | 3 | 0.6667 | 0.6667 | 0.6667 |
| 31578 | 64 | 1 | 2 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 67 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 31578 | 72 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 |
| 31578 | 73 | 1 | 2 | 0.5000 | 1.0000 | 1.0000 |
| 31578 | 75 | 1 | 2 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 78 | 1 | 2 | 0.5000 | 0.5000 | 0.5000 |
| 31578 | 79 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 81 | 1 | 4 | 0.2500 | 0.7500 | 0.7500 |
| 31578 | 82 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 85 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 89 | 1 | 1 | 0.0000 | 1.0000 | 1.0000 |
| 31578 | 93 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31578 | 94 | 1 | 3 | 0.6667 | 0.6667 | 0.6667 |
| 31578 | 96 | 1 | 2 | 0.0000 | 0.5000 | 0.5000 |
| 31578 | 97 | 1 | 1 | 0.0000 | 0.0000 | 0.0000 |
| 31955 | 12 | 10 | 83 | 0.5223 | 0.5891 | 0.5891 |

#### Pooled

| Job | n_cases | S_top1 | S_top10 | S_oracle |
|----:|--------:|-------:|--------:|---------:|
| 2180 | 77 | 0.4545 | 0.5065 | 0.5065 |
| 3771 | 80 | 0.3250 | 0.3875 | 0.3875 |
| 9156 | 26 | 0.4615 | 0.5769 | 0.5769 |
| 18852 | 83 | 0.8313 | 0.8675 | 0.8675 |
| 21046 | 84 | 0.8571 | 0.8571 | 0.8571 |
| 21822 | 84 | 0.6548 | 0.9643 | 0.9643 |
| 28081 | 84 | 0.5595 | 0.9643 | 0.9643 |
| 29332 | 84 | 0.3452 | 0.7857 | 0.7857 |
| 31578 | 24 | 0.4167 | 0.5417 | 0.5417 |
| 31955 | 84 | 0.8333 | 0.8690 | 0.8690 |

### 5.3 `FlexAID_CF+FO_results.dat`

- Poses: **19889**; jobs: 2180, 21822, 30281, 31578, 3771, 8572, 9156, 9763

| Job | n_cases | S_top1 | S_top10 | S_oracle |
|----:|--------:|-------:|--------:|---------:|
| 2180 | 77 | 0.4545 | 0.5065 | 0.5065 |
| 3771 | 80 | 0.3250 | 0.3875 | 0.3875 |
| 8572 | 84 | 0.0119 | 0.2024 | 0.2024 |
| 9156 | 26 | 0.4615 | 0.5769 | 0.5769 |
| 9763 | 84 | 0.1190 | 0.3214 | 0.3214 |
| 21822 | 84 | 0.6548 | 0.9643 | 0.9643 |
| 30281 | 34 | 0.5882 | 0.6471 | 0.6471 |
| 31578 | 24 | 0.4167 | 0.5417 | 0.5417 |

### 5.4 `BEST_result_per_case*.dat`

#### `BEST_result_per_case.dat` (feeds bootstrap oracle)

Invalid placeholder rows (empty rank/minPts, RMSD 0.0) excluded.

\* perl-style: best-RMSD pose has rank < K and RMSD ≤ 2.0 Å.

| Job | n | S_oracle | S_top1* | S_top10* | minPts | median best RMSD |
|----:|--:|---------:|--------:|---------:|:-------|-----------------:|
| 18852 | 67 | 0.8657 | 0.0746 | 0.8657 | 10,15,22 | 0.765 |
| 21046 | 84 | 0.8571 | 0.0595 | 0.8571 | 12 | 0.666 |
| 21822 | 84 | 0.9643 | 0.0833 | 0.9643 | REF | 0.632 |
| 28081 | 84 | 0.9643 | 0.0476 | 0.9643 | REF | 0.665 |
| 29332 | 81 | 0.7778 | 0.0370 | 0.7778 | REF | 0.826 |
| 31955 | 84 | 0.8690 | 0.0952 | 0.8690 | 12 | 0.574 |

#### `BEST_result_per_case_RANKED.dat`

| Job | n | S_oracle | S_top1* | S_top10* |
|----:|--:|---------:|--------:|---------:|
| 21046 | 84 | 0.8571 | 0.0000 | 0.7381 |
| 21822 | 84 | 0.9643 | 0.0000 | 0.9405 |
| 31955 | 84 | 0.8690 | 0.0000 | 0.7857 |

#### `BEST_result_per_case_per_batch.dat`

| Job | minPts | n | S_oracle |
|----:|:-------|--:|---------:|
| 18852 | 10 | 36 | 0.8889 |
| 18852 | 15 | 27 | 0.8889 |
| 18852 | 22 | 20 | 0.8000 |
| 21046 | 12 | 84 | 0.8571 |
| 21822 | REF | 84 | 0.9643 |
| 31955 | 12 | 84 | 0.8690 |

### 5.5 Other compiled tables

| File | n_poses | Jobs |
|------|--------:|------|
| `ASTEX.dat` | 45719 | 18852, 21046, 2180, 21822, 31955, 3771, 9763 |
| `with_CFtot.dat` | 35317 | 18852, 2180, 21822, 31955, 3771, 9763 |
| `with_CFtot_minus31955.dat` | 33397 | 18852, 2180, 21822, 3771, 9763 |
| `scores_compiled.out` | 26106 | 18852, 2180, 21822, 3771, 9763 |

### 5.6 Per-job pose `.dat`

| File | n_poses | Notes / pooled success |
|------|--------:|------------------------|
| `8572.dat` | 840 | 8572: S1=0.012 S10=0.202 orc=0.202 |
| `21822.dat` | 8390 | 21822: S1=0.655 S10=0.964 orc=0.964 |
| `31955+21046.dat` | 12581 | no usable RMSD |

---

## 6. `notes_05-05017.dat` (FLRP / FLFP)

```
HAP2NN

FlexAID.FLRP 0.177777777777778   
FlexAIDdS.FLRP 0.222222222222222 
Vina.FLRP 0.133333333333333     
FlexX.FLRP 0.111111111111111     
rDock.FLRP 0.222222222222222     
FlexAID.FLFP 0.288888888888889  
FlexAIDdS.FLFP 0.355555555555556 
Vina.FLFP 0.333333333333333      
rDock.FLFP 0.2

ASTEXNN

FlexAID.FLRP 0.46   
FlexAIDdS.FLRP 0.52 
Vina.FLRP 0.44      
FlexX.FLRP 0.5     
rDock.FLRP 0.64     
FlexAID.FLFP 0.5    
FlexAIDdS.FLFP 0.54 
Vina.FLFP 0.5      
rDock.FLFP 0.74
```

---

## 7. Root inventory

Total root entries: **194**.

### 7.1 Job directories (shallow)

| Dir | Children | Notes |
|-----|---------:|-------|
| `14656` | 1108 | |
| `17826` | 76 | |
| `18851_00078_3` | 1 | |
| `19786` | 1108 | |
| `20087` | 1108 | |
| `2109` | 84 | |
| `21822` | 84 | |
| `22589` | 84 | |
| `24091` | 84 | |
| `26896` | 1108 | |
| `27836` | 674 | |
| `9312` | 76 | |

Broken legacy symlinks (not followed): `23178`, `24348`, `31470`, `cluster_NRG_docking_results.pl` → `/Users/lmorency/...`.

### 7.2 Root numeric job IDs

Count: **81** — 2109, 2110, 2111, 2180, 2427, 2428, 2429, 3771, 5148, 5296, 5297, 5298, 8137, 8138, 8139, 8435, 8572, 8643, 9156, 9312, 9313, 9314, 9566, 9710, 9897, 9898, 9899, 10902, 10903, 12109, 12717, 13796, 14656, 14657, 14658, 14992, 15407, 15658, 16477, 17826, 18417, 18851, 18852, 18949, 19786, 19788, 19789, 19790, 20087, 20120, 21046, 21243, 21244, 21245, 21763, 21822, 22589, 23178, 23979, 24091, 24348, 24508, 25414, 25593, 25594, 25595, 26896, 27836, 28081, 28518, 28519, 28520, 29221, 29222, 29332, 30417, 31470, 31578, 31955, 95407, 292211

Not in `jobsID.info` (70): 2109, 2110, 2111, 2427, 2428, 2429, 5148, 5296, 5297, 5298, 8137, 8138, 8139, 8435, 8643, 9312, 9313, 9314, 9566, 9710, 9897, 9898, 9899, 10902, 10903, 12109, 12717, 13796, 14656, 14657, 14658, 14992, 15407, 15658, 16477, 17826, 18417, 18851, 18949, 19786, 19788, 19789, 19790, 20087, 20120, 21243, 21244, 21245, 21763, 22589, 23178, 23979, 24091, 24348, 24508, 25414, 25593, 25594, 25595, 26896, 27836, 28518, 28519, 28520, 29221, 29222, 30417, 31470, 95407, 292211

These extras are largely later BOINC/NRG docking batches (2016–2018) with `*.out` / `*_RMSD.lst` (often triplets e.g. 14656/14657/14658).

### 7.3 Root analysis / non-job files (selected)

| Name | Kind | Size |
|------|------|-----:|
| `.RData` | file | 3.5 MB |
| `.Rhistory` | file | 2.8 KB |
| `10903_RMSD.lst` | file | 504.9 KB |
| `12109_RMSD.lst` | file | 502.3 KB |
| `12717_RMSD.lst` | file | 530.0 KB |
| `13796_RMSD.lst` | file | 494.2 KB |
| `14656_RMSD.lst` | file | 13.0 MB |
| `14657_RMSD.lst` | file | 12.7 MB |
| `14658_RMSD.lst` | file | 12.6 MB |
| `14992_RMSD.lst` | file | 534.1 KB |
| `15407_RMSD.lst` | file | 468.2 KB |
| `15658_RMSD.lst` | file | 462.6 KB |
| `16477_RMSD.lst` | file | 507.9 KB |
| `17826_RMSD.lst` | file | 543.6 KB |
| `18417_RMSD.lst` | file | 526.6 KB |
| `18852_RMSD.lst` | file | 975.6 KB |
| `18949_RMSD.lst` | file | 494.1 KB |
| `19786_RMSD.lst` | file | 6.5 MB |
| `19788_RMSD.lst` | file | 13.0 MB |
| `19789_RMSD.lst` | file | 12.7 MB |
| `19790_RMSD.lst` | file | 12.7 MB |
| `1stp.cif` | file | 133.4 KB |
| `1stp.pse` | file | 270.3 KB |
| `20087_RMSD.lst` | file | 6.5 MB |
| `20120_RMSD.lst` | file | 402.8 KB |
| `21046_RMSD.lst` | file | 343.2 KB |
| `2109_RMSD.lst` | file | 1010.8 KB |
| `2110_RMSD.lst` | file | 983.0 KB |
| `2111_RMSD.lst` | file | 979.9 KB |
| `21243_RMSD.lst` | file | 12.9 MB |
| `21244_RMSD.lst` | file | 12.6 MB |
| `21245_RMSD.lst` | file | 12.6 MB |
| `21763_RMSD.lst` | file | 448.5 KB |
| `2180_RMSD.lst` | file | 225.0 KB |
| `21822_RMSD.lst` | file | 532.2 KB |
| `22589_RMSD.lst` | file | 505.4 KB |
| `23178_RMSD.lst` | file | 541.8 KB |
| `23979_RMSD.lst` | file | 476.4 KB |
| `24091_RMSD.lst` | file | 534.9 KB |
| `2427_RMSD.lst` | file | 13.0 MB |
| `2428_RMSD.lst` | file | 12.7 MB |
| `2429_RMSD.lst` | file | 12.7 MB |
| `24348_RMSD.lst` | file | 544.4 KB |
| `24508_RMSD.lst` | file | 466.6 KB |
| `25414_RMSD.lst` | symlink | 33 B |
| `25593_RMSD.lst` | file | 6.5 MB |
| `25594_RMSD.lst` | file | 6.4 MB |
| `25595_RMSD.lst` | file | 6.4 MB |
| `26896_RMSD.lst` | file | 6.5 MB |
| `27836_RMSD` | file | 499.7 KB |
| `27836_RMSD.lst` | file | 499.7 KB |
| `28081_RMSD.lst` | file | 495.5 KB |
| `28518_RMSD.lst` | file | 13.0 MB |
| `28519_RMSD.lst` | file | 12.7 MB |
| `28520_RMSD.lst` | file | 12.6 MB |
| `29221_RMSD.lst` | file | 427.4 KB |
| `29222_RMSD.lst` | file | 782.1 KB |
| `29332_RMSD.lst` | file | 497.9 KB |
| `30417_RMSD.lst` | file | 456.7 KB |
| `31470_RMSD.lst` | file | 6.4 MB |
| `31578_RMSD.lst` | file | 32.7 KB |
| `3771_RMSD.lst` | file | 440.4 KB |
| `5148_RMSD.lst` | file | 535.1 KB |
| `5296_RMSD.lst` | file | 12.9 MB |
| `5297_RMSD.lst` | file | 12.7 MB |
| `5298_RMSD.lst` | file | 12.6 MB |
| `8137_RMSD.lst` | file | 1.1 MB |
| `8138_RMSD.lst` | file | 1.0 MB |
| `8139_RMSD.lst` | file | 1.0 MB |
| `8435_RMSD.lst` | file | 499.2 KB |
| `8643_RMSD.lst` | file | 480.1 KB |
| `9156_RMSD.lst` | file | 36.2 KB |
| `9312_RMSD.lst` | file | 1.1 MB |
| `9313_RMSD.lst` | file | 1.0 MB |
| `9314_RMSD.lst` | file | 1.0 MB |
| `95407_RMSD.lst` | file | 466.6 KB |
| `9566_RMSD.lst` | file | 480.6 KB |
| `9710_RMSD.lst` | file | 457.7 KB |
| `9897_RMSD.lst` | file | 1.3 MB |
| `9898_RMSD.lst` | file | 1.3 MB |
| `9899_RMSD.lst` | file | 1.2 MB |
| `ALL_astex.dat` | file | 3.6 MB |
| `ASTEX.dat` | file | 2.5 MB |
| `ASTEX_bootstrap.out` | file | 48.0 KB |
| `ASTEX_final.dat` | file | 3.1 MB |
| `BEST_result_per_case.dat` | file | 28.0 KB |
| `BEST_result_per_case_RANKED.dat` | file | 13.7 KB |
| `BEST_result_per_case_per_batch.dat` | file | 17.6 KB |
| `DockingDatasetBenchmarks.py` | file | 577 B |
| `FlexAID_CF+FO_results.dat` | file | 883.4 KB |
| `RMSD.dat` | file | 532.2 KB |
| `RMSD.lst` | file | 0 B |
| `R_packages` | dir | 1.2 KB |
| `Untitled.ipynb` | file | 72 B |
| `analyze_entropy_advantages.py` | file | 851 B |
| `analyze_poses.pl` | file | 3.5 KB |
| `analyze_poses_quickfix.pl` | file | 3.5 KB |
| `cluster_NRG_docking_results.pl` | symlink | 68 B |
| `compile_RMSD.pl` | file | 4.7 KB |
| `concatenate_results.pl` | file | 2.8 KB |
| `cys.png` | file | 42.8 KB |
| `format_dat_to_RMSDlst.pl` | file | 3.7 KB |
| `get_boinc_RMSD.pl` | file | 1.2 KB |
| `get_ligands.pl` | file | 1.0 KB |
| `jobsID.info` | file | 1.2 KB |
| `loop_compile_RMSD.pl` | file | 751 B |
| `met.png` | file | 34.8 KB |
| `notes_05-05017.dat` | file | 485 B |
| `oxy.png` | file | 159.4 KB |
| `parse_FlexAID_CF+CO_minPoints.pl` | file | 1.7 KB |
| `parse_FlexAID_CF+CO_repeat.pl` | file | 1.6 KB |
| `parse_FlexAID_CF+CO_results.pl` | file | 2.8 KB |
| `parse_NRG_dockings_results.pl` | file | 6.1 KB |
| `parse_NRG_dockings_results_mRMSD.pl` | file | 9.5 KB |
| `peptine.png` | file | 43.2 KB |
| `rerank_by_CF.pl` | file | 394 B |
| `scores_compiled.out` | file | 1.1 MB |
| `ser.png` | file | 39.0 KB |
| `with_CFtot.dat` | file | 1.9 MB |
| `with_CFtot_minus31955.dat` | file | 1.8 MB |

### 7.4 Scripts at root

- `DockingDatasetBenchmarks.py`
- `analyze_entropy_advantages.py`
- `analyze_poses.pl`
- `analyze_poses_quickfix.pl`
- `cluster_NRG_docking_results.pl`
- `compile_RMSD.pl`
- `concatenate_results.pl`
- `format_dat_to_RMSDlst.pl`
- `get_boinc_RMSD.pl`
- `get_ligands.pl`
- `loop_compile_RMSD.pl`
- `parse_FlexAID_CF+CO_minPoints.pl`
- `parse_FlexAID_CF+CO_repeat.pl`
- `parse_FlexAID_CF+CO_results.pl`
- `parse_NRG_dockings_results.pl`
- `parse_NRG_dockings_results_mRMSD.pl`
- `rerank_by_CF.pl`

---

## 8. `R_packages/`

| Entry | Kind | Detail |
|-------|------|--------|
| `.RData` | file | 7.9 MB |
| `.Rhistory` | file | 24.2 KB |
| `2xa0.pdb` | file | 262.1 KB |
| `R_crit_sc.lst` | file | 22.7 KB |
| `all.bt` | file | 32.5 KB |
| `all.lt` | file | 32.1 KB |
| `all.mat` | file | 1.1 MB |
| `all.pt` | file | 33.7 KB |
| `astex_rmsd_bb.lst` | file | 44.4 KB |
| `astex_w50_crit.lst` | file | 16.6 KB |
| `density_rescore.lst` | file | 315.8 KB |
| `diff.lst` | file | 280.9 KB |
| `hap2_clashing_residues_good5.lst` | file | 91.4 KB |
| `hap2_critical_residues_good5.lst` | file | 1.3 KB |
| `hap2_ligflex.lst` | file | 438 B |
| `hap2_nn_pdbids_good5.lst` | file | 760 B |
| `hydrophobicity_index_pocket.lst` | file | 2.9 KB |
| `matchres` | file | 645.7 KB |
| `matchrot` | file | 681.5 KB |
| `matchwall` | file | 25.0 MB |
| `nflex_bonds.lst` | file | 596 B |
| `plot_flexaid_barplot.R` | file | 5.4 KB |
| `plot_flexaid_corrflex.R` | file | 1.5 KB |
| `plot_flexaid_corrhyd.R` | file | 3.3 KB |
| `plot_flexaid_interactions.R` | file | 9.6 KB |
| `plot_flexaid_params.R` | file | 19.0 KB |
| `plot_flexaid_plotdata.R` | file | 7.1 KB |
| `plot_flexaid_plotrank.R` | file | 3.8 KB |
| `plot_flexaid_plotref.R` | file | 4.9 KB |
| `plot_flexaid_readdata.R` | file | 17.8 KB |
| `solvent_exp.lst` | file | 891 B |
| `summary/` | dir | ~1400 RMSD.lst files — not walked |
| `tiff/` | dir | TIFFs: HAP2NN_TOP1.tiff, ASTEXNN_TOP1.tiff, ASTEXNN_TOP10.tiff, ASTEXNN_TOP5.tiff, HAP2NN_TOP5.tiff, ASTEXNN_TOP3.tiff, HAP2NN_TOP3.tiff, HAP2NN_TOP10.tiff |
| `training_st0r5.2.surfaces` | file | 885.5 MB |
| `vectors.R` | file | 30.0 KB |

R plot scripts: `plot_flexaid_barplot.R`, `plot_flexaid_corrflex.R`, `plot_flexaid_corrhyd.R`, `plot_flexaid_interactions.R`, `plot_flexaid_params.R`, `plot_flexaid_plotdata.R`, `plot_flexaid_plotrank.R`, `plot_flexaid_plotref.R`, `plot_flexaid_readdata.R`, `vectors.R`.

TIFF barplots: `ASTEXNN_TOP{1,3,5,10}.tiff`, `HAP2NN_TOP{1,3,5,10}.tiff`.

Large: `training_st0r5.2.surfaces` (~886 MB) — not content-parsed.

---

## 9. Pose-line schema (from R2 perl)

```
jobID  minPoints|REF  case_rN  R_k  RMSD  CF  [totCF]  [freq]
```

- `minPoints == REF` ⇒ **CF** clustering; numeric minPoints ⇒ **FastOPTICS (FO)**
- Rank `R_0` = top-1 for that (case, repeat, minPts)
- `compile_RMSD.pl` writes `JobID PDB Ligand RunID ResID RMSD CF TCF ACF freq` into `*_RMSD.lst`
- `concatenate_results.pl` builds `BEST_result_per_case.dat` + `ASTEX_bootstrap.out`

---

## 10. Interpretation notes (FlexAIDdS / 3DSig lineage)

1. **CF jobs 21822 (TEMP=0) and 28081 (TEMP=13, MC_st0r5.2)** both reach **~96.4% oracle** and **~96.4% S_top10** — strongest CF cluster results in this tree. S_top1 is lower (0.65 / 0.56), so ranking within the top-10 list matters.
2. **FO 18852** (minPts 10/15/22, 10 reps) is the multi-parameter FastOPTICS workhorse in `ASTEX_final`; pooled S_top1 **0.831**, oracle **0.867**, bootstrap median **0.893**.
3. **29332** (`NRG_mat_BEST_13`) underperforms other CF jobs (**~78.6%** bootstrap / oracle).
4. **21046 / 31955** FO with CF < 0 filter; oracle **0.857 / 0.869**.
5. Early FO/CF jobs (3771, 2180, 9156, 8572, 9763, 30281, 31578) appear mainly in `ALL_astex` / `FlexAID_CF+FO` with lower success and often incomplete case coverage (n < 84).
6. `notes_05-05017.dat` records **FLRP/FLFP** for FlexAID vs FlexAIDdS vs Vina/FlexX/rDock on HAP2NN and ASTEXNN — separate from S_topK tables.
7. Post-2016 root jobs (14656, 19788, 21243, 2427, 5296, 28518, …) are bulk NRG/BOINC exports; not in `jobsID.info` or bootstrap.

---

## 11. Full root inventory

| Kind | Size (bytes) | Name | Extra |
|------|-------------:|------|-------|
| file | 6148 | `.DS_Store` |  |
| file | 6148 | `.DS_Store 2` |  |
| file | 3621916 | `.RData` |  |
| file | 2901 | `.Rhistory` |  |
| dir | 96 | `.ipynb_checkpoints` | n_children=1 |
| file | 452769 | `10902.out` |  |
| file | 456275 | `10903.out` |  |
| file | 516974 | `10903_RMSD.lst` |  |
| file | 452682 | `12109.out` |  |
| file | 453744 | `12109_CE.out` |  |
| file | 452682 | `12109_CF.out` |  |
| file | 514357 | `12109_RMSD.lst` |  |
| file | 510551 | `12109_RNSD.lst` |  |
| file | 489422 | `12717.out` |  |
| file | 542694 | `12717_RMSD.lst` |  |
| file | 448074 | `13796.out` |  |
| file | 506086 | `13796_RMSD.lst` |  |
| dir | 35520 | `14656` | n_children=1108 |
| file | 14963454 | `14656.out` |  |
| file | 13595319 | `14656_RMSD.lst` |  |
| file | 13298617 | `14657_RMSD.lst` |  |
| file | 13244740 | `14658_RMSD.lst` |  |
| file | 493663 | `14992.out` |  |
| file | 546881 | `14992_RMSD.lst` |  |
| file | 423253 | `15407.out` |  |
| file | 479435 | `15407_RMSD.lst` |  |
| file | 418290 | `15658.out` |  |
| file | 473683 | `15658_RMSD.lst` |  |
| file | 459551 | `16477.out` |  |
| file | 520112 | `16477_RMSD.lst` |  |
| dir | 2496 | `17826` | n_children=76 |
| file | 410581 | `17826.out` |  |
| file | 556658 | `17826_RMSD.lst` |  |
| file | 486253 | `18417.out` |  |
| file | 539206 | `18417_RMSD.lst` |  |
| dir | 96 | `18851_00078_3` | n_children=1 |
| file | 119345 | `18852.out` |  |
| file | 999042 | `18852_RMSD.lst` |  |
| file | 445805 | `18949.out` |  |
| file | 505913 | `18949_RMSD.lst` |  |
| dir | 35520 | `19786` | n_children=1108 |
| file | 5984275 | `19786.out` |  |
| file | 6784620 | `19786_RMSD.lst` |  |
| file | 15366829 | `19788.out` |  |
| file | 13611231 | `19788_RMSD.lst` |  |
| file | 13338705 | `19789_RMSD.lst` |  |
| file | 13333640 | `19790_RMSD.lst` |  |
| file | 136568 | `1stp.cif` |  |
| file | 276815 | `1stp.pse` |  |
| dir | 35520 | `20087` | n_children=1108 |
| file | 5992893 | `20087.out` |  |
| file | 6793384 | `20087_RMSD.lst` |  |
| file | 364220 | `20120.out` |  |
| file | 412431 | `20120_RMSD.lst` |  |
| file | 191622 | `20120new.out` |  |
| file | 351415 | `21046_RMSD.lst` |  |
| dir | 2752 | `2109` | n_children=84 |
| file | 1115941 | `2109.out` |  |
| file | 1035054 | `2109_RMSD.lst` |  |
| file | 1006602 | `2110_RMSD.lst` |  |
| file | 1003402 | `2111_RMSD.lst` |  |
| file | 14591648 | `21243.out` |  |
| file | 13577042 | `21243_RMSD.lst` |  |
| file | 13242028 | `21244_RMSD.lst` |  |
| file | 13180142 | `21245_RMSD.lst` |  |
| file | 395084 | `21763.out` |  |
| file | 459276 | `21763_RMSD.lst` |  |
| file | 230386 | `2180_RMSD.lst` |  |
| dir | 2752 | `21822` | n_children=84 |
| file | 533752 | `21822.dat` |  |
| file | 544995 | `21822_RMSD.lst` |  |
| dir | 2752 | `22589` | n_children=84 |
| file | 356097 | `22589.dat` |  |
| file | 457502 | `22589.out` |  |
| file | 517491 | `22589_RMSD.lst` |  |
| symlink | 57 | `23178` | /Users/lmorency/Documents/PhD/FlexAID_DevFiles/jobs/23178 |
| file | 408687 | `23178.out` |  |
| file | 554821 | `23178_RMSD.lst` |  |
| file | 487800 | `23979_RMSD.lst` |  |
| dir | 2752 | `24091` | n_children=84 |
| file | 335289 | `24091.dat` |  |
| file | 494531 | `24091.out` |  |
| file | 547714 | `24091_RMSD.lst` |  |
| file | 14881312 | `2427.out` |  |
| file | 13597856 | `2427_RMSD.lst` |  |
| file | 13313689 | `2428_RMSD.lst` |  |
| file | 13271451 | `2429_RMSD.lst` |  |
| symlink | 57 | `24348` | /Users/lmorency/Documents/PhD/FlexAID_DevFiles/jobs/24348 |
| file | 411409 | `24348.out` |  |
| file | 557505 | `24348_RMSD.lst` |  |
| file | 422423 | `24508.out` |  |
| file | 477778 | `24508_RMSD.lst` |  |
| symlink | 33 | `25414_RMSD.lst` | R_packages/summary/25414_RMSD.lst |
| file | 7471568 | `25593.out` |  |
| file | 6820975 | `25593_RMSD.lst` |  |
| file | 6668871 | `25594_RMSD.lst` |  |
| file | 6661394 | `25595_RMSD.lst` |  |
| dir | 35520 | `26896` | n_children=1108 |
| file | 6007865 | `26896.out` |  |
| file | 6808131 | `26896_RMSD.lst` |  |
| dir | 21632 | `27836` | n_children=674 |
| file | 335137 | `27836.dat` |  |
| file | 451116 | `27836.out` |  |
| file | 511691 | `27836_RMSD` |  |
| file | 511741 | `27836_RMSD.lst` |  |
| file | 507379 | `28081_RMSD.lst` |  |
| file | 14773211 | `28518.out` |  |
| file | 13585176 | `28518_RMSD.lst` |  |
| file | 13271415 | `28519_RMSD.lst` |  |
| file | 13209127 | `28520_RMSD.lst` |  |
| file | 386431 | `29221.out` |  |
| file | 317929 | `292211.out` |  |
| file | 437673 | `29221_RMSD.lst` |  |
| file | 707223 | `29222.out` |  |
| file | 800883 | `29222_RMSD.lst` |  |
| file | 509823 | `29332_RMSD.lst` |  |
| file | 412973 | `30417.out` |  |
| file | 467625 | `30417_RMSD.lst` |  |
| symlink | 57 | `31470` | /Users/lmorency/Documents/PhD/FlexAID_DevFiles/jobs/31470 |
| file | 5929624 | `31470.out` |  |
| file | 6687353 | `31470_RMSD.lst` |  |
| file | 33495 | `31578_RMSD.lst` |  |
| file | 663723 | `31955+21046.dat` |  |
| file | 389806 | `3771.out` |  |
| file | 450934 | `3771_RMSD.lst` |  |
| file | 486312 | `5148.out` |  |
| file | 547985 | `5148_RMSD.lst` |  |
| file | 14668691 | `5296.out` |  |
| file | 13566787 | `5296_RMSD.lst` |  |
| file | 13266918 | `5297_RMSD.lst` |  |
| file | 13211677 | `5298_RMSD.lst` |  |
| file | 1012499 | `8137.out` |  |
| file | 1111453 | `8137_RMSD.lst` |  |
| file | 1091753 | `8138_RMSD.lst` |  |
| file | 1088733 | `8139_RMSD.lst` |  |
| file | 444456 | `8435.out` |  |
| file | 511195 | `8435_RMSD.lst` |  |
| file | 41111 | `8572.dat` |  |
| file | 426456 | `8643.out` |  |
| file | 491592 | `8643_RMSD.lst` |  |
| file | 37029 | `9156_RMSD.lst` |  |
| dir | 2496 | `9312` | n_children=76 |
| file | 991913 | `9312.out` |  |
| file | 1110840 | `9312_RMSD.lst` |  |
| file | 1089358 | `9313_RMSD.lst` |  |
| file | 1085161 | `9314_RMSD.lst` |  |
| file | 477749 | `95407_RMSD.lst` |  |
| file | 428358 | `9566.out` |  |
| file | 492161 | `9566_RMSD.lst` |  |
| file | 404631 | `9710.out` |  |
| file | 468666 | `9710_RMSD.lst` |  |
| file | 1419195 | `9897.out` |  |
| file | 1346448 | `9897_RMSD.lst` |  |
| file | 1313134 | `9898_RMSD.lst` |  |
| file | 1307379 | `9899_RMSD.lst` |  |
| file | 3739140 | `ALL_astex.dat` |  |
| file | 2606041 | `ASTEX.dat` |  |
| file | 49152 | `ASTEX_bootstrap.out` |  |
| file | 3229130 | `ASTEX_final.dat` |  |
| file | 28688 | `BEST_result_per_case.dat` |  |
| file | 14000 | `BEST_result_per_case_RANKED.dat` |  |
| file | 17989 | `BEST_result_per_case_per_batch.dat` |  |
| file | 577 | `DockingDatasetBenchmarks.py` |  |
| file | 904575 | `FlexAID_CF+FO_results.dat` |  |
| file | 544945 | `RMSD.dat` |  |
| file | 0 | `RMSD.lst` |  |
| dir | 1184 | `R_packages` | n_children=35 |
| file | 72 | `Untitled.ipynb` |  |
| file | 851 | `analyze_entropy_advantages.py` |  |
| file | 3608 | `analyze_poses.pl` |  |
| file | 3626 | `analyze_poses_quickfix.pl` |  |
| symlink | 68 | `cluster_NRG_docking_results.pl` | /Users/lmorency/Documents/PhD/Scripts/cluster_NRG_docking_results.pl |
| file | 4770 | `compile_RMSD.pl` |  |
| file | 2817 | `concatenate_results.pl` |  |
| file | 43806 | `cys.png` |  |
| file | 3799 | `format_dat_to_RMSDlst.pl` |  |
| file | 1198 | `get_boinc_RMSD.pl` |  |
| file | 1067 | `get_ligands.pl` |  |
| file | 1247 | `jobsID.info` |  |
| file | 751 | `loop_compile_RMSD.pl` |  |
| file | 35650 | `met.png` |  |
| file | 485 | `notes_05-05017.dat` |  |
| file | 163209 | `oxy.png` |  |
| file | 1738 | `parse_FlexAID_CF+CO_minPoints.pl` |  |
| file | 1651 | `parse_FlexAID_CF+CO_repeat.pl` |  |
| file | 2842 | `parse_FlexAID_CF+CO_results.pl` |  |
| file | 6200 | `parse_NRG_dockings_results.pl` |  |
| file | 9745 | `parse_NRG_dockings_results_mRMSD.pl` |  |
| file | 44280 | `peptine.png` |  |
| file | 394 | `rerank_by_CF.pl` |  |
| file | 1189918 | `scores_compiled.out` |  |
| file | 39931 | `ser.png` |  |
| file | 1984751 | `with_CFtot.dat` |  |
| file | 1871406 | `with_CFtot_minus31955.dat` |  |

*End of catalog. R2 originals untouched.*
