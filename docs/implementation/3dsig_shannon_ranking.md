# 3Dsig / Morency 2017 Shannon ranking in DatasetRunner

**Reference:** L.-P. Morency, *The Impact of Conformational Entropy on the Accuracy of FlexAID in Binding Mode Prediction*, ISMB/ECCB 2017 — 3Dsig.

## Formula (must match poster)

\[
Z = \sum_i e^{-\mathrm{CF}_i / T},\quad
p_i = \frac{e^{-\mathrm{CF}_i / T}}{Z}
\]

\[
\tilde H = \sum_i p_i\,\mathrm{CF}_i,\quad
\tilde S = -\sum_i p_i\ln p_i,\quad
\tilde G = \tilde H - T\,\tilde S
\]

Elect **lowest** \(\tilde G\). CF is a scoring **proxy** (a.u.); soft-β uses \(T\) in K (\(\beta=1/T\)), not \(k_B\).

## Where it is implemented

| Layer | Behavior |
|--------|----------|
| `cluster.cpp` ACF | \(\mathrm{ACF}=E_{\min}-T\ln Z_{\mathrm{loc}}\) — emission order when \(T>0\) |
| `BindingMode::compute_energy` | \(\tilde G=H-TS\) for FO mode sort when classic entropy ranking |
| **`DatasetRunner` S1 election** | **Same \(\tilde G\)** over cluster heads (+ `.mcf` members). Log: `[3DSIG-RANK]` |

## Env

| Variable | Default | Meaning |
|----------|---------|---------|
| `FLEXAIDDS_ELECTION_SHANNON_F` | **1 (ON)** | Use \(\tilde G=H-TS\) for S1 |
| `FLEXAIDDS_ELECTION_LEGACY_ZH` | 0 | Rollback to Z+H composite (≈ min-CF) |
| `FLEXAIDDS_ELECTION_SOFT_T` | 0 → **298** | Soft-β \(T\) in K |

Live claim binaries must be **rebuilt and restaged** to pick this up; running claim PIDs keep the old election until restart with new binary.
